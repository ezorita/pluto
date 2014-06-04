#include "plutoaligner.h"

void SIGSEGV_handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(EXIT_FAILURE);
}

int
main
(
  int argc,
  char *argv[]
)
{
   // Get a sequence and align it, return the loci list.
   if (argc != 4) {
      fprintf(stderr, "usage: pluto <tau> <threads> <queryfile> <.pif file>\n");
      exit(1);
   }

  // Backtrace handler
   signal(SIGSEGV, SIGSEGV_handler); 

   // Check seq len.
   if(strlen(argv[1]) != 14) {
      fprintf(stderr, "Error: Only seqs of len 14 are supported.\nusage: query <sequence> <.pif file>\n");
      exit(EXIT_FAILURE);
   }
   
   // *************
   // Input parsing
   // *************

   int    tau = atoi(argv[1]);
   int    maxthreads = atoi(argv[2]);
   char * queryfile = argv[3];
   char * indexfile = argv[4];

   // DEBUG:
   fprintf(stderr, "load params: tau=%d\tseq=%s\tfilename=%s\n", tau, seq, indexfile);
   
   // *********************
   // Load index/data files
   // *********************

   int fdi = open(indexfile, O_RDONLY);

   // Get file size.
   unsigned long isize = lseek(fdi, 0, SEEK_END);
   lseek(fdi, 0, SEEK_SET);
   
   // Map LUT and index to memory.
   loc_t *lut = (uint *) mmap(NULL, isize, PROT_READ, MAP_SHARED, fdi, 0);
   loc_t *index = lut + NSEQ;

   if (lut == MAP_FAILED) {
      fprintf(stderr, "error loading index (mmap): %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   // Read file.
   uint    num_query;
   char ** all_query = read_file(queryfile, &num_query);
   mergesort((void **) all_query, num_query, &ualpha, maxthreads);


   // ******************
   // Allocate variables
   // ******************

   // Alloc mstack and lstack buffers based on the length of the first seq.
   int nnodes = 2*strlen(all_query[0]/SEQLEN)-1;
   lstackbuf_t * lstacks = malloc(sizeof(lstackbuf_t));
   mstackbuf_t * mstacks = malloc(sizeof(mstackbuf_t));

   lstacks->stack = (lstack_t **) malloc(nnodes * sizeof(lstack_t *));
   mstacks->stack = (mstack_t **) malloc(nnodes * sizeof(mstack_t *));
   
   for (int i = 0; i < nnodes; i++) {
      lstacks->stack[i] = new_lstack(LOCI_STACKSIZE);
      mstacks->stack[i] = new_mstack(MISM_STACKSIZE);
   }
   
   lstacks->count = nnodes;
   mstacks->count = nnodes;

   // Alloc an empty merge-tree with one node:
   tree_t * tree = malloc(sizeof(tree_t));
   tree->node   = malloc(sizeof(node_t));
   tree->node[0].data   = (lstack_t **) malloc((tau+1)*sizeof(lstack_t *));
   tree->node[0].mstack = (mstack_t **) malloc((tau+1)*sizeof(mstack_t *));
   tree->nnodes = 1;

   // Struct arg_t:
   struct arg_t args = {
      .maxtau = tau,
      .lut    = lut,
      .index  = index,
   };


   // *****************
   // Mapping algorithm
   // *****************
   
   for (int q = 0; q < num_query; q++) {
      if (all_query[q] == NULL) continue;

      // Split in chunks.
      char * query    = all_query[q];
      int    querylen = strlen(query);
      int    nleaves  = querylen/SEQLEN;
      seq_t  seqid[nchunks];

      for (int c = 0; c < nchunks; c++)
         seqid[c]  = seqtoid(query + c*SEQLEN, SEQLEN);
      
      int   a            = 0;
      loc_t merge_output = MERGE_EMPTY;
      
      // Build tree.
      node_t * root = build_tree(nleaves, seqid, tree, lstacks, mstacks, tau);

      // Set up arg struct.
      args.nleaves = nleaves;
      args.query   = query;

      // Iterate merge algorithm.
      while (merge_output == MERGE_EMPTY && a <= tau) {
         merge_output = merge_node(root, a, &args);
         a++;
      }
      
      // Restore the last used tau.
      a--;
      
      // Vomit loci:
      fprintf(stdout, ">%s\n", all_query[q]);
      if (merge_output == MERGE_EMPTY) {
         fprintf(stdout,"\tsequence not found.\n");
      } else {
         lstack_t * loci = root->data[a];
         for (int i = 0; i < loci->pos; i++) {
            fprintf(stdout, "\t%lu\n", loci->u[i]);
         }
      }
   }
}



tnode_t *
build_tree
(
 int             nleaves,
 seq_t         * leaves,
 tree_t        * tree,
 lstackbuf_t   * lstacks,
 mstackbuf_t   * mstacks,
 char            tau
)
// First reshape the tree. Compute the number of inner nodes and leaves
//  # of node for n leaves: 2n - 1
//  # of inner nodes: n - 1
// and make the new structure creating/freeing as many nodes as needed.
//
// TODO: This algorithm is quite crappy and super-exhaustive. May be improved but the 
//       overall performance gain will most likely be zero.
// TODO (improve): During reshaping check if any previously merged stack can be reused...
//
// Then assign stacks to each leave:
// 0. Resize the stack buffer according to the num of nodes and alloc as many new stacks
//    as necessary. (Only grow!)
// 1. Go through stack buf and change MERGE_PENDING -> MERGE_AVAILABLE (NOT TO SEQS!!).
// 2. For each seq look at the stack node if it has been precomputed and assign pointer.
//    If not, just assign the first MERGE_AVAILABLE and switch to MERGE_PENDING.
// 3. For the rest of the inner nodes also assign the first available buffer.
{
   // Alloc nodes.
   int nnodes = 2*nleaves - 1;
   if (nnodes > tree->nnodes) {
      tree->node = realloc(tree->node, nnodes * sizeof(tnode_t));
      for (int i = tree->nnodes; i < nnodes; i++) {
         tree->node[i].data   = (lstack_t **) malloc((tau+1)*sizeof(lstack_t *));
         tree->node[i].mstack = (mstack_t **) malloc((tau+1)*sizeof(mstack_t *));
      }
      tree->nnodes = nnodes;
   }

   // Connect nodes.
   int firstidx = 0;
   int nextidx  = 1;
   int level    = 1;
   for (int i = 0; i < nnodes; i++) {
      tnode_t * node = tree->node + i;
      node->mintau = 0;

      // Node.
      if (nextidx + 2*(i - firstidx) < nnodes) {
         node->leaf   = -1;
         node->status = NODE_OPEN; // NODE_OPEN, NODE_SET.
         node->lchild = tree->node + nextidx + 2*(i - firstidx);
         node->rchild = node->lchild + 1;
      } 
      // Leaf.
      else {
         node->leaf   = 1;
         node->status = NODE_OPEN;
         node->lchild = node->rchild = NULL;
      }

      // Next level.
      if (i+1 < nnodes && i+1 == nextidx) {
         firstidx = nextidx;
         nextidx += (1 << level++);
      }
   }

   // Realloc stack buffers.
   // Loci stacks and Mismatch stacks share the same indices (they are taken by pairs).
   if ((tau+1) * nnodes > lstacks->count) {
      // Realloc loci and mismatch stack buffers.
      lstacks->stack = (lstack_t **) realloc(lstacks->stack, nnodes * (tau+1) * sizeof(lstack_t *));
      mstacks->stack = (mstack_t **) realloc(mstacks->stack, nnodes * (tau+1) * sizeof(mstack_t *));
      if (lstacks->stack == NULL || mstacks->stack == NULL) {
         fprintf(stderr, "Buffer stack realloc error (build_tree): %s\n", strerror(errno));
         exit(EXIT_FAILURE);
      }
      
      // Initialize new loci stacks.
      for (int i = lstacks->count; i < (tau+1)*nnodes; i++) {
         lstacks->stack[i] = new_lstack(LOCI_STACKSIZE);
         mstacks->stack[i] = new_mstack(MISM_STACKSIZE);
      }

      // Update buffer counts.
      lstacks->count = (tau+1)*nnodes;
      mstacks->count = (tau+1)*nnodes;

   }

   // Find repeats or reset stack otherwise.
   for (int i = 0; i < (tau+1)*nnodes; i++) {
      if (lstacks->stack[i]->seq > SEQMASK) {
         lstacks->stack[i]->seq = MERGE_DONE;
         lstacks->stack[i]->pos = 0;
      }
      else {
         int repeat = 0;
          for (int j = 0; j < nleaves; j++) {
             if (lstacks->stack[i]->seq == leaves[j]) {
                repeat = 1;
                break;
            }
         }
         if (!repeat) {
            lstacks->stack[i]->seq = MERGE_DONE;
            lstacks->stack[i]->pos = 0;
         }
      }
   }

   // Go through all leaf nodes by order. Assign stacks.
   int leafcnt = 0;
   int i = firstidx;
   while (leafcnt < nleaves) {
      tnode_t * leaf = tree->node + i++;
      if (leaf->leaf < 0) continue;
      leaf->leaf   = leafcnt;
      leaf->status = leaves[leafcnt++];

      // Look for precomputed loci and mismatch stacks.
      for (int a = 0; a <= tau; a++) {
         int found = 0;
         for (int j = 0; j < (tau+1)*nnodes; j++) {
            lstack_t * lstck = lstacks->stack[j];
            if (lstck->seq == leaf->status && lstck->tau == a) {
               leaf->data[a] = lstck;
               leaf->mstack[a] = mstacks->stack[j];
               found = 1;
               break;
            }
         }

         // If not found, just assign the first available.
         if (found == 0) {
            for (int j = 0; j < (tau+1)*nnodes; j++) {
               if (lstacks->stack[j]->seq == MERGE_DONE) {
                  leaf->data[a] = lstacks->stack[j];
                  leaf->data[a]->seq = leaf->status;
                  leaf->data[a]->tau = a;
                  leaf->mstack[a] = mstacks->stack[j];
                  break;
               }
            }
         }
      }
      
      // When done whith the last level, continue at the previous.
      if (i == nnodes) i = firstidx - (1 << level);
   }

   // Assign stacks to the inner nodes.
   for (int k = 0, j = 0; k < nnodes; k++) {
      tnode_t * node = tree->node + k;
      if (node->leaf >= 0) continue;
      for (int a = 0; a <= tau; a++) {
         while (j < (tau+1)*nnodes) {
            lstack_t * lstck = lstacks->stack[j++];
            if (lstck->seq != MERGE_DONE) continue;
            node->data[a] = lstck;
            node->data[a]->seq = MERGE_PENDING;
            node->data[a]->tau = a;
         }
      }
   }

   return tree->node;
}

int
merge_node
(
 tnode_t * parent,
 int       targettau,
 arg_t   * args
)

// Arg_t contains:
// - Char pointer to the sequence. This will be used for the extras.
// - Number of chunks (nleaves).
// - Pointer to lookup table.
// - Pointer to index.
{
   // Do not compute if unnecessary.
   if (targettau < parent->mintau) return MERGE_EMPTY;

   // **********
   // Leaf nodes
   // **********

   if (parent->leaf >= 0) {

      loc_t nloc;

      // Avoid recomputing.
      if (parent->data[targettau]->seq == parent->status) {
         nloc = parent->data[targettau]->pos;
      }
      else {
         // Compute extras.
         int slen   = strlen(args->query);
         int offset = parent->leaf * SEQLEN;
         char extras[tau];
         for (int i = 0; i < tau; i++) extras[i] = -1;
         for (int i = 0; (i < tau) && (offset + i < slen); i++) extras[i] = args->query[offset + i];
      
         // Generate mismatches.
         sma(parent->mstack, parent->status, SEQLEN, extras, targettau);

         // Lookup loci list.
         nloc = lookup(targettau, parent->mstack[targettau], args->lut, args->index, parent->data + targettau);
      }

      // Check output.
      if (nloc > 0) return MERGE_OK;
      else {
         // If all previous taus were empty, increase mintau.
         int empty = 1;
         for (int i = 0; i < targettau; i++)
            if (parent->data[i]->pos > 0) {
               empty = 0;
               break;
            }
         if (empty) parent->mintau = targettau + 1;

         return MERGE_EMPTY;
      }
   }

   // ***********
   // Inner nodes
   // ***********

   // Do not repeat the same merging again!
   if (parent->data[targettau]->seq == MERGE_DONE) {
      if (parent->data[targettau]->pos == 0) return MERGE_EMPTY;
      else return MERGE_OK;
   }

   // Collect children loci and merge.
   for (int a = 0; a <= tau; a++) {
      tnode_t *na, *nb;
      if (a < tau - a) {
         na = parent->lchild;
         nb = parent->rchild;
      } else {
         na = parent->rchild;
         nb = parent->lchild;
      }
      
      // If anything goes wrong, update mintau and continue (Take into account that there may be untracked insertions from previous taus).
      if (merge_node(na, min(a, tau-a)) != MERGE_OK && na->data[targettau]->pos == 0) {
            parent->mintau = parent->lchild->mintau + parent->rchild->mintau;
            if (targettau < parent->mintau) return MERGE_EMPTY;
            continue;
      }

      // This repetition looks crappy, but it may avoid computing higher order taus.
      if (merge_node(nb, max(a, tau-a)) != MERGE_OK && nb->data[targettau]->pos == 0) {
         parent->mintau = parent->lchild->mintau + parent->rchild->mintau;
         if (targettau < parent->mintau) return MERGE_EMPTY;
         continue;
      }

      // Sequence distance between nodes.
      int seqdist = seqstart(parent->rchild) - seqstart(parent->lchild);
      merge_lstack(parent->data, parent->lchild->data[a], parent->rchild->data[tau-a], seqdist, targettau, args->maxtau);
   }
   parent->data[targettau]->seq = MERGE_DONE;
   if (parent->data[targettau]->pos == 0) {
      if (parent->status == NODE_OPEN) parent->mintau = targettau + 1;
      return MERGE_EMPTY;
   }

   parent->status = NODE_SET;
   
   // Sort.
   qsort(parent->data[targettau]->u, parent->data[targettau]->pos, sizeof(loc_t), loccomp);

   return MERGE_OK;
}

int
seqstart
(
 tnode_t * node
)
{
   if (node->leaf >= 0)
      return node->leaf * SEQLEN;
   else
      return seqstart(node->lchild);
}
 
 
void
merge_lstack
(
 lstack_t ** dest,
 lstack_t  * lstack,
 lstack_t  * rstack,
 int         seqdist,
 int         targettau,
 int         maxtau
)
{
   size_t size = min(lstack->pos, rstack->pos);

   // Realloc.
   if (dest[targettau]->lim < size) {
      dest[targettau] = realloc(dest[targettau], sizeof(seq_t) + 3*sizeof(int) + size*sizeof(loc_t));
      if (dest[targettau] == NULL) {
         fprintf(stderr, "Realloc error (merge_lstack): %s.\n", strerror(errno));
         exit(EXIT_FAILURE);
      }
      dest[targettau]->lim = size;
   }

   lstack_t * stack = dest[targettau];

   // Merge.
   int offset = maxtau - targettau;
   loc_t * llist = lstack->u;
   loc_t * rlist = rstack->u;
   int l = 0, r = 0;
   while (l < lstack->pos && r < rstack->pos) {
      // Match condition. Anotate all rights that match.
      int ref = llist[l] + seqdist;
      int nextr = 0;
      while (r + nextr < rstack->pos) {
         int dist = rlist[r + (nextr++)] - ref;
         if (dist >= -offset && dist <= offset) {
            if (dist == 0) stack->u[stack->pos++] = llist[l];
            else lstack_add(dest + abs(dist), llist[l] + dist);
         } else break;
      }
      if (rlist[r] - ref > 0 || nextr > 0) l++;
      else r++;
   }
}

char **
read_file
(
   FILE *inputf,
   int *nlines
)
{
   size_t size = 1024;
   char **all_seq = malloc(size * sizeof(char *));
   if (all_seq == NULL) exit(EXIT_FAILURE);

   *nlines = 0;
   ssize_t nread;
   size_t nchar = MAX_QUERYLEN + 1;
   char *seq = malloc(nchar * sizeof(char));
   // Read sequences from input file and store in an array. Assume
   // that it contains one sequence per line and nothing else. 
   while ((nread = getline(&seq, &nchar, inputf)) != -1) {
      // Strip end of line character and copy.
      if (seq[nread-1] == '\n') seq[nread-1] = '\0';
      char *new = malloc(nread);
      if (new == NULL) exit(EXIT_FAILURE);
      strncpy(new, seq, nread);
      // Grow 'all_seq' if needed.
      if (*nlines >= size) {
         size *= 2;
         char **ptr = realloc(all_seq, size * sizeof(char *));
         if (ptr == NULL) exit(EXIT_FAILURE);
         all_seq = ptr;
      }
      all_seq[(*nlines)++] = new;
   }
   all_seq = realloc(all_seq, *nlines * sizeof(char *));

   free(seq);
   return all_seq;
}
