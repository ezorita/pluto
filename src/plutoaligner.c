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
      fprintf(stderr, "usage: pluto <sequence> <tau> <.pif file>\n");
      exit(1);
   }

  // Backtrace handler
   signal(SIGSEGV, SIGSEGV_handler); 

   // Check seq len.
   if(strlen(argv[1]) != 14) {
      fprintf(stderr, "Error: Only seqs of len 14 are supported.\nusage: query <sequence> <.pif file>\n");
      exit(EXIT_FAILURE);
   }

   // Parse input params.
   char * seq = argv[1];
   int    tau = atoi(argv[2]);
   char * filename = argv[3];

   // DEBUG:
   fprintf(stderr, "load params: tau=%d\tseq=%s\tfilename=%s\n", tau, seq, filename);

   // Open index file.
   int fdi = open(filename, O_RDONLY);

   // Get file size.
   unsigned long isize = lseek(fdi, 0, SEEK_END);
   lseek(fdi, 0, SEEK_SET);
   
   // Map LUT and index to memory.
   uint *lut = (uint *) mmap(NULL, isize, PROT_READ, MAP_SHARED, fdi, 0);
   uint *index = lut + NSEQ;

   if (lut == MAP_FAILED) {
      fprintf(stderr, "error loading index (mmap): %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   // Open tree file.
   filename[strlen(filename)-2] = 't';
   int fdt = open(filename, O_RDONLY);

   // Get file size.
   unsigned long tsize = lseek(fdt, 0, SEEK_END);
   lseek(fdt, 0, SEEK_SET);

   // Map tree.
   uchar *tree = (uchar *) mmap(NULL, tsize, PROT_READ, MAP_SHARED, fdt, 0);

   if (tree == MAP_FAILED) {
      fprintf(stderr, "error loading tree (mmap): %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   // Hit bins:
   // - Allocate a hit bin for each kind of hit.  (...2 ins, 1 ins, 0 ins/del, 1 del, 2 del...)
   // - For a given maxtau:
   //    Nbins = [tau=0] 1 + [tau = 1] (1 + 2) + ... + [tau = k] (1 + 2k) + ... + (1 + 2maxtau)
   //    General formula : aggregated Nbins = 1 + maxtau*(maxtau + 2)

   usnap_t * hits[MAXTREEQUERY];
   for (int l = 0; l < MAXTREEQUERY; l++) {
      hits[l] = (usnap_t *) malloc(tau*sizeof(usnap_t));
      for (int a = 0; a < tau; a++) {
         hits[l][a].ustack = new_uarray(2*a + 1, 4);
         hits[l][a].ustack += a;         // Center hit bins.
         hits[l][a].lastid = 0xF0000000; // Flag uninitialized.
      }
   }

   // Only the last tree will be trailed, so there's no need to allocate one cache array for each tree.
   cstack_t ** cache[tau];
   usnap_t  *  miles[tau];
   for (int a = 1; a < tau; a++) {
      cache[a] = new_carray((SEQLEN-1), 4*(2*a + 1));
      miles[a] = (usnap_t *) malloc((SEQLEN-1)*sizeof(usnap_t));
      for (int h = 0; h < SEQLEN-1; h++) {
         miles[a][h].ustack = new_uarray(1, 4);
         miles[a][h].lastid = 0xF0000000;
      }
   }

   // Algorithm args.
   struct searcharg_t arg = {
      .tau        = 0,
      .trail      = 0,
      .query      = 0,
      .tree       = tree,
      .lut        = lut,
      .index      = index,
      .hits       = NULL,
      .milestones = NULL,
      .cstack     = NULL,
   };


   // Parent cache.
   uchar rootcache[2*tau + 1];
   uchar * cachep = &rootcache[0] + tau;
   for (int i = -tau; i <= tau; i++) cachep[i] = (i < 0 ? -i : i);
   

   // TODO:
   // - Define 'nseqs' as the number of sorted input seqs.
   // - Define 'useq' as the char* vector of sorted unique input seqs.

   // Alignment snapshot:
   // - HITS:
   //   After querying the tree, a snapshot of the hit stack will be saved in hits_snapshot.
   //   The sequence associated with the snapshot will be kept and used by future queries to
   //   see whether the snapshot can be reused.
   // - MILESTONES:
   //   The local milestone stack array has an stack for each tau and height combination. After saving the
   //   milestones on a concrete stack (tau, height), the sequence id is saved at this same position on 
   //   'mile_last_sid'. These ids are later used by other sequences to see whether the milestones saved at
   //   a certain height can be reused.

   for (int i=0; i < nseqs; i++) { 
      // Query length and no. of times the tree will be queried.
      int  qlen    = strlen(useq[i]);
      int  ntrees  = qlen / SEQLEN;
      int  trail_miles, trail_trees; 
      // Alignment snapshot.
      // Make the snapshot based on the next sequence.
      // - trail_trees: if prefix > SEQLEN there's no need to query the tree again,
      //   reuse the hits of the last query.
      // - trail_miles: resume the tree query from the milestones at this height.
      if (i < nseqs-1) {
         int prefix = 0;
         while (useq[i][prefix] == useq[i+1][prefix] && prefix < qlen) prefix++;
         trail_trees = prefix / SEQLEN;
         trail_miles = prefix % SEQLEN;
      } else {
         trail_trees = 0;
         trail_miles = 0;
      }

      for (int a=0; a<=tau; a++) {
         // Set args for the current distance.
         arg.tau        = a;
         arg.milestones = miles[a];
         
         //******************************************************
         //***                   TREE QUERY                   ***
         //******************************************************
         for (int t = 0; t < min(ntrees, MAXTREEQUERY); t++) {
            // Set args for current tree.
            arg.hits = hits[t] + a;

            // Split the query in 14mers.
            uint  nids;
            uint *seqids = seqtoid(useq[i]+SEQLEN*t, &nids);
            if (nids != 1) {
               fprintf("error: BAD sequence: %s\n", useq[i]);
               break;
            }
            arg.query = seqids[0];

            // Check for available snapshots of this query.
            int snapshot = 0;
            for (int k = 0; k < MAXTREEQUERY; k++) {
               if (snapshot = (arg.query == hits[k][a].lastid)) {
                  // Copy hit stack if necessary.
                  if (k != t) {
                     for (int l = -a; l <=a; l++)
                        copy_ustack(hits[t][a].ustack + l, hits[k][a].ustack + l);
                  }
                  break;
               }
            }
            
            // No need to compute if snapshot is available.
            if (snapshot) continue;


            // ...SEARCH FUNCTION...

            // Clear hits.
            arg.hits->lastid = arg.query;
            for (int l = -arg.tau; l <= arg.tau; l++)
               arg.hits->ustack[l]->pos = 0;

            // Start the search.
            if (arg.tau == 0) {
               // Go straight to the bottom of the tree.
               // Either check the tree if node=0 OR check if the entry in the LUT is 0.
               if (arg.lut[arg.query & SEQMASK])
                  addloci(arg.query, arg.lut, arg.index, arg.hits->ustack);
            }
            else {
               // Start from the most favorable milestone snapshot.
               int hstart = 0;
               for (int h = SEQLEN-1; h > 0; h--) {
                  // Continue if snapshot is empty.
                  if (arg.milestones[h-1].lastid & ~SEQMASK) continue; 
                  if (get_prefixlen(arg.query, arg.milestones[h-1].lastid) >= h) {
                     hstart = h;
                     break;
                  }
               }

               // Clear the milestones that will be overwritten.
               for (int h = trail_miles; h > hstart; h--) {
                  arg.milestones[h-1].ustack[0]->pos = 0;
                  arg.milestones[h-1].lastid      = arg.query & SEQMASK;
               }

               // Start search.
               if (hstart == 0) {
                  _search(0, cachep - arg.tau, &arg, t);
               }
               else {
                  for (int m=0; m < arg.milestones[hstart-1]->pos; m++) {
                     _search(arg.milestones[hstart-1].ustack[0]->u[m], arg.cstack[hstart-1]->c + m*(2*arg.tau+1), &arg);
                  }
               }
            }
            // TODO:
            // - Copy the content of the hits stack if t < trail_trees.
            //   This will allow the next seq to recover the hits snapshot and continue the query.

         }

         // So far we have a list of candidate loci for each tree (BEFORE the merging filter).
         // The temptative distance of these loci is arg.tau, but this still needs
         // to be verified since we may have the following exceptions:
         // 1. Distance between loci is SEQLEN + k.
         //    There are k untracked deletions in the query. So the final tau is arg.tau + k.
         // 2. Distance between loci is SEQLEN - k.
         //    There are k untracked insertions in the query. The final tau is arg.tau + k.

         // Check whether there is any remaining hit after tree search.
         int nhits = 0;
         for (int l = -arg.tau; l <= arg.tau; l++)
            nhits += arg.hits[l]->pos;

         if (nhits > 0) {
            //******************************************************
            //***            LOCAL GENOME ALIGNMENT              ***
            //******************************************************

            // Run the local genome alignment only for the spare nucleotides.
            int  sparent = qlen % SEQLEN;
            // If there has been a match, add it to the match list and break the tau loop!
            
         }
      }
   }

   // Unmap memory.
   munmap(lut, isize);
   munmap(tree, tsize);
   
   // Close fds.
   close(fdi);
   close(fdt);

   return 0;
}

void
_search
(
 uint     nodeid,
 uchar  * pcache,
 sarg_t * arg
)
// filter: whether to fill or to zig-zag-filter the hits stack.
//         filter = 0 for the first 14mer, 1 for the next ones of the same query.
{
   // Root node has height bits = '1111', so get_height will return 0.
   uchar depth = get_height(nodeid) + 1;
   // Point to the middle of parent cache, ie the angle of the L.
   // This makes it easier to distinguish the part that goes upward,
   // with positive index and requiring the path, from the part that
   // goes horizontally, with negative index and requiring previous
   // characters of the query.
   pcache += arg->tau;
   // Risk of overflow at depth lower than 'tau'.
   int maxa = min((depth-1), arg->tau);

   // Remove the outer diagonals whose score > tau.
   int ctop, cleft;
   ctop = cleft = arg->tau;
   
   while (pcache[ctop]  > arg->tau) ctop--;
   while (pcache[cleft] > arg->tau) cleft--;

   // Penalty for match/mismatch and insertion/deletion resepectively.
   unsigned char mmatch;
   unsigned char shift;

   // Part of the cache that is shared between all the children.
   uchar cache[2*arg->tau + 3];
   uchar * common = &cache[0] + arg->tau + 1;

   // Fill outer diagonals with tau + 1
   for (int i = -arg->tau - 1; i < -cleft; i++) common[i] = arg->tau + 1;
   for (int i =  arg->tau + 1; i > ctop; i--)   common[i] = arg->tau + 1;

   // The branch of the L that is identical among all children
   // is computed separately. It will be copied later.
   for (int a = min(maxa, ctop) ; a > max(0, -cleft-1) ; a--) {
      // Upper arm of the L (need the path).
      mmatch = pcache[a] + (get_nt(nodeid, depth - a) != get_nt(arg->query, depth));
      shift = min(pcache[a-1], common[a+1]) + 1;
      common[a] = min(mmatch, shift);
   }

   // Continue with children.
   for (int i=0; i < 4; i++) {
      // Skip if current node has no child at this position.
      uint childid = get_child(nodeid, i);
      if(arg->tree[nodeaddr(childid)] == 0) continue;
      
      // Horizontal part of the L.
      for(int a = min(maxa, cleft); a > max(0, -ctop-1); a--) {
         mmatch = pcache[-a] + (get_nt(nodeid, depth) != get_nt(arg->query, depth - a));
         shift  = min(pcache[1 - a], common[-a - 1]) + 1;
         common[-a] = min(mmatch, shift);
      }

      // Center cell if center diagonal is alive.
      if (ctop >= 0 && cleft >= 0) {
         mmatch = pcache[0] + (i != get_nt(arg->query, depth));
         shift  = min(common[1], common[-1]) + 1;
         common[0] = min(mmatch, shift);
      }

      // Stop searching if tau is exceeded at all alive diagonals.
      int dels = 0;
      for (int a = maxa; a >= -maxa; a--) {
         if (common[a] > arg->tau) dels++;
         else break;
      }
      if (dels == 2*maxa + 1) continue;

      // Reached the height, check for hits.
      if (depth == SEQLEN) {
         // Store instertions and deletions in different bins.
         // Insertions:
         for (int a = maxa; a > 0; a--) {
            if (common[a] == arg->tau) {
               // Prune the tree and return all loci below this node.
               uint tmpid = childid & (0xFFFFFFFF << 2*a);
               // TODO:
               // - Profile: this func may be slow (too many reallocs inside..).
               for (int k = 0; k < (1 << 2*a); k++)
                  if (arg->tree[nodeaddr(tmpid + k)])
                     addloci(tmpid + k, arg->lut, arg->index, arg->hits->ustack + a);
            }
         }

         // Deletions and diagonal.
         for (int a = -maxa; a <= 0; a++) {
            if (common[a] == arg->tau)
               addloci(childid, arg->lut, arg->index, arg->hits->ustack + a);
         }

         continue;
      }

      // Cache node in milestones when trailing.
      if (depth <= arg->trail) save_milestone(arg->query, childid, arg->tau, common - arg->tau, arg->milestones, arg->cstack);

      _search(childid, common - arg->tau, arg);
   }
}

void
save_milestone
(
  uint        query,
  uint        nodeid,
  uint        tau,
  uchar     * cache,
  usnap_t  ** milestones,
  cstack_t ** cachestack
)
// Uses tau, nodeid, depth (inferred from nodeid) and the cache to save the milestone.
{
   uint       height = get_height(nodeid);
   ustack_t * mstack = milestones[height-1].ustack;
   cstack_t * cstack = cachestack[height-1];

   // Check data alignment.
   if (mstack->pos != cstack->pos)
      // WTF happened here?
      mstack->pos = cstack->pos = min(mstack->pos, cstack->pos);

   // Save nodeid and cache.
   ustack_add(&(milestones[height-1].ustack), nodeid);
   cstack_add(cachestack + height - 1, cache, tau);

   // Flag milestone snapshot.
   milestones[height].lastid = query;
}
