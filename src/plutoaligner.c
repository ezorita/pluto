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

   // Allocate stacks.   
   ustack_t ** miles = new_uarray((SEQLEN-1)*tau, 4);
   cstack_t ** cache = new_carray((SEQLEN-1)*tau, 4*(2*tau + 1));

   // Compute sequence IDs.
   uint nids;
   uint * sid = seqtoid(seq, &nids);

   // Algorithm vars.
   uint start = 0;
   uint trail = 0;
   struct searcharg_t arg = {
      .tau        = 0,
      .trail      = 0,
      .query      = 0,
      .tree       = tree,
      .lut        = lut,
      .index      = index,
      .hits       = new_ustack(8),
      .milestones = NULL,
      .cstack     = cache,
   };

   uchar rootcache[2*tau + 1];
   uchar * cachep = &rootcache[0] + tau;
   for (int i = -tau; i <= tau; i++) cachep[i] = (i < 0 ? -i : i);

   // DEBUG
   fprintf(stderr, "common cache:");
   for (int i = -tau; i <= tau; i++) fprintf(stderr, " %d", cachep[i]);
   fprintf(stderr, "\n");

   // For each sequence:
   for (int a=0; a<=tau; a++) {
      start = 0;
      arg.tau = a;
      arg.milestones = miles + (arg.tau-1)*(SEQLEN-1);
      
      // Output.
      fprintf(stdout, "---- DISTANCE = %d ----\n", a);
      
      for (int i=0; i<nids; i++) {
         // Add query to arguments.
         arg.query = sid[i];

         // Clear hits.
         arg.hits->pos = 0;

         // Search.
         if (arg.tau == 0) {
            // Go straight to the bottom of the tree.
            // Either check the tree if node=0 OR check if the entry in the LUT is 0.
            if (arg.lut[arg.query & SEQMASK])
               addloci(arg.query, arg.lut, arg.index, &(arg.hits));
         }
         else {
            // Find prefix len.
            arg.trail = 0;
            if (i < nids - 1 && arg.tau > 0)
               arg.trail = get_prefixlen(sid[i], sid[i+1]);

            // Reset milestones for for height > start.
            for (int h = start; h < SEQLEN-1; h++)
               arg.milestones[h]->pos = 0;

            if (start == 0)
               _search(0, cachep - arg.tau, &arg);
            else
               for (int m=0; m < arg.milestones[start-1]->pos; m++)
                  _search(arg.milestones[start-1]->u[m], arg.cstack[start-1]->c + m*(2*arg.tau+1), &arg);
            
         }

         // Output results.
         char * nseq = idtoseq(sid[i]);
         fprintf(stdout, "\tsequence: %s (%#08x)\tresults: %u\n", nseq, sid[i], arg.hits->pos);
         free(nseq);
         for (int h = 0; h < min(arg.hits->pos, 10); h++)
            fprintf(stdout, "\t\t%u\n", arg.hits->u[h]);
         
         // Starting point for next query.
         start = trail;
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

   // Penalty for match/mismatch and insertion/deletion resepectively.
   unsigned char mmatch;
   unsigned char shift;

   // Part of the cache that is shared between all the children.
   uchar cache[2*arg->tau + 3];
   uchar * common = &cache[0] + arg->tau + 1;
   for (int i = -arg->tau - 1; i <=  arg->tau + 1; i++) common[i] = (i < 0 ? -i : i);;

   // The branch of the L that is identical among all children
   // is computed separately. It will be copied later.
   for (int a = maxa ; a > 0 ; a--) {
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
      for(int a = maxa; a > 0; a--) {
         mmatch = pcache[-a] + (get_nt(nodeid, depth) != get_nt(arg->query, depth - a));
         shift  = min(pcache[1 - a], common[-a - 1]) + 1;
         common[-a] = min(mmatch, shift);
      }

      // Center cell.
      mmatch = pcache[0] + (i != get_nt(arg->query, depth));
      shift  = min(common[1], common[-1]) + 1;
      common[0] = min(mmatch, shift);

      // Stop searching if tau is exceeded.
      if (common[0] > arg->tau) continue;

      // Reached the height, it's a hit!
      if (depth == SEQLEN && common[0] == arg->tau) {
         addloci(childid, arg->lut, arg->index, &(arg->hits));
         continue;
      }

      // Dash if trail is over and no more mismatches allowed.
      if ((common[0] == arg->tau) && (depth > arg->trail)) {
         uint matchid;
         // Go straight to the bottom of the tree.
         matchid = add_suffix(childid, arg->query);

         if (matchid == 202752348)
            fprintf(stderr, "debug\n");

         // Either check the tree if node=0 OR check if the entry in the LUT is 0.
         if (arg->lut[matchid & SEQMASK]) {
            addloci(matchid, arg->lut, arg->index, &(arg->hits));
            continue;
         }
      }
      
      // Cache node in milestones when trailing.
      if (depth <= arg->trail) save_milestone(childid, arg->tau, common - arg->tau, arg->milestones, arg->cstack);

      _search(childid, common - arg->tau, arg);
   }
}

void
save_milestone
(
  uint        nodeid,
  uint        tau,
  uchar     * cache,
  ustack_t ** milestones,
  cstack_t ** cachestack
)
// Uses tau, nodeid, depth (inferred from nodeid) and the cache to save the milestone.
{
   uint       height = get_height(nodeid);
   ustack_t * mstack = milestones[height];
   cstack_t * cstack = cachestack[height];

   // Check data alignment.
   if (mstack->pos != cstack->pos)
      // WTF happened here?
      mstack->pos = cstack->pos = min(mstack->pos, cstack->pos);

   // Save nodeid and cache.
   ustack_add(milestones + height - 1, nodeid);
   cstack_add(cachestack + height - 1, cache, tau);
}
