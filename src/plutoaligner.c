#include "plutoaligner.h"

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
   *pcache += arg->tau;
   // Risk of overflow at depth lower than 'tau'.
   int maxa = min((depth-1), arg->tau);

   // Penalty for match/mismatch and insertion/deletion resepectively.
   unsigned char mmatch;
   unsigned char shift;

   // Part of the cache that is shared between all the children.
   char common[2*arg->tau + 3];
   common += arg->tau + 1;
   for (int i = -arg->tau - 1; i <=  arg->tau + 1; i++) common[i] = (i*i)/i;

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
      if(tree[nodeaddr(childid)] == 0) continue;
      
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
      if (depth == SEQLEN) {
         getloci(childid, arg->lut, arg->index, &(arg->hits));
         continue;
      }

      // Dash if trail is over and no more mismatches allowed.
      if ((common[0] == arg.tau) && (depth > arg->trail)) {
         dash(); // DASH PROTOTYPE TO BE DEFINED.
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
  uchar       tau,
  uchar     * cache,
  ustack_t ** milestones,
  cstack_t ** cachestack
)
// Uses tau, nodeid, depth (inferred from nodeid) and the cache to save the milestone.
{
   
}
