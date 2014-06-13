#include "plutocore.h"
#include "mergesort.h"


seq_t
seqtoid
(
 char * seq,
 int   slen
)
{
  seq_t seqid = 0;
   for (int j = 0; j < slen ; j++) {
      if (seq[j] == 'A' || seq[j] == 'a') { }
      else if (seq[j] == 'C' || seq[j] == 'c') seqid += 1;
      else if (seq[j] == 'G' || seq[j] == 'g') seqid += 2;
      else if (seq[j] == 'T' || seq[j] == 't') seqid += 3;
      else return BAD_SEQ;
      if (j < slen-1) seqid <<= 2;
   }

   return seqid;   
}

seq_t*
seqtoid_N
(
 char * seq,
 int  * nids,
 int    slen
 )
// SYNOPSIS:                                                              
//   Generates the ID of the sequence pointed by seq. If the sequence
//   contains up to 'MAXUNKN' unknowns ('N'), seqid generates all the
//   combinations of the sequence obtained by replacing 'N' by 'ACTG'.
//   The function returns NULL otherwise.
//                                                                        
// PARAMETERS:                                                            
//   seq: String of the sequence (null-termination is not necessary).
//                                                                        
// RETURN:                                                                
//   Pointer to an array of uints containing all the ids or
//   NULL if the number of unknowns is found to be greater than MAXUNKN.
//   *nids: number of sequence ids generated by substituting unknowns.
//                                                                        
// SIDE EFFECTS:                                                          
//   None.
//   
{
   char nunkn = 0;

   // Count number of unknowns.
   for (int j = 0; j < slen; j++)
      nunkn += (seq[j] == 'N' || seq[j] == 'n');

   // Return NULL if max amount of 'N' has been reached.
   if (nunkn > MAXUNKN) return NULL;
      
   // Initialize seq ids.
   int nseqs = 1 << (2*nunkn);
   seq_t * seqid = (seq_t *) malloc(nseqs * sizeof(seq_t));
   for (int k = 0; k < nseqs; k++) seqid[k] = 0;

   // Generate sequence ids.
   for (int j = 0, u = nseqs/4; j < slen ; j++) {
      switch(seq[j]) {
      case 'C':
         for (int k = 0; k < nseqs; k++) seqid[k] += 1;
         break;
      case 'G':
         for (int k = 0; k < nseqs; k++) seqid[k] += 2;
         break;
      case 'T':
         for (int k = 0; k < nseqs; k++) seqid[k] += 3;
         break;
      case 'N':
         for (int k = 0; k < nseqs; k++) seqid[k] += (k/u) % 4;
         u /= 4;
         break;
      }
      if (j < slen - 1)
         for (int k = 0; k < nseqs; k++) seqid[k] <<= 2;
   }

   *nids = nseqs;
   return seqid;
}

anchor_t*
anchorid_N
(
 char * seq,
 int  * nids,
 int    seqlen,
 int    anclen
 )
// SYNOPSIS:                                                              
//   Generates the ID of the anchor pointed by seq. If the sequence
//   contains up to 'MAXUNKN' unknowns ('N'), seqid generates all the
//   combinations of the sequence obtained by replacing 'N' by 'ACTG'.
//   The function returns NULL otherwise.
//                                                                        
// PARAMETERS:                                                            
//   seq: String of the sequence (null-termination is not necessary).
//   nids: int pointers where the number of returned seqs will be placed.
//   seqlen: length of the index k-mers.
//   anclen: length of each individual anchor.
//                                                                        
// RETURN:                                                                
//   Pointer to an array of uints containing all the ids or
//   NULL if the number of unknowns is found to be greater than MAXUNKN.
//   *nids: number of sequence ids generated by substituting unknowns.
//                                                                        
// SIDE EFFECTS:                                                          
//   None.
//   
{
   char nunkn = 0;

   // Count number of unknowns.
   for (int j = 0; j < anclen; j++)
      nunkn += (seq[j] == 'N' || seq[j] == 'n');
   for (int j = seqlen; j < seqlen+anclen; j++)
      nunkn += (seq[j] == 'N' || seq[j] == 'n');

   // Return NULL if max amount of 'N' has been reached.
   if (nunkn > MAXUNKN) return NULL;
      
   // Initialize seq ids.
   int nseqs = 1 << (2*nunkn);
   anchor_t * ancid = (anchor_t *) malloc(nseqs * sizeof(anchor_t));
   for (int k = 0; k < nseqs; k++) ancid[k] = 0;

   // Generate sequence ids.
   int j = 0, u = nseqs/4;
   while (j < seqlen + anclen) {
      switch(seq[j]) {
      case 'C':
         for (int k = 0; k < nseqs; k++) ancid[k] += 1;
         break;
      case 'G':
         for (int k = 0; k < nseqs; k++) ancid[k] += 2;
         break;
      case 'T':
         for (int k = 0; k < nseqs; k++) ancid[k] += 3;
         break;
      case 'N':
         for (int k = 0; k < nseqs; k++) ancid[k] += (k/u) % 4;
         u /= 4;
         break;
      }
      if (j < seqlen + anclen - 1)
         for (int k = 0; k < nseqs; k++) ancid[k] <<= 2;
      
      j++;
      if (j == anclen) j = seqlen;
   }

   *nids = nseqs;
   return ancid;
}

char*
idtoseq
(
 seq_t seqid,
 int  slen
)
// SYNOPSIS:                                                              
//   Reverts the sequence associated with the given ID.
//                                                                        
// PARAMETERS:                                                            
//   seqid: the sequence ID.
//                                                                        
// RETURN:                                                                
//   A string containing the sequence.
//                                                                        
// SIDE EFFECTS:                                                          
//   Allocates 15 bytes to fit the sequence that must be manually freed.
//   
{
   char bases[4] = {'A','C','G','T'};
   char * seq = malloc(slen+1);
   for (int i= 0; i<slen; i++)
      seq[slen-1-i] = bases[(seqid >> 2*i) & 3];

   return seq;
}

loc_t
getloci
(
 seq_t    seq_id,
 loc_t    nloci,
 loc_t  * lut,
 loc_t  * index,
 loc_t ** list
)
// SYNOPSIS:                                                              
//   Given a genome LUT and its index, retrieves the address of the loci list
//   that correspond to the specified sequence id.
//                                                                        
// PARAMETERS:                                                            
//   seq_id: the node id containing the queried sequence.
//   nloci:  size of the index. Can be read at index[0].
//   lut:    an mmaped pointer to the genome lookup table.
//   index:  an mmaped pointer to the genome index file.
//   list:   pointer where the starting address of the list will be placed.
//                                                                        
// RETURN:                                                                
//   Returns the size of the list if any match has been found for the specified
//   sequence or 0 otherwise.
//   'list' will contain either the address of the first element in the list or 
//   NULL if the sequence does not match any existing loci.
//                                                                        
// SIDE EFFECTS:                                                          
//   None.
//   
{
   // Get the address offset of the 1st locus.
   loc_t start = lut[seq_id];
   if (start == 0) {
      *list = NULL;
      return 0;
   }

   // Find the next node to infer the offset.
   loc_t end = 0;
   seq_id++;
   while (seq_id < NSEQ) {
      if (lut[seq_id] != 0) {
         end = lut[seq_id];
         break;
      }
      else seq_id++;
   }

   if (seq_id >= NSEQ) end = nloci;

   // Point to the beginning of the list and return the list size.
   *list = index + start;
   return end - start;
}

loc_t
addloci
(
 seq_t       seq_id,
 loc_t     * lut,
 loc_t     * index,
 lstack_t ** lstackp
)
// SYNOPSIS:                                                              
//   Same as getloci, but directly appends the content of the loci list to a
//   stack of uints, instead of returning the address of the list.
//                                                                        
// PARAMETERS:                                                            
//   seq_id: the node id containing the queried sequence.
//   lut:    an mmaped pointer to the genome lookup table.
//   index:  an mmaped pointer to the genome index file.
//   lstack: uint stack where the loci will be appended.
//                                                                        
// RETURN:                                                                
//   The number of inserted loci.
//                                                                        
// SIDE EFFECTS:                                                          
//   May realloc the stack.
//   
{
   lstack_t * lstack = *lstackp;
   // Get the loci list.
   loc_t * list;
   loc_t   nloc = getloci(seq_id, index[0], lut, index, &list);
   if (nloc == 0) return 0;

   // TODO: CHECK HERE FOR POSSIBLE HIGHLY REPEATED SEQUENCES.
   
   // Realloc the stack if needed.
   if (lstack->pos + nloc >= lstack->lim) {
      loc_t newsize = lstack->pos + nloc;
      lstack_t * p = realloc(lstack, sizeof(seq_t) + 3*sizeof(int) + newsize * sizeof(loc_t));
      if (p == NULL) {
         fprintf(stderr, "error extending lstack (realloc): %s\n", strerror(errno));
         exit(EXIT_FAILURE);
      }
      *lstackp = lstack = p;
      lstack->lim = newsize;
   }
   
   // Copy data and update index.
   memcpy(lstack->u + lstack->pos, list, nloc*sizeof(loc_t));
   lstack->pos += nloc;

   return nloc;
}


loc_t
lookup
(
 int         tau,
 mstack_t  * mstack,
 loc_t     * lut,
 loc_t     * index,
 lstack_t ** lstackp
)
{
   lstack_t * lstack = *lstackp;
   // Avoid recomputing the same thing.
   if (mstack->seq == (*lstackp)->seq) return mstack->pos;
   
   // Get the number of loci in the genome. (Stored at index[0]).
   loc_t nloci = index[0];

   for (int i = 0; i < mstack->pos; i++) {
      // Get the loci list.
      loc_t * list;
      loc_t   nloc = getloci(mstack->m[i].seq, nloci, lut, index, &list);
      if (nloc == 0) continue;

         // Realloc the stack if needed.
      if (lstack->pos + nloc >= lstack->lim) {
         loc_t newsize = 2*lstack->lim;
         while (newsize < lstack->pos + nloc) newsize *= 2;
         lstack_t * p = realloc(lstack, sizeof(seq_t) + 3*sizeof(int) + newsize*sizeof(loc_t));
         if (p == NULL) {
            fprintf(stderr, "error extending lstack (realloc): %s\n", strerror(errno));
            exit(EXIT_FAILURE);
         }
         *lstackp = lstack = p;
         lstack->lim = newsize;
      }

      // Copy data and update index.
      memcpy(lstack->u + lstack->pos, list, nloc*sizeof(loc_t));
      lstack->pos += nloc;

      // Add offset now and forget about it! (Less painful option)
      char offset;
      if ((offset = mstack->m[i].offset) != 0) {
         for (int i = lstack->pos - nloc; i < lstack->pos; i++) {
            lstack->u[i] += offset;
         }
      }
   
   }

   // Sort loci if different loci were mixed.
   if (mstack->pos > 1 && lstack->pos > 1) {
      //mergesort_loc_int(lstack->u, intervals, nints-1, lstack->pos);
      radix_sort(lstack->u, 0, lstack->pos, 24);
   }
      // mergesort_loc(lstack->u, lstack->pos);
      //qsort(lstack->u, lstack->pos, sizeof(loc_t), loccomp);

   return (*lstackp)->pos;
}


lstack_t *
new_lstack
(
  int size
)
{
   lstack_t * ustack = malloc(sizeof(seq_t) + 3*sizeof(int) + size*sizeof(loc_t));
   if (ustack == NULL) {
      fprintf(stderr, "error allocating lstack_t in 'new_lstack' (malloc): %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   ustack->pos = 0;
   ustack->lim = size;
   ustack->seq = BAD_SEQ;
   ustack->tau = 0;

   return ustack;
}


mstack_t *
new_mstack
(
  int size
)
{
   mstack_t * mstack = malloc(sizeof(seq_t) + 3*sizeof(int) + size*sizeof(mismatch_t));
   if (mstack == NULL) {
      fprintf(stderr, "error allocating mstack_t in 'new_mstack' (malloc): %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   mstack->pos = 0;
   mstack->lim = size;
   mstack->seq = BAD_SEQ;
   mstack->tau = 0;

   return mstack;
}


void
lstack_add
(
 lstack_t ** lstackp,
 loc_t       value
)

{
   lstack_t * lstack = *lstackp;
   
   // Realloc the stack if needed.
   if (lstack->pos >= lstack->lim) {
      uint newsize = 2 * lstack->lim;
      lstack_t * p = realloc(lstack, sizeof(seq_t) + 3*sizeof(int) + newsize*sizeof(loc_t));
      if (p == NULL) {
         fprintf(stderr, "error while extending lstack (lstack_add/realloc): %s\n", strerror(errno));
         exit(EXIT_FAILURE);
      }
      *lstackp = lstack = p;
      lstack->lim = newsize;
   }
   
   // Add value.
   lstack->u[lstack->pos++] = value;
}

void
copy_lstack
(
 lstack_t ** dstp,
 lstack_t ** srcp
)
{
   lstack_t * src = *srcp;
   lstack_t * dst;
   
   size_t dstsz = src->pos*sizeof(loc_t) + 3*sizeof(int) + sizeof(seq_t);
   dst = *dstp = realloc(*dstp, dstsz);

   memcpy(dst, src, dstsz);
   dst->lim = src->pos;
}

void
copy_mstack
(
 mstack_t ** dstp,
 mstack_t ** srcp
)
{
   mstack_t * src = *srcp;
   mstack_t * dst;
   
   size_t dstsz = src->pos*sizeof(mismatch_t) + 3*sizeof(int) + sizeof(seq_t);
   dst = *dstp = realloc(*dstp, dstsz);
   if (dst == NULL) {
      fprintf(stderr, "error while copying mstack (copy_mstack/realloc): %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   memcpy(dst, src, dstsz);
   dst->lim = src->lim;
}


int
get_prefixlen
(
 seq_t seqa,
 seq_t seqb,
 int  slen
)
{
   int len = slen;
   // TODO:
   // - Seqmask depends on slen. (Maybe a macro SEQMASK(slen))
   seqa &= SEQMASK;
   seqb &= SEQMASK;
   while ((((seqa >> 2*(slen - len))&3) != ((seqb >> 2*(slen - len))&3)) && len > 0) len--;
   
   return len;
}

int
loccomp
(
 const void * a,
 const void * b
)
{
   loc_t * la = (loc_t *) a;
   loc_t * lb = (loc_t *) b;
   if (*la > *lb)
      return 1;
   else if (*la < *lb)
      return -1;
   else return 0;
}
