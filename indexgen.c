#include "indexgen.h"

int
main
(
 int argc,
 char *argv[]
 )
{
   if (argc < 4) {
      fprintf(stderr,"Usage: trieseqs <seqlen> <offset> <genomefile>.fasta\n");
      exit(0);
   }

   // Open genome file.
   FILE * input = fopen(argv[3], "r");
   if (input == NULL) {
      fprintf(stderr, "fopen(INPUT) failed: %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   // Allocate tree.
   char * tree = (char *) calloc(sizeof(char), TREESZ);
   if (tree == NULL) {
      fprintf(stderr, "malloc(tree) failed: %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   // Allocate Loci-list pointers.
   loclst_t ** loclist = (loclst_t **) calloc(sizeof(loclst_t *), NSEQ);
   if (loclist == NULL) {
      fprintf(stderr, "malloc(loclist) failed: %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   // Sequence name stack.
   seqname_t * seqname = new_namestack(10);//(seqname_t *) malloc(sizeof(seqname_t));

   // Read genome file.   
   ssize_t nread;
   char *chunk = malloc(CHUNKSZ*sizeof(char));
   size_t nchar = INITBSZ * sizeof(char);
   char *seq  = malloc(nchar);

   int end = 0, i = 0;
   long loc = 0;
   while(end == 0) {
      nread = getline(&seq, &nchar, input);
      // Check end of chromosome or EOF.
      if (nread == -1 || seq[0] == '>') {
         // Process chunk and update local and locus index till the end of the chromosome.
         if (i > 0) i = procseqs(chunk, i, seqlen, offset, 1);
         if (seq[0] == '>') {
            // Link sequence name with the current locus.
            addseqname();
         }
         else if (nread == -1) {
            // Add 'end' sequence at the current locus.

            end = 1;
         }
         continue;
      }
      // Process data if chunk is full.
      if (i + nread - 1 > CHUNKSZ)
         i = procseqs(chunk, i, seqlen, offset, 0);

      // Add seq to chunk.
      memcpy(chunk + i, seq, nread - 1);
      i += nread - 1;
   }

   return 0;
} 


int
procseqs
(
 char      * tree,
 loclst_t ** loclst,
 int         genpos,
 int         nbases,
 char      * chunk,
 int         crop
 )
{
   // Number of seqs.
   int nseqs = nbases - SEQLEN + 1;
   int remainder = nbases - nseqs;
   char * seq = malloc((SEQLEN+1) * sizeof(char));

   // Uppercase bases.
   for (int i = 0; i < nbases; i++) {
      switch (chunk[i]) {
      case 'a':
         chunk[i] = 'A';
         break;
      case 'c':
         chunk[i] = 'C';
         break;
      case 't':
         chunk[i] = 'T';
         break;
      case 'g':
         chunk[i] = 'G';
         break;
      case 'n':
         chunk[i] = 'N';
      default:
         break;
      }
   }
   // Process seqs.
   for (int i = 0; i < nseqs; i++) {
      char * seq = chunk + i;
      int    nids;
      int  * seqids = seqid(seq, &nids);

      // Continue if there are too many 'N'.
      if(seqids == NULL) continue;

      // For all sequence ids...
      for (int j = 0; j < nids; j++) {
         int seqid = seqids[j];

         // Insert nodes in the tree.
         for (int h = 1; h <= SEQLEN; h++) {
            // Increase height.
            seq += 0x100000000;
            // ENABLE node bit.
            int bitpos = nodeaddr(seq);
            tree[bitpos/8] |= 1 << (bitpos%8);
         }

         // Add locus to loclst.
         seq &= SEQMASK;
         addlocus(loclst + seq, genpos + i);
      }

      free(seqids);
   }

   //TODO: 02/05/2014

   if (crop) {
      if (remainder+offset > seqlen) {
         strncpy(seq, chunk + nseqs*offset, remainder);
         seq[remainder] = 0;
         fprintf(stdout, "%s\n", seq);
      }
      free(seq);
      return 0;
   }

   // Save remainders.
   memmove(chunk, chunk + offset*nseqs, remainder);
   free(seq);
   return remainder;
}


int*
seqid
(
 char * seq,
 int  * nids
 )
{
   char nunkn = 0;
   // Count number of unknowns.
   for (int j = 0; j < SEQLEN; j++)
      nunkn += (seq[j] == 'N' || seq[j] == 'n');

   // Return NULL if max amount of 'N' has been reached.
   if (nunkn > MAXUNKN) return NULL;
      
   // Initialize seq ids.
   int nseqs = 1 << (2*nunkn);
   unsigned int * seqid = (unsigned int *) malloc(nseqs * sizeof(unisigned int));
   for (int k = 0; k < nseqs; k++) seqid[k] = 0;

   // Generate sequence ids.
   for (int j = SEQLEN-1, u = nseqs/4; j >= 0 ; j--) {
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
      if (j > 0)
         for (int k = 0; k < nseqs; k++) seqid[k] <<= 2;
   }

   *nids = nseqs;
   return seqid;
}


void
addlocus
(
 loclst_t ** loclst,
 int         locus
)
{
   loclst_t * list = *loclst;

   // Initialize stack.
   if (list == NULL) {
      // Initial stack size: 2 ints.
      list = (loclst_t *) malloc(4*sizeof(int));
      list->pos = 0;
   }

   // Stack grows exponentially. Resize if next pos is power of 2.
   if (!((list->pos + 1) & list->pos)) {
      int newsize = (list->pos + 1) << 1;
      list = realloc(list, (2+newsize)*sizeof(int));
      if (list == NULL) {
         fprintf(stderr, "realloc(list) failed at locus=%d: %s\n", locus, strerror(errno));
         exit(EXIT_FAILURE);
      }
   }

   // Add locus to list.
   list->l[list->pos++] = locus;
}
