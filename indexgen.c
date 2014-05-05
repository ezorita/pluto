#include "indexgen.h"

int
main
(
 int argc,
 char *argv[]
 )
{
   if (argc != 2) {
      fprintf(stderr,"Usage: plutoindex <genomefile>.fasta\n");
      exit(0);
   }

   // Open genome file.
   FILE * input = fopen(argv[1], "r");
   if (input == NULL) {
      fprintf(stderr, "fopen(INPUT) failed: %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }
   char * filename = strtok(argv[1], ".");

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
   chrstack_t * chrstack = new_chrstack(10);

   // Read genome file.   
   ssize_t nread;
   char *chunk = malloc(CHUNKSZ*sizeof(char));
   size_t nchar = INITBSZ * sizeof(char);
   char *seq  = malloc(nchar);

   int end = 0, i = 0;
   long loc = 0;

   // Read file, fill tree and create a loci list for each sequence.
   while (end == 0) {
      nread = getline(&seq, &nchar, input);
      // Check end of chromosome or EOF.
      if (nread == -1 || seq[0] == '>') {
         // Process chunk and update local and locus index till the end of the chromosome.
         if (i > 0) {
            procseqs(tree, loclist, loc, i, chunk);
            loc += i;
            i = 0;
         }

         // Link sequence name with the current locus.
         if (seq[0] == '>')
            addchrom(&chrstack, seq+1, loc);
         else if (nread == -1) {
            // Add 'end' sequence at the current locus.
            addchrom(&chrstack, "End of index.", loc);
            end = 1;
         }
         continue;
      }
      // Process data if chunk is full.
      if (i + nread - 1 > CHUNKSZ) {
         int offset = procseqs(tree, loclist, loc, i, chunk);
         loc += i;
         i = offset;
      }

      // Add seq to chunk.
      memcpy(chunk + i, seq, nread - 1);
      i += nread - 1;
   }

   // Write tree.
   char * outfile = malloc(strlen(filename) + 5);
   strcpy(outfile, filename);
   outfile = strcat(outfile, ".ptf");

   int fd = open(outfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   unsigned long wbytes = write(fd, tree, TREESZ*sizeof(char));
   if (wbytes < TREESZ)
      fprintf(stderr, "I/O error: tree (commited %lu bytes, %lu bytes written).\n", (unsigned long)TREESZ, wbytes);
   close(fd);

   free(tree);

   // Write chromosome index.
   outfile[strlen(outfile)-2] = 'c';

   FILE * chrout = fopen(outfile, "w");
   for (int i = 0; i < chrstack->pos-1; i++)
      fprintf(chrout, "%d\n%s\n", chrstack->c[i]->loc, chrstack->c[i]->name);
   fprintf(chrout, "%d\n%s", chrstack->c[chrstack->pos-1]->loc, chrstack->c[chrstack->pos-1]->name);
   fclose(chrout);

   // Allocate index.
   int * index = (int *) calloc(sizeof(int), NSEQ);
   if (index == NULL) {
      fprintf(stderr, "calloc(index) failed: %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   // Write loci list and generate index.
   outfile[strlen(outfile)-2] = 'l';
   
   fd = open(outfile, O_CREAT | O_WRONLY | O_TRUNC, 0644);
   // First int will be ignored, write some crap.
   wbytes = write(fd, &fd, sizeof(int));

   for (int i = 0, idx = 1; i < NSEQ; i++) {
      if (loclist[i] == NULL) continue;
      // Save memory offset to index.
      index[i] = idx;
      // Write loci list.
      int bytes = loclist[i]->pos * sizeof(int);
      wbytes = write(fd, &(loclist[i]->l[0]), bytes);
      if (wbytes < loclist[i]->pos * sizeof(char))
         fprintf(stderr, "I/O error: loci list %d (commited %d bytes, %lu bytes written).\n", i, bytes, wbytes);

      // Increase memory offset.
      idx += loclist[i]->pos;
   }
   close(fd);
    

   // Write index.
   outfile[strlen(outfile)-2] = 'i';
   
   fd = open(outfile, O_CREAT | O_WRONLY | O_TRUNC, 0644);
   wbytes = write(fd, index, NSEQ*sizeof(int));
   if (wbytes < NSEQ*sizeof(int))
      fprintf(stderr, "I/O error: index (commited %lu bytes, %lu bytes written).\n", NSEQ*sizeof(int), wbytes);
   close(fd);

   return 0;
} 


int
procseqs
(
 char      * tree,
 loclst_t ** loclst,
 int         genpos,
 int         nbases,
 char      * chunk
)
// SYNOPSIS:                                                              
//   Inserts nodes in the tree and updates the loci lists from a given set of bases.
//                                                                        
// PARAMETERS:                                                            
//   tree:   1-bit binary tree.
//   loclst: memory region containing pointers to a list of loci for each useq.
//   genpos: absolute position in the genome of the first base of 'chunk'.
//   nbases: number of valid bases present in 'chunk'.
//   chunk:  buffer containing the set of bases to be processed.
//                                                                        
// RETURN:                                                                
//   Returns the number of valid bases left at the buffer.
//                                                                        
// SIDE EFFECTS:                                                          
//   Modifies the contents of 'tree' and some *loclst may be allocated and/or updated.
//   
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
      unsigned int  * seqids = seqid(seq, &nids);

      // Continue if there are too many 'N'.
      if(seqids == NULL) continue;

      // For all sequence ids...
      for (int j = 0; j < nids; j++) {
         unsigned int seqid = seqids[j];

         // Insert nodes in the tree.
         for (int h = 0; h < SEQLEN; h++) {
            // Increase height.
            seq += 0x100000000;
            // ENABLE node bit.
            unsigned int bitpos = nodeaddr(seqid);
            tree[bitpos/8] |= 1 << (bitpos%8);
         }

         // Add locus to loclst.
         seqid &= SEQMASK;
         addlocus(loclst + seqid, genpos + i);
      }

      free(seqids);
   }

   free(seq);

   // Save remainders.
   memmove(chunk, chunk + nseqs, remainder);
   return remainder;
}


unsigned int*
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
   unsigned int * seqid = (unsigned int *) malloc(nseqs * sizeof(unsigned int));
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


chrstack_t *
new_chrstack
(
  int size
)
// SYNOPSIS:                                                              
//   Initializes a chromosome stack.
//                                                                        
// PARAMETERS:                                                            
//   size: initial size of the stack.
//                                                                        
// RETURN:                                                                
//   Returns a pointer to the allocated stack.
//                                                                        
// SIDE EFFECTS:
//   Allocates a memory region that must be manually freed.
//   
{
   chrstack_t * newstack = malloc(2*sizeof(int) + size*sizeof(chrom_t *));
   if (newstack == NULL) {
      fprintf(stderr, "malloc(chrstack_t) failed: %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }
   newstack->lim = size;
   newstack->pos = 0;

   return newstack;
}

void
addchrom
(
 chrstack_t ** chrstack,
 char       *  name,
 int           locus
)
// SYNOPSIS:                                                              
//   Inserts a new chromosome in a given chromosome stack.
//                                                                        
// PARAMETERS:                                                            
//   chrstack: pointer to an initialized chromosome stack pointer.
//   name:     chromosome name.
//   locus:    locus of the first nucleotide of the chromosome.
//                                                                        
// RETURN:                                                                
//   void.
//                                                                        
// SIDE EFFECTS:                                                          
//   Modifies the contents of the stack and may change its address if realloc is performed.
//   
{
   chrstack_t * stack = *chrstack;
   if (stack->pos == stack->lim) {
      int newsize = stack->lim * 2;
      *chrstack = stack  = realloc(stack, 2*sizeof(int) + newsize*sizeof(chrom_t *));
      if (stack == NULL) {
         fprintf(stderr, "realloc(chrstack_t) failed: %s\n", strerror(errno));
         exit(EXIT_FAILURE);
      }
      stack->lim = newsize;
   }

   // Alloc new chromosome.
   char * chrname = malloc(strlen(name) + 1);
   if (chrname == NULL) {
      fprintf(stderr, "malloc(char * chrname) failed: %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }
   strcpy(chrname, name);

   chrom_t * chr = malloc(sizeof(chrom_t));
   if (chr == NULL) {
      fprintf(stderr, "malloc(chrom_t) failed: %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }
   chr->loc  = locus;
   chr->name = chrname;

   //Add chromosome to the stack.
   stack->c[stack->pos++] = chr;
}

void
addlocus
(
 loclst_t ** loclst,
 int         locus
)
// SYNOPSIS:                                                              
//   Inserts a new locus in a given loci list.
//                                                                        
// PARAMETERS:                                                            
//   loclst: pointer to a loci-list pointer.
//   locus:  absolute position in the index to be added.
//                                                                        
// RETURN:                                                                
//   void.
//                                                                        
// SIDE EFFECTS:                                                          
//   Modifies the contents of the list and may change its address if realloc is performed.
//   
{
   loclst_t * list = *loclst;

   // Initialize stack.
   if (list == NULL) {
      // Initial stack size: 2 ints.
      *loclst = list = (loclst_t *) malloc(4*sizeof(int));
      list->pos = 0;
   }

   // Stack grows exponentially. Resize if next pos is power of 2.
   if (!((list->pos - 1) & list->pos) && list->pos > 1) {
      int newsize = list->pos << 1;
      *loclst = list = realloc(list, (2+newsize)*sizeof(int));
      if (list == NULL) {
         fprintf(stderr, "realloc(list) failed at locus=%d: %s\n", locus, strerror(errno));
         exit(EXIT_FAILURE);
      }
   }

   // Add locus to list.
   list->l[list->pos++] = locus;
}

unsigned int
nodeaddr
(
  unsigned int nodeid
)
{
   return (HOFFSET & (SEQMASK >> (2*(SEQLEN-((nodeid>>28) & 0x0000000F))))) + (nodeid & (SEQMASK >> (2*(SEQLEN-((nodeid>>28) & 0x0000000F)))));
}
