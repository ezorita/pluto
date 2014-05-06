#include "plutoindex.h"

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
   char * outfile = malloc(strlen(filename) + 5);
   strcpy(outfile, filename);
   outfile = strcat(outfile, ".pgf");
   int fd = open(outfile,O_CREAT | O_TRUNC | O_WRONLY, 0644);
   
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
   unsigned char last = 0;

   // Read file, fill tree and create a loci list for each sequence.
   while (end == 0) {
      // Read line and remove \n.
      nread = getline(&seq, &nchar, input);
      seq[nread - 1] = 0;

      // Check end of chromosome or EOF.
      if (nread == -1 || seq[0] == '>') {
         // Process chunk and update local and locus index till the end of the chromosome.
         if (i > 0) {
            // First write genome, then process index.
            last = writegen(fd, chunk, i, last);
            procseqs(tree, loclist, loc, i, chunk);
            loc += i;
            i = 0;
         }

         // Link sequence name with the current locus.
         if (seq[0] == '>')
            addchrom(&chrstack, seq+1, loc);
         else if (nread == -1) {
            // Add 'end' sequence at the current locus.
            addchrom(&chrstack, "END", loc);
            end = 1;
         }
         continue;
      }
      // Process data if chunk is full.
      if (i + nread - 1 > CHUNKSZ) {
         // Write genome and process index.
         last = writegen(fd, chunk, i, last);
         int offset = procseqs(tree, loclist, loc, i, chunk);
         loc += i;
         i = offset;
      }

      // Add seq to chunk.
      memcpy(chunk + i, seq, nread - 1);
      i += nread - 1;
   }
   outfile[strlen(outfile)-2] = 't';

   fd = open(outfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   unsigned long wbytes = write(fd, tree, TREESZ*sizeof(char));
   if (wbytes < TREESZ)
      fprintf(stderr, "I/O error: tree (commited %lu bytes, %lu bytes written).\n", (unsigned long)TREESZ, wbytes);
   close(fd);

   free(tree);

   // Write chromosome index.
   outfile[strlen(outfile)-2] = 'c';

   FILE * chrout = fopen(outfile, "w");
   for (int i = 0; i < chrstack->pos-1; i++)
      fprintf(chrout, "%u\n%s\n", chrstack->c[i]->loc, chrstack->c[i]->name);
   fprintf(chrout, "%u\n%s", chrstack->c[chrstack->pos-1]->loc, chrstack->c[chrstack->pos-1]->name);
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
   // First int will contain the total number of loci (Including nulls).
   wbytes = write(fd, &fd, sizeof(int));

   int idx = 1;
   for (int i = 0; i < NSEQ; i++) {
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

   // Write the size of the loci list at the first int.
   lseek(fd, SEEK_SET, 0);
   wbytes = write(fd, &idx, sizeof(int));
   if (wbytes < 1) {
      fprintf(stderr, "I/O error: loci list size (commited 4 bytes, %lu bytes written).\n", wbytes);
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

   // Process seqs.
   for (int i = 0; i < nseqs; i++) {
      char * seq = chunk + i;
      int    nids;
      uint  * seqids = seqtoid(seq, &nids);

      // Continue if there are too many 'N'.
      if(seqids == NULL) continue;

      // For all sequence ids...
      for (int j = 0; j < nids; j++) {
         uint sid = seqids[j];

         // Insert nodes in the tree.
         for (int h = 0; h < SEQLEN; h++) {
            // Increase height.
            sid += 0x10000000;
            // ENABLE node bit.
            uint bitpos = nodeaddr(sid);

            tree[bitpos/8] |= 1 << (bitpos%8);
         }

         // Add locus to loclst.
         sid &= SEQMASK;
         addlocus(loclst + sid, genpos + i);
      }

      free(seqids);
   }

   free(seq);

   // Save remainders.
   memmove(chunk, chunk + nseqs, remainder);
   return remainder;
}



unsigned char
writegen
(
  int             fd,
  char          * chunk,
  int             nbases,
  unsigned char   last
)
// SYNOPSIS:                                                              
//   Writes to the genome file coding each base with 4 bits:
//     '0000' = A
//     '0001' = C
//     '0010' = G
//     '0011' = T
//     '1111' = N
//                                                                        
// PARAMETERS:                                                            
//   fd:     open file descriptor pointing to the output genome file.
//   chunk:  buffer containing the set of bases to be processed.
//   nbases: number of valid bases present in 'chunk'.
//   last:   last 4 spare bits returned by the last writegen call.
//                                                                        
// RETURN:                                                                
//   A char containing the last 4 bits if [(nbases+(last>0))%2 == 1]. 0 otherwise.
//                                                                        
// SIDE EFFECTS:                                                          
//   None.
//   
{
   // First add remainder.
   nbases += (last != 0);
   int nbytes = nbases / 2;
   int remain = nbases % 2;
   unsigned char * values = calloc(sizeof(char), nbytes + remain);

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
         break;
      }
   }

   int i = 0;
   // Do the first base.
   if (last) {
      values[0] = last & 0X0F;
      i++;
   }

   for (int j = 0; i < 2*nbytes + remain; i++, j++) {
      switch (chunk[j]) {
      case 'C':
         values[i/2] |= 0x01 << (4*(i%2));
         break;
      case 'G':
         values[i/2] |= 0x02 << (4*(i%2));
         break;
      case 'T':
         values[i/2] |= 0x03 << (4*(i%2));
         break;
      case 'N':
         values[i/2] |= 0x0F << (4*(i%2));
         break;
      }
   }
   
   int wbytes = write(fd, values, nbytes);
   if (wbytes < nbytes) {
      fprintf(stderr, "I/O error: writegen (commited %d bytes, %d bytes written).\n", nbytes, wbytes);
      exit(EXIT_FAILURE);
   }

   if (remain) {
      unsigned char bits = values[nbytes];
      free(values);
      return bits | 0xF0;
   }

   free(values);
   return 0;
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
