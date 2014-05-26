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

   // Parse input params.
   int    tau = atoi(argv[1]);
   int    maxthreads = atoi(argv[2]);
   char * queryfile = argv[3];
   char * indexfile = argv[4];

   // DEBUG:
   fprintf(stderr, "load params: tau=%d\tseq=%s\tfilename=%s\n", tau, seq, indexfile);

   // Open index file.
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

   // Mismatch bins:
   mstack_t ** mstack[MAX_CHUNKS];
   for (int c = 0; c < MAX_CHUNKS; c++) {
      mstack[c] = malloc((tau+1)*sizeof(mstack_t *));
      for (int i = 0; i <= tau; i++) {
         int ssize = 1 << (4*i);
         mstack[c][i] = malloc(3*sizeof(seq_t) + ssize*sizeof(mismatch_t));
         mstack[c][i]->lim = ssize;
         mstack[c][i]->pos = 0;
         mstack[c][i]->seq = 0xFFFFFFFF;
      }
   }


   // Loci bins:
   // - Allocate a loci bin for each kind of mismatch.  (...2 ins, 1 ins, 0 ins/del, 1 del, 2 del...)
   // - For a given maxtau:
   //    Nbins = [tau=0] 1 + [tau = 1] (1 + 2) + ... + [tau = k] (1 + 2k) + ... + (1 + 2maxtau)
   //    General formula : aggregated Nbins = 1 + maxtau*(maxtau + 2)
   lstack_t ** lstack[MAX_CHUNKS];
   int nstacks = 1+maxtau*(maxtau+2);
   for (int c = 0; c < MAX_CHUNKS; c++) {
      lstack[c] = (lstack_t **) malloc(nstacks*sizeof(lstack_t *));
      for (int i = 0; i < nstacks; i++) lstack[c][i] = new_lstack(HITSTACK_SIZE);
   }


   // Read file:
   uint    num_query;
   char ** all_query = read_file(queryfile, &num_query);
   mergesort((void **) all_query, num_query, &ualpha, maxthreads);

   // Align sequences.
   seq_t chunk[MAX_CHUNKS];
   for (int q = 0; q < num_query; q++) {
      if (all_query[q] == NULL) continue;

      // Divide in chunks.
      char * query    = all_query[q];
      int    querylen = strlen(query);
      int    nchunks  = querylen/SEQLEN;
      int    repeat[nchunks];
      seq_t  seqid [nchunks];

      for (int c = 0; c < nchunks; c++) {
         repeat[c] = 0;
         seqid[c]  = seqtoid(query, SEQLEN);
      }
      
      int match = 0;
      while (match == 0) {
         // Get the Loci list.
         for (int c = 0; c < nchunks; c++) {
            
         }
      }
      
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
