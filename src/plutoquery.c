#include "plutocore.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <execinfo.h>
#include <signal.h>

#define min(a,b) (a < b ? a : b)

void SIGSEGV_handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}


int main(int argc, char * argv[]) {
   if (argc != 3) {
      fprintf(stderr, "usage: query <sequence> <.pif file>\n");
      exit(1);
   }

   // Backtrace handler
   signal(SIGSEGV, SIGSEGV_handler); 


   // Check seq len.
   if(strlen(argv[1]) != 14) {
      fprintf(stderr, "Error: Only seqs of len 14 are supported.\nusage: query <sequence> <.pif file>\n");
      exit(1);
   }

   // Parse input params.
   char * seq = argv[1];
   char * indexfile = argv[2];

   char * locifile = malloc(strlen(indexfile) + 1);
   strcpy(locifile, indexfile);
   locifile[strlen(locifile)-2] = 'l';
   int fdi = open(indexfile,O_RDONLY);
   int fdl = open(locifile,O_RDONLY);

   // Get file sizes.
   unsigned long isize = lseek(fdi, 0, SEEK_END);
   unsigned long lsize = lseek(fdl, 0, SEEK_END);
   lseek(fdi, 0, SEEK_SET);
   lseek(fdl, 0, SEEK_SET);
   
   // Map index and loci list to memory.
   uint *index = (uint *) mmap(NULL, isize, PROT_READ, MAP_SHARED, fdi, 0);
   uint *llist = (uint *) mmap(NULL, lsize, PROT_READ, MAP_SHARED, fdl, 0);
   
   if (index == MAP_FAILED || llist == MAP_FAILED) {
      fprintf(stderr, "error mmaping index: %s\n", strerror(errno));
      exit(1);
   }
   // Compute sequence IDs.
   int nids;
   uint * sid = seqtoid(seq, &nids);
   
   // For each sequence:
   for (int i=0; i<nids; i++) {
      uint * locus;
      int nloc = getloci(sid[i], index, llist, &locus);
      char * nseq = idtoseq(sid[i]);
      fprintf(stdout, "%s (%#08x)\t%d results\n", nseq, sid[i], nloc);
      free(nseq);
      for (int j = 0; j<min(nloc,10); j++)
         fprintf(stdout,"%u\n",locus[j]);
      fprintf(stdout,"\n");

   }

   // Unmap and close files.
   munmap(index, isize);
   munmap(llist, lsize);
   close(fdi);
   close(fdl);
}
