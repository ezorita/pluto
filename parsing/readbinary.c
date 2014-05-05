#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

int main(int argc, char * argv[]) {
   if (argc != 2) {
      fprintf(stderr, "usage: readb <file>\n");
      exit(1);
   }

   int fd = open(argv[1], O_RDONLY);
   unsigned char b;
   unsigned long count = 0;
   while (read(fd, &b, 1) == 1) {
      int bnum = 0, shift = 10000000;
      for (int i = 0; i < 8; i++) {
         bnum = bnum + ((b>>i)&1)*shift;
         count += (b>>i)&1;
         shift /= 10;
      }
      fprintf(stdout, "%08d", bnum);
   }
   fprintf(stderr, "number of nodes: %lu\n", count);

   return 0;
}
