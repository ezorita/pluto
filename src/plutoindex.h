#include "plutocore.h"
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


// INDEXER PARAMS.
// Initial size of seq name stacks.
#define NAMESZ   10

// READ BUFFER SIZE.
#define INITBSZ 10
#define CHUNKSZ 1000000

// DATA TYPES.
typedef struct llst_t  loclst_t;
typedef struct chrstack_t chrstack_t;
typedef struct chrom_t chrom_t;

struct llst_t {
   loc_t   pos;
   loc_t   l[];
};

struct chrom_t {
   loc_t   loc;
   char * name;
};

struct chrstack_t {
          int       lim;
          int       pos;
   struct chrom_t * c[];
};

// FUNCTION HEADERS.
loclst_t     * new_loclist   (int);
chrstack_t   * new_chrstack  (int);
int            procseqs      (loclst_t **, char *, int, int, char *, int);
void           addlocus      (loclst_t **, int);
void           addchrom      (chrstack_t **, char *, int);
unsigned char  writegen      (int, char *, int, unsigned char);
