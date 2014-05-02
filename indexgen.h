#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

// INDEXER PARAMS.
// Maximum number of unknowns (sequences with count('N') > MAXUNKN will not be indexed).
#define MAXUNKN  3
// Initial size of seq name stacks.
#define NAMESZ   10

// SEQUENCE LENGTH.
#define SEQLEN  14

// NUMBER OF TREE NODES.
// N_{nodes} = \sum_{i=1}^{14}{4^i} = 4^{14} * \sum_{i=0}^{13}{(1/4)^i} = 4^{14}*\frac{1-(1/4)^14}{1-1/4}
#define NNODES  357916940

// SIZE OF TREE IN MEMORY.
// Number of bytes.
#define TREESZ  44739618

// NUMBER OF POSSIBLE SEQUENCES.
// N_{seq} = 4^{14}
#define NSEQ    268435456

// HEIGHT OFFSET.
// Apply an AND mask ~(0xFFFFFFFF << 2*height) to obtain the offset in the tree for a given height.
#define HOFFSET 0x05555554

// SEQUENCE MASK.
#define SEQMASK 0x0FFFFFFF

// READ BUFFER SIZE.
#define INITBSZ 10
#define CHUNKSZ 1000000

// TREE MACROS.
#define nodeaddr(id) (HOFFSET & (SEQMASK >> (2*(SEQLEN-((id>>28) & 0x0000000F))))) + (id & (SEQMASK >> (2*(SEQLEN-((id>>28) & 0x0000000F)))))

// DATA TYPES.
typedef struct llst_t  loclst_t;
typedef struct sname_t seqname_t;

struct llst_t
{
   int    pos;
   int    l[];
};

struct sname_t
{
   int    lim;
   int    pos;
   int    l[];
   char * n[];
};

// FUNCTION HEADERS.

loclst_t  * new_loclist   (int);
seqname_t * new_namestack (int);
int         procseqs      (char *, loclst_t *, int, int, char *, int);
int       * seqid         (char *, int *);
void        addlocus      (loclst_t **, int);
