#define _GNU_SOURCE
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>


#ifndef __PLUTO_CORE_
#define __PLUTO_CORE_

// Maximum number of unknowns (sequences with count('N') > MAXUNKN will not be indexed).
#define MAXUNKN  3
// Height offset:
// ** Apply an AND mask ~(0xFFFFFFFF << 2*height) to obtain the offset in the tree for a given height.
#define HOFFSET 0x05555554

// Sequence length:
#define SEQLEN  14

// Number of tree nodes:
// ** N_{nodes} = \sum_{i=1}^{14}{4^i} = 4^{14} * \sum_{i=0}^{13}{(1/4)^i} = 4^{14}*\frac{1-(1/4)^14}{1-1/4}
#define NNODES  357916940

// Size of the tree in memory:
// ** Number of bytes.
#define TREESZ  44739618

// Number of possible sequences:
// ** N_{seq} = 4^{14}
#define NSEQ    268435456

// Sequence mask:
#define SEQMASK 0x0FFFFFFF

// Define bad sequence.
#define BAD_SEQ       0xFFFFFFFF
#define MERGE_PENDING 0x7FFFFFFF
#define MERGE_DONE    0x3FFFFFFF

#define NODE_OPEN     0
#define NODE_SET      1

#define MERGE_EMPTY   -1
#define MERGE_OK      0

// Config params.
#define LOCI_STACKSIZE 256
#define MISM_STACKSIZE 256
#define MAX_CHUNKS     5
#define MAX_QUERYLEN   100

// Some macros:
#define min(a,b) (a < b ? a : b)
#define max(a,b) (a > b ? a : b)
#define abs(a) (a < 0 ? -a : a)
// Returns the height of the node.
#define get_height(node) ((node >> (2 * SEQLEN)) & 0x0F)
// Returns the nucleotide of the sequence at the specified height.
#define get_nt(node, height) ((node >> (2*(SEQLEN - height)) & 3))
// Returns the nodeid of the 1st generation child of the given sequence.
#define get_child(node, nt) (((node & (0xFFFFFFFF << (2*(SEQLEN - get_height(node))))) | ((nt & 3) << (2*(SEQLEN - get_height(node) - 1)))) + 0x10000000)
// Returns the bottom node of the tree that matches the remaining nucleotides of the query.
#define add_suffix(node,query) ((((node & ~(SEQMASK >> 2*get_height(node))) | (query & (SEQMASK >> 2*get_height(node)))) & SEQMASK) | 0xE0000000)

// Type definitions
typedef unsigned int seq_t;
typedef unsigned int loc_t;
typedef struct mismatch_t mismatch_t;
typedef struct lstack_t lstack_t;
typedef struct mstack_t mstack_t;

struct mismatch_t {
   seq_t seq;
   char  offset;
};

struct lstack_t {
   seq_t seq;
   int   tau;
   int   lim;
   int   pos;
   loc_t u[];
};

struct mstack_t {
   seq_t      seq;
   int        tau;
   int        pos;
   int        lim;
   mismatch_t m[];
};

// Shared functions headers.
seq_t       seqtoid       (char *, int);
seq_t     * seqtoid_N     (char *, int *, int);
char      * idtoseq       (seq_t, int);
loc_t       getloci       (seq_t, loc_t, loc_t *, loc_t *, loc_t **);
loc_t       addloci       (seq_t, loc_t *, loc_t *, lstack_t **);
loc_t       lookup        (int, mstack_t *, loc_t *, loc_t *, lstack_t **);
lstack_t  * new_lstack    (int);
mstack_t  * new_mstack    (int);
void        lstack_add    (lstack_t **, loc_t);
void        copy_lstack   (lstack_t **, lstack_t **);
void        copy_mstack   (mstack_t **, mstack_t **);
int         mcomp         (const void *, const void *);
int         loccomp       (const void *, const void *);

#endif
