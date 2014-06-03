#define _GNU_SOURCE
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>


#ifndef __PLUTO_CORE_
#define __PLUTO_CODE_

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
typedef struct tnode_t tnode_t;
typedef struct lstack_t lstack_t;
typedef struct mismatch_t mismatch_t;
typedef struct mstack_t mstack_t;
typedef struct stackbuf_t stackbuf_t;

struct tnode_t {
   char               mintau;
   char               leaf;
   seq_t              status;
   struct tnode_t   * lchild;
   struct tnode_t   * rchild;
   struct lstack_t ** data;
};

struct tree_t {
   int              nnodes;
   struct tnode_t * node;
};


struct stackbuf_t {
   int                count;
   struct lstack_t ** stack;
};

struct lstack_t {
   seq_t seq;
   int   tau;
   int   lim;
   int   pos;
   loc_t u[];
};

struct mismatch_t {
   seq_t seq;
   char  offset;
};

struct mstack_t {
   seq_t      seq;
   int        pos;
   int        lim;
   mismatch_t m[];
};

// Shared functions headers.
seq_t       seqtoid       (char *, int);
seq_t     * seqtoid_N     (char *, uint *, int);
char      * idtoseq       (uint, int);
uint        nodeaddr      (uint);
uint        getloci       (uint, uint *, uint *, uint **);
uint        addloci       (uint, uint *, uint *, ustack_t **);
ustack_t  * new_ustack    (uint);
cstack_t  * new_cstack    (uint);
ustack_t ** new_uarray    (uint, uint);
cstack_t ** new_carray    (uint, uint);
void        ustack_add    (ustack_t **, uint);
void        cstack_add    (cstack_t **, uchar *, uint);
uint        get_prefixlen (uint, uint, int);

#endif
