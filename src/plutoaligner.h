#include "plutocore.h"
#include "sma.h"
#include <execinfo.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

#ifndef _PLUTO_ALIGNER__
#define _PLUTO_ALIGNER__
// Type definitions
typedef struct tnode_t     tnode_t;
typedef struct tree_t      tree_t;
typedef struct lstack_t    lstack_t;
typedef struct lstackbuf_t lstackbuf_t;
typedef struct mstackbuf_t mstackbuf_t;
typedef struct arg_t       arg_t;

struct lstackbuf_t {
   int                count;
   struct lstack_t ** stack;
};

struct mstackbuf_t {
   int                count;
   struct mstack_t ** stack;
};


struct tnode_t {
   char               mintau;
   char               leaf;
   seq_t              status;
   struct tnode_t   * lchild;
   struct tnode_t   * rchild;
   struct lstack_t ** data;
   struct mstack_t ** mstack;
};

struct tree_t {
   int              nnodes;
   struct tnode_t * node;
};

struct arg_t {
   int     nleaves;
   int     maxtau;
   char  * query;
   loc_t * lut;
   loc_t * index;
};

tnode_t  * build_tree   (int, seq_t *, tree_t *, lstackbuf_t *, mstackbuf_t *, char);
int        seqstart     (tnode_t *);
int        merge_node   (tnode_t *, int, arg_t *);
void       merge_lstack (lstack_t **, lstack_t *, lstack_t *, int, int, int);
char    ** read_file    (FILE *, int *);

#endif
