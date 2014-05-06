#define _GNU_SOURCE
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
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

// Type definitions
typedef unsigned int uint;

// Shared functions headers.
uint * seqtoid       (char *, int *);
char * idtoseq       (uint);
uint   nodeaddr      (uint);
uint   getloci       (uint, uint *, uint *, uint **);

#endif
