#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include "plutocore.h"

#ifndef _MERGESORT_SRC__
#define _MERGESORT_SRC__

// Type definitions.
typedef struct sortargs_t sortargs_t;

struct sortargs_t {
   void ** buf0;
   void ** buf1;
   int     (*compar)(const void*, const void*);
   int     size;
   int     b;
   int     thread;
   int     repeats;
};

// Function headers.
void   insertion_sort     (loc_t *, int, int);
void   radix_sort         (loc_t *, int, int, int);
void   mergesort_lstack   (lstack_t **);
void   insert_loci        (lstack_t **, loc_t *, int, char);
void   mergesort_loc_int  (loc_t *, int *, int, int);
void   _mergesort_loc_int (loc_t *, loc_t *, int *, int, int);
void   mergesort_loc      (loc_t *, int);
void   _mergesort_loc     (loc_t *, loc_t *, int);
int    mergesort          (void **, int, int (*)(const void*, const void*), int);
void * _mergesort         (void *);
int    ualpha             (const void *, const void *);

#endif
