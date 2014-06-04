#include <pthread.h>
#include <stdlib.h>
#include <string.h>

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
int    mergesort  (void **, int, int (*)(const void*, const void*), int);
void * _mergesort (void *);
int    ualpha     (const void *, const void *);

#endif
