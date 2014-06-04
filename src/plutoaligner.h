#include "plutocore.h"
#include <execinfo.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

#define MAXLUTQUERY 2

// Data types:
typedef struct arg_t arg_t;

struct arg_t {
   int     nleaves;
   int     maxtau;
   char  * query;
   loc_t * lut;
   loc_t * index;
};

// Function headers;
void _search        (uint, uchar *, sarg_t *);
void save_milestone (uint, uint, uchar *, ustack_t **, cstack_t **);
