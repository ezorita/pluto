#include "plutocore.h"
#include <execinfo.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

#define MAXTREEQUERY 4

typedef struct searcharg_t sarg_t;
typedef struct usnap_t     usnap_t;

struct usnap_t {
   uint   lastid;
   uint   pos;
   uint   lim;
   uint * u[];
};

struct searcharg_t {
   int                tau;
   char               trail;
   uint               query;
   uchar            * tree;
   uint             * lut;
   uint             * index;
   struct ustack_t ** hits;
   struct usnap_t  ** milestones;
   struct cstack_t ** cstack;
};

// Function headers;

void _search        (uint, uchar *, sarg_t *);
void save_milestone (uint, uint, uchar *, ustack_t **, cstack_t **);
