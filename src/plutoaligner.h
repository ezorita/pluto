#include "plutocore.h"
#include <execinfo.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

typedef struct searcharg_t sarg_t;

struct searcharg_t {
   int                tau;
   uint               trail;
   uint               query;
   uchar            * tree;
   uint             * lut;
   uint             * index;
   struct ustack_t  * hits;
   struct ustack_t ** milestones;
   struct cstack_t ** cstack;
};

// Function headers;

void _search        (uint, uchar *, sarg_t *);
void save_milestone (uint, uint, uchar *, ustack_t **, cstack_t **);
