#include "plutocore.h"

typedef struct searcharg_t sarg_t;
typedef struct ustack_t ustack_t;
typedef struct cstack_t cstack_t;

struct ustack_t {
   uint lim;
   uint pos;
   uint u[];
};

struct cstack_t {
   uint lim;
   uint pos;
   char c[];
};

struct searcharg_t {
   uint               tau;
   uint               trail;
   uint               query;
   uchar            * tree;
   uint             * lut;
   uint             * index;
   struct ustack_t  * hits;
   struct ustack_t ** milestones;
   struct cstack_t ** cstack;
};

