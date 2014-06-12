#include "mergesort.h"

void *
_mergesort
(
 void * args
)
{
   sortargs_t * sortargs = (sortargs_t *) args;
   if (sortargs->size > 1) {
      // Next level params.
      sortargs_t arg1 = *sortargs, arg2 = *sortargs;
      arg1.size /= 2;
      arg2.size = arg1.size + arg2.size % 2;
      arg2.buf0 += arg1.size;
      arg2.buf1 += arg1.size;
      arg1.b = arg2.b = (arg1.b + 1) % 2;

      // Either run threads or DIY.
      if (arg1.thread) {
         // Decrease one level.
         arg1.thread = arg2.thread = arg1.thread - 1;
         // Create threads.
         pthread_t thread1, thread2;
         pthread_create(&thread1, NULL, _mergesort, (void *) &arg1);
         pthread_create(&thread2, NULL, _mergesort, (void *) &arg2);
         // Wait for threads.
         pthread_join(thread1, NULL);
         pthread_join(thread2, NULL);
      }
      else {
         _mergesort((void *) &arg1);
         _mergesort((void *) &arg2);
      }

      // Separate data and buffer (b specifies which is buffer).
      void ** l = (sortargs->b ? arg1.buf0 : arg1.buf1);
      void ** r = (sortargs->b ? arg2.buf0 : arg2.buf1);
      void ** buf = (sortargs->b ? arg1.buf1 : arg1.buf0);

      // Accumulate repeats
      sortargs->repeats = arg1.repeats + arg2.repeats;

      // Merge sets
      int i = 0, j = 0, nulli = 0, nullj = 0, cmp = 0, repeats = 0;
      for (int k = 0, idx = 0, n = 0; k < sortargs->size; k++) {
         if (j == arg2.size) {
            // Insert pending nulls, if any.
            for (n = 0; n < nulli; n++)
               buf[idx++] = NULL;
            nulli = 0;
            buf[idx++] = l[i++];
         }
         else if (i == arg1.size) {
            // Insert pending nulls, if any.
            for (n = 0; n < nullj; n++)
               buf[idx++] = NULL;
            nullj = 0;
            buf[idx++] = r[j++];
         }
         else if (l[i] == NULL) {
            nulli++;
            i++;
         }
         else if (r[j] == NULL) {
            nullj++;
            j++;
         }
         else if ((cmp = sortargs->compar(l[i],r[j])) == 0) {
            j++;
            repeats++;
            // Insert sum of repeats as NULL.
            for (n = 0; n <= nulli + nullj; n++) {
               buf[idx++] = NULL;
            }
            nulli = nullj = 0;
         } 
         else if (cmp < 0) {
            // Insert repeats as NULL.
            for (n = 0; n < nulli; n++)
               buf[idx++] = NULL;
            nulli = 0;

            buf[idx++] = l[i++];
         }
         else {
            // Insert repeats as NULL.
            for (n = 0; n < nullj; n++)
               buf[idx++] = NULL;
            nullj = 0;

            buf[idx++] = r[j++];
         }
      }
      sortargs->repeats += repeats;
   }
   
   return NULL;

}

int
mergesort
(
 void ** data,
 int numels,
 int (*compar)(const void*, const void*),
 int maxthreads
)
{
   // Copy to buffer.
   void ** buffer = (void **)malloc(numels * sizeof(void *));
   memcpy(buffer, data, numels * sizeof(void *));

   // Prepare args struct.
   sortargs_t args;
   args.buf0   = data;
   args.buf1   = buffer;
   args.size   = numels;
   args.b      = 0; // Important so that sorted elements end in data (not in buffer).
   args.thread = 0;
   args.repeats = 0;
   args.compar = compar;
   while ((maxthreads >> (args.thread + 1)) > 0) args.thread++;

   _mergesort((void *) &args);

   free(buffer);
   
   return numels - args.repeats;
}



void
_mergesort_loc_int_nonrecursive
(
 loc_t * data,
 loc_t * buf,
 int   * intervals,
 int     numints,
 int     numels
)
{
   int merges = 0;
   while ((1 << merges) < numints) merges++;
  
   if (merges%2) {  
      // Swap buffers if num of merges is odd.
      loc_t * tmp = data;
      data = buf;
      buf  = tmp;    
   }
   
   for (int m = 0; m < merges; m++) {
      int offset = 2 << m;
      int next = offset / 2;
      int sets = numints/offset;
      int nextset = offset;

      // Special case for the last merge. (Data size may not be power of 2..)
      if (sets == 0 && numints > next) {
         nextset = numints % next;
         sets++;
      }

      for (int n = 0; n < sets; n++) {
         int sizel = intervals[n*offset + next] - intervals[n*offset];
         int sizer = intervals[n*offset + nextset] - intervals[n*offset + next];
         loc_t * l = data + intervals[n*offset];
         loc_t * r = data + intervals[n*offset + next];
         loc_t * b = buf + intervals[n*offset];
         // Merge sets.
         int k = 0, i = 0, j = 0;
         while (i < sizel && j < sizer) {
            if (l[i] > r[j]) {
               b[k++] = r[j++];
            }
            else {
               b[k++] = l[i++];
            }
         }
      
         // Copy remainders.
         if (i < sizel) memcpy(b + k, l + i, (sizel - i)*sizeof(loc_t));
         else           memcpy(b + k, r + j, (sizer - j)*sizeof(loc_t));
      }

      // Swap buffers.
      loc_t * tmp = data;
      data = buf;
      buf  = tmp;
   }
}


void
_mergesort_loc_int
(
 loc_t * data,
 loc_t * buf,
 int   * intervals,
 int     nints,
 int     numels
)
{
   if (nints < 2) return;

   int intsl = nints / 2;
   int intsr = nints - intsl;
   int sizel = intervals[intsl] - intervals[0];
   int sizer = numels - sizel;

   _mergesort_loc_int(buf, data, intervals, intsl, sizel);
   _mergesort_loc_int(buf + sizel, data + sizel, intervals + intsl, intsr, sizer);

   loc_t * l = data;
   loc_t * r = data + sizel; 

   // Merge sets.
   int k = 0, i = 0, j = 0;
   while (i < sizel && j < sizer) {
      if (l[i] > r[j]) {
         buf[k++] = r[j++];
      }
      else {
         buf[k++] = l[i++];
      }
   }
      
   // Copy remainders.
   if (i < sizel) memcpy(buf + k, l + i, (sizel - i)*sizeof(loc_t));
   else           memcpy(buf + k, r + j, (sizer - j)*sizeof(loc_t));
}

void insertion_sort(loc_t *array, int offset, int end) {
   int x, y;
   loc_t temp;
   for (x=offset; x<end; ++x) {
      for (y=x; y>offset && array[y-1]>array[y]; y--) {
         temp = array[y];
         array[y] = array[y-1];
         array[y-1] = temp;
      }
   }
}

void radix_sort(loc_t *array, int offset, int end, int shift) {
   int x, y;
   loc_t value, temp;
   loc_t last[256] = { 0 }, pointer[256];

   for (x=offset; x<end; ++x) {
      ++last[(array[x] >> shift) & 0xFF];
   }

   last[0] += offset;
   pointer[0] = offset;
   for (x=1; x<256; ++x) {
      pointer[x] = last[x-1];
      last[x] += last[x-1];
   }

   for (x=0; x<256; ++x) {
      while (pointer[x] != last[x]) {
         value = array[pointer[x]];
         y = (value >> shift) & 0xFF;
         while (x != y) {
            temp = array[pointer[y]];
            array[pointer[y]++] = value;
            value = temp;
            y = (value >> shift) & 0xFF;
         }
         array[pointer[x]++] = value;
      }
   }

   if (shift > 0) {
      shift -= 8;
      for (x=0; x<256; ++x) {
         temp = x > 0 ? pointer[x] - pointer[x-1] : pointer[0] - offset;
         if (temp > 64) {
            radix_sort(array, pointer[x] - temp, pointer[x], shift);
         } else if (temp > 1) {
            insertion_sort(array, pointer[x] - temp, pointer[x]);
         }
      }
   }
}

void
mergesort_loc_int
(
 loc_t * data,
 int   * intervals,
 int     nints,
 int     numels
)
{
   // Alloc buffer.
   loc_t * buf = (loc_t *)malloc(numels * sizeof(loc_t));
   memcpy(buf, data, numels*sizeof(loc_t));

   _mergesort_loc_int (buf, data, intervals, nints, numels);

   free(buf);
}


void
_mergesort_loc
(
 loc_t * data,
 loc_t * buf,
 int     numels
)
{
   int sizel = numels / 2;
   int sizer = numels / 2 + numels % 2;

   if (sizel > 1) _mergesort_loc(buf, data, sizel);
   if (sizer > 1) _mergesort_loc(buf + sizel, data + sizel, sizer);

   loc_t * l = data;
   loc_t * r = data + sizel; 

   // Merge sets.
   int k = 0, i = 0, j = 0;
   while (i < sizel && j < sizer) {
      if (l[i] > r[j]) {
         buf[k++] = r[j++];
      }
      else {
         buf[k++] = l[i++];
      }
   }
      
   // Copy remainders.
   if (i < sizel) memcpy(buf + k, l + i, (sizel - i)*sizeof(loc_t));
   else           memcpy(buf + k, r + j, (sizer - j)*sizeof(loc_t));
}


void
mergesort_loc
(
 loc_t * data,
 int     numels
)
{
   // Alloc buffer.
   loc_t * buf = (loc_t *)malloc(numels * sizeof(loc_t));
   memcpy(buf, data, numels*sizeof(loc_t));

   _mergesort_loc(buf, data, numels);

   free(buf);
}

void
insert_loci
(
 lstack_t ** lstackp,
 loc_t     * newloci,
 int         nloci,
 char        offset
)
{
   lstack_t * lstack = *lstackp;
   int cloci = lstack->pos;

   // Realloc stack if needed.
   size_t minsize = cloci + nloci;
   if (lstack->lim < minsize) {
      size_t newsize = lstack->lim*2;
      while (newsize < minsize) newsize *= 2;
      *lstackp = lstack = realloc(lstack, 3*sizeof(int) + sizeof(seq_t) + newsize*sizeof(loc_t));
      if (lstack == NULL) {
         fprintf(stderr, "error (realloc) in 'insert_loci': %s\n", strerror(errno));
         exit(EXIT_FAILURE);
      }
      lstack->lim = newsize;
   }

   // Merge.
   //loc_t * currloci = malloc(cloci*sizeof(loc_t));
   loc_t currloci[cloci];
   memcpy(currloci, lstack->u, cloci * sizeof(loc_t));

   int i=0, j=0, k=0;
   while (i < cloci && j < nloci) {
      if (currloci[i] < newloci[i] + offset)
         lstack->u[k++] = currloci[i++];
      else 
         lstack->u[k++] = newloci[j++] + offset;
   }

   // Copy remainders.
   if (i < cloci) memcpy(lstack->u + k, currloci + i, (cloci - i)*sizeof(loc_t));
   else {
      if (offset == 0)
         memcpy(lstack->u + k, newloci  + j, (nloci - j)*sizeof(loc_t));
      else
         for (; j < nloci; j++) 
            lstack->u[k++] = newloci[j];
   }
   
   lstack->pos += nloci;
   //free(currloci);
}

int
ualpha
(
   const void *a,
   const void *b
) 
{
   char *u1 = (char *)a;
   char *u2 = (char *)b;
   return strcmp(u2, u1);
}
