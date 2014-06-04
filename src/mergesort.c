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
