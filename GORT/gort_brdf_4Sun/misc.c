#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "gort1.h"

/*
 * miscellaneous functions that didn't quite fit anywhere else
 */
fp_t           *alloc_1d();
fp_t          **alloc_2d();
fp_t            sec();

void           free_1d();
void           free_2d();
void           free_3d();

fp_t           *
alloc_1d(d)
int d;
{
   static fp_t    *r;

   r = (fp_t *) calloc((unsigned) d, sizeof(fp_t));
   if (!r) {
      fprintf(stderr, "Memory allocation failed (alloc_1d)\n");
      exit(1);
   }
   return (r);
}


fp_t          **
alloc_2d(d1, d2)
int d1, d2;
{
   int             q;
   static fp_t   **r;

   r = (fp_t **) calloc((unsigned) d1, sizeof(fp_t *));
   if (!r) {
      fprintf(stderr, "Memory allocation failed (alloc_2d)\n");
      exit(1);
   }
   for (q = 0; q < d1; q++) {
      r[q] = alloc_1d(d2);
   }
   return (r);
}


fp_t
sec(th)
fp_t th;
{
   return (1.0 / cos(th));
}


void
free_1d(a)
  void  *a;
{
   free(a);
}

void
free_2d(a)
  void  **a;
{
   free(a[0]);
   free(a);
}

void
free_3d(a)
  void  ***a;
{
   free(a[0][0]);
   free(a[0]);
   free(a);
}


