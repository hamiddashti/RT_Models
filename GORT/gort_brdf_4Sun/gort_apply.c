#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gort1.h"

#define SIZE 80

GLOBAL_T g;

void
gort_apply(double *brdfValue)
{
 
   int i,j;
   int sizeth = (int) (DEG_TO_RAD(90.0) / g.dth + 0.5) + 1;    /* array size of epgap  */
   double epgap[sizeth];       /* return 1-d epgap array from pgap function  */
   double g_egap[sizeth][3];   /* return 2-d g_epgap array from pgap function */
   double opentot = 0.0;        /* return the total value from pgap function  */
   double reValue = 0.0;       /* return value from calulating gort brdf */ 
   
   /* call pgap() fuction to calculate the gap probality and openness factors
      return two arrays to respective the gap probality. one dimension array: epgap [pgap],
      two dimension array: g_egap array[zenith_angle, within_crown_pgap, between_crown_pgap] */
   pgap(epgap, g_egap, &opentot);
   

   /* print the total gap probality   */   
   /*
   printf ("the openness factors and total pgap value:  %10.3f %10.3f %10.3f\n", g.k_open, g.k_openEP, g.k_opentot);
   */

   /* print return array, between crown pgap  */
   /*
   printf("Epgap array: \n");
   for (j = 0; j < sizeth-1; j++){
             printf(" %lf\t", epgap[j]);
       }
       printf("\n");
   */
 
  /* call gort_brdf() function to calculate the gort brdf,return one value reValue  */
       gort_brdf(g_egap, &reValue);
       *brdfValue = reValue;


  /* print the return gort brdf reValue */
  /*
  printf("the return gort_brdf value: %lf\t", reValue);
  printf("the gort_brdf value: %lf\t", brdfValue);
  printf("the g.brdf value: %lf\n", g.brdf);
  */
 
}

