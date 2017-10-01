#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gort1.h"

/*
 * gort1 -- program calculates the vertical distribution of sunlit crown
 * surface in a 3-dimensional canopy, using the Li-Strahler geometric-optical
 * canopy model.  We do lots of numerical integration over 'layers' of the
 * canopy, between the lower and upper bounds of foliage.
 */

double           gort1(double epgap[],double g_egap[][3], double *openep);

extern fp_t     crown_proj_cross_section();
extern fp_t     crown_proj_volume();   

GLOBAL_T g;

double
gort1(double epgap[],double g_egap[][3], double *opentot)
{
   int             h;	/* 'height' index */
   int             t;	/* 'theta' index */
   int		   ord;
   int             i;
   int             s_max_p;
   fp_t             tmp;

   for (t = 0; t < g.nth; t++) {
      for (h = 0; h < g.nlayers; h++) {
	 /*
	  * Find P(n=0), using Vgamma.
	  */
	 V_g[h][t] = crown_proj_volume(height_p[h], theta_p[t]);
	 /*
         tmp =0.7*cos(theta[t])*cos(theta[t]) + 0.3; 
	 */
         tmp = 1.0; 
	 P_n0[h][t] = exp(-1.0 * tmp *g.Lv_p * V_g[h][t]);
      }
   }


   /*
    * Get P(s=0) for all heights.  It's just the differential of P(n=0) for
    * decreasing height.
    */
   for (t = 0; t < g.nth; t++) {
      P_s0[g.nlayers - 1][t] = 0.0;
      for (h = g.nlayers - 2; h >= 0; h--) {
	 P_s0[h][t] = P_n0[h + 1][t] - P_n0[h][t];
      }
   }

   /*
    * Get probability distributions for s for each possible h and theta.
    */
   for (t = 0; t < g.nth; t++) {  
      for (h = 0; h < g.nlayers; h++) {  
   /*
    * figure out maximum path length at this height and this angle.
    */
         s_max_p = s_to_index( (g.z2_p - height_p[h]) / cos(theta_p[t]));
         for (i = 0; i < ( s_max_p + 3 ) ; i++) {  
             PD_s[h][t][i] = 0.0; }
	 get_PD_s(h, t);
      }
   }

   /*
    * Get gap probabilities.
    */

   get_EPgap(epgap,g_egap);
    *opentot = get_K_open();
   
}


