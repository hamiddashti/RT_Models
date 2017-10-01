#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gort1.h"

extern fp_t    *height;
extern fp_t    *theta_p; 

extern fp_t   **P_n0;
extern fp_t  ***PD_s;

 
 void     get_Vb();
 void     get_fB();
 void     get_T_open(); 

static fp_t     get_Pgap();
GLOBAL_T g;


void
get_EPgap(double epgap[],double g_egap[][3] )
{
   /*
    * Fill the EPgap array.  each element of the table holds the expected gap
    * probability for a given height and angle.
    */
   int             h;
   int             t;
   int             s;
   int             maxs_p;
   fp_t            sum_Ps; 
   fp_t            m_s;

   get_Vb();
   get_fB();
   get_T_open();
   
   for (h = 0; h < 1; h++) {
      for (t = 0; t < g.nth-1; t++) {
         EPgap[h][t] = 0.0;
         sum_Ps = 0.0;
         m_s = 0.0;
	 /*
	  * Integrate over path lengths.
	  */

	 maxs_p = s_to_index((g.z2_p - height_p[h]) / cos(theta_p[t]));  
	 for (s = 0; s <= maxs_p;  s++) {
            m_s += index_to_s(s) * PD_s[h][t][s] ; 
	    sum_Ps += PD_s[h][t][s] ;
	    EPgap[h][t] += get_Pgap(index_to_s(s), theta[t]) * PD_s[h][t][s] ;
	 }

         
	 /* if(h==0)   printf("\n %4.2f %8.4f %8.4f", 
               RAD_TO_DEG(theta[t]),P_n0[h][t],EPgap[h][t]); 
	 */

         epgap[t]= EPgap[0][t];
         g_egap[t][0]= RAD_TO_DEG(theta[t]);
         g_egap[t][1]= P_n0[h][t];
         g_egap[t][2]= EPgap[0][t];
      } 
   }
}


static          fp_t
get_Pgap(s, th)
fp_t s, th;
{
   /*
    * Get a gap probability given a pathlength and a zenith angle. At the
    * moment, the zenith angle is not used because we are assuming a
    * spherical LAD.
    */
   return (exp(-s * g.tau_p));  
/*   return (exp(-s * g.tau)); */
}



void
get_Vb()
{
   /*
    * Fill the Vb array.  Vb(h) is the volume of a sphere centered at height
    * h, intersected by h1 and h2 planes.  So we first calculate the total
    * volume, then subtract the volumes of intersected sectors, if necessary.
    */

   int             i;
   fp_t            Vol;
   fp_t            tmp;

   for (i = 0; i < g.nlayers; i++) {
      Vol = 4.0 * PI * g.R_cubed / 3.0;

      if (height_p[i] + g.R > g.h2_p ) {
	 tmp = height_p[i] + g.R - g.h2_p;
	 Vol -= PI * tmp * tmp * (3.0 * g.R - tmp) / 3.0;
      }
      if (height_p[i] - g.R < g.h1_p) {
	 tmp = g.h1_p - (height_p[i] - g.R);
	 Vol -= PI * tmp * tmp * (3.0 * g.R - tmp) / 3.0;
      }

      Vb[i] = Vol;
   }
}


 void
get_fB()
{
   int             i,t;

   for(t=0; t<g.nth; t++) {
   for (i = 0; i < g.nlayers; i++) {
      fB[i][t] = (1.0 - exp(-g.Lv_p * Vb[i])) / (1.0 - P_n0[i][t]);
   }
   }
}




void 
get_T_open()
{
   /*
    * Calculate the mean passing through crown gap
    */

   int    h,z, t, n, k;
   fp_t   T, dT;
   fp_t  temp1, ds, s, s_p, P_n;

   factorial[0] = 1;
     for (k = 1; k <= MAXCROWNS; k++) {
       factorial[k] = factorial[k - 1] * (fp_t) k;
   }

 for(z = 0; z<g.nlayers; z++) { 
  for(h = g.nlayers-1; h >= z  ; h--) {  
    T_open[h][z] = 0.0;
    dT_open[h][z] = 0.0;
    if( z != h ) {
     for(t = 0; t < g.nth; t++) { 
       T = 0.0; 
       dT = 0.0;
       s_p = fabs(height_p[z] - height_p[h])/cos(theta_p[t]) ;  

       for (n = 1; n <= MAXCROWNS; n++) {
          s = s_p *(1.0 - exp(-n*get_ES(z,t)/s_p) );
          temp1 = g.Lv_p * PI * g.R * g.R * s_p; 
          P_n = (pow(temp1, (fp_t) n) * exp(-temp1))/
	      (factorial[n] * (1.0 - exp(-temp1))); 
          T += P_n * exp(- s * g.tau_p);
          ds = (1.0 - exp(-g.Lv_p * Vb[z]))*g.dz_p/cos(theta_p[t]);
          dT += P_n * exp(-s * g.tau_p)*(1.0-exp(-g.tau_p*ds));
       }
       T_open[h][z] += sin(2.0*theta[t]) * T * g.dth;
       dT_open[h][z] += sin(2.0*theta[t]) * dT * g.dth;
       T_open[z][h] = T_open[h][z];
       dT_open[z][h] = dT_open[h][z];
     }
    }
    else { 
     for(t = 0; t < g.nth; t++) { 
       ds = 0.5*(1.0 - exp(-g.Lv_p * Vb[z]))*g.dz_p/cos(theta_p[t]);
       dT = 1.0-exp(-g.tau_p * ds) ;
       dT_open[h][z] += sin(2.0*theta[t]) * dT * g.dth;
     }
    }
/*
    printf("\n T_open(%4.1f->%4.1f): %6.4f %6.4f", height[h], 
        height[z],T_open[h][z],dT_open[h][z] );
	*/
  }
 }

}




