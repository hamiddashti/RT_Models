#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gort1.h"

/*
 * Get k_open and delta(K_open).  This involves integrating P(n=0) and P(s=0)
 * over theta.  We re-transform theta back to its original space -- that is,
 * not adjusted by crown ellipticity.  Use the trapezoidal rule.
 */ 


fp_t            tmp1, tmp1_last;
fp_t            tmp2, tmp2_last;
fp_t            tmp3, tmp3_last;      

GLOBAL_T g;

double
get_K_open()
{

   fp_t            gtau0;
   fp_t            gtau;
   int             h, t;
   double          re_openep;

   for (h = 0; h < g.nlayers; h++) {   

      K_open[h] = 0.0 ;
      K_openEP[h] = 0.0;
      dK_open[h] = 0.0 ;

      tmp1_last = P_n0[h][0] * sin(2.0*theta[0]);
      tmp2_last = EPgap[h][0] * sin(2.0*theta[0]);
      tmp3_last = P_s0[h][0] * sin(2.0*theta[0]);

      for (t = 1; t < g.nth; t++) {
	 tmp1 = P_n0[h][t] * sin(2.0*theta[t]);
	 K_open[h] += (tmp1 + tmp1_last) / 2.0 * g.dth;
	 tmp1_last = tmp1;

	 tmp2 = EPgap[h][t] * sin(2.0*theta[t]);
	 K_openEP[h] += (tmp2 + tmp2_last) / 2.0 * g.dth;
	 tmp2_last = tmp2;

	 tmp3 = P_s0[h][t] * sin(2.0*theta[t]);
	 dK_open[h] += (tmp3 + tmp3_last) / 2.0 * g.dth;
	 tmp3_last = tmp3;

#ifdef OLD
	 K_open[h] += P_n0[h][t] * sin(2.0*theta[t]) * g.dth;
	 dK_open[h] += P_s0[h][t] * sin(2.0*theta[t]) * g.dth;
#endif
      }
         gtau0 =  g.k * g.ELAI / ( 1.0-K_open[h] ) ;
         gtau0 /= g.z2 - g.z1;

         gtau =  g.k * g.ELAI ;
         gtau /= g.z2 - g.z1;
	
     /* if(h==0)   printf("\n%8.4f %8.4f",K_open[h],K_openEP[h]); */
     re_openep = K_open[0]+K_openEP[0];
     g.k_open =  K_open[0];
     g.k_openEP = K_openEP[0];
     g.k_opentot =  re_openep;
    return re_openep;
    
   }
}
