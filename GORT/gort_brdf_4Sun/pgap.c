
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gort1.h"

#define SIZE 80


/*
 * structure for getting values from parameter file
 */

extern void   gort1(double epgap[], double g_egap[][3], double *opentot);
extern void     free_1d();
extern void     free_2d();
extern void     free_3d();

static void     initialize();
static void     alloc_arrays();


GLOBAL_T g;

 fp_t           *height;	/* array of values for height */
 fp_t           *height_p;	/* array of values for height */
 fp_t           *theta;	/* array of values for theta */
 fp_t           *theta_p;	/* array of values for theta prime */
 fp_t          **P_n0;	/* array of values for P(n = 0 | h) */
 fp_t          **P_s0;	/* array of values of P(s = 0 | h) */
 fp_t          **T_open;
 fp_t          **dT_open;
 fp_t          **V_g;	/* array of values of V sub gamma */
 fp_t          **s_p;	/* array of s' values */
 fp_t          **EPgap;	/* gap probability table */
 fp_t          *K_open;        
 fp_t          *K_openEP;           
 fp_t        ***PD_s;	/* 3-dimensional array holding values of P(s
				 * | h, theta) */
 fp_t          *ES;  	/* array of values of E(S(h)) */
 fp_t          *factorial;
 fp_t          *Vb;
 fp_t          **fB;
 fp_t          *dK_open;
 fp_t          *Lk_up;
 fp_t          *Lk_down;

pgap(double epgap[], double g_egap[][3], double *opentot)
{
  int kk;

   factorial = alloc_1d(MAXCROWNS + 1);   
   /* calculate factorial */ 
   factorial[0] = 1;
   for (kk = 1; kk <= MAXCROWNS; kk++) {
      factorial[kk] = factorial[kk - 1] * (fp_t) kk;
   }

   /* initialize the global variables  */ 
   initialize();

   /* Allocate arrays  */ 
   alloc_arrays();

   /* call gort1()to calculate the gap probality and openness factors */
   gort1(epgap,g_egap, opentot);

   free_1d(Vb);
   free_1d(K_open);
   free_1d(dK_open);
   free_1d(Lk_up);
   free_1d(Lk_down);
   free_1d(K_openEP);
   free_1d(ES);

/*   free_1d(factorial);
 */
   free_1d(height);
   free_1d(height_p);
   free_1d(theta);
   free_1d(theta_p);

   free_2d(fB,g.nlayers);
   free_2d(T_open,g.nlayers);
   free_2d(dT_open,g.nlayers);
   free_2d(s_p,g.nlayers);
   free_2d(V_g,g.nlayers);
   free_2d(P_s0,g.nlayers);
   free_2d(P_n0,g.nlayers);
   free_2d(EPgap,g.nlayers);
   free_3d(PD_s, g.nlayers, g.nth);

   return (0);
}

static void
initialize()
{
    int h, t;
   /*
    * Initialize global variables.
    */

   /*
    * As described on the command-line, crowns are considered ellipsoids. The
    * model treats them as spheres.  We need to effectively alter the
    * vertical range so that the vertical axis of the ellipsoids becomes
    * equal to the horizontal axis (R).  This is accomplished by dividing h1
    * and h2 by the ellipticity parameter (b/r).  This compression in turn
    * changes the solar zenith angle.  We make all necessary adjustments
    * here.
    */

   g.dz_p = g.dz / g.ellipticity ; 

   g.R_squared = g.R * g.R;
   g.R_cubed = g.R_squared * g.R;
   g.H =  2.0 * g.R * g.ellipticity  + g.h2 - g.h1;

   /*
    * figure out tau.
    */
   /* g.ELAI = g.FAVD * (1.333333*g.density*PI*g.ellipticity*g.R_cubed); original one  */
   g.FAVD=g.ELAI/((1.333333)*PI*g.density*g.R*g.R*g.R*g.ellipticity);
   
   g.tau  = g.k * g.FAVD;
   
   /* 
    * While h1 and h2 define the volume in which we find crown centers, we
    * need to integrate from h1-r to h2+r in order to consider all heights at
    * which scattering by foliage may occur.  These will be referred to as z1
    * and z2.
    */

   g.z1 = g.h1 - g.R *  g.ellipticity;
   g.z2 = g.h2 + g.R * g.ellipticity;
   g.Lv = g.density / (g.h2 - g.h1);

  /* Corresponded variables in transformed space */

   g.FAVD_p = g.FAVD * g.ellipticity;
   g.tau_p  = g.k * g.FAVD_p;
   g.Lv_p = g.Lv * g.ellipticity;

   g.z1_p = g.z1 / g.ellipticity;
   g.z2_p = g.z2 / g.ellipticity;
   g.h1_p = g.h1 / g.ellipticity;
   g.h2_p = g.h2 / g.ellipticity;

   g.nlayers = (int) ((g.z2 - g.z1) / g.dz);
   g.nlayers++;
   
   g.nh1 =   (int) ((g.h1 - g.z1) / g.dz) + 1 ;
   g.nh2 =   (int) ((g.h2 - g.z1) / g.dz) + 1 ;
   
   /*
    * The following variable relates to the discrete approximation to P(s) at
    * a given height.  The variable 'ds' is used in the conversions between
    * real values of 's' and array indices.
    */

   /*
    * Figure out the number of angles theta (nth) for use in
    * various numerical integrations. 
   g.nth = (int) (DEG_TO_RAD(90.0) / g.dth + 0.5);
   */
   g.nth = (int) (DEG_TO_RAD(90.0) / g.dth + 0.5) + 1;


   /*
    * Initialize height and theta arrays.
    */

   height = alloc_1d(g.nlayers);
   height_p = alloc_1d(g.nlayers);
   theta = alloc_1d(g.nth);
   theta_p = alloc_1d(g.nth);

   for (h = g.nlayers - 1; h >= 0; h--) {
      height[h] = g.z2 - g.dz * (fp_t) (g.nlayers - 1 - h);
      height_p[h] = height[h]/g.ellipticity ; 
   }

   for (t = 0; t < g.nth; t++) {
      theta[t] = g.dth * (fp_t) t;
      if( theta[t] >= PI/2.0 ) { 
          theta[t] = PI/2.0  - 1.0*PI/180;
      }

   /*
    * Get theta prime -- original angles adjusted by ellipticity.
    */
      theta_p[t] = atan(tan(theta[t]) * g.ellipticity);
      if( theta_p[t] >= PI/2.0 ) { 
          theta_p[t] = PI/2.0 - 1.0*PI/180;
      }
      /* update ds', 
       * ds = (b/R) (cos(theta_p[t]/cos[theta[t]) ds'  
      g.ds_p = g.ds*cos(theta[t])/cos(theta_p[t]) ;    
      g.ds_p /= g.ellipticity ;
      */
   }
}


static void
alloc_arrays()
{
   int           h, t;
   int           s_max;

   /*
    * Allocate arrays.
    */

   Vb = alloc_1d(g.nlayers);
   K_open = alloc_1d(g.nlayers);
   dK_open = alloc_1d(g.nlayers);
   Lk_up = alloc_1d(g.nlayers);
   Lk_down = alloc_1d(g.nlayers);
   K_openEP = alloc_1d(g.nlayers);
   ES = alloc_1d(g.nlayers);
/*   factorial = alloc_1d(MAXCROWNS + 1);
 */
   fB = alloc_2d(g.nlayers,g.nth);
   T_open = alloc_2d(g.nlayers,g.nlayers);
   dT_open = alloc_2d(g.nlayers,g.nlayers);
   s_p = alloc_2d(g.nlayers, g.nth);
   V_g = alloc_2d(g.nlayers, g.nth);
   P_s0 = alloc_2d(g.nlayers, g.nth);
   P_n0 = alloc_2d(g.nlayers, g.nth);
   EPgap = alloc_2d(g.nlayers, g.nth);
   
   /*
    * Allocate the probability distribution table.
    */
   PD_s = (fp_t ***) calloc(g.nlayers, sizeof(fp_t **));
   if (!PD_s) {
      fprintf(stderr, "can't get memory for PD_s array\n");
      exit(1);
   }
   for (h = 0; h < g.nlayers; h++) {
      PD_s[h] = (fp_t **) calloc(g.nth, sizeof(fp_t *));
      if (!PD_s[h]) {
	 fprintf(stderr, "out of memory\n");
	 exit(1);
      }
      for (t = 0; t < g.nth; t++) {
   /*
    * figure out maximum path length at this height and this angle.
    *
    * Allocate enough memory to hold all possible s values.	 
    */
 
    s_max = s_to_index((g.z2 - height[h]) / cos(theta_p[t]));
    PD_s[h][t] = alloc_1d( s_max + 3);

 /*    s_max = (g.z2 - height[h]) / cos(theta_p[t]);	 
	 PD_s[h][t] = alloc_1d((int) (s_max / g.dth) + 3); */


      }
   }
}
