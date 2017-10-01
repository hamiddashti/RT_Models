#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gort1.h"


static fp_t   get_S();
static fp_t   get_Pcc();

extern fp_t   crown_proj_cross_section();


fp_t
get_ES(z, th)
int z, th;
{
   /*
    * Return the expected value of S(z, h) -- on average, the distance that a
    * beam passes through a single crown in order to reach the height z. This
    * could be seriously optimized.  But anyway, we integrate over heights at
    * which such a sphere could be centered.
    */
   int             nh = 20;	/* number of heights */
   fp_t          h;
   fp_t          dh;
   fp_t          ES = 0.0;
/*
   for (h = g.nh1; h <= g.nh2; h++) {
      ES += get_S(z, h, th) * (get_Pcc(height[h]) * g.dz);
   }
 */

   dh = (g.h2_p - g.h1_p) / (fp_t) nh;
   for (h = g.h1_p + dh / 2.0; h <= g.h2_p; h += dh) {
      ES += get_S(z, h, th) * (get_Pcc(h) * dh);
   }

   return (ES);
}


static fp_t
get_S(z, h, th)
int z, th;
fp_t h;
{
      fp_t          th_p, S;
/*      th_p = atan(tan(th) * g.ellipticity);  */
   /*
    * Return the average distance a beam must pass through a single crown
    * centered at height h to reach height z.  We calculate and return the
    * volume of a sphere centered at h above the height z, divided by its
    * projected area.
    */

/*   if (height_p[z] > height_p[h] + g.R - 0.0001) { */
   if (height_p[z] > h + g.R - 0.0001) {
      /*
       * Plane is above the sphere.  Average distance is zero.
       */
      S = 0.0;
   }
/*   else if (height_p[z] < height_p[h] - g.R + 0.0001) {  */
   else if (height_p[z] < h - g.R + 0.0001) {
      /*
       * Plane is below the sphere.  Average distance is the volume of the
       * sphere divided by its projected area.
       */

      S = 4.0 * g.R / 3.0;
/*      S *= cos(theta_p[th]) * g.ellipticity / cos(theta[th]);  */
   }
   else {
      fp_t          proj_area = 0.0;
      fp_t          r_p = 0.0;
      fp_t          V_sphere = 0.0;
      fp_t          V_slice = 0.0;
      fp_t          V_tot = 0.0;
      fp_t          ht = 0.0;
      fp_t          zdiff = 0.0;

/*      h /= g.ellipticity;
      z /= g.ellipticity;
 */

      V_sphere = 4.0 * PI * g.R_cubed / 3.0;

      zdiff = fabs(h - height_p[z]);
      r_p = sqrt(g.R_squared - zdiff * zdiff);

      ht = g.R - zdiff;
      V_slice = PI * ht * ht / 3.0 * (3.0 * g.R - ht);

      if (height_p[z] > h)
	 V_tot = V_slice;
      else
	 V_tot = V_sphere - V_slice;

      /* Adjust V_tot by the solar zenith angle */
      V_tot /= cos(theta_p[th]);

      /*
       * Get projected area. The 'crown_proj_cross_section' function assumes
       * that we are projecting _towards_ the sun, as is the case for the
       * calculation of V_gamma and other things.  But here we are projecting
       * _away_ from the sun, so we flip things around.
       */

      if (h < height_p[z]) {
	 proj_area = crown_proj_cross_section(h, (h - zdiff), theta_p[th]);
      }
      else {
	 proj_area = crown_proj_cross_section(h, (h + zdiff), theta_p[th]);
      }

      S = V_tot / proj_area;

 /*     S *= cos(theta_p[th]) * g.ellipticity / cos(theta[th]);  */

#ifdef NARRATE
      printf("z: %8.4f  h: %8.4f  V_tot: %8.4f  proj: %8.4f  S: %8.4f\n",
	     z,
	     h,
	     V_tot,
	     proj_area,
	     S);
#endif
   }
   return (S);
}


static          fp_t
get_Pcc(h)
fp_t h;
{
   /*
    * return the probability of finding a crown center at the height 'h'
    */
   return (1.0 / (g.h2_p - g.h1_p)); 
/*   return (1.0 / (g.h2 - g.h1)); */
}


