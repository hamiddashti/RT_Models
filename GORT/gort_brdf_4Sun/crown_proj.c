#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gort1.h"


fp_t            crown_proj_cross_section();
static fp_t     wierd_cross_section();
static fp_t     left_circle_area();
static fp_t     right_ellipse_area();


fp_t
crown_proj_volume(h, th)
fp_t h, th;
{
   /*
    * Calculate crown projection volume for a crown at the given height, by
    * numerical integration, using the midpoint rule.
    */

   fp_t            vol = 0.0;
   fp_t            z;

   for (z = g.h1_p + g.dz_p / 2.0; z <= g.h2_p; z += g.dz_p) {
      vol += crown_proj_cross_section(h, z, th) * (g.dz_p);
   }

   return (vol);
}


fp_t
crown_proj_cross_section(h, z, th)
fp_t h, z, th;
{
   /* all the calculation is in transformed dimension */
   /*
    * Calculate and return the cross-sectional area of the crown projection
    * for height h at height z.  The shape is very strange.  For a given
    * height h, the horizontal cross section at z < h - r sin (th) is
    * circular. For z > h + r sin (th), the cross section is elliptical.
    * In between, the shape is part circle and part ellipse.
    */

   fp_t            h_low;
   fp_t            h_high;
   fp_t            csa;		/* cross-sectional area */

#ifdef NARRATE
   printf("h: %8.4lf  z: %8.4lf\n", h, z);
#endif

   if (z < h - g.R)
      return (0.0);

   h_low = h - g.R * sin(th);
   h_high = h + g.R * sin(th);

   if (z <= h_low) {

      /*
       * Cross-section is circular.
       */
      fp_t            r_p;

      r_p = sqrt(g.R_squared - (h - z) * (h - z));
      csa = PI * r_p * r_p;
   }

   else if (z > h_low && z < h_high) {

      /*
       * Cross-section is neither a circle nor an ellipse.
       */

      csa = wierd_cross_section(h, z, th);
   }

   else {

      /*
       * Cross-section is an ellipse.
       */

      csa = PI * g.R_squared * sec(th);
   }

#ifdef NARRATE
   printf("h: %8.4lf  z: %8.4lf  csa: %8.4lf\n", h, z, csa);
#endif

   return (csa);
}


static          fp_t
wierd_cross_section(h, z, th)
fp_t h, z, th;
{
   fp_t            r_p;		/* r' - radius of the circular part of the
				 * cross section at the height 'z' */
   fp_t            x_cc;	/* x-coord of the circle center, relative to
				 * the ellipse center */
   fp_t            x_p;		/* x-coord of the line separating the circle
				 * part and the ellipse part of the cross
				 * section. */
   fp_t            A_cp;	/* area of the circle part of the cross
				 * section. */
   fp_t            A_ep;	/* area of the ellipse part of the cross
				 * section. */
   fp_t            A;		/* total area */
   fp_t            zdiff;

   zdiff = h - z;

   r_p = sqrt(g.R_squared - zdiff * zdiff);
   x_cc = zdiff * tan(th);
   x_p = x_cc / (1.0 - cos(th) * cos(th));

   /*
    * So now we know the cut-off line for each partial area. We need to find:
    * 
    * 1. The area of the circle of radius r', centered at zero, and cut off by
    * the line x = x_p - x_cc.
    * 
    * 2. The area of the ellipse to the right of the line x = x_p.
    * 
    * Of course, the terms 'left' and 'right' here are only relative.
    */

   A_cp = left_circle_area(r_p, x_p - x_cc);
   A_ep = right_ellipse_area(g.R, g.R * sec(th), x_p);
   A = A_cp + A_ep;

   return (A);
}


static          fp_t
left_circle_area(r, x_cut)
fp_t r, x_cut;
{
   fp_t            Area_tot;
   fp_t            Ang_sector;
   fp_t            Area_sector;
   fp_t            Area_triangle;
   fp_t            A;

   Area_tot = PI * r * r;
   Ang_sector = acos(fabs(x_cut) / r) * 2.0;
   Area_sector = Area_tot * Ang_sector / (2.0 * PI);
   Area_triangle = fabs(x_cut) * sqrt(r * r - x_cut * x_cut);
   if (x_cut > 0.0)
      A = Area_tot - (Area_sector - Area_triangle);
   else
      A = Area_sector - Area_triangle;

   return (A);
}


static          fp_t
right_ellipse_area(r,  b, x_cut)
fp_t r, b, x_cut;
{
   fp_t            x_cut_p;
   fp_t            A_p;

   x_cut_p = x_cut / (b / r);
   A_p = PI * r * r;

   A_p -= left_circle_area(r, x_cut_p);
   return (A_p * (b / r));
}



