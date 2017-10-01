#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gort1.h"

/*
 * */


void
get_Lk()
{
   int             i;

   /*
    * Figure out upward and downward leakage.  These quantities are the
    * proportion of scattered radiation that escapes through gaps between
    * crowns from a given height.  Since we make the assumption that crown
    * centers are distributed symetrically about the height (h1 + h2) / 2,
    * the downward array is the reverse of the upward array.
    */
   for (i = 0; i < g.nlayers - 1; i++) {
      Lk_up[i] = Lk_down[g.nlayers - i - 1] =
	 dK_open[i] / (1.0 - K_open[i]);
   }
   Lk_up[g.nlayers - 1] = Lk_down[0] = 1.0;
}
