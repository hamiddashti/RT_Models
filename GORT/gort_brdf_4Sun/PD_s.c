#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gort1.h"



/*
 * Routines associated with the calculation of P(s|h, theta).
 */

/* static void     init(); */
static fp_t Vol();
static fp_t sector();    /* the intersected spherical  volume */ 
static fp_t triang();    /* the intersected triangle   volume */ 
static fp_t cylind();    /* the intersected cylinder   volume */
static fp_t fcn();       /* function used in trian()          */
static fp_t fun1();      /* function used in cylind           */ 
static fp_t trisec();     /* sum of sector and cylind          */


void
get_PD_s(h, t)
int h, t;
{
   /*
    * Get a probability distribution for the within-crown path length 's' at
    * the given height and zenith angle. This routine fills the array
    * PD_s[h][t][...].
    */

   static int      first_time = 0; 
 
   int             sp_i;	/* index to s' array */
   int             ps0_i;	/* index to P_s0 array */
   int             i;
   int             n, maxs;		/* number-of-crowns counter */
   fp_t            P_s_p;	/* probability associated with a given s'
				 * value at a given height. */
   fp_t            P_n;		/* probability of penetrating 'n' crowns at a
				 * given s' value */
   fp_t            temp1;
   fp_t            temp2;
   fp_t            s;
   fp_t            n_mean;       /*  mean crowns */
   fp_t            tmp_Lv; 

   /*
    * Make some preparations.
    

   if (first_time == 0) {
      init();
      first_time = 1 ;
    }
   
    * The ES[] array contains values indicating the average distance that a
    * beam must pass through a single crown in order to reach a given height.
    */
     ES[h] = get_ES(h,t);
/*     printf("\n %4.1f %5.2f ES(h)=%6.3f", RAD_TO_DEG(theta[t]), height[h], ES[h]); */

     n_mean = 0.0 ; 
   /*
    * Iterate over our possible s' (s_p[]) values.  The s' values are
    * after-entering-the-crown pathlengths, which are not equivalent to
    * within-crown pathlengths since they don't account for exiting from
    * crowns.
    */
 
   for (sp_i = g.nlayers - 1; sp_i > h; sp_i--) {
	 /*
	  * The s' array (s_p[][]) holds after-entering-the-crown pathlengths
	  * associated with a given zenith angle and within-canopy height.
	  */

      s_p[sp_i][t] = (height_p[sp_i] - height_p[h]) / cos(theta_p[t]); 
      if (sp_i == g.nlayers - 1) {
	 /*
	  * Path length is zero.  Calculations wouldn't work, but we know
	  * that the appropriate index is zero and we know what the
	  * probability is.
	  */         
	 PD_s[h][t][0] += P_s0[sp_i][t];
	 continue;
      }
 
      /*
       * Get the probability associated with the current value of s' at this
       * height.
       */
 /*       P_s_p = P_s0[ps0_i][t];  */
      P_s_p = P_s0[sp_i][t];  

      /*
       * Iterate over number of crowns penetrated by a beam given the
       * pathlength s'.
       */
      for (n = 1; n <= MAXCROWNS; n++) {

	 /*
	  * Get the probability that 'n' crowns are penetrated, given the
 	  * current after-entering-crown pathlength s_p[i].
	  */
         temp1 = Vol(h, sp_i, t, g.h2_p) - Vol(h, sp_i, t, g.h1_p);

	 temp1 *= g.Lv_p;  
/*         if(n==1) 
 printf("\n tmp1=%8.3f,h=%6.2f sp_i=%5.1f theta=%5.1f Vol2=%8.2f Vol1=%8.2f h2=%5.1f h1=%5.1f ", 
 printf("\n %8.3f %6.2f %5.1f %5.1f %8.2f %8.2f %5.1f %5.1f ", 
                 temp1, height_p[h], height_p[sp_i], RAD_TO_DEG(theta_p[t]),
                  Vol(h, sp_i, t, g.h2_p), Vol(h, sp_i, t, g.h1_p), g.h2_p, g.h1_p) ;  */
	 P_n = (pow(temp1, (fp_t) n) * exp(-temp1)) /
	    (factorial[n] * (1.0 - exp(-temp1)));

	 /*
	  * Get the mean within-crown pathlength.
	  */

	 s = s_p[sp_i][t] * (1.0 - exp(-1.0 * (fp_t) n *
				       ES[h] / s_p[sp_i][t]));

         /* convert back to original elliptical dimension 
         s *= g.ellipticity * cos(theta_p[t])/cos(theta[t]); */

	 /*
	  * Accumulate the appropriate value.  We use the calculated path
	  * length 's' as an index into the array that we are filling, and
	  * increment its value by the product of the probability that 'n'
	  * crowns are penetrated and the probability of finding this
	  * particular after-entering-the-crown pathlength.
	  */
 
         n_mean   += (float) n * P_n * P_s_p;  
	 PD_s[h][t][s_to_index(s)] += P_n * P_s_p;  
/*         
         if(h==0 && (n==1 || n==2 || n==3 || n==4 || n==5) )
             printf("\n%4.1f %3d %5.2f %6.3f %5.2f  %5.2f %3d ",
             height[sp_i], n, P_n, P_s_p, s_p[sp_i][t], s, s_to_index(s));
 */
#ifdef NOT
	 printf("%4.1f  %6.2f  n: %2d  s_i: %2d  s': %4.2f  s: %4.2f  i: %2d  prob: %8.6f\n",
		height[h], RAD_TO_DEG(theta_p[t]), n,
	   sp_i, s_p[sp_i][t], s, s_to_index(s), PD_s[h][t][s_to_index(s)]);
#endif
      }
   }

/*      if ( t==0 || t == 12 )   
         if( h == 0) { 
               printf("\n %4.1f, n_mean=%4.1f", RAD_TO_DEG(theta[t]), n_mean ); 
               printf("\n PD_s");   
               maxs = s_to_index((g.z2_p - height_p[h]) / cos(theta_p[t]));
               for (i = 0; i <= maxs; i++) {
                  printf("\n %3d %4.2f %10.7f", i, index_to_s(i),PD_s[h][t][i] ); 
               }
	  }
*/
}           

/*
static void
init()
{
   int             k;
    * Fill our factorial array.
   factorial[0] = 1;
   for (k = 1; k <= MAXCROWNS; k++) {
      factorial[k] = factorial[k - 1] * (fp_t) k;
   }
}
*/
int
s_to_index( s)
fp_t s;
{
   return ((int) (s / g.ds + 0.5));
}


fp_t
index_to_s(index)
int index;
{
/*   return ((fp_t) index * g.ds);   */
   return ((fp_t) index * g.ds);
}

/* Get the overlay or shadow  volume for P_n, which is the volume intersected by 
 * h1 and h2 
 */
static fp_t Vol(h, h_s, t, h_b)
int h, h_s, t;
fp_t h_b;
{
   fp_t   V, V_0, V_sp1, V_sp2, V_cyln;
   fp_t   tmp_s;
   fp_t   h_t, h_tt;

   tmp_s = (height_p[h_s] - height_p[h]) / cos(theta_p[t]) ;
   V_0 =  PI * g.R_squared * tmp_s;
   V_0 += (4.0/3.0) * PI * g.R_cubed;    
 
      if     ( (height_p[h] - g.R) >= h_b ) {
          V = 0.0; 
/*          printf(" V1 = %5.2f", V);  */
     }

      else if( (height_p[h] - g.R*sin(theta_p[t])) >= h_b ) {
         h_t = g.R - ( height_p[h] - h_b );
         V = (PI/3.0) * h_t * h_t * (3.0*g.R - h_t);
/*         printf(" V2 = %5.2f", V);  */
      }

      else if( (height_p[h] + g.R*sin(theta_p[t])) >= h_b ) {
         V_sp1 = (2.0/3.0)*PI*g.R_cubed ;
         V_sp1 -= trisec(height_p[h], h_b, theta_p[t], g.R);
         h_tt = (h_b-(height_p[h]-g.R*sin(theta_p[t])))/cos(theta_p[t]);
/*        if( h_tt <= (height_p[h_s]-height[h])/cos(theta_p[t]) ) { */
         if( height_p[h_s]-g.R*sin(theta_p[t]) >= h_b ) { 
            fp_t hh1, hh2, hh;
            hh1 = (height_p[h]-h_b)/sin(theta_p[t]);
            hh2 = g.R;
            hh = h_tt;
            V_cyln = cylind(g.R, hh1, hh2, hh);
            V_sp2 = 0.0;
	  }
         else {
            fp_t hh1, hh2, hh;
            V_sp2 = 0.0;
            hh1 = (height_p[h]-h_b)/sin(theta_p[t]);
            hh2 =(height_p[h_s]-h_b)/sin(theta_p[t]);
            hh = (height_p[h_s]-height_p[h])/cos(theta_p[t]) ;
            V_cyln = cylind(g.R, hh1, hh2, hh);
            V_sp2 = trisec( h_b,height_p[h_s], theta_p[t], g.R);
	 }
         V = V_sp1 + V_cyln + V_sp2;
/*      printf(" V3 = %5.2f V_sp1=%5.2f V_cyln=%5.2f V_sp2=%5.2f", V,  V_sp1, V_cyln, V_sp2); */
      }  
      else if( height_p[h_s] - g.R*sin(theta_p[t]) >= h_b ) {
         fp_t tmp_h;
         tmp_h = (h_b-height_p[h])/cos(theta_p[t]);
         V_cyln = PI*g.R*g.R*tmp_h;       
         V_sp1 = (2.0/3.0)*PI*g.R_cubed ;  
         V = V_sp1 + V_cyln;
/*         printf(" V4 = %5.2f V_sp1=%5.2f V_cyln=%5.2f", V,  V_sp1, V_cyln);  */
      }
      else if( height_p[h_s] + g.R*sin(theta_p[t]) >= h_b ) {
         fp_t hh1, hh2, hh;
         fp_t tmp_h;
         h_tt = (height_p[h_s]+g.R*sin(theta_p[t]) - h_b)/cos(theta_p[t]);  
         hh1 = (h_b-height_p[h_s])/sin(theta_p[t]); 
         hh2 = g.R;
         hh = h_tt;
         tmp_h = (height_p[h_s]-height_p[h])/cos(theta_p[t]);
         V_cyln = PI*g.R*g.R*tmp_h - cylind(g.R, hh1, hh2, hh);
         V_sp2 = trisec( h_b, height_p[h_s], theta_p[t], g.R);
	 V_sp1 = (2.0/3.0)*PI*g.R_cubed ;  
         V = V_cyln + V_sp2 + V_sp1;
/*         printf(" V5 = %5.2f V_sp1=%5.2f V_cyln=%5.2f V_sp2=%5.2f", V, V_sp1, V_cyln, V_sp2);  */
      }

      else if( height_p[h_s] + g.R  >= h_b ){
         h_t = g.R - (h_b - height_p[h_s]);
         V_sp1 = (PI/3.0) * h_t * h_t * ( 3.0*g.R - h_t);
         V = V_0 - V_sp1;
/*         printf(" V6 = %5.2f", V);  */
      }
      else {      /* height[h_s] + g.R < g.h2_p */ 
         V  = V_0;
/*         printf(" V7 = %5.2f", V);   */
      }
  return (V);
}


static fp_t trisec(hh, hh_b, th, r)
fp_t hh, hh_b, th, r;
{
    fp_t x, h_0, tmp;
    fp_t b, a1, a2;
/*    int noint = 20; */
    int noint = 20; 

    tmp = (hh - hh_b);

/*    printf("\n trisec: %6.2f", r*r-tmp*tmp); */
    h_0 =  1.0 * tmp*cos(th) + sqrt(r*r- tmp*tmp)*sin(th);
    x   = -1.0 * tmp*sin(th) + sqrt(r*r- tmp*tmp)*cos(th);
 
    b = - tmp/sin(th);    /* input parm for triangle volume calculation  */
/*    b = x - h_0/tan(th);   input parm for triangle volume calculation  */
    a1 = x ;              /* input parm for sector   volume calculation  */
    a2 = r;               /* input parm for sector   volume calculation  */

    return(triang(b,r,th,noint)+sector(a1,a2,r));

}

 
/*     a1  and  a2  lie in the range (-r,r)  */
static fp_t sector(a1,a2,r)
fp_t a1,a2,r;
{
 fp_t b1, b2,volume;

      b1=r*r*a1-(a1*a1*a1)/3.0;
      b2=r*r*a2-(a2*a2*a2)/3.0;
      volume=PI*(b2-b1)/2.0;
 
      return(volume);
}
 
/*    b  lies in the range (-r,r)
      theta  in (0,pi/2)
 */
static fp_t triang(b,r,the,noint)
fp_t b,r,the;
int noint;
{
fp_t sint, cost, a1, x0, h, volume;
int m, i;
fp_t tmpb;
fp_t sum1, sum2;

      sint=sin(the);
      cost=cos(the);
/*
      tmpb = (b*b*sint*sint);
      a1=tmpb*tmpb +(r*r*cost*cost) - tmpb;
      x0=b*(sint*sint) + sqrt(a1);
 */
     
      a1 = r*r - b*b*sint*sint;
/*      printf("\n triangle:  %5.2f", a1); */
      x0=b*(sint*sint) + sqrt(a1)*cost;
 
/*     composite Simpson's rule */
      m=noint;
      h=.50*(x0-b)/(float) m;
 
      sum1=0.0;
      for(i=0; i<m; i++) {
        sum1+=fcn(b+(float)(2*i+1)*h, b, r,the);
      } 
      volume=4.0*sum1;
 
      sum2=0.0;
      for(i=0; i<m-1; i++){
         sum2 += fcn(b+(float)(2*(i+1))*h, b, r, the);
      }       
      volume += 2.0*sum2;

/*    note  that   fcn(x0,b,r,theta)=0  by definition of x0
             and   fcn(b,b,r,theta)=(r**2-x**2)*pi/2
      volume += (r*r - b*b)*PI*.50;
 */

      volume +=  fcn(x0,b,r,the);
      volume +=  fcn(b,b,r,the);
      volume *=h/3.0;
 
/*     end of composite Simpson's rule */
      return(volume);
}  

static fp_t fcn(x,b,r,the)  
fp_t x,b,r,the;
{
fp_t a1,a2,a3,func;
 
      a1=tan(the)*(x-b);
      a2=r*r-x*x ;
      a3=a2-a1*a1;
      if(fabs(a3)<0.0000000001) a3 = 0.0;
/*      printf("\n fcn: %5.2f", a3); */
      func = 2.0* a1*sqrt(a3);  
/*      func = -a1*sqrt(a3)+a2*asin(sqrt(a3/a2)) ;  */

      return(func);
}
 
 
/*     integrand of  sqrt(r**2-x**2) */

static fp_t  fun1(x,r)
fp_t x,r;
{ fp_t func;  
/*      printf("\n fun1: %5.2f", (r*r-x*x) ); */
      func =.50*x*sqrt(r*r-x*x)+.50*r*r*asin(x/r); 
      return(func);
} 
 
/*
      plane cuts both ends of the cylinder
      h1        is the bottom x-intercept, lies in (-r,r)
      h2        is the bottom x-intercept, lies in (-r,r)
      h2 > h1
      height    is the height of the cylinder
 */

static fp_t cylind(r,h1,h2,h)
fp_t r,h1,h2,h;
{ fp_t slope, volume;
  fp_t tmp1, tmp2;

      slope=h/(h2-h1);
/*      printf("\n cylinder:(r*r-h1*h1): %5.2f",(r*r-h1*h1) ); 
      printf("\n cylinder:(r*r-h2*h2): %5.2f",(r*r-h2*h2) );   */
      tmp1 = sqrt(r*r-h1*h1); 
      tmp2 = sqrt(r*r-h2*h2); 

      volume= tmp1*tmp1*tmp1 - tmp2*tmp2*tmp2; 

      volume /= 3.0;
      volume -= h1*(fun1(h2,r)-fun1(h1,r));
      volume *= 2.0*slope; 

      if(h2 < r) {
        fp_t phi, s1,s2;
        phi = acos(h2/r); 
        s1 = r*r*phi;
        s2 = r*sin(phi)*h2;
        volume += (s1-s2)*h;
      }
      return(volume);
}

