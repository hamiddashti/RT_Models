/*  A hybrid GORT-BRDF model for discontinuous plant canopies
 *
 *  Input parameters are:
 *  h2:      upper boundary of crown centers (m)
 *  h1:      lower boundary of crown center (m) 
 *  radius:  horizontal mean crown radius (m)
 *  brratio: crown spheroid ellipticity (b/r)
 *  lambda:  trees per unit area             
 *  omega:   single scattering albedo of leaf
 *  rg:      background albedo      
 *  ELAI: effective leaf area index
 *
 *  This model also needs the gap probabilities and
 *  openness factors (Ni et al., 1997) as inputs. The gap probabilities 
 *  are functions of  tree geometry,  tree density (h2,h1,radius,brratio, 
 *  lambda) and solar and viewing geometry. The openness factors 
 *  are the inregrations of gap probabilities over the hermishere.
 *  They are calculated by pgap function. 
 * 
 *  
 *  Anytime feel free to ask if you have any questions. 
 *
 *  Wenge Ni-Meister 
 *  212-772-5321
 *  email: wenge.ni-meister@hunter.cuny.edu
 * 
 *  References:
 *  
 *  Ni, W,  X. Li, C.E. Woodcock, M.R. Caetano and A. Strahler, 1999.
 *  An Analytical Hybrid GORT Model for Bidirectional Reflectance
 *  over Discontinuous Plant Canopies.  
 *  IEEE Transactions on Geoscience and Remote Sensing, 
 *  vol. 37(2):987-999
 *
 *  Ni, W,  X., Li, C.E. Woodcock, J. L. Roujean, R. Davis, 
 *  Transmission  of Solar Radiation In Boreal Conifer Forests:
 *  Measurements and Models, Journal of Geophysical Research, 1997, 
 *  102(D24):29555-29566. 
 *
 * Li, X. and A. Strahler, 1992. Geometric-Optical Bidirectional Reflectance
 * Modeling of the Discrete Crown Vegetation Canopy: Effect of Crown Shape 
 * and Mutual Shadowing. IEEE Trans. Geosci. Rem. Sen. 30(2):276-292.
 *
 * Schaaf, C.B., X. Li, A.H. Strahler, 1994. Topographic Effects on Bidirectional
 * and Hemispherical Reflectance Calculated with a Geometric-Optical Canopy Model.
 * IEEE Trans. Geosci. Rem. Sen. 32(6):11861193.
 * 
 */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "gort1.h"

#define NPT 45
#define npn 90

fp_t g_k_open1, g_k_open2;         

double g_pn[npn][3];
double g_h1, g_h2, gg_b, g_radius, g_lambda, g_brratio;
double g_rg, g_omega, g_rl, g_tl, g_ELAI, g_FAVD;
double g_SAZ,g_SZN,g_VAZ,g_VZN,g_SLOPE, g_ASPECT;
double func4();
double g_tsp,g_phsp,g_tip,g_phip,g_sints,g_costs;

GLOBAL_T g;

void gort_brdf(double g_egap[][3], double *reValue)
{
 
  int i, m, k, j;
  double fret;    /* return value from calulating gort brdf */
  double gg_pn[g.nth][3]; /* 2-d array [zenith_angle, within_crown_pgap, between_crown_pgap] */ 
 
 /* assign the g structure to the variables */ 
  g_h1 = g.h1;
  g_h2 = g.h2;
  g_radius = g.R;
  g_brratio = g.ellipticity;
  g_lambda = g.density;
  g_rl = g.rl;
  g_tl = g.tl;
  g_rg = g.rg;
  g_FAVD = g.FAVD;
  g_SZN = g.szn;
  g_VZN = g.vzn;
  g_SAZ = g.saz;
  g_VAZ = g.saz;
  g_SLOPE = g.slope;
  g_ASPECT = g.aspect;

  /* read the gap probabilities data calculated from directory /gort */ 
  g_k_open1 = g.k_open;
  g_k_open2 = g.k_openEP;


 /* read the gap probabilities array calculated from directory /gort */ 
  for (i=0; i< g.nth - 1; i++)
    for (j=0; j< 3; j++)
     g_pn[i][j] = g_egap[i][j];

  
  /* call func4() to calculate brdf  */
     fret=func4();
     *reValue = fret;
     g.brdf = fret;
}

/* External variable for slop angle g_tsp, g_phsp, input direction and
output direction g_tip, g_phip. All in radians */

slope()
{          
  float sinti,costi,cosphi,sinphi;
  float x,y,z,r;
  
   if(g_SLOPE < 5.0) return 0;
  
  g_phip -=  g_phsp;
  cosphi = cos(g_phip);
  sinphi = sin(g_phip);
  
  sinti = sin(g_tip);
  costi = cos(g_tip);
  
  x = cosphi*sinti*g_costs - costi*g_sints;
  y = sinti*sinphi;
  z = sinti*cosphi*g_sints + g_costs*costi;
  r = x*x + y*y;
  g_tip = acos(z);
  g_phip = asin(y/sqrt(r));
  if(x < 0 ) {
    if( y > 0) g_phip = PI - g_phip;
    else g_phip = -PI - g_phip;
  }
  
  return 0;
}

double func4()
{
  float sumsq,x,y,sinpv,sinpi,cospi,cospv,cost,temp;
  float g_h, g_b, beta, numtree, r2, hbratio, pixsiz;
  float sunzn, rc, rz, rt, sunaz,secti,sinti,costi,tanti;
  float vaz, vzn, mv, mi, stdh, htdif, D, Ac, projAci;
  float overpp,cosphi,sinphi,overpp0,overpp180,iv0;
  float iv180, f0, f180, kg0, kg180, kc0;
  float kc180, kz0, kz180;
  float dist0, dist180, absdist0, absdist180, t0, t180, cost0, cost180;
  float bigF,bigF0, bigF180, m0, m180, thetai, thetav, fhat0;
  float fhat180, gammac180, pm0,pm180, L, L0, L180;
  float pvmv0, pvmv180, pimi0, pimi180, gamma0, gamma180, gammac0;
  float sintv, costv, sectv, tantv, sumax, thetami, thetamv;
  float gammax, gammai, gammav, ptshotspot;
  float over, Sa, rho, xcenter, secdiff, phiprime, tanhalf, kg;
  float kc, kz, kt, f, iv,  maxz,minz, dthetadphi, albedo, lambt;
  double rc0, rc1, rz0, rz1, rt0, rt1, tmp1;
  double lambda0, vazpi;
 
  /*for C, G,Z */
  float  H;
  double th, th_p, mu_s, sin_s, rmu_s, gama, V_g, P_n0; 
  double gtau_d, gtau_f;
  double t_o, R_ff, R_df, T_ff, T_df;
  double t_df, t_ff, rho_df, rho_ff;
  
  /* for K_open */  
  int     t, nth;
  double  dth, the, the_p; 
  double f_d, f_f, t_ff_p, t_df_p, rho_df_p, rho_ff_p, t_o_p;  
  int ki, kv;
  double tau, si, sv, cosa, siv, piv, p, prop, mu_v;
  double A, B, C, DD, E, F, HH;
  double g;  

  g_h = (g_h1+g_h2)/2.0;
  g_b = g_brratio*g_radius; 
  hbratio = g_h/g_b ;
  
  htdif = (g_h2-g_h1)/g_b; 
  H = htdif*g_b;
  /*
  H = htdif*b + 2.0*b;
  */
  sunzn = DEG_TO_RAD(g_SZN);
  g_phip =  DEG_TO_RAD(g_SAZ);
  tanti = tan(sunzn);

  tanti *= g_brratio;
  if(tanti < 0) tanti = 0.;
  g_tip = atan(tanti);

  /*  Slope  and aspect are considered by Conghe
  slope effect, but here the slope is set to zero  
  g_sints = 0. ;
  g_costs = 1. ; 
  */ 

  /* slope effect considered here by Conghe Song */
  g_tsp=DEG_TO_RAD(g_SLOPE);
  g_tsp=PI/2.0-atan(g_brratio*tan(PI/2.0-g_tsp));
  g_phsp=DEG_TO_RAD(g_ASPECT);
  g_sints=sin(g_tsp);
  g_costs=cos(g_tsp);
  slope();  
 /* Conghe
  vaz = g_phip;
 */
  sunaz=g_phip;
  thetai = g_tip;
  
  sinti = sin(thetai);
  costi = cos(thetai);
  secti = 1/costi;
  
  
  /* This is for hot-spot quantit, C and beta*/
  gammai = PI* secti;
  gamma0 = gammai;
  
  lambda0 = g_lambda*g_radius*g_radius;
  kg0 = exp( -lambda0*gammai );
  
  /* Calculate the beta - Based on Li and Strahler IGARSS92 paper and 
     formulation by Abdelgadir Abuelgasim (RSE 93 paper)*/
  
  D= g_brratio*htdif*sin(thetai*.5) / cos(thetai*.5);
  temp = lambda0*gammai ;
  if(temp < 0.00001) beta = 1.0;
  else
    beta=(lambda0*gammai/(lambda0*gammai+D))*
      ((1.0-exp(-(lambda0*gammai+D)))/(1.0-exp(-lambda0*gammai)));
  
  /* This is M - eqn 14*/
  temp = lambda0*gamma0 ;
  if(temp < 0.00001) m0 = 0.0;
  else m0 = 1.0-((1.0-kg0)/ temp);
  mi = m0;
  
  /* This is theta sub Mi - eqn 18*/ 
  if(mi < 0.) mi = 0.;
  thetami = 2.0*asin(sqrt(mi));
  /*
    if(thetami < PI/2.0) thetami = atan(tan(thetami)*g_radius/b);
    else thetami = PI - atan(tan(PI-thetami)*g_radius/b);
    */
  

  /* This is the Ac -- the physical surface area of an illuminated crown*/
  vzn =  DEG_TO_RAD(g_VZN);
  tantv = tan(vzn);
  tantv *= g_brratio;
  if(tantv < 0) tantv = 0.;
  g_tip=vzn;
  g_phip =  DEG_TO_RAD(g_VAZ) ;
  
  slope();
  
  thetav = g_tip;
  vaz = g_phip - sunaz;
  /*
    if(vaz > 0.) vaz *= -1;
    */
  if(vaz < 0.) vaz *= -1;
  cosphi = cos(vaz);
  sinphi = sin(vaz);
  
  sintv = sin(thetav);
  costv = cos(thetav);
  sectv = 1/costv;
  /* if(Iter == 100) printf("%f %f %f\n ",vaz,cosphi,sinphi); */
  x = tanti*tanti + tantv*tantv - 2*tanti*tantv*cosphi;
  if (x > 0.00001) {
    x = sqrt(x);
    /* x: distance between centers of two projection */
    sinpv = tanti*sinphi/x;
    cospv = (tanti* cosphi - tantv ) / x;
    sinpi = tantv*sinphi/x;
    cospi = (tantv* cosphi - tanti ) / x;
  }
  else {
    x = 0.;
    sinpv = 0;
    cospv = 1.;
    sinpi = 0;
    cospi = 1.;
  }
  temp =sinpv*sinpv + cospv*cospv*sectv*sectv;
  if (temp < 0.) temp = 0.;
  y = sqrt(temp);
  temp = sinpi*sinpi + cospi*cospi*secti*secti;
  if (temp < 0.) temp = 0.;
  y += sqrt(temp);
  /* sum of two half axis */
  cost = hbratio*x/y;
  if(cost < 0.) cost = 0.;
  if (cost >= 1.0) temp = 0.0;
  else temp = acos(cost);
  over = (temp - sin(temp)*cos(temp)) * (secti+sectv)/PI;
  /*  accurate PP, very good approximation on PC, reasonably
      good off PP/PC unless b/r very large AND h/b very small */
  sumax = secti+ sectv;
  
  /*Calculate iv explicitly --eqn 9 -- the illuminated portion 
    of the area of a single spheriod */
  iv=costi*costv +sinti*sintv*cosphi;
  
  /*Calculate F (with gamma and gammac) explicitly--eqn 12*/
  
  /* This is gamma - pg 14, at given sunzn and vzn, vaz*/
  gammax = PI*(sumax-over);
  
  /*Use this gammax in the kg calculation*/
  kg=exp(-lambda0*gammax);
  
  /* Calculate f0, bigF0 on PP at give sunzn and vzn*/
  
  dist0 = tanti - tantv;
  
  /*Absdist is the second part of eqn 4 -- (cos phi = +1 along PP)*/
  
  absdist0 = dist0*hbratio;
  if (absdist0 < 0.0) {
    absdist0 *= -1.0;
  }
  
  /* eqn 6*/
  cost0 = absdist0 /sumax;
  if(cost0 < 0.) cost0 = 0.;
  if (cost0 >= 1.0) t0 =0.0;
  else  t0 = acos(cost0);
  /* eqn 5 -- the exact overlap for the PP*/
  overpp0 = (t0-sin(2.0*t0)/2.0)*sumax;
  
  /* Calculate f180, bigF180 */
       
  dist180 = tanti + tantv;
  
  /*Absdist is the second part of eqn 4 -- (cos phi = -1 along PP)*/
  
  absdist180 = dist180*hbratio;
  
  /* eqn 6*/
  cost180 = absdist180 /sumax;
  if(cost180 < 0.) cost180 = 0.;
  if (cost180 >= 1.0) t180 = 0.0;
  else t180 = acos(cost180);
  
  /* eqn 5 -- the exact overlap for the PP*/
  overpp180 = (t180-sin(2.0*t180)/2.0)*sumax;
  
  /*Calculate the mutual shadowing effect along the PP*/
  
  /* This is gamma - pg 14*/
  gamma0 = PI*sumax- overpp0;
  gamma180 = PI*sumax- overpp180;
  
  /* the kg computation -- eqn 3 (the overlap should be Over divided by PI)*/
  kg0 = exp( -lambda0*gamma0);
  kg180 = exp( -lambda0*gamma180);
  
  /* the illuminated portion of the area of a single spheriod - eqn 9*/
  iv0 =  costi*costv + sinti*sintv;
  iv180 =  costi*costv - sinti*sintv;
  /* This is gamma sub v - pg 14*/
  gammav = PI* sectv;
  /* This is gamma sub i - pg 14*/
  gammai = PI* secti;
  /* This is gamma sub c - top of eqn 12*/
  gammac0 = 0.5 * (1+iv0)*gammav;
  gammac180 = 0.5 * (1+iv180)*gammav;
  
  /* F - eqn 12*/
  bigF0 = gammac0/gamma0;
  bigF180 = gammac180/gamma180;
  
  
  /* This is M - eqn 14*/
  temp = lambda0*gamma0 ;
  if(temp < 0.00001) m0 = 0.0;
  else
    m0 = 1.0-((1.0-kg0)/(temp));
  temp = lambda0*gamma180 ;
  if(temp < 0.00001) m180 = 0.0;
  else
    m180 = 1.0-((1.0-kg180)/(temp));
  
  /* This is Mv - eqn 11*/
  temp = lambda0*gammav ;
  if(temp < 0.00001) mv = 0.0;
  else
    mv = 1.0 - ((1.0 - exp(-temp))/(temp));
  
  /* This is Mi - eqn 10*/
  temp = lambda0*gammai ;
  if(temp < 0.00001) mi = 0.0;
  else
    mi = 1.0 - ((1.0 - exp(-temp))/(temp));
  
  /* This is theta sub Mi - eqn 18*/ 
  if(mi < 0.) mi = 0.;
  thetami = 2.0*asin(sqrt(mi));
  
  /* This is theta sub Mv - eqn 18*/
  if(mv < 0.) mv = 0.;
  thetamv = 2.0*asin(sqrt(mv));
  
  /*pvmv and pimi are eqns 4 and 5 in the IGARSS 1992 paper*/
  
  pvmv0=mv-(1-cos(thetav-thetai))/2.0;
  pvmv180=mv-(1-cos(-thetav-thetai))/2.0;
  pimi0=(1-cos(thetami*(1 - (thetai - thetav)/PI)))/2.0;
  pimi180=(1-cos(thetami*(1 - (thetai + thetav)/PI)))/2.0;
  pm0 = pvmv0;
  if(pimi0 > pvmv0) pm0 = pimi0;
  pm180 = pvmv180;
  if(pimi180 > pm180) pm180 = pimi180;
  
  /* This f is from eqn 17 and the first part of eqn 20*/
  if (thetav > thetai) { 
    f0= (1-(gammav*pvmv0/gammac0))/(1-m0);
  }
  else  {
    f0= (1-(gammav*pm0/gammac0))/(1-m0);
         }
  
  f180= (1-(gammav*pm180/gammac180))/(1-m180);  
  fhat0=beta*f0*bigF0 + (1.0-beta)*bigF0;      
  fhat180=beta*f180*bigF180 + (1.0-beta)*bigF180;
  
  /* This is gamma sub c - top of eqn 12*/
  gammac0 = 0.5 * (1+iv)*gammav;
  
  bigF = gammac0/gammax;
  /* The following is for preparing interpretation of Mutual shadowing */
  
  /*Linearly interpolate f (the non F portion)*/
  vazpi = vaz/PI;
  if(vazpi > 1.0) vazpi = 2.0-vazpi;
  f = (1.0-vazpi)*f0*bigF0+ (vazpi)*f180*bigF180;
  fhat0 = beta*f + (1.0-beta)*bigF;
  kc = fhat0*(1.0-kg);
  kz =exp(-lambda0*gammav) - kg;
  kt = 1.0 - kg - kc - kz;

  kt = MAX(0.0,kt);
  /* A: areal proportion of sunlit background (including viewed + non-viewed)
   * B: areal proportion of shaded background (including viewed + non-viewed) (1-A) 
   * C: areal proportion of viewed background (including sunlit + shaded ) 
   * kz: areal proportion of viewed and shaded background 
   * DD: areal proportion of non-viewed background (including sunlit + shaded ) (1-C)
   * kg: areal proportion of sunlit and viewed background 
   * E: areal proportion of sunlit but non-viewed background (A-kg) 
   * HH: areal proportion of non-viewed  and shaded background (DD-E=B-kz)
   */
  A =  exp(-lambda0*gammai);
  B = 1 - A;
  C =  exp(-lambda0*gammav);
  DD = 1 - C;
  E = A - kg;  
  HH = D - E;

  /* calculation of the component spectral signature */
  /* The following two lines are no necessary. Angle transformation
     has been done earlier Conghe
  th =  DEG_TO_RAD(g_SZN);
  th_p = atan(tan(th)*g_brratio); 
  */
  g_omega = g_rl + g_tl;
  
  /* Conghe
  mu_s  = cos(th);
  */
  mu_s = costi; /*Conghe */
  mu_v  = costv;
  sin_s = sqrt ( 1.0 - mu_s * mu_s ); 
  rmu_s = 1.0 / mu_s ;
  gama = sqrt( 1.0 - g_omega );
  
  /* by Conghe  
  ki = (int) g_SZN;
  kv = (int) g_VZN;
  */
  /* The above two lines are changed to the following by Conghe */
  ki = (int) RAD_TO_DEG(thetai);
  kv = (int) RAD_TO_DEG(thetav);

  tau = 0.5* g_FAVD;
  g_ELAI = g_FAVD * (1.333333*lambda0*PI*g_b);

  if(g_pn[ki][2]>0)
    si = -log(g_pn[ki][2])/tau ;
  else si = 1000;
  if(g_pn[kv][2]>0)
    sv = -log(g_pn[kv][2])/tau ;
  else sv = 1000;

  cosa = costi*costv+sinti*sintv*cosphi;
  siv = sqrt(si*si+sv*sv-2*si*sv*cosa);
  if(siv > 0.0) 
    piv = exp( tau*sqrt(si*sv)*(1.0-exp(-siv/(1.0*g_radius)))/(siv/(1.0*g_radius)) );
  else piv = exp( tau*sqrt(si*sv) );
  p = (g_pn[ki][2])*(g_pn[kv][2])*piv;

  gtau_d = 0.5 * g_ELAI; 
  gtau_d /= H ;
  gtau_f = gtau_d;

  t_o = exp(- H * gtau_d / mu_s );   
  R_ff = ( 1.0 - gama ) / ( 1.0 + gama );
  R_df = ( 1.0 - gama ) / ( 1.0 + 2.0 * mu_s * gama );    
  T_ff = exp( - 2.0 * gama * gtau_f * H ) ;
  T_df = g_omega/2.0;
  T_df *= (T_ff - t_o);
  T_df *= (1.0 + 2.0*mu_s)/(1.0-(2.0*gama*mu_s)*(2.0*gama*mu_s));
  
  t_ff = ( 1.0 - R_ff * R_ff ) / ( 1.0 - T_ff*T_ff*R_ff*R_ff );
  t_ff *= T_ff ;
  
  rho_ff = (1.0 -  T_ff*T_ff ) / ( 1.0 - T_ff*T_ff*R_ff*R_ff );
  rho_ff *= R_ff ; 
  t_df = T_df - rho_ff*(t_o*R_df+T_df*R_ff);
  /* Conghe
  rho_df = R_df - rho_ff*(t_o*R_df+T_df*R_ff);
  */
  rho_df = R_df - t_ff*(t_o*R_df+T_df*R_ff);
  /* the souce energy is 1.0 at zenith and change as a consine function 
  f_d = 1.0 - mu_s/(0.09 + mu_s);
  */
  f_d = mu_s;
  f_d *=mu_s/(mu_s +0.09);
  f_f =mu_s - f_d;
  t_ff_p = t_ff*(1.0 - g_k_open1- g_k_open2) + g_k_open1 + g_k_open2;
  t_df_p = (1.0 - g_pn[ki][1]- g_pn[ki][2]) * t_df;
  rho_ff_p = rho_ff;
  rho_df_p = rho_df;
  t_o_p =  g_pn[ki][1]+ g_pn[ki][2] ;
  rz0 = (t_df_p+g_pn[ki][2])* g_rg;
  rz1 = (t_ff_p-g_k_open1) * g_rg;
  /* Conghe
  rz = (1.0-f_d)*rz0 + f_d*rz1; 
  */
  rz=f_d*rz0+f_d*rz1;  
  rc0 = rz*(DD-E) + g_rg*E ;

  /* rc0 /= DD; -conghe */
  rc0 *= g_k_open2;

  /* prop = 0.25*g_omega / R_ff; -conghe */
  g = -(4.0/9.0)*(g_rl-g_tl)/g_omega;
  rc0 += rho_df_p + (1-g_omega)*g_omega*p*(1-g)/(2.0*costi*costv);
  rc0 += (t_o_p+t_df_p) * g_rg * (t_ff_p - g_k_open1) /(1.0-g_rg*rho_ff_p);

  rc1 = rho_ff_p;
  rc1 += (g_k_open2+g_k_open1)*g_rg*g_k_open2 +(1-g_k_open1-g_k_open2)*rz*g_k_open2;
  rc1 +=  t_ff_p*g_rg*(t_ff_p - g_k_open1)/(1.0-g_rg*rho_ff_p);

 /* Conghe
  rc = (1.0-f_d)*rc0 + f_d*rc1;   
  */

  rc = f_d*rc0 + f_f*rc1;

  /* printf("gort_brdf g, rc0, rc1, rc  value:%lf %lf %lf  %lf\n",g, rc0, rc1, rc); */

  rt0 = (t_df_p + g_pn[ki][1] + g_pn[ki][2])* g_rg *(t_ff_p- g_k_open1) /(1.0-g_rg*rho_ff_p);
  rt1 = rho_ff_p;
  rt1 += t_ff_p * g_rg * (t_ff_p - g_k_open1) /(1.0-g_rg*rho_ff_p); 
  /*  rt0 = t_df_p + (t_df_p + g_pn[ki][2])* g_rg *(t_ff_p- g_k_open1) /(1.0-g_rg*rho_ff_p);
  rt1 = t_ff*(1.0 - g_k_open1- g_k_open2)+ (t_ff_p-g_k_open1) * g_rg * (t_ff_p - g_k_open1) /(1.0-g_rg*rho_ff_p); */
  /* Conghe
  rt = (1.0-f_d)*rt0 + f_d*rt1;   
  */
  rt = f_d*rt0 + f_f*rt1;
  
  /* printf("gort_brdf f_d, rt0, rt1, rt  value:%lf %lf %lf  %lf\n",f_d, rt0, rt1, rt); */
  /*
  if(g_SAZ>120 || g_SAZ < 0 ) {
    printf("%5.1f  %5.1f  %5.1f  %5.3f  %5.3f %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %6.3f\n",
   g_SZN,-g_VZN, g_SAZ, kc, rc, kg, g_rg, kt, rt, kz, rz, kg*g_rg + kc*rc + kz*rz + kt*rt); 
  }
  else {
    printf("%5.1f  %5.1f  %5.1f  %5.3f  %5.3f  %5.3f %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %6.3f\n",
	   g_SZN, g_VZN, g_SAZ, kc, rc, kg, g_rg, kt, rt, kz, rz, kg*g_rg + kc*rc + kz*rz + kt*rt); 
  }
  */

  return kg*g_rg + kc*rc + kz*rz + kt*rt;

}

