#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gort1.h"

#define SIZE 80


GLOBAL_T g;

int
main()
{
   FILE    *pfp;		/* parameter file pointer */
   char gortFile[] = "gort_inputs";   /* parameter for gort1,gort_brdf */
 
   int i,j;
   double brdfValue = 0.0;


   /* declare the fp_t type variables   */
  fp_t h1; 
  fp_t h2; 
  fp_t R ;    
  fp_t ellipticity; 
  fp_t density;  
  fp_t ELAI;   
  fp_t dth;
  fp_t dz;
  fp_t ds; 
  fp_t k; 
  fp_t lambda;
 
  double dth_num;  
  fp_t rl; 
  fp_t tl; 
  fp_t rg;  

  fp_t szn ; 
  fp_t vzn; 
  fp_t saz;
  fp_t vaz; 
  fp_t slope; 
  fp_t aspect; 

   /* read gort parameter file (gort_parm)to get the input value for pgap function   */
   pfp =fopen(gortFile,"r");
   if(pfp==NULL){printf("can not open file input parameter file\n");
    exit(1);}

  /* read input parameters */
  char   p_file[30];
  char   dummy[200];
  
  fscanf (pfp, "%lf", &h1);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &h2);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &R);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &ellipticity);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &density);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &ELAI);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &dz);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &ds);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &dth_num);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &k);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &rl);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &tl);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &rg);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &szn);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &vzn);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &saz);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &vaz);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &slope);
  fgets (dummy, 200, pfp);
  fscanf (pfp, "%lf", &aspect);
  fgets (dummy, 200, pfp);
  fclose (pfp);

  dth =  DEG_TO_RAD(dth_num);
  
   g.h1 = h1; 
   g.h2 = h2; 
   g.R = R;    
   g.ellipticity = ellipticity; 
   g.density = density;
   g.ELAI = ELAI;  
   g.dth = dth;
   g.dz = dz;
   g.ds = ds; 
   g.k = k; 
   g.rl = rl; 
   g.tl= tl; 
   g.rg = rg;  

   g.szn = szn ; 
   g.vzn = vzn; 
   g.saz = saz;
   g.vaz = vaz; 
   g.slope= slope; 
   g.aspect = aspect; 

    /* call gort_apply fuction to calculate the gap probality and openness factors, then get the gort brdf return */
   gort_apply(&brdfValue);

 /* print the return gort brdf reValue */
  printf("the gort_brdf value: %lf\n", brdfValue);

}

