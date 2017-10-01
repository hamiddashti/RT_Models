#ifndef GORT1_H
#define GORT1_H

#define RAD_TO_DEG(r) ((r) * 57.29577951)
#define DEG_TO_RAD(d) ((d) / 57.29577951)
#define PI 3.141592653589793
#define NORDERS 10		/* number of orders of scattering */

#define MIN(a, b) (a < b ? a : b)
#define MAX(a, b) (a > b ? a : b)
 
#define MAXCROWNS 30 


/*
 * Type of floating-point variables to use -- either float or double.
 */
typedef double   fp_t;

/*
 * Global variable structure
 */
typedef struct {
   fp_t            szn;		/* solar zenith */
   int             szni;	/* array index associated w/ solar zenith */
   fp_t            h1;		/* lower boundary of the volume containing
				 * crown centers */
   fp_t            h2;		/* upper boundary of the volume containing
				 * crown centers */
   fp_t            z1;		/* lower boundary of integration */
   fp_t            z2;		/* upper boundary of integration */
   fp_t            density;	/* tree centers per unit ground area */
   fp_t            Lv;		/* tree centers per unit canopy volume */

   fp_t	 	   FAVD;        /* foliage area volume density  (1/m) */
   fp_t            ELAI;
   fp_t		   k;		/* extinction coefficient */
   fp_t            tau;         /* projected foliage density in the direction */ 


/* begin the variable declariation in transformed dimension */
   int             szni_p;	/* array index associated w/ solar zenith */
   fp_t            h1_p;	/* lower boundary of the volume containing
				 * crown centers */
   fp_t            h2_p;        /* upper boundary of the volume containing
				 * crown centers */
   fp_t            z1_p;	/* lower boundary of integration */
   fp_t            z2_p;	/* upper boundary of integration */
   fp_t            Lv_p;	/* tree centers per unit canopy volume */
   fp_t	 	   FAVD_p;      /* foliage area volume density  (1/m) */
   fp_t            tau_p;       /* projected foliage density in the direction */ 

/* begin the variable declariation in transformed dimension */

   fp_t            ellipticity;	/* the b/r ratio, where 'b' is vertical crown
				 * radius and 'r' is horizontal crown radius */
   fp_t            R;		/* crown radius */
   fp_t            R_squared;   /* R*R  */
   fp_t            R_cubed;     /* R*R*R  */
   fp_t            H;

   fp_t            albedo;      /* single scattering albedo of leaf */
   fp_t            rs;		/* surface reflectance */
   fp_t           dth;		/* theta differential */
   fp_t            dz;		/* height differential */
   fp_t            dz_p;	/* height differential in the transformed dimension */
   fp_t            ds;		/* path length differential */
   int             nlayers;	/* in numerical integration */
   int             nth;		/* number of theta values */
   int             nh1;		/* number at h1 */
   int             nh2;		/* number at h2  */

   /* new add for read file   */
  fp_t        rl;         /* rl: leaf reflectance */
  fp_t        tl;         /* tl: leaf transmittance */
  fp_t        rg;          /* rg: background reflectance	*/
  fp_t       vzn;	      /* viewing zenith angle (vzn) */
  fp_t       saz;          /* solar azimuth angle (saz)  */
  fp_t       vaz;	      /* viewing azimuth angle (vaz) */
  fp_t     slope;          /* slope */
  fp_t     aspect;         /* aspect */   

  fp_t       k_open;      /* the first openess factor calculate from gotr1.c */
  fp_t       k_openEP;	  /* the second openess factor calculate from gotr1.c */
  fp_t       k_opentot;   /* sum of the first and second openess factors */
  fp_t       brdf;	  /* calculating gort brdf from gort_brdf.c  */
}               GLOBAL_T;

extern GLOBAL_T g;


extern fp_t          *height;	/* array of values for height */
extern fp_t          *height_p;	/* array of values for height prime  */
extern fp_t          *theta;	/* array of values for theta */
extern fp_t          *theta_p;	/* array of values for theta prime */
extern fp_t          **P_n0;	/* array of values for P(n = 0 | h) */
extern fp_t          **P_s0;	/* array of values of P(s = 0 | h) */
extern fp_t          **Js;
extern fp_t          **T_open;
extern fp_t          **dT_open;
extern fp_t          **Js_sun;
extern fp_t          **Js_sky;
extern fp_t          **V_g;	/* array of values of V sub gamma */
extern fp_t          **s_p;	/* array of s' values */
extern fp_t          **EPgap;	/* gap probability table */
extern fp_t          *K_open;        
extern fp_t          *K_openEP;           
extern fp_t        ***PD_s;	/* 3-dimensional array holding values of P(s
				 * | h, theta) */
extern fp_t          *ES;  
extern fp_t          *factorial;
extern fp_t   *Vb;
extern fp_t   **fB;
extern fp_t   *dK_open;
extern fp_t   *Lk_up;
extern fp_t   *Lk_down;

extern void     get_PD_s();
extern void     get_EPgap();
extern double   get_K_open();
extern void     get_T_open();
extern void     get_Vb();
extern void     get_fB();


extern fp_t   **alloc_2d();
extern fp_t    *alloc_1d();
extern fp_t     get_ES();
extern fp_t     index_to_s();
extern int      s_to_index();
extern fp_t     sec();

extern fp_t     get_T();

#endif












