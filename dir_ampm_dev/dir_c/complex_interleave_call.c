#ifndef _MONOLITH
#include "ampm_header.h"
#endif /* _MONOLITH */

/* %%%%%%%% */
/* complex interleaving: */
/* %%%%%%%% */
double complex * double_complex_malloc_and_interleave(int n_a,double *double_real_,double *double_imag_)
{
  int na=0;
  double complex *double_complex_=NULL;
  double_complex_ = (double complex *) malloc(n_a*sizeof(double complex));
  for (na=0;na<n_a;na++){ double_complex_[na] = double_real_[na] + _Complex_I*double_imag_[na];}
  return double_complex_;
}

/* %%%%%%%% */
/* aligned float allocation: */
/* %%%%%%%% */
float* float_m256_malloc_from_double_(int n_a,double *d_)
{
  int na=0;
  int n_a_rup = rup(n_a,8);
  int n_a_256 = n_a_rup/8; //%<-- 8 floats per __m256. ;
  float *f_=NULL;
  f_ = (float *) _mm_malloc(n_a_256*sizeof(__m256),32);
  memset(f_,0,n_a_256*8*sizeof(float));
  if (d_!=NULL){
    for (na=0;na<n_a;na++){ f_[na] = (float)(d_[na]);}
  /* if (d_!=NULL){ } */}
  return f_;
}

float* float__m256_malloc_from_double__(int n_a,int n_b,double *d__)
{
  int na=0;
  int n_a_rup = rup(n_a,8);
  int n_a_256 = n_a_rup/8; //%<-- 8 floats per __m256. ;
  int nb=0;
  int n_b_rup = rup(n_b,8);
  unsigned long long int ulli=0;
  float *f__=NULL;
  ulli = (unsigned long long int)n_a_256*(unsigned long long int)n_b_rup;
  f__ = (float *) _mm_malloc(ulli*sizeof(__m256),32);
  memset(f__,0,ulli*8*sizeof(float));
  if (d__!=NULL){
    for (nb=0;nb<n_b;nb++){
      for (na=0;na<n_a;na++){
	f__[na + nb*n_a_rup] = (float)(d__[na + nb*n_a]);
	/* for (na=0;na<n_a;na++){  } */}
      /* for (nb=0;nb<n_b;nb++){  } */}
    /* if (d__!=NULL){ } */}
  return f__;
}

float* float___m256_malloc_from_double___(int n_a,int n_b,int n_c,double *d___)
{
  int na=0;
  int n_a_rup = rup(n_a,8);
  int n_a_256 = n_a_rup/8; //%<-- 8 floats per __m256. ;
  int nb=0;
  int n_b_rup = rup(n_b,8);
  int nc=0;
  int n_c_rup = rup(n_c,8);
  unsigned long long int ulli=0;
  float *f___=NULL;
  ulli = (unsigned long long int)n_a_256*(unsigned long long int)n_b_rup*(unsigned long long int)n_c_rup;
  f___ = (float *) _mm_malloc(ulli*sizeof(__m256),32);
  memset(f___,0,ulli*8*sizeof(float));
  if (d___!=NULL){
    for (nc=0;nc<n_c;nc++){
      for (nb=0;nb<n_b;nb++){
	for (na=0;na<n_a;na++){
	  f___[na + (nb + nc*n_b_rup)*n_a_rup] = (float)(d___[na + (nb + nc*n_b)*n_a]);
	  /* for (na=0;na<n_a;na++){  } */}
	/* for (nb=0;nb<n_b;nb++){  } */}
      /* for (nc=0;nc<n_c;nc++){  } */}
    /* if (d___!=NULL){ } */}
  return f___;
}

float* float____m256_malloc_from_double____(int n_a,int n_b,int n_c,int n_d,double *d____)
{
  int na=0;
  int n_a_rup = rup(n_a,8);
  int n_a_256 = n_a_rup/8; //%<-- 8 floats per __m256. ;
  int nb=0;
  int n_b_rup = rup(n_b,8);
  int nc=0;
  int n_c_rup = rup(n_c,8);
  int nd=0;
  int n_d_rup = rup(n_d,8);
  unsigned long long int ulli=0;
  float *f____=NULL;
  ulli = (unsigned long long int)n_a_256*(unsigned long long int)n_b_rup*(unsigned long long int)n_c_rup*(unsigned long long int)n_d_rup;
  f____ = (float *) _mm_malloc(ulli*sizeof(__m256),32);
  memset(f____,0,ulli*8*sizeof(float));
  if (d____!=NULL){
    for (nd=0;nd<n_d;nd++){
      for (nc=0;nc<n_c;nc++){
	for (nb=0;nb<n_b;nb++){
	  for (na=0;na<n_a;na++){
	    f____[na + (nb + (nc + nd*n_c_rup)*n_b_rup)*n_a_rup] = (float)(d____[na + (nb + (nc + nd*n_c)*n_b)*n_a]);
	    /* for (na=0;na<n_a;na++){  } */}
	  /* for (nb=0;nb<n_b;nb++){  } */}
	/* for (nc=0;nc<n_c;nc++){  } */}
      /* for (nd=0;nd<n_d;nd++){  } */}
    /* if (d____!=NULL){ } */}
  return f____;
}

double l2_double_complex__vs_float_m256__(int n_a,int n_b,double complex *z__,float *f_R__,float *f_I__)
{
  int na=0;
  int n_a_rup = rup(n_a,8);
  int nb=0;
  int n_b_rup = rup(n_b,8);
  unsigned long long int tabA=0,tabB=0;
  double l2_R=0,l2_I=0,l2=0;
  for (nb=0;nb<n_b;nb++){
    for (na=0;na<n_a;na++){
      tabA =
	(unsigned long long int)na +
	(unsigned long long int)nb
	*(unsigned long long int)n_a;
      tabB =
	(unsigned long long int)na +
	(unsigned long long int)nb
	*(unsigned long long int)n_a_rup;
      l2_R = (creal(z__[tabA]) - (double)f_R__[tabB]);
      l2_I = (cimag(z__[tabA]) - (double)f_I__[tabB]);
      l2 += l2_R*l2_R + l2_I*l2_I;
      /* for (na=0;na<n_a;na++){  } */}
    /* for (nb=0;nb<n_b;nb++){  } */}
  return l2;
}

double l2_double_complex___vs_float_m256___(int n_a,int n_b,int n_c,double complex *z___,float *f_R___,float *f_I___)
{
  int na=0;
  int n_a_rup = rup(n_a,8);
  int nb=0;
  int n_b_rup = rup(n_b,8);
  int nc=0;
  int n_c_rup = rup(n_c,8);
  unsigned long long int tabA=0,tabB=0;
  double l2_R=0,l2_I=0,l2=0;
  for (nc=0;nc<n_c;nc++){
    for (nb=0;nb<n_b;nb++){
      for (na=0;na<n_a;na++){
	tabA =
	  (unsigned long long int)na +
	  ((unsigned long long int)nb +
	   (unsigned long long int)nc
	   *(unsigned long long int)n_b)
	  *(unsigned long long int)n_a;
	tabB =
	  (unsigned long long int)na +
	  ((unsigned long long int)nb +
	   (unsigned long long int)nc
	   *(unsigned long long int)n_b_rup)
	  *(unsigned long long int)n_a_rup;
	l2_R = (creal(z___[tabA]) - (double)f_R___[tabB]);
	l2_I = (cimag(z___[tabA]) - (double)f_I___[tabB]);
	l2 += l2_R*l2_R + l2_I*l2_I;
	/* for (na=0;na<n_a;na++){  } */}
      /* for (nb=0;nb<n_b;nb++){  } */}
    /* for (nc=0;nc<n_c;nc++){  } */}
  return l2;
}

double l2_double_complex____vs_float_m256____(int n_a,int n_b,int n_c,int n_d,double complex *z____,float *f_R____,float *f_I____)
{
  int na=0;
  int n_a_rup = rup(n_a,8);
  int nb=0;
  int n_b_rup = rup(n_b,8);
  int nc=0;
  int n_c_rup = rup(n_c,8);
  int nd=0;
  int n_d_rup = rup(n_d,8);
  unsigned long long int tabA=0,tabB=0;
  double l2_R=0,l2_I=0,l2=0;
  for (nd=0;nd<n_d;nd++){
    for (nc=0;nc<n_c;nc++){
      for (nb=0;nb<n_b;nb++){
	for (na=0;na<n_a;na++){
	  tabA =
	    (unsigned long long int)na +
	    ((unsigned long long int)nb +
	     ((unsigned long long int)nc +
	      (unsigned long long int)nd
	      *(unsigned long long int)n_c)
	     *(unsigned long long int)n_b)
	    *(unsigned long long int)n_a;
	  tabB =
	    (unsigned long long int)na +
	    ((unsigned long long int)nb +
	     ((unsigned long long int)nc +
	      (unsigned long long int)nd
	      *(unsigned long long int)n_c_rup)
	     *(unsigned long long int)n_b_rup)
	    *(unsigned long long int)n_a_rup;
	  l2_R = (creal(z____[tabA]) - (double)f_R____[tabB]);
	  l2_I = (cimag(z____[tabA]) - (double)f_I____[tabB]);
	  l2 += l2_R*l2_R + l2_I*l2_I;
	  /* for (na=0;na<n_a;na++){  } */}
	/* for (nb=0;nb<n_b;nb++){  } */}
      /* for (nc=0;nc<n_c;nc++){  } */}
    /* for (nd=0;nd<n_d;nd++){  } */}
  return l2;
}

