/*==========================================================
 * mex_ampmh_X_wSM___10.c ;
 * compile from Matlab with: ;
 * mex CFLAGS='$CFLAGS -fPIC -w -O2 -D_AVX -mavx -D_CBLAS -D_FMA -mfma' mex_ampmh_X_wSM___10.c -lm -lfftw3 -lgslcblas -output /data/rangan/dir_cryoem/dir_rangan_playroom/mex_ampmh_X_wSM___10 ;
 * compile from shell with: ;
 * gcc -fPIC -w -O2 -DWITHOUT_MEX -D_AVX -mavx -D_CBLAS -D_FMA -mfma mex_ampmh_X_wSM___10.c -lm -lfftw3 -lgslcblas ; ./a.out;
 * test with: ampmh_X_wSM_mex___8.m ;
 *========================================================*/
 
/*========================================================== 
 * trying to build a faster ampmh_X_wSM___8; not there yet! ;
 * still need to actually calculate the innerproducts (loop through nS). ;
 * also need to implement: ;
 * aligned memory ;
 * intrinsic multiplication for segregated complex ;
 * fast transpose ;
 * correct configuration for fftw_plan_many ;
 * another fast transpose ;
 *========================================================*/

/*========================================================== 
 * testing out different kinds of transposes ;
 *========================================================*/

#ifndef WITHOUT_MEX
#include "mex.h"
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <cblas.h>
#include <string.h>
#include <complex.h> /* <-- must include before fftw3.h to ensure that fftw_complex is compatible */
#include "fftw3.h"  /* <-- must include after complex.h to ensure that fftw_complex is compatible */
/* #include <omp.h> */
#include <emmintrin.h>
#include <immintrin.h>

#define PI (3.141592653589793)
#define rup(A,B) ((A) + !!((A)%(B))*((B) - ((A)%(B))))
#define maximum(A,B) ((A) > (B) ? (A) : (B))
#define minimum(A,B) ((A) < (B) ? (A) : (B))
#define periodize(A,B,C) ((A) < (B) ? (A) + (C) - (B) : ((A) >= (C) ? (A) - (C) + (B) : (A)))
#define crop(A,B,C) ((A) < (B) ? (B) : ((A) > (C) ? (C) : (A)))

/* %%%%%%%% */
/* array_printf: */
/* %%%%%%%% */
void array_sub_printf(void *v_,const char *type,int n_row,int n_row_sub,int n_col,int n_col_sub,const char *prefix)
{
  /* prints out arrays of varying types */
  int nr=0,nc=0,tmp=0;
  float *f_=NULL;
  double *d_=NULL; int d_l_flag=0;
  int *i_=NULL;
  unsigned int *ui_=NULL;
  long *l_=NULL;
  long long *ll_=NULL;
  unsigned long int *ul_=NULL;
  unsigned long long int *ull_=NULL;
  char *char_=NULL;
  unsigned char *uchar_=NULL;
  float complex *c_=NULL;
  double complex *z_=NULL;
  double printftol=0.000000001;
  if (strcmp(type,"float")==0){
    f_ = (float *) v_;
    for (nr=0;nr<n_row_sub;nr++){ 
      printf("%s",prefix); 
      for (nc=0;nc<n_col_sub;nc++){ 
	if (fabs(f_[nr+nc*n_row]-(int)f_[nr+nc*n_row])<printftol){ printf(" %s%d",(int)f_[nr+nc*n_row]>0 ? "+" : ((int)f_[nr+nc*n_row]<0 ? "" : " "),(int)f_[nr+nc*n_row]);}
	else{ printf(" %s%f",f_[nr+nc*n_row]>0 ? "+" : (f_[nr+nc*n_row]<0 ? "" : " "),f_[nr+nc*n_row]);}}
      printf("\n");}}
  else if (strcmp(type,"double")==0){
    if (n_row>1 && n_col==1){ tmp=n_row; n_row=n_col; n_col=tmp;}
    d_ = (double *) v_;
    if (0){
      d_l_flag=1; for (nr=0;nr<n_row_sub;nr++){ for (nc=0;nc<n_col_sub;nc++){ if (fabs(d_[nr+nc*n_row]-(int)d_[nr+nc*n_row])>printftol){ d_l_flag=0;}}}
      for (nr=0;nr<n_row_sub;nr++){ 
	printf("%s",prefix); 
	for (nc=0;nc<n_col_sub;nc++){ 
	  if (d_l_flag && fabs(d_[nr+nc*n_row]-(int)d_[nr+nc*n_row])<printftol){ printf(" %s%d",(int)d_[nr+nc*n_row]>0 ? "+" : ((int)d_[nr+nc*n_row]<0 ? "" : " "),(int)d_[nr+nc*n_row]);}
	  else{ /* printf(" %s%.3f",d_[nr+nc*n_row]>0 ? "+" : (d_[nr+nc*n_row]<0 ? "" : " "),d_[nr+nc*n_row]); */ 
	    if (fabs(d_[nr+nc*n_row])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",d_[nr+nc*n_row]);}
	    /* if !d_l_flag */}
	  /* for (nc=0;nc<n_col_sub;nc++){ } */}
	printf("\n");} /* if 0 */}
    if (1){
      for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); 
	for (nc=0;nc<n_col_sub;nc++){ if (fabs(d_[nr+nc*n_row])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",d_[nr+nc*n_row]);} /* for (nc=0;nc<n_col_sub;nc++){ } */}
	printf("\n"); /* for (nr=0;nr<n_row_sub;nr++){ } */} /* if 1 */}
    /* if double */}
  else if (strcmp(type,"double_trn")==0 && n_row>1){
    d_ = (double *) v_;
    for (nc=0;nc<n_col_sub;nc++){ printf("%s",prefix); 
      for (nr=0;nr<n_row_sub;nr++){ if (fabs(d_[nr+nc*n_row])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",d_[nr+nc*n_row]);} /* for (nr=0;nr<n_row_sub;nr++){ } */}
      printf("\n"); /* for (nc=0;nc<n_col_sub;nc++){ } */}
    /* if double_trn */}
  else if (strcmp(type,"float complex")==0){
    if (n_row>1 && n_col==1){ tmp=n_row; n_row=n_col; n_col=tmp;}
    c_ = (float complex *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); 
      for (nc=0;nc<n_col_sub;nc++){
	if (cabsf(c_[nr+nc*n_row])<printftol){ printf("    .        .   ");} else{ printf(" [%+07.3f %+07.3fi]",crealf(c_[nr+nc*n_row]),cimagf(c_[nr+nc*n_row]));} 
	/* for (nc=0;nc<n_col_sub;nc++){ } */}
      printf("\n"); /* for (nr=0;nr<n_row_sub;nr++){ } */}
    /* else if (strcmp(type,"float complex")==0){ } */}
  else if (strcmp(type,"double complex")==0){
    if (n_row>1 && n_col==1){ tmp=n_row; n_row=n_col; n_col=tmp;}
    z_ = (double complex *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); 
      for (nc=0;nc<n_col_sub;nc++){
	if (cabs(z_[nr+nc*n_row])<printftol){ printf("    .        .   ");} else{ printf(" [%+07.3f %+07.3fi]",creal(z_[nr+nc*n_row]),cimag(z_[nr+nc*n_row]));} 
	/* for (nc=0;nc<n_col_sub;nc++){ } */}
      printf("\n"); /* for (nr=0;nr<n_row_sub;nr++){ } */}
    /* else if (strcmp(type,"double complex")==0){ } */}
  else if (strcmp(type,"long")==0){
    l_ = (long *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %ld",l_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"long long int")==0){
    ll_ = (long long *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %lld",ll_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"int")==0){
    i_ = (int *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %d",i_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned int")==0){
    ui_ = (unsigned int *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %d",(int)ui_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned long int")==0){
    ul_ = (unsigned long int *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %ld",ul_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned long long int")==0){
    ull_ = (unsigned long long int *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %lld",ull_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"alpha")==0){
    char_ = (char *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %c",(int)char_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"char")==0){
    char_ = (char *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %d",(int)char_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned char")==0){
    uchar_ = (unsigned char *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %d",(int)uchar_[nr+nc*n_row]);} printf("\n");}}
  else{ printf(" warning, poor type %s in array_printf\n",type);}
}
void array_printf(void *v_,const char *type,int n_row,int n_col,const char *prefix)
{
  /* prints out arrays of varying types */
  int nr=0,nc=0,tmp=0;
  float *f_=NULL;
  double *d_=NULL; int d_l_flag=0;
  int *i_=NULL;
  unsigned int *ui_=NULL;
  long *l_=NULL;
  long long *ll_=NULL;
  unsigned long int *ul_=NULL;
  unsigned long long int *ull_=NULL;
  char *char_=NULL;
  unsigned char *uchar_=NULL;
  float complex *c_=NULL;
  double complex *z_=NULL;
  double printftol=0.000000001;
  if (strcmp(type,"float")==0){
    f_ = (float *) v_;
    for (nr=0;nr<n_row;nr++){ 
      printf("%s",prefix); 
      for (nc=0;nc<n_col;nc++){ 
	if (fabs(f_[nr+nc*n_row]-(int)f_[nr+nc*n_row])<printftol){ printf(" %s%d",(int)f_[nr+nc*n_row]>0 ? "+" : ((int)f_[nr+nc*n_row]<0 ? "" : " "),(int)f_[nr+nc*n_row]);}
	else{ printf(" %s%f",f_[nr+nc*n_row]>0 ? "+" : (f_[nr+nc*n_row]<0 ? "" : " "),f_[nr+nc*n_row]);}}
      printf("\n");}}
  else if (strcmp(type,"double")==0){
    if (n_row>1 && n_col==1){ tmp=n_row; n_row=n_col; n_col=tmp;}
    d_ = (double *) v_;
    if (0){
      d_l_flag=1; for (nr=0;nr<n_row;nr++){ for (nc=0;nc<n_col;nc++){ if (fabs(d_[nr+nc*n_row]-(int)d_[nr+nc*n_row])>printftol){ d_l_flag=0;}}}
      for (nr=0;nr<n_row;nr++){ 
	printf("%s",prefix); 
	for (nc=0;nc<n_col;nc++){ 
	  if (d_l_flag && fabs(d_[nr+nc*n_row]-(int)d_[nr+nc*n_row])<printftol){ printf(" %s%d",(int)d_[nr+nc*n_row]>0 ? "+" : ((int)d_[nr+nc*n_row]<0 ? "" : " "),(int)d_[nr+nc*n_row]);}
	  else{ /* printf(" %s%.3f",d_[nr+nc*n_row]>0 ? "+" : (d_[nr+nc*n_row]<0 ? "" : " "),d_[nr+nc*n_row]); */ 
	    if (fabs(d_[nr+nc*n_row])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",d_[nr+nc*n_row]);}
	    /* if !d_l_flag */}
	  /* for (nc=0;nc<n_col;nc++){ } */}
	printf("\n");} /* if 0 */}
    if (1){
      for (nr=0;nr<n_row;nr++){ printf("%s",prefix); 
	for (nc=0;nc<n_col;nc++){ if (fabs(d_[nr+nc*n_row])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",d_[nr+nc*n_row]);} /* for (nc=0;nc<n_col;nc++){ } */}
	printf("\n"); /* for (nr=0;nr<n_row;nr++){ } */} /* if 1 */}
    /* if double */}
  else if (strcmp(type,"double_trn")==0 && n_row>1){
    d_ = (double *) v_;
    for (nc=0;nc<n_col;nc++){ printf("%s",prefix); 
      for (nr=0;nr<n_row;nr++){ if (fabs(d_[nr+nc*n_row])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",d_[nr+nc*n_row]);} /* for (nr=0;nr<n_row;nr++){ } */}
      printf("\n"); /* for (nc=0;nc<n_col;nc++){ } */}
    /* if double_trn */}
  else if (strcmp(type,"float complex")==0){
    if (n_row>1 && n_col==1){ tmp=n_row; n_row=n_col; n_col=tmp;}
    c_ = (float complex *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); 
      for (nc=0;nc<n_col;nc++){
	if (cabsf(c_[nr+nc*n_row])<printftol){ printf("    .        .   ");} else{ printf(" [%+07.3f %+07.3fi]",crealf(c_[nr+nc*n_row]),cimagf(c_[nr+nc*n_row]));} 
	/* for (nc=0;nc<n_col;nc++){ } */}
      printf("\n"); /* for (nr=0;nr<n_row;nr++){ } */}
    /* else if (strcmp(type,"float complex")==0){ } */}
  else if (strcmp(type,"double complex")==0){
    if (n_row>1 && n_col==1){ tmp=n_row; n_row=n_col; n_col=tmp;}
    z_ = (double complex *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); 
      for (nc=0;nc<n_col;nc++){
	if (cabs(z_[nr+nc*n_row])<printftol){ printf("    .        .   ");} else{ printf(" [%+07.3f %+07.3fi]",creal(z_[nr+nc*n_row]),cimag(z_[nr+nc*n_row]));} 
	/* for (nc=0;nc<n_col;nc++){ } */}
      printf("\n"); /* for (nr=0;nr<n_row;nr++){ } */}
    /* else if (strcmp(type,"double complex")==0){ } */}
  else if (strcmp(type,"long")==0){
    l_ = (long *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %ld",l_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"long long int")==0){
    ll_ = (long long *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %lld",ll_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"int")==0){
    i_ = (int *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %d",i_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned int")==0){
    ui_ = (unsigned int *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %d",(int)ui_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned long int")==0){
    ul_ = (unsigned long int *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %ld",ul_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned long long int")==0){
    ull_ = (unsigned long long int *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %lld",ull_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"alpha")==0){
    char_ = (char *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %c",(int)char_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"char")==0){
    char_ = (char *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %d",(int)char_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned char")==0){
    uchar_ = (unsigned char *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %d",(int)uchar_[nr+nc*n_row]);} printf("\n");}}
  else{ printf(" warning, poor type %s in array_printf\n",type);}
}

/* %%%%%%%% */
/* frobenius norms: */
/* %%%%%%%% */
double ifnormn(unsigned long long int ulli_length,int *i_0_,int *i_1_)
{
  unsigned long long int ulli=0;
  int i_01 = 0;
  int i_0 = 0;
  double d_f01 = 0;
  double d_f0 = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ i_01 = i_0_[ulli]-i_1_[ulli]; d_f01 += i_01*i_01; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  d_f01 = sqrt(d_f01);
  for (ulli=0;ulli<ulli_length;ulli++){ i_0 = i_0_[ulli]; d_f0 += i_0*i_0; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  d_f0 = sqrt(d_f0);
  return d_f01/maximum(1e-12,d_f0);
}
double ullifnormn(unsigned long long int ulli_length,unsigned long long int *ulli_0_,unsigned long long int *ulli_1_)
{
  unsigned long long int ulli=0;
  unsigned long long int ulli_01 = 0;
  unsigned long long int ulli_0 = 0;
  double d_f01 = 0;
  double d_f0 = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ ulli_01 = ulli_0_[ulli]-ulli_1_[ulli]; d_f01 += ulli_01*ulli_01; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  d_f01 = sqrt(d_f01);
  for (ulli=0;ulli<ulli_length;ulli++){ ulli_0 = ulli_0_[ulli]; d_f0 += ulli_0*ulli_0; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  d_f0 = sqrt(d_f0);
  return d_f01/maximum(1e-12,d_f0);
}
float ffnorm(unsigned long long int ulli_length,float *f_0_,float *f_1_)
{
  float output = 0;
  unsigned long long int ulli=0;
  float f_diff = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ f_diff = f_0_[ulli]-f_1_[ulli]; output += f_diff*f_diff; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  return sqrtf(output);
}
float ffnormn(unsigned long long int ulli_length,float *f_0_,float *f_1_)
{
  unsigned long long int ulli=0;
  float f_01 = 0;
  float f_0 = 0;
  float f_f01 = 0;
  float f_f0 = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ f_01 = f_0_[ulli]-f_1_[ulli]; f_f01 += f_01*f_01; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  f_f01 = sqrt(f_f01);
  for (ulli=0;ulli<ulli_length;ulli++){ f_0 = f_0_[ulli]; f_f0 += f_0*f_0; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  f_f0 = sqrtf(f_f0);
  return f_f01/maximum(1e-12,f_f0);
}
double dfnormn(unsigned long long int ulli_length,double *d_0_,double *d_1_)
{
  unsigned long long int ulli=0;
  double d_01 = 0;
  double d_0 = 0;
  double d_f01 = 0;
  double d_f0 = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ d_01 = d_0_[ulli]-d_1_[ulli]; d_f01 += d_01*d_01; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  d_f01 = sqrt(d_f01);
  for (ulli=0;ulli<ulli_length;ulli++){ d_0 = d_0_[ulli]; d_f0 += d_0*d_0; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  d_f0 = sqrt(d_f0);
  return d_f01/maximum(1e-12,d_f0);
}
double dfnorm(unsigned long long int ulli_length,double *d_0_,double *d_1_)
{
  double output = 0;
  unsigned long long int ulli=0;
  double d_diff = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ d_diff = d_0_[ulli]-d_1_[ulli]; output += d_diff*d_diff; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  return sqrt(output);
}
double cfnorm(unsigned long long int ulli_length,float complex *c_0_,float complex *c_1_)
{
  double output = 0;
  unsigned long long int ulli=0;
  float complex c_diff = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ c_diff = c_0_[ulli]-c_1_[ulli]; output += crealf(c_diff*conjf(c_diff)); /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  return sqrt(output);
}
double zfnorm(unsigned long long int ulli_length,double complex *z_0_,double complex *z_1_)
{
  double output = 0;
  unsigned long long int ulli=0;
  double complex z_diff = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ z_diff = z_0_[ulli]-z_1_[ulli]; output += creal(z_diff*conj(z_diff)); /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  return sqrt(output);
}
double zfnormn(unsigned long long int ulli_length,double complex *z_0_,double complex *z_1_)
{
  double z_f01 = 0;
  double z_f0 = 0;
  unsigned long long int ulli=0;
  double complex z_01=0,z_0=0;
  for (ulli=0;ulli<ulli_length;ulli++){ z_01 = z_0_[ulli]-z_1_[ulli]; z_f01 += creal(z_01*conj(z_01)); /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  z_f01 = sqrt(z_f01);
  for (ulli=0;ulli<ulli_length;ulli++){ z_0 = z_0_[ulli]; z_f0 += creal(z_0*conj(z_0)); /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  z_f0 = sqrt(z_f0);
  return z_f01/maximum(1e-12,z_f0);
}

/* %%%%%%%% */
/* MDA io call */
/* %%%%%%%% */
void MDA_write_i4(int n_dim,int *dim_,int *i4_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(dim_,sizeof(int),n_dim,fp); if (s!=n_dim){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  for (ndim=0;ndim<n_dim;ndim++){ n_i*=dim_[ndim];}
  s=fwrite(i4_,sizeof(int),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  fclose(fp);fp=NULL;
}
void MDA_read_i4(int *n_dim_p,int **dim_p_,int **i4_p_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  int n_dim=0;
  int *dim_=NULL;
  int *i4_=NULL;
  if (verbose){ printf(" %% reading %s\n",fname);}
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fread(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
  if (n_dim_p!=NULL){ *n_dim_p = n_dim;}
  dim_=NULL;
  if (dim_p_!=NULL){
    if ( (*dim_p_)==NULL ){ (*dim_p_) = (int *) malloc(n_dim*sizeof(int)); }
    dim_ = *dim_p_;
    /* if (dim_p_!=NULL){ } */}
  if (dim_!=NULL){
    s=fread(dim_,sizeof(int),n_dim,fp);
    if (s!=n_dim){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
    n_i = 1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)(dim_[ndim]);}
    i4_=NULL;
    if (i4_p_!=NULL){
      if ( (*i4_p_)==NULL ){ (*i4_p_) = (int *) malloc(n_i*sizeof(int)); }
      i4_ = *i4_p_;
      /* if (i4_p_!=NULL){ } */}
    if (i4_!=NULL){
      s=fread(i4_,sizeof(int),n_i,fp);
      if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}      
      /* if (i4_!=NULL){ } */}
    /* if (dim_!=NULL){ } */}
  fclose(fp);fp=NULL;
}
void MDA_write_r4(int n_dim,int *dim_,float *r4_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(dim_,sizeof(int),n_dim,fp); if (s!=n_dim){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  for (ndim=0;ndim<n_dim;ndim++){ n_i*=dim_[ndim];}
  s=fwrite(r4_,sizeof(float),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  fclose(fp);fp=NULL;
}
void MDA_read_r4(int *n_dim_p,int **dim_p_,float **r4_p_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  int n_dim=0;
  int *dim_=NULL;
  float *r4_=NULL;
  if (verbose){ printf(" %% reading %s\n",fname);}
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fread(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
  if (n_dim_p!=NULL){ *n_dim_p = n_dim;}
  dim_=NULL;
  if (dim_p_!=NULL){
    if ( (*dim_p_)==NULL ){ (*dim_p_) = (int *) malloc(n_dim*sizeof(int)); }
    dim_ = *dim_p_;
    /* if (dim_p_!=NULL){ } */}
  if (dim_!=NULL){
    s=fread(dim_,sizeof(int),n_dim,fp);
    if (s!=n_dim){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
    n_i = 1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)(dim_[ndim]);}
    r4_=NULL;
    if (r4_p_!=NULL){
      if ( (*r4_p_)==NULL ){ (*r4_p_) = (float *) malloc(n_i*sizeof(float)); }
      r4_ = *r4_p_;
      /* if (r4_p_!=NULL){ } */}
    if (r4_!=NULL){
      s=fread(r4_,sizeof(float),n_i,fp);
      if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}      
      /* if (r4_!=NULL){ } */}
    /* if (dim_!=NULL){ } */}
  fclose(fp);fp=NULL;
}
void MDA_write_r8(int n_dim,int *dim_,double *r8_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(dim_,sizeof(int),n_dim,fp); if (s!=n_dim){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  for (ndim=0;ndim<n_dim;ndim++){ n_i*=dim_[ndim];}
  s=fwrite(r8_,sizeof(double),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  fclose(fp);fp=NULL;
}
void MDA_read_r8(int *n_dim_p,int **dim_p_,double **r8_p_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  int n_dim=0;
  int *dim_=NULL;
  double *r8_=NULL;
  if (verbose){ printf(" %% reading %s\n",fname);}
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fread(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
  if (n_dim_p!=NULL){ *n_dim_p = n_dim;}
  dim_=NULL;
  if (dim_p_!=NULL){
    if ( (*dim_p_)==NULL ){ (*dim_p_) = (int *) malloc(n_dim*sizeof(int)); }
    dim_ = *dim_p_;
    /* if (dim_p_!=NULL){ } */}
  if (dim_!=NULL){
    s=fread(dim_,sizeof(int),n_dim,fp);
    if (s!=n_dim){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
    n_i = 1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)(dim_[ndim]);}
    r8_=NULL;
    if (r8_p_!=NULL){
      if ( (*r8_p_)==NULL ){ (*r8_p_) = (double *) malloc(n_i*sizeof(double)); }
      r8_ = *r8_p_;
      /* if (r8_p_!=NULL){ } */}
    if (r8_!=NULL){
      s=fread(r8_,sizeof(double),n_i,fp);
      if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}      
      /* if (r8_!=NULL){ } */}
    /* if (dim_!=NULL){ } */}
  fclose(fp);fp=NULL;
}
void MDA_write_c16(int n_dim,int *dim_,double complex *c16_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(dim_,sizeof(int),n_dim,fp); if (s!=n_dim){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  for (ndim=0;ndim<n_dim;ndim++){ n_i*=dim_[ndim];}
  s=fwrite(c16_,sizeof(double complex),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  fclose(fp);fp=NULL;
}
void MDA_read_c16(int *n_dim_p,int **dim_p_,double complex **c16_p_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  int n_dim=0;
  int *dim_=NULL;
  double complex *c16_=NULL;
  if (verbose){ printf(" %% reading %s\n",fname);}
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fread(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
  if (n_dim_p!=NULL){ *n_dim_p = n_dim;}
  dim_=NULL;
  if (dim_p_!=NULL){
    if ( (*dim_p_)==NULL ){ (*dim_p_) = (int *) malloc(n_dim*sizeof(int)); }
    dim_ = *dim_p_;
    /* if (dim_p_!=NULL){ } */}
  if (dim_!=NULL){
    s=fread(dim_,sizeof(int),n_dim,fp);
    if (s!=n_dim){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
    n_i = 1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)(dim_[ndim]);}
    c16_=NULL;
    if (c16_p_!=NULL){
      if ( (*c16_p_)==NULL ){ (*c16_p_) = (double complex *) malloc(n_i*sizeof(double complex)); }
      c16_ = *c16_p_;
      /* if (c16_p_!=NULL){ } */}
    if (c16_!=NULL){
      s=fread(c16_,sizeof(double complex),n_i,fp);
      if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}      
      /* if (c16_!=NULL){ } */}
    /* if (dim_!=NULL){ } */}
  fclose(fp);fp=NULL;
}
void MDA_write_ulli(int n_dim,int *dim_,unsigned long long int *ulli_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(dim_,sizeof(int),n_dim,fp); if (s!=n_dim){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  for (ndim=0;ndim<n_dim;ndim++){ n_i*=dim_[ndim];}
  s=fwrite(ulli_,sizeof(unsigned long long int),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  fclose(fp);fp=NULL;
}
void MDA_read_ulli(int *n_dim_p,int **dim_p_,unsigned long long int **ulli_p_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  int n_dim=0;
  int *dim_=NULL;
  unsigned long long int *ulli_=NULL;
  if (verbose){ printf(" %% reading %s\n",fname);}
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fread(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
  if (n_dim_p!=NULL){ *n_dim_p = n_dim;}
  dim_=NULL;
  if (dim_p_!=NULL){
    if ( (*dim_p_)==NULL ){ (*dim_p_) = (int *) malloc(n_dim*sizeof(int)); }
    dim_ = *dim_p_;
    /* if (dim_p_!=NULL){ } */}
  if (dim_!=NULL){
    s=fread(dim_,sizeof(int),n_dim,fp);
    if (s!=n_dim){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
    n_i = 1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)(dim_[ndim]);}
    ulli_=NULL;
    if (ulli_p_!=NULL){
      if ( (*ulli_p_)==NULL ){ (*ulli_p_) = (unsigned long long int *) malloc(n_i*sizeof(unsigned long long int)); }
      ulli_ = *ulli_p_;
      /* if (ulli_p_!=NULL){ } */}
    if (ulli_!=NULL){
      s=fread(ulli_,sizeof(unsigned long long int),n_i,fp);
      if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}      
      /* if (ulli_!=NULL){ } */}
    /* if (dim_!=NULL){ } */}
  fclose(fp);fp=NULL;
}
void MDA_io_test()
{
  int verbose=1;
  int n_dim_0 = 2;
  int *dim_0_=NULL;
  int n_i=0;
  int ndim=0;
  int *i4_0_=NULL;
  int *i4_1_=NULL;
  int *dim_1_=NULL;
  int n_dim_1=NULL;
  char fname[32];
  double *r8_0_;
  int n_dim_2=NULL;
  int *r8_2_=NULL;
  int *dim_2_=NULL;
  float *r4_0_;
  int n_dim_3=NULL;
  int *r4_3_=NULL;
  int *dim_3_=NULL;
  unsigned long long int *ulli_0_;
  int n_dim_4=NULL;
  int *ulli_4_=NULL;
  int *dim_4_=NULL;
  double complex *c16_0_;
  int n_dim_5=NULL;
  int *c16_5_=NULL;
  int *dim_5_=NULL;
  dim_0_ = (int *) malloc(n_dim_0*sizeof(int));
  dim_0_[0] = 2; dim_0_[1] = 3;
  n_i = 1; for (ndim=0;ndim<n_dim_0;ndim++){ n_i*=dim_0_[ndim];}
  /* %%%%%%%% */
  i4_0_ = (int *) malloc(n_i*sizeof(int));
  i4_0_[0] = 0; i4_0_[2] = 10; i4_0_[4] = 100;
  i4_0_[1] = 1; i4_0_[3] = 20; i4_0_[5] = 200;
  sprintf(fname,"MDA_io.test");
  array_printf(i4_0_,"int",dim_0_[0],dim_0_[1]," %% i4_0_: ");
  MDA_write_i4(n_dim_0,dim_0_,i4_0_,fname);
  MDA_read_i4(&n_dim_1,&dim_1_,&i4_1_,fname);
  array_printf(i4_1_,"int",dim_1_[0],dim_1_[1]," %% i4_1_: ");
  printf(" %% i4_0_ vs i4_1_: relative error %0.16f\n",ifnormn(n_i,i4_0_,i4_1_));
  /* %%%%%%%% */
  r8_0_ = (double *) malloc(n_i*sizeof(double));
  r8_0_[0] = 0.1; r8_0_[2] = 10.1; r8_0_[4] = 100.1;
  r8_0_[1] = 1.1; r8_0_[3] = 20.1; r8_0_[5] = 200.1;
  array_printf(r8_0_,"double",dim_0_[0],dim_0_[1]," %% r8_0_: ");
  MDA_write_r8(n_dim_0,dim_0_,r8_0_,fname);
  MDA_read_r8(&n_dim_2,&dim_2_,&r8_2_,fname);
  array_printf(r8_2_,"double",dim_2_[0],dim_2_[1]," %% r8_2_: ");
  printf(" %% r8_0_ vs r8_2_: relative error %0.16f\n",dfnormn(n_i,r8_0_,r8_2_));
  /* %%%%%%%% */
  r4_0_ = (float *) malloc(n_i*sizeof(float));
  r4_0_[0] = 0.1; r4_0_[2] = 10.1; r4_0_[4] = 100.1;
  r4_0_[1] = 1.1; r4_0_[3] = 20.1; r4_0_[5] = 200.1;
  array_printf(r4_0_,"float",dim_0_[0],dim_0_[1]," %% r4_0_: ");
  MDA_write_r4(n_dim_0,dim_0_,r4_0_,fname);
  MDA_read_r4(&n_dim_3,&dim_3_,&r4_3_,fname);
  array_printf(r4_3_,"float",dim_3_[0],dim_3_[1]," %% r4_3_: ");
  printf(" %% r4_0_ vs r4_3_: relative error %0.16f\n",ffnormn(n_i,r4_0_,r4_3_));
  /* %%%%%%%% */
  ulli_0_ = (unsigned long long int *) malloc(n_i*sizeof(unsigned long long int));
  ulli_0_[0] = 100l; ulli_0_[2] = 1010l; ulli_0_[4] = 10100l;
  ulli_0_[1] = 101l; ulli_0_[3] = 1020l; ulli_0_[5] = 10200l;
  array_printf(ulli_0_,"unsigned long long int",dim_0_[0],dim_0_[1]," %% ulli_0_: ");
  MDA_write_ulli(n_dim_0,dim_0_,ulli_0_,fname);
  MDA_read_ulli(&n_dim_4,&dim_4_,&ulli_4_,fname);
  array_printf(ulli_4_,"unsigned long long int",dim_4_[0],dim_4_[1]," %% ulli_4_: ");
  printf(" %% ulli_0_ vs ulli_4_: relative error %0.16f\n",ullifnormn(n_i,ulli_0_,ulli_4_));
  /* %%%%%%%% */
  c16_0_ = (double complex *) malloc(n_i*sizeof(double complex));
  c16_0_[0] = (double complex)100.0 + (double complex)10.0*_Complex_I;
  c16_0_[2] = (double complex)1010.0 + (double complex)12.0*_Complex_I;
  c16_0_[4] = (double complex)10100.0 + (double complex)14.0*_Complex_I;
  c16_0_[1] = (double complex)101.0 + (double complex)11.0*_Complex_I;
  c16_0_[3] = (double complex)1020.0 + (double complex)13.0*_Complex_I;
  c16_0_[5] = (double complex)10200.0 + (double complex)15.0*_Complex_I;
  array_printf(c16_0_,"double complex",dim_0_[0],dim_0_[1]," %% c16_0_: ");
  MDA_write_c16(n_dim_0,dim_0_,c16_0_,fname);
  MDA_read_c16(&n_dim_5,&dim_5_,&c16_5_,fname);
  array_printf(c16_5_,"double complex",dim_5_[0],dim_5_[1]," %% c16_5_: ");
  printf(" %% c16_0_ vs c16_5_: relative error %0.16f\n",zfnormn(n_i,c16_0_,c16_5_));
  /* %%%%%%%% */
  free(dim_0_); dim_0_=NULL;
  free(i4_0_); i4_0_=NULL;
  free(r8_0_); r8_0_=NULL;
  free(r4_0_); r4_0_=NULL;
  free(ulli_0_); ulli_0_=NULL;
  free(ulli_4_); ulli_4_=NULL;
  free(dim_5_); dim_5_=NULL;
  free(dim_4_); dim_4_=NULL;
  free(dim_3_); dim_3_=NULL;
  free(dim_2_); dim_2_=NULL;
  free(dim_1_); dim_1_=NULL;
  free(c16_0_); c16_0_=NULL;
  free(c16_5_); c16_5_=NULL;
  free(r4_3_); r4_3_=NULL;
  free(r8_2_); r8_2_=NULL;
  free(i4_1_); i4_1_=NULL;
}
/* %%%%%%%% */

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
/* index of max: */
/* %%%%%%%% */
int dmax_index(int n_a,double *a_)
{
  int na=0;
  int index=0;
  double a=0;
  if (n_a>0){
    a = a_[index];
    for (na=0;na<n_a;na++){ if (a_[na]>a){ a=a_[na]; index=na;}}
    /* if (n_a>0){ } */}
  return index;
}

/* %%%%%%%% */
/* local timing: */
/* %%%%%%%% */
void local_tic(int ntick,clock_t *t_start_,struct timeval *d_start_)
{
  t_start_[ntick] = clock(); gettimeofday(&d_start_[ntick],NULL);
}
void local_toc(int ntick,clock_t *t_start_,clock_t *t_final_,struct timeval *d_start_,struct timeval *d_final_,long *l_msec_,long *l_ssec_,long *l_usec_,double *elct_,double *elrt_,double n_op,int verbose,const char *prefix)
{
  double r=0,s=0;
  t_final_[ntick] = clock(); gettimeofday(&d_final_[ntick],NULL);
  l_ssec_[ntick] =  d_final_[ntick].tv_sec -  d_start_[ntick].tv_sec;
  l_usec_[ntick] = d_final_[ntick].tv_usec - d_start_[ntick].tv_usec;
  l_msec_[ntick] = ((l_ssec_[ntick]*1000) + l_usec_[ntick]/1000.0) + 0.5;
  elct_[ntick] = (double)(1000*(t_final_[ntick]-t_start_[ntick])/CLOCKS_PER_SEC)/(double)1000;
  elrt_[ntick] = (double)l_msec_[ntick]/(double)1000;
  r = elct_[ntick]/elrt_[ntick];
  s = n_op/elrt_[ntick];
  if (verbose>=1){ 
    if (finite(r)){ printf("%sct/rt %0.3f/%0.3f = %.1f <-- %0.2f Mhz %0.2f Ghz\n",prefix,elct_[ntick],elrt_[ntick],r,s/1e6,s/1e9);}
    else{ printf("%sct/rt %0.3f/%0.3f = 0 <-- %0.2f Mhz %0.2f Ghz\n",prefix,elct_[ntick],elrt_[ntick],s/1e6,s/1e9);}
    /* if (verbose>=1){ } */}
}

/* %%%%%%%% */
/* testing transposes found at: */
/* https://stackoverflow.com/questions/16737298/what-is-the-fastest-way-to-transpose-a-matrix-in-c */
/* %%%%%%%% */
inline void transpose_float_rc_bruteforce(float *A__, float *B__, const int n_row, const int n_col)
{
  int nr=0,nc=0,na=0,nb=0;
  for (nr=0;nr<n_row;nr++){
    nb=nr;
    for (nc=0;nc<n_col;nc++){
      B__[na++] = A__[nb]; nb+=n_row;
      /* for (nc=0;nc<n_col;nc++){ } */}
    /* for (nr=0;nr<n_row;nr++){ } */}
}
inline void transpose_float_cr_bruteforce(float *A__, float *B__, const int n_row, const int n_col)
{
  int nr=0,nc=0,na=0,nb=0;
  for (nc=0;nc<n_col;nc++){
    nb=nc;
    for (nr=0;nr<n_row;nr++){
      B__[nb] = A__[na++]; nb+=n_col;
      /* for (nr=0;nr<n_row;nr++){ } */}
    /* for (nc=0;nc<n_col;nc++){ } */}
}
inline void transpose_float_block0(float *A__, float *B__, const int lda, const int ldb, const int block_size)
{
  int i=0,j=0;
  for(i=0; i<block_size; i++) {
    for(j=0; j<block_size; j++) {
      B__[j*ldb + i] = A__[i*lda +j];
      /* for(j=0; j<block_size; j++) { } */}
    /* for(i=0; i<block_size; i++) { } */}
}
inline void transpose_float_block1(float *A__, float *B__, const int n, const int m, const int lda, const int ldb, const int block_size)
{
  int i=0,j=0;
  for(i=0; i<n; i+=block_size) {
    for(j=0; j<m; j+=block_size) {
      transpose_float_block0(&A__[i*lda +j], &B__[j*ldb + i], lda, ldb, block_size);
      /* for(j=0; j<m; j+=block_size) { } */}
    /* for(i=0; i<n; i+=block_size) { } */}
}
inline void transpose_ps_block0(float *A, float *B, const int lda, const int ldb)
{
  __m128 row1 = _mm_load_ps(&A[0*lda]);
  __m128 row2 = _mm_load_ps(&A[1*lda]);
  __m128 row3 = _mm_load_ps(&A[2*lda]);
  __m128 row4 = _mm_load_ps(&A[3*lda]);
  _MM_TRANSPOSE4_PS(row1, row2, row3, row4);
  _mm_store_ps(&B[0*ldb], row1);
  _mm_store_ps(&B[1*ldb], row2);
  _mm_store_ps(&B[2*ldb], row3);
  _mm_store_ps(&B[3*ldb], row4);
}
inline void transpose_ps_block1(float *A, float *B, const int n, const int m, const int lda, const int ldb ,const int block_size)
{
  int i=0,j=0,max_i2=0,max_j2=0,i2=0,j2=0;
  for(i=0; i<n; i+=block_size) {
    for(j=0; j<m; j+=block_size) {
      max_i2 = i+block_size < n ? i + block_size : n;
      max_j2 = j+block_size < m ? j + block_size : m;
      for(i2=i; i2<max_i2; i2+=4) {
	for(j2=j; j2<max_j2; j2+=4) {
	  transpose_ps_block0(&A[i2*lda +j2], &B[j2*ldb + i2], lda, ldb);
	  /* for(j2=j; j2<max_j2; j2+=4) { } */}
	/* for(i2=i; i2<max_i2; i2+=4) { } */}
      /* for(j=0; j<m; j+=block_size) { } */}
    /* for(i=0; i<n; i+=block_size) { } */}
}
void test_transpose_float()
{
  int verbose=1;
  const int n_row = 17001;
  const int n_col = 15001;
  int n_row_rup = rup(n_row, 16);
  int n_col_rup = rup(n_col, 16);
  unsigned long long int ulli_total = (unsigned long long int)n_row_rup*(unsigned long long int)n_col_rup;
  int n_row_sub = 5;
  int n_col_sub = 6;
  unsigned long long int tab=0;
  float *A__ = NULL;
  float *B_bf__ = NULL;
  float *B_al__ = NULL;
  /* %%%% */
  int nr=0,nc=0;
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  /* %%%% */
  if (verbose){ printf(" %% [entering test_transpose_float]\n");}
  /* %%%% */
  A__ = (float*)_mm_malloc(sizeof(float)*ulli_total, 64);
  for (nr=0;nr<n_row;nr++){
    for (nc=0;nc<n_col;nc++){
      A__[nr + nc*n_row_rup] = nr*100 + nc;
      /* for (nc=0;nc<n_col;nc++){ } */}
    /* for (nr=0;nr<n_row;nr++){ } */}
  if (verbose>1){ array_sub_printf(A__,"float",n_row_rup,n_row_sub,n_col_rup,n_col_sub," %% A_sub__: ");}
  B_bf__ = (float*)_mm_malloc(sizeof(float)*ulli_total, 64);
  B_al__ = (float*)_mm_malloc(sizeof(float)*ulli_total, 64);
  /* %%%% */
  memset(B_bf__,0,ulli_total*sizeof(float));
  local_tic(0,t_start_,d_start_);
  tab = ulli_total;
  transpose_float_rc_bruteforce(A__,B_bf__,n_row_rup,n_col_rup);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"transpose_float_rc_bruteforce: ");
  if (verbose>1){ array_sub_printf(B_bf__,"float",n_col_rup,n_col_sub,n_row_rup,n_row_sub," %% B_sub_bf__: ");}
  /* %%%% */
  memset(B_al__,0,ulli_total*sizeof(float));
  local_tic(0,t_start_,d_start_);
  tab = ulli_total;
  transpose_float_cr_bruteforce(A__,B_al__,n_row_rup,n_col_rup);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"transpose_float_cr_bruteforce: ");
  if (verbose){ printf(" %% error: %0.16f\n",ffnormn(ulli_total,B_bf__,B_al__));}
  /* %%%% */
  memset(B_al__,0,ulli_total*sizeof(float));
  local_tic(0,t_start_,d_start_);
  tab = ulli_total;
  transpose_float_block1(A__,B_al__,n_col_rup,n_row_rup,n_row_rup,n_col_rup,4);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"transpose_float_block1(4): ");
  if (verbose){ printf(" %% error: %0.16f\n",ffnormn(ulli_total,B_bf__,B_al__));}
  /* %%%% */
  memset(B_al__,0,ulli_total*sizeof(float));
  local_tic(0,t_start_,d_start_);
  tab = ulli_total;
  transpose_float_block1(A__,B_al__,n_col_rup,n_row_rup,n_row_rup,n_col_rup,8);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"transpose_float_block1(8): ");
  if (verbose){ printf(" %% error: %0.16f\n",ffnormn(ulli_total,B_bf__,B_al__));}
  /* %%%% */
  memset(B_al__,0,ulli_total*sizeof(float));
  local_tic(0,t_start_,d_start_);
  tab = ulli_total;
  transpose_float_block1(A__,B_al__,n_col_rup,n_row_rup,n_row_rup,n_col_rup,16);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"transpose_float_block1(16): ");
  if (verbose){ printf(" %% error: %0.16f\n",ffnormn(ulli_total,B_bf__,B_al__));}
  /* %%%% */
  memset(B_al__,0,ulli_total*sizeof(float));
  local_tic(0,t_start_,d_start_);
  tab = ulli_total;
  transpose_ps_block1(A__,B_al__,n_col_rup,n_row_rup,n_row_rup,n_col_rup,16);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"transpose_ps_block1(16): ");
  if (verbose){ printf(" %% error: %0.16f\n",ffnormn(ulli_total,B_bf__,B_al__));}
  /* %%%% */
  _mm_free(A__); A__=NULL;
  _mm_free(B_bf__); B_bf__=NULL;
  _mm_free(B_al__); B_al__=NULL;
  /* %%%% */
  if (verbose){ printf(" %% [finished test_transpose_float]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* testing matrix multiplication */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void dp_ps_bruteforce(int n_col_X,float *f_A_,float *f_B_,float *f_C_)
{
  int ncol_X=0;
  float *f_A_point0_=NULL,*f_B_point0_=NULL;
  float f_ac0=0;
  f_A_point0_ = f_A_;
  f_B_point0_ = f_B_;
  for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
    f_ac0 += *(f_A_point0_++) * *(f_B_point0_++);
    /* for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ } */}
  *f_C_ = f_ac0;
}
void dp_ps_mult_bruteforce(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
	  f_C__[nrow_A + nrow_B*n_row_A] += f_A_trn__[ncol_X + nrow_A*n_col_X_rup] * f_B_trn__[ncol_X + nrow_B*n_col_X_rup];
	  /* for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ } */}
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}
void dp_ps_mult_cblas_sgemm(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
#ifdef _CBLAS
    cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,n_row_A,n_row_B,n_col_X,1.0,f_A_trn__,n_col_X_rup,f_B_trn__,n_col_X_rup,0.0,f_C__,n_row_A);
#endif /* _CBLAS */
#ifndef _CBLAS
    printf(" %% Warning, _CBLAS not defined in dp_ps_mult_cblas_sgemm, using bruteforce instead.\n");
    dp_ps_mult_bruteforce(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,f_C_p_);
#endif /* _CBLAS */
  /* if (f_C__!=NULL){ } */}
}
void dp_ps_mult_cblas_sdot(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
#ifdef _CBLAS
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	f_C__[nrow_A + nrow_B*n_row_A] = cblas_sdot(n_col_X,f_A_trn__+nrow_A*n_col_X_rup,1,f_B_trn__+nrow_B*n_col_X_rup,1);
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
#endif /* _CBLAS */
#ifndef _CBLAS
    printf(" %% Warning, _CBLAS not defined in dp_ps_mult_cblas_sdot, using bruteforce instead.\n");
    dp_ps_mult_bruteforce(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,f_C_p_);
#endif /* _CBLAS */
    /* if (f_C__!=NULL){ } */}
}
void dp_ps_immintrin_loadu_avx(int n_col_X,float *f_A_,float *f_B_,float *f_C_)
{
  /* does not assume alignment */
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int ncol_X=0,ncol_X_256=0;
  float *f_A_point0_=NULL,*f_B_point0_=NULL;
  __m256 ps_A0,ps_B0;
  __m256 ps_ac0;
  float *f_ac0_ = (float *) &ps_ac0;
  float f_ac0=0;
  f_A_point0_ = f_A_;
  f_B_point0_ = f_B_;
  ps_ac0  = _mm256_set1_ps((float)0.0);
  for (ncol_X_256=0;ncol_X_256<n_col_X_256-1;ncol_X_256++){
    ps_A0 = _mm256_loadu_ps(f_A_point0_);
    ps_B0 = _mm256_loadu_ps(f_B_point0_);
    ps_ac0 = _mm256_dp_ps(ps_A0,ps_B0,0xF1);
    f_ac0 += f_ac0_[0] + f_ac0_[4];
    f_A_point0_+=8;f_B_point0_+=8;
    /* for (ncol_X_256=0;ncol_X_256<n_col_X_256-1;ncol_X_256++){ } */}
  f_A_point0_ = &(f_A_[8*ncol_X_256]);
  f_B_point0_ = &(f_B_[8*ncol_X_256]);
  for (ncol_X=8*ncol_X_256;ncol_X<n_col_X;ncol_X++){
    f_ac0 += *(f_A_point0_++) * *(f_B_point0_++);
    /* for (ncol_X=8*ncol_X_256;ncol_X<n_col_X;ncol_X++){ } */}
  *f_C_ = f_ac0;
}
void dp_ps_mult_immintrin_loadu_avx(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  /* does not assume alignment */
  int nrow_A=0,nrow_B=0;
  float *tmp_f_A_point0_=NULL;
  float *tmp_f_B_point0_=NULL;
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_A_point0_ = &(f_A_trn__[nrow_A*n_col_X]);
	tmp_f_B_point0_ = &(f_B_trn__[nrow_B*n_col_X]);
	dp_ps_immintrin_loadu_avx(n_col_X,tmp_f_A_point0_,tmp_f_B_point0_,&(f_C__[nrow_A + nrow_B*n_row_A]));
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}
void dp_ps_immintrin_loadu_fma(int n_col_X,float *f_A_,float *f_B_,float *f_C_)
{
  /* does not assume alignment */
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int ncol_X=0,ncol_X_256=0;
  float *f_A_point0_=NULL,*f_B_point0_=NULL;
  __m256 ps_A0,ps_B0;
  __m256 ps_ac0;
  float *f_ac0_ = (float *) &ps_ac0;
  float f_ac0=0;
  f_A_point0_ = f_A_;
  f_B_point0_ = f_B_;
  ps_ac0  = _mm256_set1_ps((float)0.0);
  for (ncol_X_256=0;ncol_X_256<n_col_X_256-1;ncol_X_256++){
    ps_A0 = _mm256_loadu_ps(f_A_point0_);
    ps_B0 = _mm256_loadu_ps(f_B_point0_);
#ifdef _FMA
    ps_ac0 = _mm256_fmadd_ps(ps_A0,ps_B0,ps_ac0);
#endif /* _FMA */
    f_A_point0_+=8;f_B_point0_+=8;
    /* for (ncol_X_256=0;ncol_X_256<n_col_X_256-1;ncol_X_256++){ } */}
  f_ac0 = f_ac0_[0] + f_ac0_[1] + f_ac0_[2] + f_ac0_[3] + f_ac0_[4] + f_ac0_[5] + f_ac0_[6] + f_ac0_[7];
  f_A_point0_ = &(f_A_[8*ncol_X_256]);
  f_B_point0_ = &(f_B_[8*ncol_X_256]);
  for (ncol_X=8*ncol_X_256;ncol_X<n_col_X;ncol_X++){
    f_ac0 += *(f_A_point0_++) * *(f_B_point0_++);
    /* for (ncol_X=8*ncol_X_256;ncol_X<n_col_X;ncol_X++){ } */}
  *f_C_ = f_ac0;
}
void dp_ps_mult_immintrin_loadu_fma(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  /* does not assume alignment */
  int nrow_A=0,nrow_B=0;
  float *tmp_f_A_point0_=NULL;
  float *tmp_f_B_point0_=NULL;
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_A_point0_ = &(f_A_trn__[nrow_A*n_col_X]);
	tmp_f_B_point0_ = &(f_B_trn__[nrow_B*n_col_X]);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_A_point0_,tmp_f_B_point0_,&(f_C__[nrow_A + nrow_B*n_row_A]));
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}
void dp_ps_mult_immintrin_load1_fma(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float *tmp_f_A_point0_=NULL;
  float *tmp_f_B_point0_=NULL;
  __m256 ps_A0,ps_B0;
  __m256 ps_ac0;
  float *tmp_f_ac0_ = (float *) &ps_ac0;
  float tmp_f_ac0=0;
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_A_point0_ = &(f_A_trn__[nrow_A*n_col_X_rup]);
	tmp_f_B_point0_ = &(f_B_trn__[nrow_B*n_col_X_rup]);
	ps_ac0 = _mm256_set1_ps((float)0.0);
	for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){
	  ps_A0 = _mm256_load_ps(tmp_f_A_point0_);
	  ps_B0 = _mm256_load_ps(tmp_f_B_point0_);
#ifdef _FMA
	  ps_ac0 = _mm256_fmadd_ps(ps_A0,ps_B0,ps_ac0);
#endif /* _FMA */
#ifdef _FMA
	  tmp_f_A_point0_+=8;tmp_f_B_point0_+=8;
#endif /* _FMA */
	  /* for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){ } */}
	tmp_f_ac0 = tmp_f_ac0_[0] + tmp_f_ac0_[1] + tmp_f_ac0_[2] + tmp_f_ac0_[3] + tmp_f_ac0_[4] + tmp_f_ac0_[5] + tmp_f_ac0_[6] + tmp_f_ac0_[7];
	f_C__[nrow_A + nrow_B*n_row_A] = tmp_f_ac0;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}
void dp_ps_mult_immintrin_fma2(int n_row_A,int n_col_X,__m256 *ps_A_trn__,int n_row_B,__m256 *ps_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int nrow_A=0,nrow_B=0,ncol_X=0;
  __m256 *tmp_ps_A_trn__=NULL;
  __m256 *tmp_ps_B_trn__=NULL;
  __m256 *tmp_ps_A_trn_=NULL;
  __m256 *tmp_ps_B_trn_=NULL;
  __m256 tmp_ps_A_trn,tmp_ps_B_trn;
  __m256 ps_ac0,ps_ac1;
  float *tmp_f_ac0_ = (float *) &ps_ac0;
  float *tmp_f_ac1_ = (float *) &ps_ac1;
  float tmp_f_ac0=0,tmp_f_ac1=0; 
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    printf(" n_col_X %d n_col_X_rup %d n_col_X_256 %d \n",n_col_X,n_col_X_rup,n_col_X_256);
    tmp_ps_B_trn__ = ps_B_trn__;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      tmp_ps_A_trn__ = ps_A_trn__;
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_ps_A_trn_ = tmp_ps_A_trn__;
	tmp_ps_B_trn_ = tmp_ps_B_trn__;
	ps_ac0 = _mm256_set1_ps((float)0.0);
	ps_ac1 = _mm256_set1_ps((float)0.0);
	for (ncol_X=0;ncol_X<n_col_X_256-1;ncol_X+=2){
#ifdef _FMA
	  ps_ac0 = _mm256_fmadd_ps(*(tmp_ps_A_trn_),*(tmp_ps_B_trn_),ps_ac0);
#endif /* _FMA */
#ifdef _FMA
	  ps_ac1 = _mm256_fmadd_ps(*(tmp_ps_A_trn_+1),*(tmp_ps_B_trn_+1),ps_ac1);
#endif /* _FMA */
	  tmp_ps_A_trn_+=2;
	  tmp_ps_B_trn_+=2;
	  /* for (ncol_X=0;ncol_X<n_col_X_256-1;ncol_X+=2){ } */}
	if ( (nrow_A==0) && (nrow_B==0) ){ printf(" ncol_X %d\n",ncol_X);}
#ifdef _FMA
	if (n_col_X_256%2==1){ ps_ac0 = _mm256_fmadd_ps(*(tmp_ps_A_trn_++),*(tmp_ps_B_trn_++),ps_ac0);}
#endif /* _FMA */
	tmp_f_ac0 = tmp_f_ac0_[0] + tmp_f_ac0_[1] + tmp_f_ac0_[2] + tmp_f_ac0_[3] + tmp_f_ac0_[4] + tmp_f_ac0_[5] + tmp_f_ac0_[6] + tmp_f_ac0_[7];
	tmp_f_ac1 = tmp_f_ac1_[0] + tmp_f_ac1_[1] + tmp_f_ac1_[2] + tmp_f_ac1_[3] + tmp_f_ac1_[4] + tmp_f_ac1_[5] + tmp_f_ac1_[6] + tmp_f_ac1_[7];
	f_C__[nrow_A + nrow_B*n_row_A] = tmp_f_ac0 + tmp_f_ac1;
	tmp_ps_A_trn__ += n_col_X_256;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      tmp_ps_B_trn__ += n_col_X_256;
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}
void dp_ps_mult_immintrin_fma(int n_row_A,int n_col_X,__m256 *ps_A_trn__,int n_row_B,__m256 *ps_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int nrow_A=0,nrow_B=0,ncol_X=0;
  __m256 *tmp_ps_A_trn__=NULL;
  __m256 *tmp_ps_B_trn__=NULL;
  __m256 *tmp_ps_A_trn_=NULL;
  __m256 *tmp_ps_B_trn_=NULL;
  __m256 tmp_ps_A_trn,tmp_ps_B_trn;
  __m256 ps_acc;
  float *tmp_f_ = (float *) &ps_acc;
  float tmp_f=0; 
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    tmp_ps_B_trn__ = ps_B_trn__;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      tmp_ps_A_trn__ = ps_A_trn__;
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	ps_acc = _mm256_set1_ps((float)0.0);
	tmp_ps_A_trn_ = tmp_ps_A_trn__;
	tmp_ps_B_trn_ = tmp_ps_B_trn__;
	for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){
#ifdef _FMA
	  ps_acc = _mm256_fmadd_ps(*(tmp_ps_A_trn_++),*(tmp_ps_B_trn_++),ps_acc);
#endif /* _FMA */
	  /* for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){ } */}
	tmp_f = tmp_f_[0] + tmp_f_[1] + tmp_f_[2] + tmp_f_[3] + tmp_f_[4] + tmp_f_[5] + tmp_f_[6] + tmp_f_[7];
	f_C__[nrow_A + nrow_B*n_row_A] = tmp_f;
	tmp_ps_A_trn__ += n_col_X_256;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      tmp_ps_B_trn__ += n_col_X_256;
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}
void dp_ps_mult_immintrin_avx(int n_row_A,int n_col_X,__m256 *ps_A_trn__,int n_row_B,__m256 *ps_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int nrow_A=0,nrow_B=0,ncol_X=0;
  __m256 *tmp_ps_A_trn__=NULL;
  __m256 *tmp_ps_B_trn__=NULL;
  __m256 *tmp_ps_A_trn_=NULL;
  __m256 *tmp_ps_B_trn_=NULL;
  __m256 tmp_ps,tmp_ps_A_trn,tmp_ps_B_trn;
  float *tmp_f_ = (float *) &tmp_ps;
  float tmp_f=0; 
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    tmp_ps_B_trn__ = ps_B_trn__;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      tmp_ps_A_trn__ = ps_A_trn__;
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f=0;
	tmp_ps_A_trn_ = tmp_ps_A_trn__;
	tmp_ps_B_trn_ = tmp_ps_B_trn__;
	for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){
	  tmp_ps = _mm256_dp_ps(*(tmp_ps_A_trn_++),*(tmp_ps_B_trn_++),0xF1);
	  tmp_f += tmp_f_[0] + tmp_f_[4];
	  /* for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){ } */}
	f_C__[nrow_A + nrow_B*n_row_A] = tmp_f;
	tmp_ps_A_trn__ += n_col_X_256;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      tmp_ps_B_trn__ += n_col_X_256;
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}
void dp_ps_immintrin_loadu_wrap(int n_col_X,float *f_A_,float *f_B_,float *f_C_)
{
  /* does not assume alignment */
#ifdef _FMA
  dp_ps_immintrin_loadu_fma(n_col_X,f_A_,f_B_,f_C_);
#endif /* _FMA */
#ifndef _FMA
  dp_ps_immintrin_loadu_avx(n_col_X,f_A_,f_B_,f_C_);
#endif /* _FMA */
}
void dp_ps_mult_immintrin_loadu_wrap(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  /* does not assume alignment */
#ifdef _FMA
  dp_ps_mult_immintrin_loadu_fma(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,f_C_p_);
#endif /* _FMA */
#ifndef _FMA
  dp_ps_mult_immintrin_loadu_avx(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,f_C_p_);
#endif /* _FMA */
}
void dp_ps_mult_immintrin_test()
{
  int verbose=1;
  /* int n_row_A = 1024*1 + 723; */
  /* int n_col_X = 1024*6 + 817; */
  /* int n_row_B = 1024*1 + 511; */
  int n_row_A = 531;
  int n_col_X = 15039;
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int n_row_B = 430;
  int n_row_A_sub = 5;
  int n_col_X_sub = 15;
  int n_row_B_sub = 4;
  int nrow_A=0,ncol_X=0,nrow_B=0;
  unsigned long long int ulli_A_total = (unsigned long long int) n_row_A * (unsigned long long int) n_col_X_rup ;
  unsigned long long int ulli_B_total = (unsigned long long int) n_row_B * (unsigned long long int) n_col_X_rup ;
  unsigned long long int ulli_C_total = (unsigned long long int) n_row_A * (unsigned long long int) n_row_B ;
  unsigned long long int ulli=0;
  unsigned long long int tab=0;
  __m256 *ps_A_trn__ = NULL;
  __m256 *ps_B_trn__ = NULL;
  float *f_A_trn__ = NULL;
  float *f_A_u_trn__ = NULL;
  float *f_B_trn__ = NULL;
  float *f_B_u_trn__ = NULL;
  float *f_C_bf__ = NULL;
  float *f_C_ps__ = NULL;
  float *f_A_sub__ = NULL;
  float *f_B_sub__ = NULL;
  float *f_C_sub__ = NULL;
  float ferror=0;
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  /* %%%% */
  ps_A_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  ps_B_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  f_A_trn__ = (float *) ps_A_trn__;
  f_B_trn__ = (float *) ps_B_trn__;
  f_A_u_trn__ = (float *) malloc(n_col_X*n_row_A*sizeof(float));
  f_B_u_trn__ = (float *) malloc(n_col_X*n_row_A*sizeof(float));
  f_C_bf__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_C_ps__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_A_sub__ = (float *) malloc(n_row_A_sub*n_col_X_sub*sizeof(float));
  f_B_sub__ = (float *) malloc(n_row_B_sub*n_col_X_sub*sizeof(float));
  f_C_sub__ = (float *) malloc(n_row_A_sub*n_row_B_sub*sizeof(float));
  tab = (unsigned long long int) n_row_A * (unsigned long long int) n_row_B * (unsigned long long int) n_col_X_rup ;
  local_tic(0,t_start_,d_start_);
  ulli=0;
  for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      f_A_trn__[ncol_X + nrow_A*n_col_X_rup] = (float)(((int)ulli%7)-3); 
      f_A_u_trn__[ncol_X + nrow_A*n_col_X] = (float)(((int)ulli%7)-3); 
      ulli++;
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ f_A_sub__[nrow_A+ncol_X*n_row_A_sub] = f_A_trn__[ncol_X + nrow_A*n_col_X_rup];}}
  if (verbose>1){
    printf(" %% upper corner of f_A_trn__: \n");
    array_printf(f_A_sub__,"float",n_row_A_sub,n_col_X_sub," % f_A_sub__: ");
    /* if (verbose>1){ } */}
  ulli=0;
  for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      f_B_trn__[ncol_X + nrow_B*n_col_X_rup] = (float)(((int)ulli%5)-2);
      f_B_u_trn__[ncol_X + nrow_B*n_col_X] = (float)(((int)ulli%5)-2);
      ulli++;
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}
  for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ f_B_sub__[nrow_B+ncol_X*n_row_B_sub] = f_B_trn__[ncol_X + nrow_B*n_col_X_rup];}}
  if (verbose>1){
    printf(" %% upper corner of f_B_trn__: \n");
    array_printf(f_B_sub__,"float",n_row_B_sub,n_col_X_sub," % f_B_sub__: ");
    /* if (verbose>1){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," initialize: ");
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_bf__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_bruteforce(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,&f_C_bf__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_bruteforce: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_bf__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_bf__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_bf__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_bf__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
#ifdef _CBLAS
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_cblas_sgemm(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_cblas_sgemm: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
#endif /* _CBLAS */
  /* %%%%%%%% */
#ifdef _CBLAS
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_cblas_sdot(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_cblas_sdot: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
#endif /* _CBLAS */
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_avx(n_row_A,n_col_X,ps_A_trn__,n_row_B,ps_B_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_immintrin_avx: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_fma(n_row_A,n_col_X,ps_A_trn__,n_row_B,ps_B_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_immintrin_fma: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_fma2(n_row_A,n_col_X,ps_A_trn__,n_row_B,ps_B_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_immintrin_fma2: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_load1_fma(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_immintrin_load1_fma: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_loadu_fma(n_row_A,n_col_X,f_A_u_trn__,n_row_B,f_B_u_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_immintrin_loadu_fma: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_loadu_avx(n_row_A,n_col_X,f_A_u_trn__,n_row_B,f_B_u_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_immintrin_loadu_avx: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_loadu_wrap(n_row_A,n_col_X,f_A_u_trn__,n_row_B,f_B_u_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_immintrin_loadu_wrap: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  _mm_free(ps_A_trn__); ps_A_trn__ = NULL;
  _mm_free(ps_B_trn__); ps_B_trn__ = NULL;
  free(f_A_u_trn__); f_A_u_trn__=NULL;
  free(f_B_u_trn__); f_B_u_trn__=NULL;
  free(f_C_bf__); f_C_bf__=NULL;
  free(f_C_ps__); f_C_ps__=NULL;
  free(f_A_sub__); f_A_sub__=NULL;
  free(f_B_sub__); f_B_sub__=NULL;
  free(f_C_sub__); f_C_sub__=NULL;
  //wkspace_printf();
}
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void hp_interleave_mult_bruteforce(int n_row_A,int n_col_X,float complex *c_A_trn__,int n_row_B,float complex *c_B_trn__,float complex **c_C_p_)
{
  /* unaligned */
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float complex *c_C__=NULL;
  c_C__=NULL;
  if (c_C_p_!=NULL){
    if ( (*c_C_p_)==NULL ){ (*c_C_p_) = (float complex *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float complex));}
    c_C__ = *c_C_p_;
    /* if (c_C_p_!=NULL){ } */}
  if (c_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
	  c_C__[nrow_A + nrow_B*n_row_A] += conjf(c_A_trn__[ncol_X + nrow_A*n_col_X]) * c_B_trn__[ncol_X + nrow_B*n_col_X];
	  /* for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ } */}
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (c_C__!=NULL){ } */}
}
void hp_interleave_mult_cblas_cgemm(int n_row_A,int n_col_X,float complex *c_A_trn__,int n_row_B,float complex *c_B_trn__,float complex **c_C_p_)
{
  /* unaligned */
  float complex *c_C__=NULL;
  float complex calpha = (float complex)1.0;
  float complex cbeta = (float complex)0.0;
  c_C__=NULL;
  if (c_C_p_!=NULL){
    if ( (*c_C_p_)==NULL ){ (*c_C_p_) = (float complex *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float complex));}
    c_C__ = *c_C_p_;
    /* if (c_C_p_!=NULL){ } */}
  if (c_C__!=NULL){
#ifdef _CBLAS
    cblas_cgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,n_row_A,n_row_B,n_col_X,&calpha,c_A_trn__,n_col_X,c_B_trn__,n_col_X,&cbeta,c_C__,n_row_A);
#endif /* _CBLAS */
#ifndef _CBLAS
    printf(" %% Warning, _CBLAS not defined in hp_interleave_mult_cblas_cgemm, using bruteforce instead.\n");
    hp_interleave_mult_bruteforce(n_row_A,n_col_X,c_A_trn__,n_row_B,c_B_trn__,c_C_p_);
#endif /* _CBLAS */
  /* if (c_C__!=NULL){ } */}
}
void hp_segregated_mult_bruteforce(int n_row_A,int n_col_X,float *f_AR_trn__,float *f_AI_trn__,int n_row_B,float *f_BR_trn__,float *f_BI_trn__,float complex **c_C_p_)
{
  /* unaligned */
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float complex *c_C__=NULL;
  c_C__=NULL;
  if (c_C_p_!=NULL){
    if ( (*c_C_p_)==NULL ){ (*c_C_p_) = (float complex *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float complex));}
    c_C__ = *c_C_p_;
    /* if (c_C_p_!=NULL){ } */}
  if (c_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	c_C__[nrow_A + nrow_B*n_row_A] = (float complex) 0.0;
	for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
	  c_C__[nrow_A + nrow_B*n_row_A] +=
	    ((float complex)f_AR_trn__[ncol_X + nrow_A*n_col_X] - _Complex_I * (float complex)f_AI_trn__[ncol_X + nrow_A*n_col_X]) *
	    ((float complex)f_BR_trn__[ncol_X + nrow_B*n_col_X] + _Complex_I * (float complex)f_BI_trn__[ncol_X + nrow_B*n_col_X])
	    ;
	  /* for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ } */}
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (c_C__!=NULL){ } */}
}
void hp_segregated_to_interleaved_mult_immintrin_loadu_fma(int n_row_A,int n_col_X,float *f_AR_trn__,float *f_AI_trn__,int n_row_B,float *f_BR_trn__,float *f_BI_trn__,float complex **c_C_p_)
{
  /* does not assume alignment */
  int nrow_A=0,nrow_B=0;
  float *tmp_f_AR_point0_=NULL;
  float *tmp_f_AI_point0_=NULL;
  float *tmp_f_BR_point0_=NULL;
  float *tmp_f_BI_point0_=NULL;
  float complex *c_C__=NULL;
  float f_ARBR=0.0,f_ARBI=0.0,f_AIBR=0.0,f_AIBI=0.0;
  int na=0;
  c_C__=NULL;
  if (c_C_p_!=NULL){
    if ( (*c_C_p_)==NULL ){ (*c_C_p_) = (float complex *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float complex));}
    c_C__ = *c_C_p_;
    /* if (c_C_p_!=NULL){ } */}
  if (c_C__!=NULL){
    na=0;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_AR_point0_ = &(f_AR_trn__[nrow_A*n_col_X]);
	tmp_f_AI_point0_ = &(f_AI_trn__[nrow_A*n_col_X]);
	tmp_f_BR_point0_ = &(f_BR_trn__[nrow_B*n_col_X]);
	tmp_f_BI_point0_ = &(f_BI_trn__[nrow_B*n_col_X]);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AR_point0_,tmp_f_BR_point0_,&f_ARBR);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AR_point0_,tmp_f_BI_point0_,&f_ARBI);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AI_point0_,tmp_f_BI_point0_,&f_AIBI);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AI_point0_,tmp_f_BR_point0_,&f_AIBR);
	c_C__[na] = (float complex)(f_ARBR + f_AIBI) + _Complex_I * (float complex)(f_ARBI - f_AIBR);
	na += 1;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (c_C__!=NULL){ } */}
}
void hp_segregated_to_segregated_mult_immintrin_loadu_fma(int n_row_A,int n_col_X,float *f_AR_trn__,float *f_AI_trn__,int n_row_B,float *f_BR_trn__,float *f_BI_trn__,float **f_CR_p_,float **f_CI_p_)
{
  /* does not assume alignment */
  int nrow_A=0,nrow_B=0;
  float *tmp_f_AR_point0_=NULL;
  float *tmp_f_AI_point0_=NULL;
  float *tmp_f_BR_point0_=NULL;
  float *tmp_f_BI_point0_=NULL;
  float *f_CR__=NULL;
  float *f_CI__=NULL;
  float f_ARBR=0.0,f_ARBI=0.0,f_AIBR=0.0,f_AIBI=0.0;
  int na=0;
  f_CR__=NULL;
  if (f_CR_p_!=NULL){
    if ( (*f_CR_p_)==NULL ){ (*f_CR_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_CR__ = *f_CR_p_;
    /* if (f_CR_p_!=NULL){ } */}
  f_CI__=NULL;
  if (f_CI_p_!=NULL){
    if ( (*f_CI_p_)==NULL ){ (*f_CI_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_CI__ = *f_CI_p_;
    /* if (f_CI_p_!=NULL){ } */}
  if ((f_CR__!=NULL) && (f_CI__!=NULL)){
    na=0;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_AR_point0_ = &(f_AR_trn__[nrow_A*n_col_X]);
	tmp_f_AI_point0_ = &(f_AI_trn__[nrow_A*n_col_X]);
	tmp_f_BR_point0_ = &(f_BR_trn__[nrow_B*n_col_X]);
	tmp_f_BI_point0_ = &(f_BI_trn__[nrow_B*n_col_X]);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AR_point0_,tmp_f_BR_point0_,&f_ARBR);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AR_point0_,tmp_f_BI_point0_,&f_ARBI);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AI_point0_,tmp_f_BI_point0_,&f_AIBI);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AI_point0_,tmp_f_BR_point0_,&f_AIBR);
	f_CR__[na] = (f_ARBR + f_AIBI);
	f_CI__[na] = (f_ARBI - f_AIBR);
	na += 1;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if ((f_CR__!=NULL) && (f_CI__!=NULL)){ } */}
}
void hp_segregated_to_segregated_mult_immintrin_load1_fma(int n_row_A,int n_col_X,float *f_AR_trn__,float *f_AI_trn__,int n_row_B,float *f_BR_trn__,float *f_BI_trn__,float **f_CR_p_,float **f_CI_p_)
{
  /* assumes alignment */
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float *tmp_f_AR_point0_=NULL;
  float *tmp_f_AI_point0_=NULL;
  float *tmp_f_BR_point0_=NULL;
  float *tmp_f_BI_point0_=NULL;
  __m256 ps_AR0,ps_BR0;
  __m256 ps_AI0,ps_BI0;
  __m256 ps_acARBR0;
  __m256 ps_acARBI0;
  __m256 ps_acAIBR0;
  __m256 ps_acAIBI0;
  float *tmp_f_acARBR0_ = (float *) &ps_acARBR0;
  float *tmp_f_acARBI0_ = (float *) &ps_acARBI0;
  float *tmp_f_acAIBR0_ = (float *) &ps_acAIBR0;
  float *tmp_f_acAIBI0_ = (float *) &ps_acAIBI0;
  float tmp_f_acARBR0=0;
  float tmp_f_acARBI0=0;
  float tmp_f_acAIBR0=0;
  float tmp_f_acAIBI0=0;
  float *f_CR__=NULL;
  float *f_CI__=NULL;
  int na=0;
  f_CR__=NULL;
  if (f_CR_p_!=NULL){
    if ( (*f_CR_p_)==NULL ){ (*f_CR_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_CR__ = *f_CR_p_;
    /* if (f_CR_p_!=NULL){ } */}
  f_CI__=NULL;
  if (f_CI_p_!=NULL){
    if ( (*f_CI_p_)==NULL ){ (*f_CI_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_CI__ = *f_CI_p_;
    /* if (f_CI_p_!=NULL){ } */}
  if ((f_CR__!=NULL) && (f_CI__!=NULL)){
    na=0;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_AR_point0_ = &(f_AR_trn__[nrow_A*n_col_X_rup]);
	tmp_f_AI_point0_ = &(f_AI_trn__[nrow_A*n_col_X_rup]);
	tmp_f_BR_point0_ = &(f_BR_trn__[nrow_B*n_col_X_rup]);
	tmp_f_BI_point0_ = &(f_BI_trn__[nrow_B*n_col_X_rup]);
	ps_acARBR0 = _mm256_set1_ps((float)0.0);
	ps_acARBI0 = _mm256_set1_ps((float)0.0);
	ps_acAIBR0 = _mm256_set1_ps((float)0.0);
	ps_acAIBI0 = _mm256_set1_ps((float)0.0);
	for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){
	  ps_AR0 = _mm256_load_ps(tmp_f_AR_point0_);
	  ps_BR0 = _mm256_load_ps(tmp_f_BR_point0_);
	  ps_AI0 = _mm256_load_ps(tmp_f_AI_point0_);
	  ps_BI0 = _mm256_load_ps(tmp_f_BI_point0_);
#ifdef _FMA
	  ps_acARBR0 = _mm256_fmadd_ps(ps_AR0,ps_BR0,ps_acARBR0);
	  ps_acARBI0 = _mm256_fmadd_ps(ps_AR0,ps_BI0,ps_acARBI0);
	  ps_acAIBR0 = _mm256_fmadd_ps(ps_AI0,ps_BR0,ps_acAIBR0);
	  ps_acAIBI0 = _mm256_fmadd_ps(ps_AI0,ps_BI0,ps_acAIBI0);
#endif /* _FMA */
#ifdef _FMA
	  tmp_f_AR_point0_+=8;tmp_f_BR_point0_+=8;
	  tmp_f_AI_point0_+=8;tmp_f_BI_point0_+=8;
#endif /* _FMA */
	  /* for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){ } */}
	tmp_f_acARBR0 = tmp_f_acARBR0_[0] + tmp_f_acARBR0_[1] + tmp_f_acARBR0_[2] + tmp_f_acARBR0_[3] + tmp_f_acARBR0_[4] + tmp_f_acARBR0_[5] + tmp_f_acARBR0_[6] + tmp_f_acARBR0_[7];
	tmp_f_acARBI0 = tmp_f_acARBI0_[0] + tmp_f_acARBI0_[1] + tmp_f_acARBI0_[2] + tmp_f_acARBI0_[3] + tmp_f_acARBI0_[4] + tmp_f_acARBI0_[5] + tmp_f_acARBI0_[6] + tmp_f_acARBI0_[7];
	tmp_f_acAIBR0 = tmp_f_acAIBR0_[0] + tmp_f_acAIBR0_[1] + tmp_f_acAIBR0_[2] + tmp_f_acAIBR0_[3] + tmp_f_acAIBR0_[4] + tmp_f_acAIBR0_[5] + tmp_f_acAIBR0_[6] + tmp_f_acAIBR0_[7];
	tmp_f_acAIBI0 = tmp_f_acAIBI0_[0] + tmp_f_acAIBI0_[1] + tmp_f_acAIBI0_[2] + tmp_f_acAIBI0_[3] + tmp_f_acAIBI0_[4] + tmp_f_acAIBI0_[5] + tmp_f_acAIBI0_[6] + tmp_f_acAIBI0_[7];
	f_CR__[na] = tmp_f_acARBR0 + tmp_f_acAIBI0;
	f_CI__[na] = tmp_f_acARBI0 - tmp_f_acAIBR0;
	na += 1;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if ((f_CR__!=NULL) && (f_CI__!=NULL)){ } */}
}
void hp_ps_mult_immintrin_test()
{
  int verbose = 1;
  int n_row_A = 153;
  int n_col_X = 1152;
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int n_row_B = 147;
  int n_row_A_sub = 5;
  int n_col_X_sub = 15;
  int n_row_B_sub = 4;
  int nrow_A=0,ncol_X=0,nrow_B=0;
  unsigned long long int ulli_A_total = (unsigned long long int) n_row_A * (unsigned long long int) n_col_X_rup ;
  unsigned long long int ulli_B_total = (unsigned long long int) n_row_B * (unsigned long long int) n_col_X_rup ;
  unsigned long long int ulli_C_total = (unsigned long long int) n_row_A * (unsigned long long int) n_row_B ;
  unsigned long long int ulli=0;
  unsigned long long int tab=0;
  __m256 *ps_AR_trn__ = NULL;
  __m256 *ps_AI_trn__ = NULL;
  __m256 *ps_BR_trn__ = NULL;
  __m256 *ps_BI_trn__ = NULL;
  float *f_AR_trn__ = NULL;
  float *f_AI_trn__ = NULL;
  float *f_AR_u_trn__ = NULL;
  float *f_AI_u_trn__ = NULL;
  float *f_BR_trn__ = NULL;
  float *f_BI_trn__ = NULL;
  float *f_BR_u_trn__ = NULL;
  float *f_BI_u_trn__ = NULL;
  float *f_CR_bf__ = NULL;
  float *f_CI_bf__ = NULL;
  float *f_CR_al__ = NULL;
  float *f_CI_al__ = NULL;
  float *f_AR_sub__ = NULL;
  float *f_AI_sub__ = NULL;
  float *f_BR_sub__ = NULL;
  float *f_BI_sub__ = NULL;
  float *f_CR_sub__ = NULL;
  float *f_CI_sub__ = NULL;
  float complex *c_A_u_trn__ = NULL;
  float complex *c_B_u_trn__ = NULL;
  float complex *c_C_bf__ = NULL;
  float complex *c_C_al__ = NULL;
  float complex *c_A_sub__ = NULL;
  float complex *c_B_sub__ = NULL;
  float complex *c_C_sub__ = NULL;
  float ferror=0;
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  /* %%%% */
  ps_AR_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  ps_AI_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  ps_BR_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  ps_BI_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  f_AR_trn__ = (float *) ps_AR_trn__;
  f_AI_trn__ = (float *) ps_AI_trn__;
  f_BR_trn__ = (float *) ps_BR_trn__;
  f_BI_trn__ = (float *) ps_BI_trn__;
  f_AR_u_trn__ = (float *) malloc(n_col_X*n_row_A*sizeof(float));
  f_AI_u_trn__ = (float *) malloc(n_col_X*n_row_A*sizeof(float));
  f_BR_u_trn__ = (float *) malloc(n_col_X*n_row_A*sizeof(float));
  f_BI_u_trn__ = (float *) malloc(n_col_X*n_row_A*sizeof(float));
  f_CR_bf__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_CI_bf__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_CR_al__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_CI_al__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_AR_sub__ = (float *) malloc(n_row_A_sub*n_col_X_sub*sizeof(float));
  f_AI_sub__ = (float *) malloc(n_row_A_sub*n_col_X_sub*sizeof(float));
  f_BR_sub__ = (float *) malloc(n_row_B_sub*n_col_X_sub*sizeof(float));
  f_BI_sub__ = (float *) malloc(n_row_B_sub*n_col_X_sub*sizeof(float));
  f_CR_sub__ = (float *) malloc(n_row_A_sub*n_row_B_sub*sizeof(float));
  f_CI_sub__ = (float *) malloc(n_row_A_sub*n_row_B_sub*sizeof(float));
  c_A_u_trn__ = (float complex *) malloc(n_col_X*n_row_A*sizeof(float complex));
  c_B_u_trn__ = (float complex *) malloc(n_col_X*n_row_A*sizeof(float complex));
  c_C_bf__ = (float complex *) malloc(n_row_A*n_row_B*sizeof(float complex));
  c_C_al__ = (float complex *) malloc(n_row_A*n_row_B*sizeof(float complex));
  c_A_sub__ = (float complex *) malloc(n_row_A_sub*n_col_X_sub*sizeof(float complex));
  c_B_sub__ = (float complex *) malloc(n_row_B_sub*n_col_X_sub*sizeof(float complex));
  c_C_sub__ = (float complex *) malloc(n_row_A_sub*n_row_B_sub*sizeof(float complex));
  tab = (unsigned long long int) n_row_A * (unsigned long long int) n_row_B * (unsigned long long int) n_col_X_rup;
  local_tic(0,t_start_,d_start_);
  ulli=0;
  for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      f_AR_trn__[ncol_X + nrow_A*n_col_X_rup] = (float)(((int)ulli%7)-3);
      f_AI_trn__[ncol_X + nrow_A*n_col_X_rup] = (float)(((int)ulli%7)-2);
      f_AR_u_trn__[ncol_X + nrow_A*n_col_X] = (float)(((int)ulli%7)-3);
      f_AI_u_trn__[ncol_X + nrow_A*n_col_X] = (float)(((int)ulli%7)-2);
      c_A_u_trn__[ncol_X + nrow_A*n_col_X] = (float complex) f_AR_u_trn__[ncol_X + nrow_A*n_col_X] + _Complex_I * (float complex) f_AI_u_trn__[ncol_X + nrow_A*n_col_X];
      ulli++;
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}  
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ 
      f_AR_sub__[nrow_A+ncol_X*n_row_A_sub] = f_AR_trn__[ncol_X + nrow_A*n_col_X_rup];
      f_AI_sub__[nrow_A+ncol_X*n_row_A_sub] = f_AI_trn__[ncol_X + nrow_A*n_col_X_rup];
      c_A_sub__[nrow_A+ncol_X*n_row_A_sub] = (float complex) f_AR_sub__[nrow_A+ncol_X*n_row_A_sub] + _Complex_I * (float complex) f_AI_sub__[nrow_A+ncol_X*n_row_A_sub];
      /* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ }} */}}
  if (verbose>1){
    printf(" %% upper corner of f_AR_trn__: \n");
    array_printf(f_AR_sub__,"float",n_row_A_sub,n_col_X_sub," % f_AR_sub__: ");
    printf(" %% upper corner of f_AI_trn__: \n");
    array_printf(f_AI_sub__,"float",n_row_A_sub,n_col_X_sub," % f_AI_sub__: ");
    printf(" %% upper corner of c_A_trn__: \n");
    array_printf(c_A_sub__,"float complex",n_row_A_sub,n_col_X_sub," % c_A_sub__: ");
    /* if (verbose>1){ } */}
  ulli=0;
  for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      f_BR_trn__[ncol_X + nrow_B*n_col_X_rup] = (float)(((int)ulli%11)-2);
      f_BI_trn__[ncol_X + nrow_B*n_col_X_rup] = (float)(((int)ulli%11)-1);
      f_BR_u_trn__[ncol_X + nrow_B*n_col_X] = (float)(((int)ulli%11)-2);
      f_BI_u_trn__[ncol_X + nrow_B*n_col_X] = (float)(((int)ulli%11)-1);
      c_B_u_trn__[ncol_X + nrow_B*n_col_X] = (float complex) f_BR_u_trn__[ncol_X + nrow_B*n_col_X] + _Complex_I * (float complex) f_BI_u_trn__[ncol_X + nrow_B*n_col_X];
      ulli++;
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}
  if (verbose>1){
    for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){
	f_BR_sub__[nrow_B+ncol_X*n_row_B_sub] = f_BR_trn__[ncol_X + nrow_B*n_col_X_rup];
	f_BI_sub__[nrow_B+ncol_X*n_row_B_sub] = f_BI_trn__[ncol_X + nrow_B*n_col_X_rup];
	c_B_sub__[nrow_B+ncol_X*n_row_B_sub] = (float complex) f_BR_sub__[nrow_B+ncol_X*n_row_B_sub] + _Complex_I * (float complex) f_BI_sub__[nrow_B+ncol_X*n_row_B_sub];
	/* for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ }} */}}
    printf(" %% upper corner of f_BR_trn__: \n");
    array_printf(f_BR_sub__,"float",n_row_B_sub,n_col_X_sub," % f_BR_sub__: ");
    printf(" %% upper corner of f_BI_trn__: \n");
    array_printf(f_BI_sub__,"float",n_row_B_sub,n_col_X_sub," % f_BI_sub__: ");
    printf(" %% upper corner of c_B_trn__: \n");
    array_printf(c_B_sub__,"float complex",n_row_B_sub,n_col_X_sub," % c_B_sub__: ");
    /* if (verbose>1){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," initialize: ");
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(c_C_bf__,0,ulli_C_total*sizeof(float complex));
  hp_interleave_mult_bruteforce(n_row_A,n_col_X,c_A_u_trn__,n_row_B,c_B_u_trn__,&c_C_bf__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_ps_mult_bruteforce: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_bf__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_bf__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_bf__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_bf__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(c_C_al__,0,ulli_C_total*sizeof(float complex));
  hp_segregated_mult_bruteforce(n_row_A,n_col_X,f_AR_u_trn__,f_AI_u_trn__,n_row_B,f_BR_u_trn__,f_BI_u_trn__,&c_C_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_segregated_mult_bruteforce: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(c_C_al__,0,ulli_C_total*sizeof(float complex));
  hp_segregated_to_interleaved_mult_immintrin_loadu_fma(n_row_A,n_col_X,f_AR_u_trn__,f_AI_u_trn__,n_row_B,f_BR_u_trn__,f_BI_u_trn__,&c_C_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_segregated_to_interleaved_mult_immintrin_loadu_fma: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_CR_al__,0,ulli_C_total*sizeof(float));
  memset(f_CI_al__,0,ulli_C_total*sizeof(float));
  hp_segregated_to_segregated_mult_immintrin_loadu_fma(n_row_A,n_col_X,f_AR_u_trn__,f_AI_u_trn__,n_row_B,f_BR_u_trn__,f_BI_u_trn__,&f_CR_al__,&f_CI_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_segregated_to_segregated_mult_immintrin_loadu_fma: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      c_C_al__[nrow_A+nrow_B*n_row_A] = f_CR_al__[nrow_A+nrow_B*n_row_A] + _Complex_I * (float complex) f_CI_al__[nrow_A+nrow_B*n_row_A];
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ }} */}}
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_CR_al__,0,ulli_C_total*sizeof(float));
  memset(f_CI_al__,0,ulli_C_total*sizeof(float));
  hp_segregated_to_segregated_mult_immintrin_load1_fma(n_row_A,n_col_X,f_AR_trn__,f_AI_trn__,n_row_B,f_BR_trn__,f_BI_trn__,&f_CR_al__,&f_CI_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_segregated_to_segregated_mult_immintrin_load1_fma: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      c_C_al__[nrow_A+nrow_B*n_row_A] = f_CR_al__[nrow_A+nrow_B*n_row_A] + _Complex_I * (float complex) f_CI_al__[nrow_A+nrow_B*n_row_A];
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ }} */}}
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
#ifdef _CBLAS
  local_tic(0,t_start_,d_start_);
  memset(c_C_al__,0,ulli_C_total*sizeof(float complex));
  hp_interleave_mult_cblas_cgemm(n_row_A,n_col_X,c_A_u_trn__,n_row_B,c_B_u_trn__,&c_C_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_interleave_mult_cblas_cgemm: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
#endif /* _CBLAS */
  /* %%%%%%%% */
  _mm_free(ps_AR_trn__); ps_AR_trn__ = NULL;
  _mm_free(ps_AI_trn__); ps_AI_trn__ = NULL;
  _mm_free(ps_BR_trn__); ps_BR_trn__ = NULL;
  _mm_free(ps_BI_trn__); ps_BI_trn__ = NULL;
  free(f_AR_u_trn__); f_AR_u_trn__=NULL;
  free(f_AI_u_trn__); f_AI_u_trn__=NULL;
  free(f_BR_u_trn__); f_BR_u_trn__=NULL;
  free(f_BI_u_trn__); f_BI_u_trn__=NULL;
  free(f_CR_bf__); f_CR_bf__=NULL;
  free(f_CI_bf__); f_CI_bf__=NULL;
  free(f_CR_al__); f_CR_al__=NULL;
  free(f_CI_al__); f_CI_al__=NULL;
  free(f_AR_sub__); f_AR_sub__=NULL;
  free(f_AI_sub__); f_AI_sub__=NULL;
  free(f_BR_sub__); f_BR_sub__=NULL;
  free(f_BI_sub__); f_BI_sub__=NULL;
  free(f_CR_sub__); f_CR_sub__=NULL;
  free(f_CI_sub__); f_CI_sub__=NULL;
  free(c_A_u_trn__); c_A_u_trn__=NULL;
  free(c_B_u_trn__); c_B_u_trn__=NULL;
  free(c_C_bf__); c_C_bf__=NULL;
  free(c_C_al__); c_C_al__=NULL;
  free(c_A_sub__); c_A_sub__=NULL;
  free(c_B_sub__); c_B_sub__=NULL;
  free(c_C_sub__); c_C_sub__=NULL;
  //wkspace_printf();
}
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void test_cblas()
{
  int verbose=1;
  int M=0,N=0,K=0,T=0,nT=0;
  float *sA__=NULL;
  float *sB__=NULL;
  float *sC__=NULL;
  double *dA__=NULL;
  double *dB__=NULL;
  double *dC__=NULL;
  complex *cA__=NULL;
  complex *cB__=NULL;
  complex *cC__=NULL;
  double complex *zA__=NULL;
  double complex *zB__=NULL;
  double complex *zC__=NULL;
  unsigned long long int tab=0;
  /* %%%% */
  float cblas_salpha = (float) 1.0;
  float cblas_sbeta  = (float) 0.0;
  double cblas_dalpha = (double) 1.0;
  double cblas_dbeta  = (double) 0.0;
  complex cblas_calpha = (complex) 1.0;
  complex cblas_cbeta  = (complex) 0.0;
  double complex cblas_zalpha = (double complex) 1.0;
  double complex cblas_zbeta  = (double complex) 0.0;
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  if (verbose){ printf(" %% [entering test_cblas]\n");}
  /* %%%% */
  M = 130; N = 310; K = 160; T = 200;
  sA__ = (float *) calloc(M*K,sizeof(float)); sB__ = (float *) calloc(N*K,sizeof(float)); sC__ = (float *) calloc(M*N,sizeof(float));
  dA__ = (double *) calloc(M*K,sizeof(double)); dB__ = (double *) calloc(N*K,sizeof(double)); dC__ = (double *) calloc(M*N,sizeof(double));
  cA__ = (complex *) calloc(M*K,sizeof(complex)); cB__ = (complex *) calloc(N*K,sizeof(complex)); cC__ = (complex *) calloc(M*N,sizeof(complex));
  zA__ = (double complex *) calloc(M*K,sizeof(double complex)); zB__ = (double complex *) calloc(N*K,sizeof(double complex)); zC__ = (double complex *) calloc(M*N,sizeof(double complex));
  /* %%%% */
  local_tic(0,t_start_,d_start_);
  tab = (unsigned long long int)(M*N*K)*(unsigned long long int)T;
  for (nT=0;nT<T;nT++){
    cblas_sgemm(
		CblasColMajor
		,CblasTrans
		,CblasNoTrans
		,M
		,N
		,K
		,cblas_salpha
		,sA__
		,K
		,sB__
		,K
		,cblas_sbeta
		,sC__
		,M
		);
    /* for (nT=0;nT<T;nT++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"sgemm: ");  
  /* %%%% */
  local_tic(0,t_start_,d_start_);
  tab = (unsigned long long int)(M*N*K)*(unsigned long long int)T;
  for (nT=0;nT<T;nT++){
    cblas_dgemm(
		CblasColMajor
		,CblasTrans
		,CblasNoTrans
		,M
		,N
		,K
		,cblas_dalpha
		,dA__
		,K
		,dB__
		,K
		,cblas_dbeta
		,dC__
		,M
		);
    /* for (nT=0;nT<T;nT++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"dgemm: ");  
  /* %%%% */
  local_tic(0,t_start_,d_start_);
  tab = (unsigned long long int)(M*N*K)*(unsigned long long int)T;
  for (nT=0;nT<T;nT++){
    cblas_cgemm(
		CblasColMajor
		,CblasConjTrans
		,CblasNoTrans
		,M
		,N
		,K
		,&cblas_calpha
		,cA__
		,K
		,cB__
		,K
		,&cblas_cbeta
		,cC__
		,M
		);
    /* for (nT=0;nT<T;nT++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"cgemm: ");  
  /* %%%% */
  local_tic(0,t_start_,d_start_);
  tab = (unsigned long long int)(M*N*K)*(unsigned long long int)T;
  for (nT=0;nT<T;nT++){
    cblas_zgemm(
		CblasColMajor
		,CblasConjTrans
		,CblasNoTrans
		,M
		,N
		,K
		,&cblas_zalpha
		,zA__
		,K
		,zB__
		,K
		,&cblas_zbeta
		,zC__
		,M
		);
    /* for (nT=0;nT<T;nT++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"zgemm: ");  
  /* %%%% */
  free(sA__);free(sB__);free(sC__);
  free(dA__);free(dB__);free(dC__);
  free(cA__);free(cB__);free(cC__);
  free(zA__);free(zB__);free(zC__);
  /* %%%% */
  if (verbose){ printf(" %% [finished test_cblas]\n");}  
}
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* The computational routine */
void mex_ampmh_X_wSM___10
(
  int n_M_per_Mbatch
 ,int n_S_per_Sbatch
 ,int flag_optimize_over_gamma_z
 ,int flag_compute_I_value
 ,double tolerance_master
 ,int FTK_n_svd_l
 ,int FTK_n_delta_v
 ,double *FTK_svd_U_d_expiw_s_real__
 ,double *FTK_svd_U_d_expiw_s_imag__
 ,double *FTK_delta_x_
 ,double *FTK_delta_y_
 ,int n_w_max
 ,int pm_n_UX_rank
 ,int n_S
 ,double *CTF_UX_S_k_q_wnS_real__
 ,double *CTF_UX_S_k_q_wnS_imag__
 ,double *CTF_UX_S_l2_
 ,int n_M
 ,double *svd_VUXM_lwnM_real____
 ,double *svd_VUXM_lwnM_imag____
 ,double *UX_M_l2_dM__
 ,double *X_wSM___
 ,double *delta_x_wSM___
 ,double *delta_y_wSM___
 ,double *gamma_z_wSM___
 ,double *I_value_wSM___
 )
{
  int verbose=2;
  int flag_disp=0;
  int flag_dump=0;
  double complex *FTK_svd_U_d_expiw_s__=NULL;
  double complex *CTF_UX_S_k_q_wnS__=NULL;
  double complex *svd_VUXM_lwnM____=NULL;
  unsigned long long int na=0,n_a=0;
  /* %%%% */
  int nMbatch=0,n_Mbatch=0,nSbatch=0,n_Sbatch=0;
  int *index_M_in_Mbatch_=NULL;
  int *index_S_in_Sbatch_=NULL;
  int n_M_sub=0,nM_sub=0;
  int n_S_sub=0,nS_sub=0;
  int nl=0,ndelta_v=0,nw=0,ndw=0,nS=0,nM=0;
  int pm_nUX_rank=0;
  unsigned long long int tab=0,tabA=0,tabB=0,tabC=0;
  double complex *CTF_UX_S_k_q_nSw___=NULL;
  double complex *svd_VUXM_nMwl____=NULL;
  double complex *svd_SVUXM_SMwl____=NULL;
  fftw_complex *svd_SVUXM_0in_wlSM____;
  fftw_complex *svd_SVUXM_out_wlSM____;
  double complex *svd_SVUXM_lwSM____=NULL;
  double complex *svd_SVUXM_lwsM____=NULL;
  double complex *svd_USESVUXM_dwSM____=NULL;
  double complex cblas_alpha = (double complex) 1.0;
  double complex cblas_beta  = (double complex) 0.0;
  double *l2_dSM___=NULL;
  double *n2_dSM___=NULL;
  double *f2_dSM___=NULL;
  double *X_sub_dwSM____=NULL;
  double *I_value_use_dwSM____=NULL;
  double *gamma_z_=NULL;
  /* %%%% */
  fftw_plan fftw_plan_many_plan;
  int fftw_plan_many_rank=0;
  int fftw_plan_many_n[1]={0};
  int fftw_plan_many_howmany=0;
  int fftw_plan_many_istride=0;
  int fftw_plan_many_idist=0;
  int fftw_plan_many_ostride=0;
  int fftw_plan_many_odist=0;
  int fftw_plan_many_sign=0;
  unsigned int fftw_plan_many_flags=0;
  /* %%%% */
  int MDA_n_dim=0,MDA_ndim = 0;
  int MDA_dim_[8];
  char MDA_fname[32];
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  /* %%%% */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___10]\n");}
  /* %%%%%%%%%%%%%%%% */
  n_a = FTK_n_delta_v*FTK_n_svd_l;
  FTK_svd_U_d_expiw_s__ = double_complex_malloc_and_interleave(n_a,FTK_svd_U_d_expiw_s_real__,FTK_svd_U_d_expiw_s_imag__);
  n_a = n_w_max*pm_n_UX_rank*n_S;
  CTF_UX_S_k_q_wnS__ = double_complex_malloc_and_interleave(n_a,CTF_UX_S_k_q_wnS_real__,CTF_UX_S_k_q_wnS_imag__);
  n_a = FTK_n_svd_l*n_w_max*pm_n_UX_rank*n_M;
  svd_VUXM_lwnM____ = double_complex_malloc_and_interleave(n_a,svd_VUXM_lwnM_real____,svd_VUXM_lwnM_imag____);
  /* %%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% n_M_per_Mbatch %d\n",n_M_per_Mbatch);}
  if (verbose){ printf(" %% n_S_per_Sbatch %d\n",n_S_per_Sbatch);}
  if (verbose){ printf(" %% flag_optimize_over_gamma_z %d\n",flag_optimize_over_gamma_z);}
  if (verbose){ printf(" %% flag_compute_I_value %d\n",flag_compute_I_value);}
  if (verbose){ printf(" %% tolerance_master %0.16f\n",tolerance_master);}
  if (verbose){ printf(" %% FTK_n_svd_l %d\n",FTK_n_svd_l);}
  if (verbose){ printf(" %% FTK_n_delta_v %d\n",FTK_n_delta_v);}
  n_a = FTK_n_delta_v*FTK_n_svd_l;
  if (verbose){ printf(" %% FTK_svd_U_d_expiw_s__ %+0.16f %+0.16f*i --> %+0.16f %+0.16f*i\n",creal(FTK_svd_U_d_expiw_s__[0]),cimag(FTK_svd_U_d_expiw_s__[0]),creal(FTK_svd_U_d_expiw_s__[n_a-1]),cimag(FTK_svd_U_d_expiw_s__[n_a-1]));}
  if (verbose){ printf(" %% FTK_delta_x_ %0.16f --> %0.16f\n",FTK_delta_x_[0],FTK_delta_x_[FTK_n_delta_v-1]);}
  if (verbose){ printf(" %% FTK_delta_y_ %0.16f --> %0.16f\n",FTK_delta_y_[0],FTK_delta_y_[FTK_n_delta_v-1]);}
  if (verbose){ printf(" %% n_w_max %d\n",n_w_max);}
  if (verbose){ printf(" %% pm_n_UX_rank %d\n",pm_n_UX_rank);}
  if (verbose){ printf(" %% n_S %d\n",n_S);}
  n_a = n_w_max*pm_n_UX_rank*n_S;
  if (verbose){ printf(" %% CTF_UX_S_k_q_wnS__ %+0.16f %+0.16f*i --> %+0.16f %+0.16f*i\n",creal(CTF_UX_S_k_q_wnS__[0]),cimag(CTF_UX_S_k_q_wnS__[0]),creal(CTF_UX_S_k_q_wnS__[n_a-1]),cimag(CTF_UX_S_k_q_wnS__[n_a-1]));}
  if (verbose){ printf(" %% CTF_UX_S_l2_ %0.16f --> %0.16f\n",CTF_UX_S_l2_[0],CTF_UX_S_l2_[n_S-1]);}
  if (verbose){ printf(" %% n_M %d\n",n_M);}
  n_a = FTK_n_svd_l*n_w_max*pm_n_UX_rank*n_M;
  if (verbose){ printf(" %% svd_VUXM_lwnM____ %+0.16f %+0.16f*i --> %+0.16f %+0.16f*i\n",creal(svd_VUXM_lwnM____[0]),cimag(svd_VUXM_lwnM____[0]),creal(svd_VUXM_lwnM____[n_a-1]),cimag(svd_VUXM_lwnM____[n_a-1]));}
  if (verbose){ printf(" %% UX_M_l2_dM__ %0.16f --> %0.16f\n",UX_M_l2_dM__[0],UX_M_l2_dM__[FTK_n_delta_v*n_M-1]);}
  /* %%%% */
  /* if ( (nargout>5) & (isempty(I_value_wSM___)) ); I_value_wSM___ = ones(size(X_wSM___)); end; */
  /* %%%% */
  if (flag_compute_I_value){
    if (flag_optimize_over_gamma_z==0){ n_a = (unsigned long long int)n_M*(unsigned long long int)n_S*(unsigned long long int)n_w_max;}
    if (flag_optimize_over_gamma_z==1){ n_a = (unsigned long long int)n_M*(unsigned long long int)n_S;}
    for (na=0;na<n_a;na++){ I_value_wSM___[na] = 1.0;}
    /* if (flag_compute_I_value){ } */}
  /* %%%% */
  /* CTF_UX_S_k_q_nSw___ = permute(CTF_UX_S_k_q_wnS___,[2,3,1]); */
  /* %%%% */
  local_tic(0,t_start_,d_start_);
  tab = (unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_S*(unsigned long long int)n_w_max;
  CTF_UX_S_k_q_nSw___ = (double complex *) malloc(tab*sizeof(double complex));
  na=0;
  for (nS=0;nS<n_S;nS++){
    for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){
      for (nw=0;nw<n_w_max;nw++){
	/* tabA = */
	/*   (unsigned long long int)pm_nUX_rank */
	/*   + (unsigned long long int)nS*(unsigned long long int)pm_n_UX_rank */
	/*   + (unsigned long long int)nw*(unsigned long long int)n_S*(unsigned long long int)pm_n_UX_rank; */
	tabA = (unsigned long long int)pm_nUX_rank + ((unsigned long long int)nS + (unsigned long long int)nw*(unsigned long long int)n_S)*(unsigned long long int)pm_n_UX_rank;
	CTF_UX_S_k_q_nSw___[tabA] = CTF_UX_S_k_q_wnS__[na]; na++;
	/* for (nw=0;nw<n_w_max;nw++){ } */}
      /* for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){ } */}
    /* for (nS=0;nS<n_S;nS++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"CTF_UX_S_k_q_nSw___: ");
  /* %%%%%%%%%%%%%%%% */
  tab = (unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
  svd_VUXM_nMwl____ = (double complex *) malloc(tab*sizeof(double complex));
  tab = (unsigned long long int)n_S*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
  svd_SVUXM_SMwl____ = (double complex *) malloc(tab*sizeof(double complex));
  svd_SVUXM_0in_wlSM____ = (fftw_complex *) fftw_malloc(tab*sizeof(fftw_complex));
  svd_SVUXM_out_wlSM____ = (fftw_complex *) fftw_malloc(tab*sizeof(fftw_complex));
  fftw_plan_many_rank = 1;
  fftw_plan_many_n[0] = n_w_max;
  fftw_plan_many_howmany = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_S*(unsigned long long int)n_M_per_Mbatch;
  fftw_plan_many_istride = 1;
  fftw_plan_many_idist = n_w_max;
  fftw_plan_many_ostride = 1;
  fftw_plan_many_odist = n_w_max;
  fftw_plan_many_sign = FFTW_BACKWARD;
  fftw_plan_many_flags = FFTW_ESTIMATE;
  fftw_plan_many_plan = fftw_plan_many_dft(
					   fftw_plan_many_rank
					   ,fftw_plan_many_n
					   ,fftw_plan_many_howmany
					   ,svd_SVUXM_0in_wlSM____
					   ,NULL
					   ,fftw_plan_many_istride
					   ,fftw_plan_many_idist
					   ,svd_SVUXM_out_wlSM____
					   ,NULL
					   ,fftw_plan_many_ostride
					   ,fftw_plan_many_odist
					   ,fftw_plan_many_sign
					   ,fftw_plan_many_flags
					   );
  svd_SVUXM_lwSM____ = (double complex *) malloc(tab*sizeof(double complex));
  tab = (unsigned long long int)n_S_per_Sbatch*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
  svd_SVUXM_lwsM____ = (double complex *) malloc(tab*sizeof(double complex));
  tab = (unsigned long long int)n_S_per_Sbatch*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
  svd_USESVUXM_dwSM____ = (double complex *) malloc(tab*sizeof(double complex));
  tab = (unsigned long long int)n_S_per_Sbatch*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)FTK_n_delta_v;
  l2_dSM___ = (double *) malloc(tab*sizeof(double));
  n2_dSM___ = (double *) malloc(tab*sizeof(double));
  f2_dSM___ = (double *) malloc(tab*sizeof(double));
  tab = (unsigned long long int)n_S_per_Sbatch*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
  X_sub_dwSM____ = (double *) malloc(tab*sizeof(double));
  if (flag_compute_I_value){ I_value_use_dwSM____ = (double *) malloc(tab*sizeof(double));}
  tab = (unsigned long long int)n_w_max;
  gamma_z_ = (double *) malloc(tab*sizeof(double));
  for (nw=0;nw<n_w_max;nw++){ gamma_z_[nw] = 2*PI*(double)nw/(double)n_w_max;}
  /* %%%%%%%%%%%%%%%% */
  n_Mbatch = ceil((double)n_M/(double)n_M_per_Mbatch);
  if (verbose>1){ printf(" %% n_Mbatch %d\n",n_Mbatch);}
  index_M_in_Mbatch_ = (int *) malloc(n_M_per_Mbatch*sizeof(int));
  n_Sbatch = ceil((double)n_S/(double)n_S_per_Sbatch);
  if (verbose>1){ printf(" %% n_Sbatch %d\n",n_Sbatch);}
  index_S_in_Sbatch_ = (int *) malloc(n_S_per_Sbatch*sizeof(int));
  for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){
    n_M_sub = n_M_per_Mbatch; if (nMbatch==n_Mbatch-1){ n_M_sub = n_M - n_M_per_Mbatch*nMbatch;}
    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ index_M_in_Mbatch_[nM_sub] = n_M_per_Mbatch*nMbatch + nM_sub;}
    if (verbose>1){ printf(" %% nMbatch %d/%d: index_M_in_Mbatch_ %d --> %d\n",nMbatch,n_Mbatch,index_M_in_Mbatch_[0],index_M_in_Mbatch_[n_M_sub-1]);}
    if (n_M_sub>0){
      /* %%%% */
      /* svd_VUXM_nMwl____ = permute(svd_VUXM_lwnM____(:,:,:,1+index_M_in_Mbatch_),[3,4,2,1]); */
      /* %%%% */      
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)n_M_sub*(unsigned long long int)pm_n_UX_rank;
      na = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)(n_M_per_Mbatch*nMbatch);
      for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){
	  for (nw=0;nw<n_w_max;nw++){
	    for (nl=0;nl<FTK_n_svd_l;nl++){
	      /* tabA = */
	      /* 	(unsigned long long int)pm_nUX_rank */
	      /* 	+ (unsigned long long int)nM_sub*(unsigned long long int)pm_n_UX_rank */
	      /* 	+ (unsigned long long int)nw*(unsigned long long int)n_M_sub*(unsigned long long int)pm_n_UX_rank */
	      /* 	+ (unsigned long long int)nl*(unsigned long long int)n_w_max*(unsigned long long int)n_M_sub*(unsigned long long int)pm_n_UX_rank; */
	      tabA = (unsigned long long int)pm_nUX_rank + ((unsigned long long int)nM_sub + ((unsigned long long int)nw + (unsigned long long int)nl*(unsigned long long int)n_w_max)*(unsigned long long int)n_M_sub)*(unsigned long long int)pm_n_UX_rank;
	      svd_VUXM_nMwl____[tabA] = svd_VUXM_lwnM____[na]; na++;
	      /* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
	    /* for (nw=0;nw<n_w_max;nw++){ } */}
	  /* for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){ } */}
	/* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"svd_VUXM_nMwl____: ");
      /* %%%% */
      /* 
	 svd_SVUXM_SMwl____ = zeros(n_S,n_M_sub,n_w_max,FTK.n_svd_l);
	 for nl=0:FTK.n_svd_l-1;
	 for nw=0:n_w_max-1;
	 svd_SVUXM_SMwl____(:,:,1+nw,1+nl) = ctranspose(CTF_UX_S_k_q_nSw___(:,:,1+nw))*svd_VUXM_nMwl____(:,:,1+nw,1+nl);
	 end;%for nw=0:n_w_max-1;
	 end;%for nl=0:FTK.n_svd_l-1;
      */
      /* %%%% */
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)(FTK_n_svd_l*n_w_max)*(unsigned long long int)(n_S*n_M_sub)*(unsigned long long int)pm_n_UX_rank;
      for (nl=0;nl<FTK_n_svd_l;nl++){
	for (nw=0;nw<n_w_max;nw++){
	  tabA = (unsigned long long int)nw*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_S;
	  tabB = (unsigned long long int)(nw+nl*n_w_max)*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_M_sub;
	  tabC = (unsigned long long int)(nw+nl*n_w_max)*(unsigned long long int)n_S*(unsigned long long int)n_M_sub;
	  cblas_zgemm(
		      CblasColMajor
		      ,CblasConjTrans
		      ,CblasNoTrans
		      ,n_S
		      ,n_M_sub
		      ,pm_n_UX_rank
		      ,&cblas_alpha
		      ,CTF_UX_S_k_q_nSw___ + tabA
		      ,pm_n_UX_rank
		      ,svd_VUXM_nMwl____ + tabB
		      ,pm_n_UX_rank
		      ,&cblas_beta
		      ,svd_SVUXM_SMwl____ + tabC
		      ,n_S
		      );
	  /* for (nw=0;nw<n_w_max;nw++){ } */}
	/* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"svd_SVUXM_SMwl____: ");
      /* %%%% */
      /* svd_SVUXM_lwSM____ = permute(ifft(permute(svd_SVUXM_SMwl____,[3,4,1,2]),[],1)*n_w_max,[2,1,3,4]); */
      /* %%%% */      
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)n_S*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
      memset(svd_SVUXM_0in_wlSM____,0,tab*sizeof(fftw_complex));
      na=0; tabA=0;
      for (nl=0;nl<FTK_n_svd_l;nl++){
	for (nw=0;nw<n_w_max;nw++){
	  for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	    for (nS=0;nS<n_S;nS++){
	      /* tabA = */
	      /* 	(unsigned long long int)nw */
	      /* 	+ (unsigned long long int)nl*(unsigned long long int)n_w_max */
	      /* 	+ (unsigned long long int)nS*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l */
	      /* 	+ (unsigned long long int)nM_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_S; */
	      tabA = (unsigned long long int)nw + ((unsigned long long int)nl + ((unsigned long long int)nS + (unsigned long long int)nM_sub*(unsigned long long int)n_S)*(unsigned long long int)FTK_n_svd_l)*(unsigned long long int)n_w_max;
	      svd_SVUXM_0in_wlSM____[tabA] = (fftw_complex) (svd_SVUXM_SMwl____[na]); na++;
	      /* for (nS=0;nS<n_S;nS++){ } */}
	    /* for (nM_sub=0;nM_sub<n_M_sub;nS++){ } */}
	  /* for (nw=0;nw<n_w_max;nw++){ } */}
	/* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"svd_SVUXM_0in_wlSM____: ");
      local_tic(0,t_start_,d_start_);
      fftw_execute(fftw_plan_many_plan);
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"svd_SVUXM_out_wlSM____: ");
      tab = (unsigned long long int)n_S*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
      memset(svd_SVUXM_lwSM____,0,tab*sizeof(double complex));
      na=0; tabA=0;
      for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	for (nS=0;nS<n_S;nS++){
	  for (nl=0;nl<FTK_n_svd_l;nl++){
	    for (nw=0;nw<n_w_max;nw++){
	      /* tabA = */
	      /* 	(unsigned long long int)nl */
	      /* 	+ (unsigned long long int)nw*(unsigned long long int)FTK_n_svd_l */
	      /* 	+ (unsigned long long int)nS*(unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max */
	      /* 	+ (unsigned long long int)nM_sub*(unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)n_S; */
	      tabA = (unsigned long long int)nl + ((unsigned long long int)nw + ((unsigned long long int)nS + (unsigned long long int)nM_sub*(unsigned long long int)n_S)*(unsigned long long int)n_w_max)*(unsigned long long int)FTK_n_svd_l;
	      svd_SVUXM_lwSM____[tabA] = (double complex) (svd_SVUXM_out_wlSM____[na]); na++;
	      /* for (nw=0;nw<n_w_max;nw++){ } */}
	    /* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
	  /* for (nS=0;nS<n_S;nS++){ } */}
	/* for (nM_sub=0;nM_sub<n_M_sub;nS++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"svd_SVUXM_lwSM____: ");
      /* %%%% */
      /* display output */
      /* %%%% */
      if (flag_disp){
	printf(" %% nMbatch %d/%d\n",nMbatch,n_Mbatch);
	tab = minimum(8,(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_S*(unsigned long long int)n_M_sub);
	array_printf(svd_SVUXM_SMwl____,"double complex",1,tab," %% svd_SVUXM_SMwl____: ");
	array_printf(svd_SVUXM_lwSM____,"double complex",1,tab," %% svd_SVUXM_lwSM____: ");
	/* if (flag_disp){ } */}
      if (flag_dump){
	MDA_ndim=0;
	MDA_dim_[MDA_ndim++] = n_S;
	MDA_dim_[MDA_ndim++] = n_M_sub;
	MDA_dim_[MDA_ndim++] = n_w_max;
	MDA_dim_[MDA_ndim++] = FTK_n_svd_l;
	MDA_n_dim = MDA_ndim;
	sprintf(MDA_fname,"mex_svd_SVUXM_SMwl____.mat");
	MDA_write_c16(MDA_n_dim,MDA_dim_,svd_SVUXM_SMwl____,MDA_fname);
	MDA_ndim=0;
	MDA_dim_[MDA_ndim++] = FTK_n_svd_l;
	MDA_dim_[MDA_ndim++] = n_w_max;
	MDA_dim_[MDA_ndim++] = n_S;
	MDA_dim_[MDA_ndim++] = n_M_sub;
	MDA_n_dim = MDA_ndim;
	sprintf(MDA_fname,"mex_svd_SVUXM_lwSM____.mat");
	MDA_write_c16(MDA_n_dim,MDA_dim_,svd_SVUXM_lwSM____,MDA_fname);
	/* if (flag_dump){ } */}
      for (nSbatch=0;nSbatch<n_Sbatch;nSbatch++){
	n_S_sub = n_S_per_Sbatch; if (nSbatch==n_Sbatch-1){ n_S_sub = n_S - n_S_per_Sbatch*nSbatch;}
	for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ index_S_in_Sbatch_[nS_sub] = n_S_per_Sbatch*nSbatch + nS_sub;}
	if (verbose>1){ printf(" %% nSbatch %d/%d: index_S_in_Sbatch_ %d --> %d\n",nSbatch,n_Sbatch,index_S_in_Sbatch_[0],index_S_in_Sbatch_[n_S_sub-1]);}
	if (n_S_sub>0){
	  /* %%%% */
	  /* svd_SVUXM_lwsM____ = svd_SVUXM_lwSM____(:,:,1+index_S_in_Sbatch_,:); */
	  /* %%%% */	  
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
	  memset(svd_SVUXM_lwsM____,0,tab*sizeof(double complex));
	  for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	    for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
	      nS = index_S_in_Sbatch_[nS_sub];
	      tabA = (((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*(unsigned long long int)n_w_max)*(unsigned long long int)FTK_n_svd_l;
	      tabB = (((unsigned long long int)nS     + (unsigned long long int)nM_sub*(unsigned long long int)n_S    )*(unsigned long long int)n_w_max)*(unsigned long long int)FTK_n_svd_l;
	      tabC = (unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
	      memcpy( svd_SVUXM_lwsM____ + tabA , svd_SVUXM_lwSM____ + tabB , tabC*sizeof(double complex));
	      /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"svd_SVUXM_lwsM____: ");
	  /* %%%% */
	  if (flag_dump){
	    MDA_ndim=0;
	    MDA_dim_[MDA_ndim++] = FTK_n_svd_l;
	    MDA_dim_[MDA_ndim++] = n_w_max;
	    MDA_dim_[MDA_ndim++] = n_S_sub;
	    MDA_dim_[MDA_ndim++] = n_M_sub;
	    MDA_n_dim = MDA_ndim;
	    sprintf(MDA_fname,"mex_svd_SVUXM_lwsM____.mat");
	    MDA_write_c16(MDA_n_dim,MDA_dim_,svd_SVUXM_lwsM____,MDA_fname);
	    /* if (flag_dump){ } */}
	  /* %%%% */
	  /* svd_USESVUXM_dwSM____ = real(reshape(FTK.svd_U_d_expiw_s__*reshape(svd_SVUXM_lwsM____,[FTK.n_svd_l,n_w_max*n_S_sub*n_M_sub]),[FTK.n_delta_v,n_w_max,n_S_sub,n_M_sub])); */
	  /* %%%% */
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	  memset(svd_USESVUXM_dwSM____,0,tab*sizeof(double complex));
	  tabC = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max;
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v*(unsigned long long int)FTK_n_svd_l;
	  cblas_zgemm(
		      CblasColMajor
		      ,CblasNoTrans
		      ,CblasNoTrans
		      ,FTK_n_delta_v
		      ,tabC
		      ,FTK_n_svd_l
		      ,&cblas_alpha
		      ,FTK_svd_U_d_expiw_s__
		      ,FTK_n_delta_v
		      ,svd_SVUXM_lwsM____
		      ,FTK_n_svd_l
		      ,&cblas_beta
		      ,svd_USESVUXM_dwSM____
		      ,FTK_n_delta_v
		      );
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"svd_USESVUXM_dwSM____: ");
	  /* %%%% */
	  if (flag_dump){
	    MDA_ndim=0;
	    MDA_dim_[MDA_ndim++] = FTK_n_delta_v;
	    MDA_dim_[MDA_ndim++] = n_w_max;
	    MDA_dim_[MDA_ndim++] = n_S_sub;
	    MDA_dim_[MDA_ndim++] = n_M_sub;
	    MDA_n_dim = MDA_ndim;
	    sprintf(MDA_fname,"mex_svd_USESVUXM_dwSM____.mat");
	    MDA_write_c16(MDA_n_dim,MDA_dim_,svd_USESVUXM_dwSM____,MDA_fname);
	    /* if (flag_dump){ } */}
	  /* %%%% */
	  /* 
	     l2_dSM___ = permute(reshape(reshape(sqrt(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_)),[n_S_sub,1])*reshape(sqrt(UX_M_l2_dM_use__(:,1+index_M_in_Mbatch_)),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
	     n2_dSM___ = 1./max(1e-14,l2_dSM___);
	     f2_dSM___ = permute(reshape(reshape(sqrt(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_)),[n_S_sub,1])*reshape(1./max(1e-14,sqrt(UX_M_l2_dM_use__(:,1+index_M_in_Mbatch_))),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
	     ss_S_ = reshape(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_),[n_S_sub,1]);
	  */
	  /* %%%% */
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)FTK_n_delta_v;
	  memset(l2_dSM___,0,tab*sizeof(double));
	  memset(n2_dSM___,0,tab*sizeof(double));
	  memset(f2_dSM___,0,tab*sizeof(double));
	  for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
	    nS = index_S_in_Sbatch_[nS_sub];
	    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	      nM = index_M_in_Mbatch_[nM_sub];
	      for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){
		tabA = (unsigned long long int)ndelta_v + (unsigned long long int)nM*(unsigned long long int)FTK_n_delta_v;
		tabB = (unsigned long long int)ndelta_v + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*(unsigned long long int)FTK_n_delta_v;
		l2_dSM___[tabB] = sqrt(CTF_UX_S_l2_[nS]) * sqrt(UX_M_l2_dM__[tabA]);
		n2_dSM___[tabB] = 1.0 / maximum(1.0e-14,l2_dSM___[tabB]);
		f2_dSM___[tabB] = sqrt(CTF_UX_S_l2_[nS]) / maximum(1.0e-14,sqrt(UX_M_l2_dM__[tabA]));
		/* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */}
	      /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	    /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"l2_dSM___: ");
	  /* %%%% */
	  if (flag_dump){
	    MDA_ndim=0;
	    MDA_dim_[MDA_ndim++] = FTK_n_delta_v;
	    MDA_dim_[MDA_ndim++] = n_S_sub;
	    MDA_dim_[MDA_ndim++] = n_M_sub;
	    MDA_n_dim = MDA_ndim;
	    sprintf(MDA_fname,"mex_l2_dSM___.mat");
	    MDA_write_r8(MDA_n_dim,MDA_dim_,l2_dSM___,MDA_fname);
	    sprintf(MDA_fname,"mex_n2_dSM___.mat");
	    MDA_write_r8(MDA_n_dim,MDA_dim_,n2_dSM___,MDA_fname);
	    sprintf(MDA_fname,"mex_f2_dSM___.mat");
	    MDA_write_r8(MDA_n_dim,MDA_dim_,f2_dSM___,MDA_fname);
	    /* if (flag_dump){ } */}
	  /* %%%% */
	  /* X_sub_dwSM____ = repmat(reshape(n2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* svd_USESVUXM_dwSM____; %<-- correlation. ; */
	  /* %%%% */
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	  memset(X_sub_dwSM____,0,tab*sizeof(double));
	  na=0;
	  for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	    for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
	      for (nw=0;nw<n_w_max;nw++){
		for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){
		  tabA = (unsigned long long int)ndelta_v + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*FTK_n_delta_v;
		  X_sub_dwSM____[na] = n2_dSM___[tabA] * creal(svd_USESVUXM_dwSM____[na]) ;
		  na += 1;
		  /* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */}
		/* for (nw=0;nw<n_w_max;nw++){ } */}
	      /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	  /* %%%% */
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"X_sub_dwSM____: ");
	  /* %%%% */
	  if (flag_dump){
	    MDA_ndim=0;
	    MDA_dim_[MDA_ndim++] = FTK_n_delta_v;
	    MDA_dim_[MDA_ndim++] = n_w_max;
	    MDA_dim_[MDA_ndim++] = n_S_sub;
	    MDA_dim_[MDA_ndim++] = n_M_sub;
	    MDA_n_dim = MDA_ndim;
	    sprintf(MDA_fname,"mex_X_sub_dwSM____.mat");
	    MDA_write_r8(MDA_n_dim,MDA_dim_,X_sub_dwSM____,MDA_fname);
	    /* if (flag_dump){ } */}
	  /* %%%% */
	  /* 
	     I_value_sub_dwSM____ = repmat(reshape(f2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* X_sub_dwSM____; %<-- I_value. ;
	     I_value_use_dwSM____ = max(0,real(I_value_sub_dwSM____));
	  */
	  /* %%%% */
	  if (flag_compute_I_value){
	    local_tic(0,t_start_,d_start_);
	    tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	    memset(I_value_use_dwSM____,0,tab*sizeof(double));
	    na=0;
	    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	      for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
		for (nw=0;nw<n_w_max;nw++){
		  for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){
		    tabA = (unsigned long long int)ndelta_v + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*FTK_n_delta_v;
		    I_value_use_dwSM____[na] = maximum( 0.0 , f2_dSM___[tabA] * X_sub_dwSM____[na] ) ;
		    na += 1;
		    /* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */}
		  /* for (nw=0;nw<n_w_max;nw++){ } */}
	      /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	    local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"I_value_use_dwSM____: ");
	    /* %%%% */
	    if (flag_dump){
	      MDA_ndim=0;
	      MDA_dim_[MDA_ndim++] = FTK_n_delta_v;
	      MDA_dim_[MDA_ndim++] = n_w_max;
	      MDA_dim_[MDA_ndim++] = n_S_sub;
	      MDA_dim_[MDA_ndim++] = n_M_sub;
	      MDA_n_dim = MDA_ndim;
	      sprintf(MDA_fname,"mex_I_value_use_dwSM____.mat");
	      MDA_write_r8(MDA_n_dim,MDA_dim_,I_value_use_dwSM____,MDA_fname);
	      /* if (flag_dump){ } */}
	    /* if (flag_compute_I_value){ } */}
	  /* %%%% */
	  /*
	    if (flag_optimize_over_gamma_z == 0);
	    tmp_t = tic(); nop=0;
	    [tmp_X_wSM___,tmp_delta_ij___] = max(reshape(X_sub_dwSM____,[FTK.n_delta_v,n_w_max*n_S_sub*n_M_sub]),[],1); %<-- maximize correlation. ;
	    assert(min(tmp_delta_ij___)>=1); assert(max(tmp_delta_ij___)<=FTK.n_delta_v);
	    tmp_X_wSM___ = reshape(tmp_X_wSM___,[n_w_max,n_S_sub,n_M_sub]);
	    tmp_delta_ij___ = reshape(tmp_delta_ij___,[n_w_max,n_S_sub,n_M_sub]);
	    tmp_delta_x___ = FTK.delta_x_(tmp_delta_ij___);
	    tmp_delta_y___ = FTK.delta_y_(tmp_delta_ij___);
	    tmp_gamma_z___ = 2*pi*(0:n_w_max-1)/n_w_max;
	    X_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = tmp_X_wSM___;
	    delta_x_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_x___,[n_w_max,n_S_sub,n_M_sub]);
	    delta_y_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_y___,[n_w_max,n_S_sub,n_M_sub]);
	    gamma_z_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = repmat(tmp_gamma_z___(:),[1,n_S_sub,n_M_sub]);
	    if (flag_compute_I_value);
	    tmp_I_value_use_dwSM__ = reshape(I_value_use_dwSM____,[FTK.n_delta_v,n_w_max*n_S_sub*n_M_sub]);
	    tmp_I_value_use_wSM_ = zeros(n_w_max*n_S_sub*n_M_sub,1);
	    tmp_t2=tic();
	    for nl=0:n_w_max*n_S_sub*n_M_sub-1;
	    tmp_I_value_use_wSM_(1+nl) = tmp_I_value_use_dwSM__(tmp_delta_ij___(1+nl),1+nl);
	    end;%for nl=0:n_w_max*n_S_sub*n_M_sub-1;
	    I_value_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_I_value_use_wSM_,[n_w_max,n_S_sub,n_M_sub]);
	    tmp_t2 = toc(tmp_t2); if (verbose>1); disp(sprintf(' %% I_value_wSM___ %0.6fs',tmp_t2)); end;
	    end;%if (flag_compute_I_value);
	    tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wSM___: %0.6f',tmp_t)); end;
	    parameter = parameter_timing_update(parameter,'ampmh_X_wSM_mex___8: X_wSM___',tmp_t,1,nop);
	    end;%if (flag_optimize_over_gamma_z == 0);
	  */
	  /* %%%% */
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	  if (flag_optimize_over_gamma_z==0){
	    na=0;
	    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	      nM = index_M_in_Mbatch_[nM_sub];
	      for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
		nS = index_S_in_Sbatch_[nS_sub];
		for (nw=0;nw<n_w_max;nw++){
		  tabA = (unsigned long long int)nw + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*(unsigned long long int)n_w_max;
		  ndelta_v = dmax_index(FTK_n_delta_v,X_sub_dwSM____ + tabA*(unsigned long long int)FTK_n_delta_v);
		  tabB = (unsigned long long int)nw + ((unsigned long long int)nS     + (unsigned long long int)nM    *(unsigned long long int)n_S    )*(unsigned long long int)n_w_max;
		  tabC = ndelta_v + tabA*(unsigned long long int)FTK_n_delta_v;
		  X_wSM___[tabB] = X_sub_dwSM____[tabC];
		  delta_x_wSM___[tabB] = FTK_delta_x_[ndelta_v];
		  delta_y_wSM___[tabB] = FTK_delta_y_[ndelta_v];
		  gamma_z_wSM___[tabB] = gamma_z_[nw];
		  if (flag_compute_I_value){ I_value_wSM___[tabB] = I_value_use_dwSM____[tabC];}
		  na += FTK_n_delta_v;
		  /* for (nw=0;nw<n_w_max;nw++){ } */}
		/* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	      /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	    /* if (flag_optimize_over_gamma_z==0){ } */}
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"X_wSM___: ");
	  /* %%%% */
	  /*
	    if (flag_optimize_over_gamma_z == 1);
	    tmp_t = tic(); nop=0;
	    [tmp_X_SM__,tmp_dw_ij__] = max(reshape(X_sub_dwSM____,[FTK.n_delta_v*n_w_max,n_S_sub*n_M_sub]),[],1); %<-- maximize correlation. ;
	    [tmp_delta_ij__,tmp_gamma_ij__] = ind2sub([FTK.n_delta_v,n_w_max],tmp_dw_ij__);
	    assert(min(tmp_delta_ij__)>=1); assert(max(tmp_delta_ij__)<=FTK.n_delta_v);
	    assert(min(tmp_gamma_ij__)>=1); assert(max(tmp_gamma_ij__)<=n_w_max);
	    tmp_X_SM__ = reshape(tmp_X_SM__,[n_S_sub,n_M_sub]);
	    tmp_delta_ij__ = reshape(tmp_delta_ij__,[n_S_sub,n_M_sub]);
	    tmp_gamma_ij__ = reshape(tmp_gamma_ij__,[n_S_sub,n_M_sub]);
	    tmp_delta_x__ = FTK.delta_x_(tmp_delta_ij__);
	    tmp_delta_y__ = FTK.delta_y_(tmp_delta_ij__);
	    tmp_gamma_z__ = 2*pi*(tmp_gamma_ij__-1)/n_w_max;
	    X_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = tmp_X_SM__;
	    delta_x_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_x__,[n_S_sub,n_M_sub]);
	    delta_y_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_y__,[n_S_sub,n_M_sub]);
	    gamma_z_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_gamma_z__,[n_S_sub,n_M_sub]);
	    if (flag_compute_I_value);
	    tmp_I_value_use_dwSM___ = reshape(I_value_use_dwSM____,[FTK.n_delta_v,n_w_max,n_S_sub*n_M_sub]);
	    tmp_I_value_use_SM_ = zeros(n_S_sub*n_M_sub,1);
	    tmp_t2=tic();
	    for nl=0:n_S_sub*n_M_sub-1;
	    tmp_I_value_use_SM_(1+nl) = tmp_I_value_use_dwSM___(tmp_delta_ij__(1+nl),tmp_gamma_ij__(1+nl),1+nl);
	    end;%for nl=0:n_S_sub*n_M_sub-1;
	    I_value_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_I_value_use_SM_,[n_S_sub,n_M_sub]);
	    tmp_t2 = toc(tmp_t2); if (verbose>1); disp(sprintf(' %% I_value_SM__ %0.6fs',tmp_t2)); end;
	    end;%if (flag_compute_I_value);
	    tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wSM___: %0.6f',tmp_t)); end;
	    parameter = parameter_timing_update(parameter,'ampmh_X_wSM_mex___8: X_wSM___',tmp_t,1,nop);
	    end;%if (flag_optimize_over_gamma_z == 1);	    
	  */
	  /* %%%% */
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	  if (flag_optimize_over_gamma_z==1){
	    na=0;
	    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	      nM = index_M_in_Mbatch_[nM_sub];
	      for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
		nS = index_S_in_Sbatch_[nS_sub];
		tabA = (unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub;
		ndw = dmax_index((unsigned long long int)FTK_n_delta_v*(unsigned long long int)n_w_max,X_sub_dwSM____ + tabA*(unsigned long long int)FTK_n_delta_v*(unsigned long long int)n_w_max);
		nw = ndw / maximum(1,FTK_n_delta_v);
		ndelta_v = ndw % FTK_n_delta_v;
		tabB = (unsigned long long int)nS     + (unsigned long long int)nM    *(unsigned long long int)n_S    ;
		tabC = ndw + tabA*(unsigned long long int)FTK_n_delta_v*(unsigned long long int)n_w_max;
		X_wSM___[tabB] = X_sub_dwSM____[tabC];
		delta_x_wSM___[tabB] = FTK_delta_x_[ndelta_v];
		delta_y_wSM___[tabB] = FTK_delta_y_[ndelta_v];
		gamma_z_wSM___[tabB] = gamma_z_[nw];
		if (flag_compute_I_value){ I_value_wSM___[tabB] = I_value_use_dwSM____[tabC];}
		na += (unsigned long long int)FTK_n_delta_v*(unsigned long long int)n_w_max;
		/* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	      /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	    /* if (flag_optimize_over_gamma_z==1){ } */}
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"X_wSM___: ");
	  /* %%%% */
	  if (flag_disp){
	    printf(" %% nMbatch %d/%d nSbatch %d/%d\n",nMbatch,n_Mbatch,nSbatch,n_Sbatch);
	    if (flag_optimize_over_gamma_z==0){ tab = minimum(8,(unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max);}
	    if (flag_optimize_over_gamma_z==1){ tab = minimum(8,(unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub);}
	    array_printf(X_wSM___,"double",1,tab," %% X_wSM___: ");
	    array_printf(delta_x_wSM___,"double",1,tab," %% delta_x_wSM___: ");
	    array_printf(delta_y_wSM___,"double",1,tab," %% delta_y_wSM___: ");
	    array_printf(gamma_z_wSM___,"double",1,tab," %% gamma_z_wSM___: ");
	    if (flag_compute_I_value){ array_printf(I_value_wSM___,"double",1,tab," %% I_value_wSM___: ");}
	    /* if (flag_disp){ } */}
	  if (flag_dump){
	    MDA_ndim=0;
	    if (flag_optimize_over_gamma_z==0){ MDA_dim_[MDA_ndim++] = n_w_max;}
	    MDA_dim_[MDA_ndim++] = n_S;
	    MDA_dim_[MDA_ndim++] = n_M;
	    MDA_n_dim = MDA_ndim;
	    sprintf(MDA_fname,"mex_X_wSM___.mat");MDA_write_r8(MDA_n_dim,MDA_dim_,X_wSM___,MDA_fname);
	    sprintf(MDA_fname,"mex_delta_x_wSM___.mat");MDA_write_r8(MDA_n_dim,MDA_dim_,delta_x_wSM___,MDA_fname);
	    sprintf(MDA_fname,"mex_delta_y_wSM___.mat");MDA_write_r8(MDA_n_dim,MDA_dim_,delta_y_wSM___,MDA_fname);
	    sprintf(MDA_fname,"mex_gamma_z_wSM___.mat");MDA_write_r8(MDA_n_dim,MDA_dim_,gamma_z_wSM___,MDA_fname);
	    if (flag_compute_I_value){ sprintf(MDA_fname,"mex_I_value_wSM___.mat");MDA_write_r8(MDA_n_dim,MDA_dim_,I_value_wSM___,MDA_fname);}
	    /* if (flag_dump){ } */}
	  /* if (n_S_sub>0){ } */}
	/* for (nSbatch=0;nSbatch<n_Sbatch;nSbatch++){ } */}
      /* if (n_M_sub>0){ } */}
    /* for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){ } */}
  /* %%%%%%%%%%%%%%%% */
  free(gamma_z_); gamma_z_=NULL;
  if (flag_compute_I_value){ free(I_value_use_dwSM____); I_value_use_dwSM____=NULL;}
  free(X_sub_dwSM____); X_sub_dwSM____=NULL;
  free(l2_dSM___); l2_dSM___=NULL;
  free(n2_dSM___); n2_dSM___=NULL;
  free(f2_dSM___); f2_dSM___=NULL;
  free(svd_USESVUXM_dwSM____); svd_USESVUXM_dwSM____=NULL;
  free(svd_SVUXM_lwsM____); svd_SVUXM_lwsM____=NULL;
  free(svd_SVUXM_lwSM____); svd_SVUXM_lwSM____=NULL;
  fftw_free(svd_SVUXM_0in_wlSM____); svd_SVUXM_0in_wlSM____=NULL;
  fftw_free(svd_SVUXM_out_wlSM____); svd_SVUXM_out_wlSM____=NULL;
  fftw_destroy_plan(fftw_plan_many_plan);
  free(svd_SVUXM_SMwl____); svd_SVUXM_SMwl____=NULL;
  free(svd_VUXM_nMwl____); svd_VUXM_nMwl____=NULL;
  free(index_S_in_Sbatch_); index_S_in_Sbatch_=NULL;
  free(index_M_in_Mbatch_); index_M_in_Mbatch_=NULL;
  free(CTF_UX_S_k_q_nSw___); CTF_UX_S_k_q_nSw___=NULL;
  free(FTK_svd_U_d_expiw_s__); FTK_svd_U_d_expiw_s__=NULL;
  free(CTF_UX_S_k_q_wnS__); CTF_UX_S_k_q_wnS__=NULL;
  free(svd_VUXM_lwnM____); svd_VUXM_lwnM____=NULL;
  /* %%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___10]\n");}
}

#ifdef WITHOUT_MEX
int main()
{
  int verbose=0;
  int flag_test=1;
  /* 0in */
  int n_M_per_Mbatch=24;
  int n_S_per_Sbatch=24;
  int flag_optimize_over_gamma_z=0;
  int flag_compute_I_value=1;
  double tolerance_master=0.01;
  int FTK_n_svd_l=103;
  int FTK_n_delta_v=81;
  double *FTK_svd_U_d_expiw_s_real__=NULL,*FTK_svd_U_d_expiw_s_imag__=NULL;
  double *FTK_delta_x_=NULL;
  double *FTK_delta_y_=NULL;
  int n_w_max=192;
  int pm_n_UX_rank=60;
  int n_S=30;
  double *CTF_UX_S_k_q_wnS_real__=NULL,*CTF_UX_S_k_q_wnS_imag__=NULL;
  double *CTF_UX_S_l2_=NULL;
  int n_M=30;
  double *svd_VUXM_lwnM_real____=NULL,*svd_VUXM_lwnM_imag____=NULL;
  double *UX_M_l2_dM__=NULL;
  /* out */
  double *X_wSM___=NULL;
  double *delta_x_wSM___=NULL;
  double *delta_y_wSM___=NULL;
  double *gamma_z_wSM___=NULL;
  double *I_value_wSM___=NULL;
  /* %%%% */
  unsigned long long int na=0,n_a=0;
  /* %%%% */
  if (verbose){ printf(" %% [entering main]\n");}
  /* %%%%; */
  n_a = (unsigned long long int)FTK_n_delta_v*(unsigned long long int)FTK_n_svd_l;
  FTK_svd_U_d_expiw_s_real__ = (double *) calloc(n_a,sizeof(double));
  for (na=0;na<n_a;na++){ FTK_svd_U_d_expiw_s_real__[na] = (double) (na%3)-1;}
  FTK_svd_U_d_expiw_s_imag__ = (double *) calloc(n_a,sizeof(double));
  for (na=0;na<n_a;na++){ FTK_svd_U_d_expiw_s_imag__[na] = (double) (na%5)-2;}
  n_a = FTK_n_delta_v;
  FTK_delta_x_ = (double *) calloc(n_a,sizeof(double));
  for (na=0;na<n_a;na++){ FTK_delta_x_[na] = (double) (na%7)-3;}
  FTK_delta_y_ = (double *) calloc(n_a,sizeof(double));
  for (na=0;na<n_a;na++){ FTK_delta_y_[na] = (double) (na%9)-5;}
  n_a = (unsigned long long int)n_w_max*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_S;
  CTF_UX_S_k_q_wnS_real__ = (double *) calloc(n_a,sizeof(double));
  for (na=0;na<n_a;na++){ CTF_UX_S_k_q_wnS_real__[na] = (double) (na%3)-1;}
  CTF_UX_S_k_q_wnS_imag__ = (double *) calloc(n_a,sizeof(double));
  for (na=0;na<n_a;na++){ CTF_UX_S_k_q_wnS_imag__[na] = (double) (na%5)-2;}
  n_a = n_S;
  CTF_UX_S_l2_ = (double *) calloc(n_a,sizeof(double));
  for (na=0;na<n_a;na++){ CTF_UX_S_l2_[na] = (double) (1+na%3);}
  n_a = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_M;
  svd_VUXM_lwnM_real____ = (double *) calloc(n_a,sizeof(double));
  for (na=0;na<n_a;na++){ svd_VUXM_lwnM_real____[na] = (double) (na%3)-1;}
  svd_VUXM_lwnM_imag____ = (double *) calloc(n_a,sizeof(double));
  for (na=0;na<n_a;na++){ svd_VUXM_lwnM_imag____[na] = (double) (na%5)-2;}
  n_a = FTK_n_delta_v*n_M;
  UX_M_l2_dM__ = (double *) calloc(n_a,sizeof(double));
  for (na=0;na<n_a;na++){ UX_M_l2_dM__[na] = (double) (1+na%3);}
  /* %%%%; */
  n_a = (unsigned long long int)n_w_max*(unsigned long long int)n_S*(unsigned long long int)n_M; 
  if (flag_optimize_over_gamma_z){ n_a = (unsigned long long int)n_S*(unsigned long long int)n_M;}
  X_wSM___ = (double *) calloc(n_a,sizeof(double));
  delta_x_wSM___ = (double *) calloc(n_a,sizeof(double));
  delta_y_wSM___ = (double *) calloc(n_a,sizeof(double));
  gamma_z_wSM___ = (double *) calloc(n_a,sizeof(double));
  I_value_wSM___ = (double *) calloc(n_a,sizeof(double));  
  /* %%%%; */
  if (flag_test){
    printf(" %% calling MDA_io_test()\n"); MDA_io_test(); exit(0);
    printf(" %% calling test_cblas()\n"); test_cblas();
    printf(" %% calling test_transpose_float()\n"); test_transpose_float();
    printf(" %% calling dp_ps_mult_immintrin_test()\n"); dp_ps_mult_immintrin_test();
    printf(" %% calling hp_ps_mult_immintrin_test()\n"); hp_ps_mult_immintrin_test();
    /* if (flag_test){ } */}
  /* %%%%; */
  mex_ampmh_X_wSM___10
    (
      n_M_per_Mbatch
     ,n_S_per_Sbatch
     ,flag_optimize_over_gamma_z
     ,flag_compute_I_value
     ,tolerance_master
     ,FTK_n_svd_l
     ,FTK_n_delta_v
     ,FTK_svd_U_d_expiw_s_real__
     ,FTK_svd_U_d_expiw_s_imag__
     ,FTK_delta_x_
     ,FTK_delta_y_
     ,n_w_max
     ,pm_n_UX_rank
     ,n_S
     ,CTF_UX_S_k_q_wnS_real__
     ,CTF_UX_S_k_q_wnS_imag__
     ,CTF_UX_S_l2_
     ,n_M
     ,svd_VUXM_lwnM_real____
     ,svd_VUXM_lwnM_imag____
     ,UX_M_l2_dM__
     ,X_wSM___
     ,delta_x_wSM___
     ,delta_y_wSM___
     ,gamma_z_wSM___
     ,I_value_wSM___
     );
  /* %%%%; */
  free(FTK_svd_U_d_expiw_s_real__);FTK_svd_U_d_expiw_s_real__=NULL;
  free(FTK_svd_U_d_expiw_s_imag__);FTK_svd_U_d_expiw_s_imag__=NULL;
  free(FTK_delta_x_);FTK_delta_x_=NULL;
  free(FTK_delta_y_);FTK_delta_y_=NULL;
  free(CTF_UX_S_k_q_wnS_real__);CTF_UX_S_k_q_wnS_real__=NULL;
  free(CTF_UX_S_k_q_wnS_imag__);CTF_UX_S_k_q_wnS_imag__=NULL;
  free(CTF_UX_S_l2_);CTF_UX_S_l2_=NULL;
  free(svd_VUXM_lwnM_real____);svd_VUXM_lwnM_real____=NULL;
  free(svd_VUXM_lwnM_imag____);svd_VUXM_lwnM_imag____=NULL;
  free(UX_M_l2_dM__);UX_M_l2_dM__=NULL;
  free(X_wSM___);X_wSM___=NULL;
  free(delta_x_wSM___);delta_x_wSM___=NULL;
  free(delta_y_wSM___);delta_y_wSM___=NULL;
  free(gamma_z_wSM___);gamma_z_wSM___=NULL;
  free(I_value_wSM___);I_value_wSM___=NULL;
  /* %%%%; */
  if (verbose){ printf(" %% [finished main]\n");}
  return 0;
}
#endif

#ifndef WITHOUT_MEX
/* The gateway function */
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int verbose=0;
  int flag_test=0;
  /* 0in */
  int n_M_per_Mbatch=0;
  int n_S_per_Sbatch=0;
  int flag_optimize_over_gamma_z=0;
  int flag_compute_I_value=0;
  double tolerance_master=0;
  int FTK_n_svd_l=0;
  int FTK_n_delta_v=0;
  double *FTK_svd_U_d_expiw_s_real__=NULL,*FTK_svd_U_d_expiw_s_imag__=NULL;
  double *FTK_delta_x_=NULL;
  double *FTK_delta_y_=NULL;
  int n_w_max=0;
  int pm_n_UX_rank=0;
  int n_S=0;
  double *CTF_UX_S_k_q_wnS_real__=NULL,*CTF_UX_S_k_q_wnS_imag__=NULL;
  double *CTF_UX_S_l2_=NULL;
  int n_M=0;
  double *svd_VUXM_lwnM_real____=NULL,*svd_VUXM_lwnM_imag____=NULL;
  double *UX_M_l2_dM__=NULL;
  /* out */
  double *X_wSM___=NULL;
  double *delta_x_wSM___=NULL;
  double *delta_y_wSM___=NULL;
  double *gamma_z_wSM___=NULL;
  double *I_value_wSM___=NULL;
  /* tmp */
  int na=0,n_a;
  mwSize n_out=0;
  if (verbose){ printf(" %% [entering gateway mex_ampmh_X_wSM___10]\n");}
  /* %%%%%%%%%%%%%%%% */
  if(nrhs!=18) { mexErrMsgIdAndTxt("MyToolbox:mex_ampmh_X_wSM___10:nrhs","18 0in required.");}
  if(nlhs<  4) { mexErrMsgIdAndTxt("MyToolbox:mex_ampmh_X_wSM___10:nlhs"," 4 out required.");}
  /* %%%%%%%%%%%%%%%% */
  na=0;
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  n_M_per_Mbatch = (int)*mxGetPr(prhs[na]); na++;
  if (verbose){ printf(" %% n_M_per_Mbatch %d\n",n_M_per_Mbatch);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  n_S_per_Sbatch = (int)*mxGetPr(prhs[na]); na++;
  if (verbose){ printf(" %% n_S_per_Sbatch %d\n",n_S_per_Sbatch);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  flag_optimize_over_gamma_z = (int)*mxGetPr(prhs[na]); na++;
  if (verbose){ printf(" %% flag_optimize_over_gamma_z %d\n",flag_optimize_over_gamma_z);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  flag_compute_I_value = (int)*mxGetPr(prhs[na]); na++;
  if (verbose){ printf(" %% flag_compute_I_value %d\n",flag_compute_I_value);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  tolerance_master = (double)*mxGetPr(prhs[na]); na++;
  if (verbose){ printf(" %% tolerance_master %0.16f\n",tolerance_master);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  FTK_n_svd_l = (int)*mxGetPr(prhs[na]); na++;
  if (verbose){ printf(" %% FTK_n_svd_l %d\n",FTK_n_svd_l);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  FTK_n_delta_v = (int)*mxGetPr(prhs[na]); na++;
  if (verbose){ printf(" %% FTK_n_delta_v %d\n",FTK_n_delta_v);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  FTK_svd_U_d_expiw_s_real__ = (double *)mxGetPr(prhs[na]); FTK_svd_U_d_expiw_s_imag__ = (double *)mxGetPi(prhs[na]); na++;
  n_a = FTK_n_delta_v*FTK_n_svd_l;
  if (verbose){ printf(" %% FTK_svd_U_d_expiw_s__ %+0.16f %+0.16f*i --> %+0.16f %+0.16f*i\n",FTK_svd_U_d_expiw_s_real__[0],FTK_svd_U_d_expiw_s_imag__[0],FTK_svd_U_d_expiw_s_real__[n_a-1],FTK_svd_U_d_expiw_s_imag__[n_a-1]);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  FTK_delta_x_ = (double *)mxGetPr(prhs[na]); na++;
  if (verbose){ printf(" %% FTK_delta_x_ %0.16f --> %0.16f\n",FTK_delta_x_[0],FTK_delta_x_[FTK_n_delta_v-1]);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  FTK_delta_y_ = (double *)mxGetPr(prhs[na]); na++;
  if (verbose){ printf(" %% FTK_delta_y_ %0.16f --> %0.16f\n",FTK_delta_y_[0],FTK_delta_y_[FTK_n_delta_v-1]);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  n_w_max = (int)*mxGetPr(prhs[na]); na++;
  if (verbose){ printf(" %% n_w_max %d\n",n_w_max);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  pm_n_UX_rank = (int)*mxGetPr(prhs[na]); na++;
  if (verbose){ printf(" %% pm_n_UX_rank %d\n",pm_n_UX_rank);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  n_S = (int)*mxGetPr(prhs[na]); na++;
  if (verbose){ printf(" %% n_S %d\n",n_S);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  CTF_UX_S_k_q_wnS_real__ = (double *)mxGetPr(prhs[na]); CTF_UX_S_k_q_wnS_imag__ = (double *)mxGetPi(prhs[na]); na++;
  n_a = n_w_max*pm_n_UX_rank*n_S;
  if (verbose){ printf(" %% CTF_UX_S_k_q_wnS__ %+0.16f %+0.16f*i --> %+0.16f %+0.16f*i\n",CTF_UX_S_k_q_wnS_real__[0],CTF_UX_S_k_q_wnS_imag__[0],CTF_UX_S_k_q_wnS_real__[n_a-1],CTF_UX_S_k_q_wnS_imag__[n_a-1]);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  CTF_UX_S_l2_ = (double *)mxGetPr(prhs[na]); na++;
  if (verbose){ printf(" %% CTF_UX_S_l2_ %0.16f --> %0.16f\n",CTF_UX_S_l2_[0],CTF_UX_S_l2_[n_S-1]);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  n_M = (int)*mxGetPr(prhs[na]); na++;
  if (verbose){ printf(" %% n_M %d\n",n_M);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  svd_VUXM_lwnM_real____ = (double *)mxGetPr(prhs[na]); svd_VUXM_lwnM_imag____ = (double *)mxGetPi(prhs[na]); na++;
  n_a = FTK_n_svd_l*n_w_max*pm_n_UX_rank*n_M;
  if (verbose){ printf(" %% svd_VUXM_lwnM____ %+0.16f %+0.16f*i --> %+0.16f %+0.16f*i\n",svd_VUXM_lwnM_real____[0],svd_VUXM_lwnM_imag____[0],svd_VUXM_lwnM_real____[n_a-1],svd_VUXM_lwnM_imag____[n_a-1]);}
  /* %%%% */
  if (verbose){ printf(" %% input %d: numel (%d,%d) = %d\n",na,mxGetM(prhs[na]),mxGetN(prhs[na]),mxGetNumberOfElements(prhs[na]));}
  UX_M_l2_dM__ = (double *)mxGetPr(prhs[na]); na++;
  n_a = FTK_n_delta_v*n_M;
  if (verbose){ printf(" %% UX_M_l2_dM__ %0.16f --> %0.16f\n",UX_M_l2_dM__[0],UX_M_l2_dM__[n_a-1]);}
  /* %%%%%%%%%%%%%%%% */
  na=0;
  n_out = n_w_max*n_S*n_M; if (flag_optimize_over_gamma_z){ n_out = n_S*n_M;}
  plhs[na] = mxCreateDoubleMatrix(n_out,1,mxREAL); na++;
  plhs[na] = mxCreateDoubleMatrix(n_out,1,mxREAL); na++;
  plhs[na] = mxCreateDoubleMatrix(n_out,1,mxREAL); na++;
  plhs[na] = mxCreateDoubleMatrix(n_out,1,mxREAL); na++;
  if (nlhs==5){ plhs[na] = mxCreateDoubleMatrix(n_out,1,mxREAL); na++;}
  na=0;
  X_wSM___ = mxGetPr(plhs[na]); na++;
  delta_x_wSM___ = mxGetPr(plhs[na]); na++;
  delta_y_wSM___ = mxGetPr(plhs[na]); na++;
  gamma_z_wSM___ = mxGetPr(plhs[na]); na++;
  if (nlhs==5){ I_value_wSM___ = mxGetPr(plhs[na]); na++;}
  /* %%%%%%%%%%%%%%%% */
  if (flag_test){
    printf(" %% calling test_cblas()\n"); test_cblas();
    printf(" %% calling test_transpose_float()\n"); test_transpose_float();
    printf(" %% calling dp_ps_mult_immintrin_test()\n"); dp_ps_mult_immintrin_test();
    printf(" %% calling hp_ps_mult_immintrin_test()\n"); hp_ps_mult_immintrin_test();
    /* if (flag_test){ } */}
  /* %%%%%%%%%%%%%%%% */
  mex_ampmh_X_wSM___10
    (
      n_M_per_Mbatch
     ,n_S_per_Sbatch
     ,flag_optimize_over_gamma_z
     ,flag_compute_I_value
     ,tolerance_master
     ,FTK_n_svd_l
     ,FTK_n_delta_v
     ,FTK_svd_U_d_expiw_s_real__
     ,FTK_svd_U_d_expiw_s_imag__
     ,FTK_delta_x_
     ,FTK_delta_y_
     ,n_w_max
     ,pm_n_UX_rank
     ,n_S
     ,CTF_UX_S_k_q_wnS_real__
     ,CTF_UX_S_k_q_wnS_imag__
     ,CTF_UX_S_l2_
     ,n_M
     ,svd_VUXM_lwnM_real____
     ,svd_VUXM_lwnM_imag____
     ,UX_M_l2_dM__
     ,X_wSM___
     ,delta_x_wSM___
     ,delta_y_wSM___
     ,gamma_z_wSM___
     ,I_value_wSM___
     );
  /* %%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished gateway mex_ampmh_X_wSM___10]\n");}
}
#endif
