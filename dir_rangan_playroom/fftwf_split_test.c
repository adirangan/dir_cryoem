/* 
   gcc fftwf_split_test.c -lm -lfftw3 -lfftw3f ; ./a.out;
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <complex.h> /* <-- must include before fftw3.h to ensure that fftw_complex is compatible */
#include "fftw3.h"  /* <-- must include after complex.h to ensure that fftw_complex is compatible */

void FFTW1DSplit(const int N, double *XReal, double *XImag, double *YReal, double *YImag)
{
  fftw_plan Plan;
  fftw_iodim Dim;
  Dim.n = N;
  Dim.is = 1;
  Dim.os = 1;
  Plan = fftw_plan_guru_split_dft(1,&Dim,0,NULL,XReal,XImag,YReal,YImag,FFTW_ESTIMATE);
  fftw_execute_split_dft(Plan,XReal,XImag,YReal,YImag);
  fftw_execute_split_dft(Plan,XReal+N+7,XImag+N+7,YReal+N+7,YImag+N+7);
  fftw_destroy_plan(Plan);
}

void FFTWF1DSplit(const int N, float *XReal, float *XImag, float *YReal, float *YImag)
{
  fftwf_plan Plan;
  fftwf_iodim Dim;
  Dim.n = N;
  Dim.is = 1;
  Dim.os = 1;
  Plan = fftwf_plan_guru_split_dft(1,&Dim,0,NULL,XReal,XImag,YReal,YImag,FFTW_ESTIMATE);
  fftwf_execute_split_dft(Plan,XReal,XImag,YReal,YImag);
  fftwf_execute_split_dft(Plan,XReal+N+7,XImag+N+7,YReal+N+7,YImag+N+7);
  fftwf_destroy_plan(Plan);
}

int main()
{
  /* 
     compare output with: ;
     N = 8; n_ = [0:N-1]; x_ = (mod(n_,13)-4) + i*(mod(n_,15)-7); y_ = fft(x_); disp(y_);
     --> ;
  -4.0000 -28.0000i
 -13.6569 + 5.6569i
  -8.0000 + 0.0000i
  -5.6569 - 2.3431i
  -4.0000 - 4.0000i
  -2.3431 - 5.6569i
   0.0000 - 8.0000i
   5.6569 -13.6569i
  */
  int N = 8,n=0;
  double *d_xr=malloc(2*(N+7)*sizeof(double));
  double *d_xi=malloc(2*(N+7)*sizeof(double));
  double *d_yr=malloc(2*(N+7)*sizeof(double));
  double *d_yi=malloc(2*(N+7)*sizeof(double));
  memset(d_xr,0,2*(N+7)*sizeof(double));
  memset(d_xi,0,2*(N+7)*sizeof(double));
  memset(d_yr,0,2*(N+7)*sizeof(double));
  memset(d_yi,0,2*(N+7)*sizeof(double));
  for (n=0;n<N;n++){ d_xr[n] = n%13-4;}
  for (n=0;n<N;n++){ d_xi[n] = n%15-7;}
  for (n=0;n<N;n++){ d_xr[N+7+n] = n%13-4;}
  for (n=0;n<N;n++){ d_xi[N+7+n] = n%15-7;}
  FFTW1DSplit(N,d_xr,d_xi,d_yr,d_yi);
  printf(" d_xr: "); for (n=0;n<8;n++){ printf(" %+0.6f",d_xr[n]);} printf("\n");
  printf(" d_xi: "); for (n=0;n<8;n++){ printf(" %+0.6f",d_xi[n]);} printf("\n");
  printf(" d_yr: "); for (n=0;n<8;n++){ printf(" %+0.6f",d_yr[n]);} printf("\n");
  printf(" d_yi: "); for (n=0;n<8;n++){ printf(" %+0.6f",d_yi[n]);} printf("\n");
  printf(" d_xr: "); for (n=0;n<8;n++){ printf(" %+0.6f",d_xr[N+7+n]);} printf("\n");
  printf(" d_xi: "); for (n=0;n<8;n++){ printf(" %+0.6f",d_xi[N+7+n]);} printf("\n");
  printf(" d_yr: "); for (n=0;n<8;n++){ printf(" %+0.6f",d_yr[N+7+n]);} printf("\n");
  printf(" d_yi: "); for (n=0;n<8;n++){ printf(" %+0.6f",d_yi[N+7+n]);} printf("\n");
  float *f_xr=malloc(2*(N+7)*sizeof(float));
  float *f_xi=malloc(2*(N+7)*sizeof(float));
  float *f_yr=malloc(2*(N+7)*sizeof(float));
  float *f_yi=malloc(2*(N+7)*sizeof(float));
  memset(f_xr,0,2*(N+7)*sizeof(float));
  memset(f_xi,0,2*(N+7)*sizeof(float));
  memset(f_yr,0,2*(N+7)*sizeof(float));
  memset(f_yi,0,2*(N+7)*sizeof(float));
  for (n=0;n<N;n++){ f_xr[n] = n%13-4;}
  for (n=0;n<N;n++){ f_xi[n] = n%15-7;}
  for (n=0;n<N;n++){ f_xr[N+7+n] = n%13-4;}
  for (n=0;n<N;n++){ f_xi[N+7+n] = n%15-7;}
  FFTWF1DSplit(N,f_xr,f_xi,f_yr,f_yi);
  printf(" f_xr: "); for (n=0;n<8;n++){ printf(" %+0.6f",f_xr[n]);} printf("\n");
  printf(" f_xi: "); for (n=0;n<8;n++){ printf(" %+0.6f",f_xi[n]);} printf("\n");
  printf(" f_yr: "); for (n=0;n<8;n++){ printf(" %+0.6f",f_yr[n]);} printf("\n");
  printf(" f_yi: "); for (n=0;n<8;n++){ printf(" %+0.6f",f_yi[n]);} printf("\n");
  printf(" f_xr: "); for (n=0;n<8;n++){ printf(" %+0.6f",f_xr[N+7+n]);} printf("\n");
  printf(" f_xi: "); for (n=0;n<8;n++){ printf(" %+0.6f",f_xi[N+7+n]);} printf("\n");
  printf(" f_yr: "); for (n=0;n<8;n++){ printf(" %+0.6f",f_yr[N+7+n]);} printf("\n");
  printf(" f_yi: "); for (n=0;n<8;n++){ printf(" %+0.6f",f_yi[N+7+n]);} printf("\n");
  return 0;
}
