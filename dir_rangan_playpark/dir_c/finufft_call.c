#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

#ifndef SINGLE
int finufft_1d_test()
{
/* Simple example of calling the FINUFFT library from C, using C complex type,
   with a math test. Double-precision. Barnett 3/10/17. Opts control 6/19/18.

   Compile with:
   gcc -fopenmp example1d1c.c ../lib-static/libfinufft.a -o example1d1c -lfftw3 -lfftw3_omp -lm -lstdc++
   or if you have built a single-core version:
   gcc example1d1c.c ../lib-static/libfinufft.a -o example1d1c -lfftw3 -lm -lstdc++

   Usage: ./example1d1c
*/
  int M = 1e6;            // number of nonuniform points
  int N = 1e6;            // number of modes
  double tol = 1e-9;      // desired accuracy
  int j=0,ier=0,k=0,m=0,kout=0;
  double aF=0;

  // generate some random nonuniform points (x) and complex strengths (c):
  double* x = (double *)malloc1(sizeof(double)*M);
  double complex* c = (double complex*)malloc1(sizeof(double complex)*M);
  for (j=0; j<M; ++j) {
    x[j] = M_PI*(2*((double)rand()/RAND_MAX)-1);  // uniform random in [-pi,pi)
    c[j] = 2*((double)rand()/RAND_MAX)-1 + I*(2*((double)rand()/RAND_MAX)-1);
  }
  // allocate complex output array for the Fourier modes
  double complex* F = (double complex*)malloc1(sizeof(double complex)*N);

  nufft_opts opts;                      // opts struct (not ptr)
  finufft_default_opts(&opts);          // set default opts (must do this)
  opts.debug = 2;                       // show how to override a default
  //opts.upsampfac = 1.25;              // other opts...
  
  // call the NUFFT (with iflag=+1), passing pointers...
  ier = finufft1d1(M,x,c,+1,tol,N,F,&opts);

  k = 142519;            // check the answer just for this mode...
  assert(k>=-(double)N/2 && k<(double)N/2);
  double complex Ftest = 0.0 + 0.0*I;   // defined in complex.h (I too)
  for (j=0; j<M; ++j)
    Ftest += c[j] * cexp(I*(double)k*x[j]);
  double Fmax = 0.0;         // compute inf norm of F
  for (m=0; m<N; ++m) {
    aF = cabs(F[m]);
    if (aF>Fmax) Fmax=aF;
  }
  kout = k+N/2;          // index in output array for freq mode k
  double err = cabs(F[kout] - Ftest)/Fmax;
  printf("1D type 1 NUFFT done. ier=%d, err in F[%d] rel to max(F) is %.3g\n",ier,k,err);

  free1(&x);
  return ier;
 }
#endif /* SINGLE */

#ifdef SINGLE
int finufft_1d_test()
{
/* Simple example of calling the FINUFFT library from C, using C complex type,
   with a math test.
   Single-precision version (must be linked with single-precision libfinufft.a)
   Barnett 4/5/17. opts ctrl, t1 prefac convention, smaller prob size 9/14/18

   Compile with:
   gcc -fopenmp example1d1cf.c ../lib-static/libfinufft.a -o example1d1cf -lfftw3f -lfftw3f_omp -lm -lstdc++
   or if you have built a single-core version:
   gcc example1d1cf.c ../lib-static/libfinufft.a -o example1d1cf -lfftw3f -lm -lstdc++

   Usage: ./example1d1cf
*/
  int M = 1e5;            // number of nonuniform points
  int N = 1e5;            // number of modes (NB if too large lose acc in 1d)
  float tol = 1e-3;       // desired accuracy
  int j=0,ier=0,k=0,m=0,kout=0;
  float aF=0;

  // generate some random nonuniform points (x) and complex strengths (c):
  float* x = (float *)malloc1(sizeof(float)*M);
  float complex* c = (float complex*)malloc1(sizeof(float complex)*M);
  for (j=0; j<M; ++j) {
    x[j] = M_PI*(2*((float)rand()/RAND_MAX)-1);  // uniform random in [-pi,pi)
    c[j] = 2*((float)rand()/RAND_MAX)-1 + I*(2*((float)rand()/RAND_MAX)-1);
  }
  // allocate complex output array for the Fourier modes
  float complex* F = (float complex*)malloc1(sizeof(float complex)*N);

  nufft_opts opts;                         // opts struct (not ptr)
  finufftf_default_opts(&opts);            // set default opts (must do this)
  opts.debug = 2;                          // show how to override a default
  //opts.upsampfac = 1.25;                 // other opts...
  
  // call the NUFFT (with iflag=+1), passing pointers...
  ier = finufftf1d1(M,x,c,+1,tol,N,F,&opts);

  k = 14251;          // check the answer just for this mode...
  assert(k>=-(double)N/2 && k<(double)N/2);
  float complex Ftest = 0.0f + 0.0f*I;    // defined in complex.h (I too)
  for (j=0; j<M; ++j)
    Ftest += c[j] * cexpf(I*(float)k*x[j]);
  float Fmax = 0.0;       // compute inf norm of F
  for (m=0; m<N; ++m) {
    aF = cabsf(F[m]);
    if (aF>Fmax) Fmax=aF;
  }
  kout = k+N/2;       // index in output array for freq mode k
  float err = cabsf(F[kout] - Ftest)/Fmax;
  printf("1D type 1 NUFFT, single-prec. ier=%d, err in F[%d] rel to max(F) is %.3g\n",ier,k,err);

  free1(&x);
  return ier;
 }
#endif /* SINGLE */
