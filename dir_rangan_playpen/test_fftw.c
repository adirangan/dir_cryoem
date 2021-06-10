// gcc -o test_fftw.out test_fftw.c -I/usr/local/include/ -L/usr/local/lib -lfftw3 -lm ; ./test_fftw.out ;

#include <fftw3.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>

int main()
{
  int verbose=2;
  int N=8,n=0;
  double pi=0,theta=0;
  fftw_complex *in, *out;
  fftw_plan p;
  if (verbose){ printf(" %% [entering test_fftw], N %d\n",N);}
  pi = 4.0*atan(1.0);
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  for (n=0;n<N;n++){ theta = 2*pi*n/((double) N); in[n][0] = cos(theta); in[n][1] = sin(theta) ;}
  if (verbose>1){ printf(" %%%%  in: "); for (n=0;n<N;n++){ printf(" %+0.3f+i%+0.3f ",in[n][0],in[n][1]);} printf("\n");}
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p); /* repeat as needed */
  if (verbose>1){ printf(" %%%% out: "); for (n=0;n<N;n++){ printf(" %+0.3f+i%+0.3f ",out[n][0],out[n][1]);} printf("\n");}
  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
  if (verbose){ printf(" %% [finished test_fftw], N %d\n",N);}
  return 0;
}
