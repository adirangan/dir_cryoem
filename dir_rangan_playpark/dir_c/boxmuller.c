#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

unsigned long int lrand()
{
  /* generates a unsigned long int (64 bit) */
  int verbose=0;
  int sr=(int)round(log(RAND_MAX)/log(2)),sl=8*sizeof(unsigned long int),sld=8*sizeof(long double);
  int nshift=sl/sr + (sl - sr*(sl/sr) > 0),ns=0;
  unsigned long int ld=0;
  if (verbose){ printf(" %% rand() size %d unsigned long int size %d, factor %d, long double size %d\n",sr,sl,nshift,sld);}
  ld=0;
  for (ns=0;ns<nshift;ns++){ 
    ld <<= sr; ld+=abs(rand()); 
    if (verbose){ printf(" %% ld set to %ld 2^(%0.1f) out of %d\n",ld,(double)logl((long double)ld)/log(2),sl);}
    /* for (ns=0;ns<nshift;ns++){ } */}
  return ld;
}

double randn()
{
  /* box muller polar form */
  double u=0,v=0,s=0;
  while (s==0 || s>1){ u=2*rand01-1;v=2*rand01-1;s=u*u+v*v;}
  return u*sqrt(-2*log(s)/s);
}

unsigned long int RGET(unsigned long int *rseed_p)
{
  /* basically:
     RNEXT = (RPREV*((unsigned long int)(pow(2,RPOW)+RADD))%((unsigned long int) pow(2,2*RPOW-1)));
     return RNEXT */
  *rseed_p = (*rseed_p*POW2RPOWPLUSRADD%POW22RPOWMINUSONE);
  return *rseed_p;
}

double R01GET(unsigned long int *rseed_p)
{
 /* basically:
    RNEXT = (RPREV*((unsigned long int)(pow(2,RPOW)+RADD))%((unsigned long int) pow(2,2*RPOW-1)));
    return = (double)RNEXT/(double)pow(2,2*RPOW-1); */
  int verbose=0;
  *rseed_p = (*rseed_p*POW2RPOWPLUSRADD%POW22RPOWMINUSONE);
  if (verbose){ printf(" %% %ld, %ld/%ld, %lf\n",*rseed_p,*rseed_p,POW22RPOWMINUSONE,(double)((long double)*rseed_p/(long double)POW22RPOWMINUSONE));}
  return (double)((long double)*rseed_p/(long double)POW22RPOWMINUSONE);
}

void RSEED_adv8(unsigned long int *rseed_p)
{
  int n_iteration = 8,niteration=0;
  for (niteration=0;niteration<n_iteration;niteration++){ RGET(rseed_p);}
}

double RNGET(unsigned long int *rseed_p)
{
  /* box muller polar form */
  double u=0,v=0,s=0;
  while (s==0 || s>1){ u=2*R01GET(rseed_p)-1;v=2*R01GET(rseed_p)-1;s=u*u+v*v;}
  return u*sqrt(-2*log(s)/s);
}

double RISIGET(unsigned long int *rseed_p,double rate)
{
  double r = R01GET(rseed_p);
  if (r==0){ r=1;}
  return -log(r)/rate;
}

void R01GET_test()
{
  unsigned long long int rseed=0;
  int n_iteration = 1024*8,niteration=0;
  double *d_out_=NULL;
  d_out_ = (double *) malloc1(n_iteration*sizeof(double));
  GLOBAL_tic(0);  
  rseed=3; RSEED_adv8(&rseed); for (niteration=0;niteration<n_iteration;niteration++){ d_out_[niteration] = R01GET(&rseed);}
  GLOBAL_toc(0,1," R01GET: ");
  array_printf(d_out_,"double",1,minimum(16,n_iteration)," % d_out_ start: ");
  printf(" %% rseed %lld\n",rseed);
  free1(&d_out_);
}
