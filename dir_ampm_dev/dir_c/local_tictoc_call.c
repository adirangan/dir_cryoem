#ifndef _MONOLITH
#include "ampm_header.h"
#endif /* _MONOLITH */

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
  if (elrt_[ntick]<=1e-16){ r = INFINITY; s = INFINITY;}
  if (verbose>=1){ 
    if (finite(r)){ printf("%sct/rt %0.3f/%0.3f = %.1f <-- %0.2f Mhz %0.2f Ghz\n",prefix,elct_[ntick],elrt_[ntick],r,s/1e6,s/1e9);}
    else{ printf("%sct/rt %0.3f/%0.3f = 0 <-- %0.2f Mhz %0.2f Ghz\n",prefix,elct_[ntick],elrt_[ntick],s/1e6,s/1e9);}
    /* if (verbose>=1){ } */}
}
