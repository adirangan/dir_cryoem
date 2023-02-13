#ifndef _MONOLITH
#include "ampm_header.h"
#endif /* _MONOLITH */

void GLOBAL_tic(int nx)
{
  nx = maximum(0,minimum(GLOBAL_NTICKS-1,nx));
  GLOBAL_t_start[nx] = clock();
  gettimeofday(&GLOBAL_d_start[nx],NULL);
  GLOBAL_n_malloc1_[nx] = GLOBAL_n_malloc1;
}

void GLOBAL_toc(int nx,int verbose,const char *prefix)
{ 
  double r=0;
  nx = maximum(0,minimum(GLOBAL_NTICKS-1,nx)); 
  GLOBAL_t_final[nx] = clock(); gettimeofday(&GLOBAL_d_final[nx],NULL); 
  GLOBAL_l_ssec[nx] =  GLOBAL_d_final[nx].tv_sec -  GLOBAL_d_start[nx].tv_sec; 
  GLOBAL_l_usec[nx] = GLOBAL_d_final[nx].tv_usec - GLOBAL_d_start[nx].tv_usec; 
  GLOBAL_l_msec[nx] = ((GLOBAL_l_ssec[nx]*1000) + GLOBAL_l_usec[nx]/1000.0) + 0.5; 
  GLOBAL_elct[nx] = (double)(1000*(GLOBAL_t_final[nx]-GLOBAL_t_start[nx])/CLOCKS_PER_SEC)/(double)1000; 
  GLOBAL_elrt[nx] = (double)GLOBAL_l_msec[nx]/(double)1000; 
  r = GLOBAL_elct[nx]/GLOBAL_elrt[nx];
  if (verbose>=1){ 
    if (finite(r)){ printf("%sct/rt %0.3f/%0.3f = %.1f, diff_n_malloc %d\n",prefix,GLOBAL_elct[nx],GLOBAL_elrt[nx],r,GLOBAL_n_malloc1 - GLOBAL_n_malloc1_[nx]);}
    else{ printf("%sct/rt %0.3f/%0.3f = 0, diff_n_malloc %d\n",prefix,GLOBAL_elct[nx],GLOBAL_elrt[nx],GLOBAL_n_malloc1 - GLOBAL_n_malloc1_[nx]);}
    /* if (verbose>=1){ } */}
}
