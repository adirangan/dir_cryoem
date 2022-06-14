/*==========================================================
 * mex_ampmh_X_wSM___9.c ;
 * compile from Matlab with: ;
 * mex mex_ampmh_X_wSM___9.c ;
 * compile from shell with: ;
 * gcc -O2 mex_ampmh_X_wSM___9.c -lm -lgslcblas -DWITHOUT_MEX; ./a.out ;
 * trying to build a faster ampmh_X_wSM___8; not there yet! ;
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
/* #include "fftw3.h" */
/* #include <omp.h> */
#include <emmintrin.h>
#include <immintrin.h>
#include <complex.h>

double complex * double_complex_malloc_and_interleave(int n_a,double *double_real_,double *double_imag_)
{
  int na=0;
  double complex *double_complex_=NULL;
  double_complex_ = (double complex *) malloc(n_a*sizeof(double complex));
  for (na=0;na<n_a;na++){ double_complex_[na] = double_real_[na] + _Complex_I*double_imag_[na];}
  return double_complex_;
}

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

/* The computational routine */
void mex_ampmh_X_wSM___9
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
  int nl=0,nw=0,nS=0,nM=0;
  int pm_nUX_rank=0;
  unsigned long long int tab=0,tabA=0,tabB=0,tabC=0;
  double complex *CTF_UX_S_k_q_nSw___=NULL;
  double complex *svd_VUXM_nMwl____=NULL;
  double complex *svd_SVUXM_SMwl____=NULL;
  double complex cblas_alpha = (double complex) 1.0;
  double complex cblas_beta  = (double complex) 0.0;
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  /* %%%% */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___9]\n");}
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
  if (verbose){ printf(" %% FTK_svd_U_d_expiw_s__ %0.16f + %0.16f*i --> %0.16f + %0.16f*i\n",creal(FTK_svd_U_d_expiw_s__[0]),cimag(FTK_svd_U_d_expiw_s__[0]),creal(FTK_svd_U_d_expiw_s__[n_a-1]),cimag(FTK_svd_U_d_expiw_s__[n_a-1]));}
  if (verbose){ printf(" %% FTK_delta_x_ %0.16f --> %0.16f\n",FTK_delta_x_[0],FTK_delta_x_[FTK_n_delta_v-1]);}
  if (verbose){ printf(" %% FTK_delta_y_ %0.16f --> %0.16f\n",FTK_delta_y_[0],FTK_delta_y_[FTK_n_delta_v-1]);}
  if (verbose){ printf(" %% n_w_max %d\n",n_w_max);}
  if (verbose){ printf(" %% pm_n_UX_rank %d\n",pm_n_UX_rank);}
  if (verbose){ printf(" %% n_S %d\n",n_S);}
  n_a = n_w_max*pm_n_UX_rank*n_S;
  if (verbose){ printf(" %% CTF_UX_S_k_q_wnS__ %0.16f + %0.16f*i --> %0.16f + %0.16f*i\n",creal(CTF_UX_S_k_q_wnS__[0]),cimag(CTF_UX_S_k_q_wnS__[0]),creal(CTF_UX_S_k_q_wnS__[n_a-1]),cimag(CTF_UX_S_k_q_wnS__[n_a-1]));}
  if (verbose){ printf(" %% CTF_UX_S_l2_ %0.16f --> %0.16f\n",CTF_UX_S_l2_[0],CTF_UX_S_l2_[n_S-1]);}
  if (verbose){ printf(" %% n_M %d\n",n_M);}
  n_a = FTK_n_svd_l*n_w_max*pm_n_UX_rank*n_M;
  if (verbose){ printf(" %% svd_VUXM_lwnM____ %0.16f + %0.16f*i --> %0.16f + %0.16f*i\n",creal(svd_VUXM_lwnM____[0]),cimag(svd_VUXM_lwnM____[0]),creal(svd_VUXM_lwnM____[n_a-1]),cimag(svd_VUXM_lwnM____[n_a-1]));}
  if (verbose){ printf(" %% UX_M_l2_dM__ %0.16f --> %0.16f\n",UX_M_l2_dM__[0],UX_M_l2_dM__[FTK_n_delta_v*n_M-1]);}
  /* %%%%%%%%%%%%%%%% */
  local_tic(0,t_start_,d_start_);
  tab = (unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_S*(unsigned long long int)n_w_max;
  CTF_UX_S_k_q_nSw___ = (double complex *) malloc(tab*sizeof(double complex));
  na=0;
  for (nS=0;nS<n_S;nS++){
    for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){
      for (nw=0;nw<n_w_max;nw++){
	tabA =
	  (unsigned long long int)pm_nUX_rank
	  + (unsigned long long int)nS*(unsigned long long int)pm_n_UX_rank
	  + (unsigned long long int)nw*(unsigned long long int)n_S*(unsigned long long int)pm_n_UX_rank;
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
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)n_M_sub*(unsigned long long int)pm_n_UX_rank;
      na = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)(n_M_per_Mbatch*nMbatch);
      for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){
	  for (nw=0;nw<n_w_max;nw++){
	    for (nl=0;nl<FTK_n_svd_l;nl++){
	      tabA =
		(unsigned long long int)pm_nUX_rank
		+ (unsigned long long int)nM_sub*(unsigned long long int)pm_n_UX_rank
		+ (unsigned long long int)nw*(unsigned long long int)n_M_sub*(unsigned long long int)pm_n_UX_rank
		+ (unsigned long long int)nl*(unsigned long long int)n_w_max*(unsigned long long int)n_M_sub*(unsigned long long int)pm_n_UX_rank;
	      svd_VUXM_nMwl____[tabA] = svd_VUXM_lwnM____[na]; na++;
	      /* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
	    /* for (nw=0;nw<n_w_max;nw++){ } */}
	  /* for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){ } */}
	/* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"svd_VUXM_nMwl____: ");
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
      
      for (nSbatch=0;nSbatch<n_Sbatch;nSbatch++){
	n_S_sub = n_S_per_Sbatch; if (nSbatch==n_Sbatch-1){ n_S_sub = n_S - n_S_per_Sbatch*nSbatch;}
	for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ index_S_in_Sbatch_[nS_sub] = n_S_per_Sbatch*nSbatch + nS_sub;}
	if (verbose>1){ printf(" %% nSbatch %d/%d: index_S_in_Sbatch_ %d --> %d\n",nSbatch,n_Sbatch,index_S_in_Sbatch_[0],index_S_in_Sbatch_[n_S_sub-1]);}
	if (n_S_sub>0){
	  
	  /* if (n_S_sub>0){ } */}
	/* for (nSbatch=0;nSbatch<n_Sbatch;nSbatch++){ } */}
      /* if (n_M_sub>0){ } */}
    /* for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){ } */}
  /* %%%%%%%%%%%%%%%% */
  free(svd_SVUXM_SMwl____); svd_SVUXM_SMwl____=NULL;
  free(svd_VUXM_nMwl____); svd_VUXM_nMwl____=NULL;
  free(index_S_in_Sbatch_); index_S_in_Sbatch_=NULL;
  free(index_M_in_Mbatch_); index_M_in_Mbatch_=NULL;
  free(CTF_UX_S_k_q_nSw___); CTF_UX_S_k_q_nSw___=NULL;
  free(FTK_svd_U_d_expiw_s__); FTK_svd_U_d_expiw_s__=NULL;
  free(CTF_UX_S_k_q_wnS__); CTF_UX_S_k_q_wnS__=NULL;
  free(svd_VUXM_lwnM____); svd_VUXM_lwnM____=NULL;
  /* %%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___9]\n");}
}

#ifdef WITHOUT_MEX
int main()
{
  int verbose=1;
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
  int n_w_max=98;
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
  FTK_svd_U_d_expiw_s_imag__ = (double *) calloc(n_a,sizeof(double));
  n_a = FTK_n_delta_v;
  FTK_delta_x_ = (double *) calloc(n_a,sizeof(double));
  FTK_delta_y_ = (double *) calloc(n_a,sizeof(double));
  n_a = (unsigned long long int)n_w_max*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_S;
  CTF_UX_S_k_q_wnS_real__ = (double *) calloc(n_a,sizeof(double));
  CTF_UX_S_k_q_wnS_imag__ = (double *) calloc(n_a,sizeof(double));
  n_a = n_S;
  CTF_UX_S_l2_ = (double *) calloc(n_a,sizeof(double));
  n_a = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_M;
  svd_VUXM_lwnM_real____ = (double *) calloc(n_a,sizeof(double));
  svd_VUXM_lwnM_imag____ = (double *) calloc(n_a,sizeof(double));
  n_a = FTK_n_delta_v*n_M;
  UX_M_l2_dM__ = (double *) calloc(n_a,sizeof(double));
  /* %%%%; */
  n_a = (unsigned long long int)n_w_max*(unsigned long long int)n_S*(unsigned long long int)n_M; 
  if (flag_optimize_over_gamma_z){ n_a = (unsigned long long int)n_S*(unsigned long long int)n_M;}
  X_wSM___ = (double *) calloc(n_a,sizeof(double));
  delta_x_wSM___ = (double *) calloc(n_a,sizeof(double));
  delta_y_wSM___ = (double *) calloc(n_a,sizeof(double));
  gamma_z_wSM___ = (double *) calloc(n_a,sizeof(double));
  I_value_wSM___ = (double *) calloc(n_a,sizeof(double));  
  /* %%%%; */
  test_cblas(); exit(0);
  /* %%%%; */
  mex_ampmh_X_wSM___9
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
  int verbose=1;
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
  if (verbose){ printf(" %% [entering gateway mex_ampmh_X_wSM___9]\n");}
  /* %%%%%%%%%%%%%%%% */
  if(nrhs!=18) { mexErrMsgIdAndTxt("MyToolbox:mex_ampmh_X_wSM___9:nrhs","18 0in required.");}
  if(nlhs<  4) { mexErrMsgIdAndTxt("MyToolbox:mex_ampmh_X_wSM___9:nlhs"," 4 out required.");}
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
  if (verbose){ printf(" %% FTK_svd_U_d_expiw_s__ %0.16f + %0.16f*i --> %0.16f + %0.16f*i\n",FTK_svd_U_d_expiw_s_real__[0],FTK_svd_U_d_expiw_s_imag__[0],FTK_svd_U_d_expiw_s_real__[n_a-1],FTK_svd_U_d_expiw_s_imag__[n_a-1]);}
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
  if (verbose){ printf(" %% CTF_UX_S_k_q_wnS__ %0.16f + %0.16f*i --> %0.16f + %0.16f*i\n",CTF_UX_S_k_q_wnS_real__[0],CTF_UX_S_k_q_wnS_imag__[0],CTF_UX_S_k_q_wnS_real__[n_a-1],CTF_UX_S_k_q_wnS_imag__[n_a-1]);}
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
  if (verbose){ printf(" %% svd_VUXM_lwnM____ %0.16f + %0.16f*i --> %0.16f + %0.16f*i\n",svd_VUXM_lwnM_real____[0],svd_VUXM_lwnM_imag____[0],svd_VUXM_lwnM_real____[n_a-1],svd_VUXM_lwnM_imag____[n_a-1]);}
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
  mex_ampmh_X_wSM___9
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
  if (verbose){ printf(" %% [finished gateway mex_ampmh_X_wSM___9]\n");}
}
#endif
