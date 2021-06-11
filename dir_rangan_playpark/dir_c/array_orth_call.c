#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

void array_gram_schmidt_inplace_f(int n_r,int n_c,float *Q_rc__)
{
  int nr0=0,nc0=0;
  int nr1=0,nc1=0;
  float *Q_c0_r_=NULL,*Q_c1_r_=NULL;
  float Q_c0c0=0,Q_c0_f=0;
  float Q_c0c1=0;
  for (nc0=0;nc0<n_c;nc0++){
    Q_c0_r_ = Q_rc__ + (unsigned long long int)nc0*(unsigned long long int)n_r;
    dp_ps_immintrin_loadu(n_r,Q_c0_r_,Q_c0_r_,&Q_c0c0); Q_c0_f = sqrtf(Q_c0c0);
    for (nr0=0;nr0<n_r;nr0++){ Q_c0_r_[nr0] /= maximum(1e-12,Q_c0_f);}
    for (nc1=nc0+1;nc1<n_c;nc1++){
      Q_c1_r_ = Q_rc__ + (unsigned long long int)nc1*(unsigned long long int)n_r;
      dp_ps_immintrin_loadu(n_r,Q_c0_r_,Q_c1_r_,&Q_c0c1);
      for (nr1=0;nr1<n_r;nr1++){ Q_c1_r_[nr1] -= Q_c0c1*Q_c0_r_[nr1];}
      /* for (nc1=nc0+1;nc1<n_c;nc1++){ } */}
    /* for (nc0=0;nc0<n_c;nc0++){ } */}
}

void array_orth_f(int n_r,int n_c,float **Q_rc_p_,unsigned long int *rseed_p)
{
  float *Q_rc__=NULL;
  unsigned long long int n_ulli = (unsigned long long int)n_r*(unsigned long long int)n_c;
  unsigned long long int nulli=0;
  Q_rc__=NULL;
  if (Q_rc_p_!=NULL){
    if ( (*Q_rc_p_)==NULL ){ (*Q_rc_p_) = (float *) malloc1(n_ulli*sizeof(float));}
    Q_rc__ = *Q_rc_p_;
    /* if (Q_rc_p_!=NULL){ } */}
  if (Q_rc__!=NULL){
    for (nulli=0;nulli<n_ulli;nulli++){ Q_rc__[nulli] = R01GET(rseed_p)-0.5;}
    array_gram_schmidt_inplace_f(n_r,n_c,Q_rc__);
    /* if (Q_rc__!=NULL){ } */}
}

void array_orth_f_test_error()
{
  int n_r=4;
  int n_c=3;
  float *Q_rc__=NULL;
  float *QQ_cc__=NULL;
  float QQ_cc_ans__[9] = { 1 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 1 };
  unsigned long int rseed=1;
  RSEED_adv8(&rseed);
  GLOBAL_tic(0);
  array_orth_f(n_r,n_c,&Q_rc__,&rseed);
  QQ_cc__ = (float *) malloc1((unsigned long long int)n_c*(unsigned long long int)n_c*sizeof(float));
  dp_ps_mult_immintrin_loadu(n_c,n_r,Q_rc__,n_c,Q_rc__,&QQ_cc__);
  array_printf(Q_rc__,"float",n_r,n_c," % Q_rc__: ");
  array_printf(QQ_cc__,"float",n_c,n_c," % QQ_cc__: ");
  printf(" %% QQ_cc_ans__ vs QQ_cc__: relative error %0.16f\n",ffnormn((unsigned long long int)n_c*(unsigned long long int)n_c,QQ_cc_ans__,QQ_cc__));
  free1(&QQ_cc__);
  free1(&Q_rc__);
  GLOBAL_toc(0,1," % array_orth_f: ");
}

void array_orth_f_test_speed()
{
  int n_r=10000;
  int n_c=90;
  int nc=0;
  float *Q_rc__=NULL;
  float *QQ_cc__=NULL;
  float *QQ_cc_ans__=NULL;
  unsigned long int rseed=1;
  RSEED_adv8(&rseed);
  GLOBAL_tic(0);
  GLOBAL_tic(1);
  array_orth_f(n_r,n_c,&Q_rc__,&rseed);
  GLOBAL_toc(1,1," % array_orth_f: ");
  printf(" %% Gops %0.6f\n",(double)n_c*(double)n_c*(double)n_r/(double)2/GLOBAL_elrt[1]/1e9);
  QQ_cc_ans__ = (float *) malloc1((unsigned long long int)n_c*(unsigned long long int)n_c*sizeof(float));
  memset(QQ_cc_ans__,0,(unsigned long long int)n_c*(unsigned long long int)n_c*sizeof(float));
  for (nc=0;nc<n_c;nc++){ QQ_cc_ans__[nc+nc*n_c]=1;}
  QQ_cc__ = (float *) malloc1((unsigned long long int)n_c*(unsigned long long int)n_c*sizeof(float));
  dp_ps_mult_immintrin_loadu(n_c,n_r,Q_rc__,n_c,Q_rc__,&QQ_cc__);
  printf(" %% QQ_cc_ans__ vs QQ_cc__: relative error %0.16f\n",ffnormn((unsigned long long int)n_c*(unsigned long long int)n_c,QQ_cc_ans__,QQ_cc__));
  free1(&QQ_cc_ans__);
  free1(&QQ_cc__);
  free1(&Q_rc__);
  GLOBAL_toc(0,1," % array_orth_f_test_speed: ");
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void array_gram_schmidt_inplace_d(int n_r,int n_c,double *Q_rc__)
{
  int nr0=0,nc0=0;
  int nr1=0,nc1=0;
  double *Q_c0_r_=NULL,*Q_c1_r_=NULL;
  double Q_c0c0=0,Q_c0_d=0;
  double Q_c0c1=0;
  for (nc0=0;nc0<n_c;nc0++){
    Q_c0_r_ = Q_rc__ + (unsigned long long int)nc0*(unsigned long long int)n_r;
    dp_pd_immintrin_loadu(n_r,Q_c0_r_,Q_c0_r_,&Q_c0c0); Q_c0_d = sqrt(Q_c0c0);
    for (nr0=0;nr0<n_r;nr0++){ Q_c0_r_[nr0] /= maximum(1e-12,Q_c0_d);}
    for (nc1=nc0+1;nc1<n_c;nc1++){
      Q_c1_r_ = Q_rc__ + (unsigned long long int)nc1*(unsigned long long int)n_r;
      dp_pd_immintrin_loadu(n_r,Q_c0_r_,Q_c1_r_,&Q_c0c1);
      for (nr1=0;nr1<n_r;nr1++){ Q_c1_r_[nr1] -= Q_c0c1*Q_c0_r_[nr1];}
      /* for (nc1=nc0+1;nc1<n_c;nc1++){ } */}
    /* for (nc0=0;nc0<n_c;nc0++){ } */}
}

void array_orth_d(int n_r,int n_c,double **Q_rc_p_,unsigned long int *rseed_p)
{
  double *Q_rc__=NULL;
  unsigned long long int n_ulli = (unsigned long long int)n_r*(unsigned long long int)n_c;
  unsigned long long int nulli=0;
  Q_rc__=NULL;
  if (Q_rc_p_!=NULL){
    if ( (*Q_rc_p_)==NULL ){ (*Q_rc_p_) = (double *) malloc1(n_ulli*sizeof(double));}
    Q_rc__ = *Q_rc_p_;
    /* if (Q_rc_p_!=NULL){ } */}
  if (Q_rc__!=NULL){
    for (nulli=0;nulli<n_ulli;nulli++){ Q_rc__[nulli] = R01GET(rseed_p)-0.5;}
    array_gram_schmidt_inplace_d(n_r,n_c,Q_rc__);
    /* if (Q_rc__!=NULL){ } */}
}

void array_orth_d_test_error()
{
  int n_r=4;
  int n_c=3;
  double *Q_rc__=NULL;
  double *QQ_cc__=NULL;
  double QQ_cc_ans__[9] = { 1 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 1 };
  unsigned long int rseed=1;
  RSEED_adv8(&rseed);
  GLOBAL_tic(0);
  array_orth_d(n_r,n_c,&Q_rc__,&rseed);
  QQ_cc__ = (double *) malloc1((unsigned long long int)n_c*(unsigned long long int)n_c*sizeof(double));
  dp_pd_mult_immintrin_loadu(n_c,n_r,Q_rc__,n_c,Q_rc__,&QQ_cc__);
  array_printf(Q_rc__,"double",n_r,n_c," % Q_rc__: ");
  array_printf(QQ_cc__,"double",n_c,n_c," % QQ_cc__: ");
  printf(" %% QQ_cc_ans__ vs QQ_cc__: relative error %0.16f\n",dfnormn((unsigned long long int)n_c*(unsigned long long int)n_c,QQ_cc_ans__,QQ_cc__));
  free1(&QQ_cc__);
  free1(&Q_rc__);
  GLOBAL_toc(0,1," % array_orth_d: ");
}

void array_orth_d_test_speed()
{
  int n_r=10000;
  int n_c=90;
  int nc=0;
  double *Q_rc__=NULL;
  double *QQ_cc__=NULL;
  double *QQ_cc_ans__=NULL;
  unsigned long int rseed=1;
  RSEED_adv8(&rseed);
  GLOBAL_tic(0);
  GLOBAL_tic(1);
  array_orth_d(n_r,n_c,&Q_rc__,&rseed);
  GLOBAL_toc(1,1," % array_orth_d: ");
  printf(" %% Gops %0.6f\n",(double)n_c*(double)n_c*(double)n_r/(double)2/GLOBAL_elrt[1]/1e9);
  QQ_cc_ans__ = (double *) malloc1((unsigned long long int)n_c*(unsigned long long int)n_c*sizeof(double));
  memset(QQ_cc_ans__,0,(unsigned long long int)n_c*(unsigned long long int)n_c*sizeof(double));
  for (nc=0;nc<n_c;nc++){ QQ_cc_ans__[nc+nc*n_c]=1;}
  QQ_cc__ = (double *) malloc1((unsigned long long int)n_c*(unsigned long long int)n_c*sizeof(double));
  dp_pd_mult_immintrin_loadu(n_c,n_r,Q_rc__,n_c,Q_rc__,&QQ_cc__);
  printf(" %% QQ_cc_ans__ vs QQ_cc__: relative error %0.16f\n",dfnormn((unsigned long long int)n_c*(unsigned long long int)n_c,QQ_cc_ans__,QQ_cc__));
  free1(&QQ_cc_ans__);
  free1(&QQ_cc__);
  free1(&Q_rc__);
  GLOBAL_toc(0,1," % array_orth_d_test_speed: ");
}
