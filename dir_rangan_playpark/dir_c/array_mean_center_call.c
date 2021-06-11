#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

void array_mean_center_row_f(int n_r,int n_c,float *f_rc_0in__,float *f_cr_0in__,float **f_rc_out_p_,float **f_cr_out_p_)
{
  int nr=0,nc=0;
  float *f_rc_out__=NULL;
  float *f_cr_out__=NULL;
  float *f_r_=NULL;
  double mean=0,stdev=0;
  f_rc_out__=NULL;
  if (f_rc_out_p_!=NULL){
    if ( (*f_rc_out_p_)==NULL ){ (*f_rc_out_p_) = (float *) malloc1((unsigned long long int)n_r*(unsigned long long int)n_c*sizeof(float));}
    f_rc_out__ = *f_rc_out_p_;
    /* if (f_rc_out_p_!=NULL){ } */}
  f_cr_out__=NULL;
  if (f_cr_out_p_!=NULL){
    if ( (*f_cr_out_p_)==NULL ){ (*f_cr_out_p_) = (float *) malloc1((unsigned long long int)n_c*(unsigned long long int)n_r*sizeof(float));}
    f_cr_out__ = *f_cr_out_p_;
    /* if (f_rc_out_p_!=NULL){ } */}
  if (f_rc_out__!=NULL){
    for (nc=0;nc<n_c;nc++){
      f_r_ = f_rc_0in__ + (unsigned long long int)nc*(unsigned long long int)n_r;
      array_stats(f_r_,"float",n_r,NULL,NULL,&mean,&stdev);
      for (nr=0;nr<n_r;nr++){
	f_rc_out__[nr+nc*n_r] = (f_r_[nr] - mean);
	if (f_cr_out__!=NULL){ f_cr_out__[nc+nr*n_c] = f_rc_out__[nr+nc*n_r];}
	/* for (nr=0;nr<n_r;nr++){ } */}
      /* for (nc=0;nc<n_c;nc++){ } */}
    /* if (f_rc_out__!=NULL){ } */}
}

void array_mean_center_row_d(int n_r,int n_c,double *d_rc_0in__,double *d_cr_0in__,double **d_rc_out_p_,double **d_cr_out_p_)
{
  int nr=0,nc=0;
  double *d_rc_out__=NULL;
  double *d_cr_out__=NULL;
  double *d_r_=NULL;
  double mean=0,stdev=0;
  d_rc_out__=NULL;
  if (d_rc_out_p_!=NULL){
    if ( (*d_rc_out_p_)==NULL ){ (*d_rc_out_p_) = (double *) malloc1((unsigned long long int)n_r*(unsigned long long int)n_c*sizeof(double));}
    d_rc_out__ = *d_rc_out_p_;
    /* if (d_rc_out_p_!=NULL){ } */}
  d_cr_out__=NULL;
  if (d_cr_out_p_!=NULL){
    if ( (*d_cr_out_p_)==NULL ){ (*d_cr_out_p_) = (double *) malloc1((unsigned long long int)n_c*(unsigned long long int)n_r*sizeof(double));}
    d_cr_out__ = *d_cr_out_p_;
    /* if (d_cr_out_p_!=NULL){ } */}
  if (d_rc_out__!=NULL){
    for (nc=0;nc<n_c;nc++){
      d_r_ = d_rc_0in__ + (unsigned long long int)nc*(unsigned long long int)n_r;
      array_stats(d_r_,"double",n_r,NULL,NULL,&mean,&stdev);
      for (nr=0;nr<n_r;nr++){
	d_rc_out__[nr+nc*n_r] = (d_r_[nr] - mean);
	if (d_cr_out__!=NULL){ d_cr_out__[nc+nr*n_c] = d_rc_out__[nr+nc*n_r];}
	/* for (nr=0;nr<n_r;nr++){ } */}
      /* for (nc=0;nc<n_c;nc++){ } */}
    /* if (d_rc_out__!=NULL){ } */}
}

void array_mean_center_row(int n_r,int n_c,void *v_rc_0in__,void *v_cr_0in__,char *type,void *v_rc_out_p_,void *v_cr_out_p_)
{
  if (strcmp(type,"float")==0){ array_mean_center_row_f(n_r,n_c,(float *)v_rc_0in__,(float *)v_cr_0in__,(float **)v_rc_out_p_,(float **)v_cr_out_p_);}
  if (strcmp(type,"double")==0){ array_mean_center_row_d(n_r,n_c,(double *)v_rc_0in__,(double *)v_cr_0in__,(double **)v_rc_out_p_,(double **)v_cr_out_p_);}
}

void array_mean_center_row_test()
{
  int n_r = 4;
  int n_c = 3;
  float *f_A_rc__=NULL;
  float *f_A_cr__=NULL;
  double *d_A_rc__=NULL;
  double *d_A_cr__=NULL;
  float *f_B_rc__=NULL;
  float *f_B_cr__=NULL;
  double *d_B_rc__=NULL;
  double *d_B_cr__=NULL;
  float f_B_rc_ans__[12] = { -2 , -2 , +2 , +2 , -2 , -2 , +2 , +2 , -2 , -2 , +2 , +2 };
  float f_B_cr_ans__[12] = { -2 , -2 , -2 , -2 , -2 , -2 , +2 , +2 , +2 , +2 , +2 , +2 };
  double d_B_rc_ans__[12] = { -2 , -2 , +2 , +2 , -2 , -2 , +2 , +2 , -2 , -2 , +2 , +2 };
  double d_B_cr_ans__[12] = { -2 , -2 , -2 , -2 , -2 , -2 , +2 , +2 , +2 , +2 , +2 , +2 };
  int nr=0,nc=0;
  f_A_rc__ = (float *) malloc1(n_r*n_c*sizeof(float));
  f_A_cr__ = (float *) malloc1(n_c*n_r*sizeof(float));
  f_A_rc__[0] = -2; f_A_rc__[1] = -2; f_A_rc__[2] = +2; f_A_rc__[3] = +2;
  for (nr=0;nr<n_r;nr++){ for (nc=1;nc<n_c;nc++){ f_A_rc__[nr+nc*n_r] = f_A_rc__[nr+0*n_r] + nc;}}
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ f_A_cr__[nc+nr*n_c] = f_A_rc__[nr+nc*n_r];}}
  array_printf(f_A_rc__,"float",n_r,n_c," % f_A_rc__: ");
  array_printf(f_A_cr__,"float",n_c,n_r," % f_A_cr__: ");
  d_A_rc__ = (double *) malloc1(n_r*n_c*sizeof(double));
  d_A_cr__ = (double *) malloc1(n_c*n_r*sizeof(double));
  d_A_rc__[0] = -2; d_A_rc__[1] = -2; d_A_rc__[2] = +2; d_A_rc__[3] = +2;
  for (nr=0;nr<n_r;nr++){ for (nc=1;nc<n_c;nc++){ d_A_rc__[nr+nc*n_r] = d_A_rc__[nr+0*n_r] + nc;}}
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ d_A_cr__[nc+nr*n_c] = d_A_rc__[nr+nc*n_r];}}
  array_printf(d_A_rc__,"double",n_r,n_c," % d_A_rc__: ");
  array_printf(d_A_cr__,"double",n_c,n_r," % d_A_cr__: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_mean_center_row(n_r,n_c,f_A_rc__,f_A_cr__,"float",&f_B_rc__,&f_B_cr__);
  array_printf(f_B_rc__,"float",n_r,n_c," % f_B_rc__: ");
  array_printf(f_B_cr__,"float",n_c,n_r," % f_B_cr__: ");
  printf(" %% f_B_rc_ans__ vs f_B_rc__: relative error %0.16f\n",ffnormn(n_r*n_c,f_B_rc_ans__,f_B_rc__));
  printf(" %% f_B_cr_ans__ vs f_B_cr__: relative error %0.16f\n",ffnormn(n_r*n_c,f_B_cr_ans__,f_B_cr__));
  free1(&f_B_rc__);
  free1(&f_B_cr__);
  GLOBAL_toc(0,1," array_extract: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_mean_center_row(n_r,n_c,d_A_rc__,d_A_cr__,"double",&d_B_rc__,&d_B_cr__);
  array_printf(d_B_rc__,"double",n_r,n_c," % d_B_rc__: ");
  array_printf(d_B_cr__,"double",n_c,n_r," % d_B_cr__: ");
  printf(" %% d_B_rc_ans__ vs d_B_rc__: relative error %0.16f\n",dfnormn(n_r*n_c,d_B_rc_ans__,d_B_rc__));
  printf(" %% d_B_cr_ans__ vs d_B_cr__: relative error %0.16f\n",dfnormn(n_r*n_c,d_B_cr_ans__,d_B_cr__));
  free1(&d_B_rc__);
  free1(&d_B_cr__);
  GLOBAL_toc(0,1," array_extract: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_mean_center_row(n_r,n_c,f_A_rc__,f_A_cr__,"float",&f_A_rc__,&f_A_cr__);
  array_printf(f_A_rc__,"float",n_r,n_c," % f_A_rc__: ");
  array_printf(f_A_cr__,"float",n_c,n_r," % f_A_cr__: ");
  printf(" %% f_B_rc_ans__ vs f_A_rc__: relative error %0.16f\n",ffnormn(n_r*n_c,f_B_rc_ans__,f_A_rc__));
  printf(" %% f_B_cr_ans__ vs f_A_cr__: relative error %0.16f\n",ffnormn(n_r*n_c,f_B_cr_ans__,f_A_cr__));
  GLOBAL_toc(0,1," array_extract: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_mean_center_row(n_r,n_c,d_A_rc__,d_A_cr__,"double",&d_A_rc__,&d_A_cr__);
  array_printf(d_A_rc__,"double",n_r,n_c," % d_A_rc__: ");
  array_printf(d_A_cr__,"double",n_c,n_r," % d_A_cr__: ");
  printf(" %% d_B_rc_ans__ vs d_A_rc__: relative error %0.16f\n",dfnormn(n_r*n_c,d_B_rc_ans__,d_A_rc__));
  printf(" %% d_B_cr_ans__ vs d_A_cr__: relative error %0.16f\n",dfnormn(n_r*n_c,d_B_cr_ans__,d_A_cr__));
  GLOBAL_toc(0,1," array_extract: ");
  /* %%%%%%%% */
  free1(&f_A_rc__);
  free1(&f_A_cr__);
  free1(&d_A_rc__);  
  free1(&d_A_cr__);  
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void array_normalize_row_f(int n_r,int n_c,float *f_rc_0in__,float *f_cr_0in__,float **f_rc_out_p_,float **f_cr_out_p_)
{
  int nr=0,nc=0;
  float *f_rc_out__=NULL;
  float *f_cr_out__=NULL;
  float *f_r_=NULL;
  double mean=0,stdev=0;
  f_rc_out__=NULL;
  if (f_rc_out_p_!=NULL){
    if ( (*f_rc_out_p_)==NULL ){ (*f_rc_out_p_) = (float *) malloc1((unsigned long long int)n_r*(unsigned long long int)n_c*sizeof(float));}
    f_rc_out__ = *f_rc_out_p_;
    /* if (f_rc_out_p_!=NULL){ } */}
  f_cr_out__=NULL;
  if (f_cr_out_p_!=NULL){
    if ( (*f_cr_out_p_)==NULL ){ (*f_cr_out_p_) = (float *) malloc1((unsigned long long int)n_c*(unsigned long long int)n_r*sizeof(float));}
    f_cr_out__ = *f_cr_out_p_;
    /* if (f_rc_out_p_!=NULL){ } */}
  if (f_rc_out__!=NULL){
    for (nc=0;nc<n_c;nc++){
      f_r_ = f_rc_0in__ + (unsigned long long int)nc*(unsigned long long int)n_r;
      array_stats(f_r_,"float",n_r,NULL,NULL,&mean,&stdev);
      for (nr=0;nr<n_r;nr++){
	f_rc_out__[nr+nc*n_r] = (f_r_[nr] - mean)/maximum(1e-12,stdev);
	if (f_cr_out__!=NULL){ f_cr_out__[nc+nr*n_c] = f_rc_out__[nr+nc*n_r];}
	/* for (nr=0;nr<n_r;nr++){ } */}
      /* for (nc=0;nc<n_c;nc++){ } */}
    /* if (f_rc_out__!=NULL){ } */}
}

void array_normalize_row_d(int n_r,int n_c,double *d_rc_0in__,double *d_cr_0in__,double **d_rc_out_p_,double **d_cr_out_p_)
{
  int nr=0,nc=0;
  double *d_rc_out__=NULL;
  double *d_cr_out__=NULL;
  double *d_r_=NULL;
  double mean=0,stdev=0;
  d_rc_out__=NULL;
  if (d_rc_out_p_!=NULL){
    if ( (*d_rc_out_p_)==NULL ){ (*d_rc_out_p_) = (double *) malloc1((unsigned long long int)n_r*(unsigned long long int)n_c*sizeof(double));}
    d_rc_out__ = *d_rc_out_p_;
    /* if (d_rc_out_p_!=NULL){ } */}
  d_cr_out__=NULL;
  if (d_cr_out_p_!=NULL){
    if ( (*d_cr_out_p_)==NULL ){ (*d_cr_out_p_) = (double *) malloc1((unsigned long long int)n_c*(unsigned long long int)n_r*sizeof(double));}
    d_cr_out__ = *d_cr_out_p_;
    /* if (d_cr_out_p_!=NULL){ } */}
  if (d_rc_out__!=NULL){
    for (nc=0;nc<n_c;nc++){
      d_r_ = d_rc_0in__ + (unsigned long long int)nc*(unsigned long long int)n_r;
      array_stats(d_r_,"double",n_r,NULL,NULL,&mean,&stdev);
      for (nr=0;nr<n_r;nr++){
	d_rc_out__[nr+nc*n_r] = (d_r_[nr] - mean)/maximum(1e-12,stdev);
	if (d_cr_out__!=NULL){ d_cr_out__[nc+nr*n_c] = d_rc_out__[nr+nc*n_r];}
	/* for (nr=0;nr<n_r;nr++){ } */}
      /* for (nc=0;nc<n_c;nc++){ } */}
    /* if (d_rc_out__!=NULL){ } */}
}

void array_normalize_row(int n_r,int n_c,void *v_rc_0in__,void *v_cr_0in__,char *type,void *v_rc_out_p_,void *v_cr_out_p_)
{
  if (strcmp(type,"float")==0){ array_normalize_row_f(n_r,n_c,(float *)v_rc_0in__,(float *)v_cr_0in__,(float **)v_rc_out_p_,(float **)v_cr_out_p_);}
  if (strcmp(type,"double")==0){ array_normalize_row_d(n_r,n_c,(double *)v_rc_0in__,(double *)v_cr_0in__,(double **)v_rc_out_p_,(double **)v_cr_out_p_);}
}

void array_normalize_row_test()
{
  int n_r = 4;
  int n_c = 3;
  float *f_A_rc__=NULL;
  float *f_A_cr__=NULL;
  double *d_A_rc__=NULL;
  double *d_A_cr__=NULL;
  float *f_B_rc__=NULL;
  float *f_B_cr__=NULL;
  double *d_B_rc__=NULL;
  double *d_B_cr__=NULL;
  float f_B_rc_ans__[12] = { -1 , -1 , +1 , +1 , -1 , -1 , +1 , +1 , -1 , -1 , +1 , +1 };
  float f_B_cr_ans__[12] = { -1 , -1 , -1 , -1 , -1 , -1 , +1 , +1 , +1 , +1 , +1 , +1 };
  double d_B_rc_ans__[12] = { -1 , -1 , +1 , +1 , -1 , -1 , +1 , +1 , -1 , -1 , +1 , +1 };
  double d_B_cr_ans__[12] = { -1 , -1 , -1 , -1 , -1 , -1 , +1 , +1 , +1 , +1 , +1 , +1 };
  int nr=0,nc=0;
  f_A_rc__ = (float *) malloc1(n_r*n_c*sizeof(float));
  f_A_cr__ = (float *) malloc1(n_c*n_r*sizeof(float));
  f_A_rc__[0] = -2; f_A_rc__[1] = -2; f_A_rc__[2] = +2; f_A_rc__[3] = +2;
  for (nr=0;nr<n_r;nr++){ for (nc=1;nc<n_c;nc++){ f_A_rc__[nr+nc*n_r] = f_A_rc__[nr+0*n_r] + nc;}}
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ f_A_cr__[nc+nr*n_c] = f_A_rc__[nr+nc*n_r];}}
  array_printf(f_A_rc__,"float",n_r,n_c," % f_A_rc__: ");
  array_printf(f_A_cr__,"float",n_c,n_r," % f_A_cr__: ");
  d_A_rc__ = (double *) malloc1(n_r*n_c*sizeof(double));
  d_A_cr__ = (double *) malloc1(n_c*n_r*sizeof(double));
  d_A_rc__[0] = -2; d_A_rc__[1] = -2; d_A_rc__[2] = +2; d_A_rc__[3] = +2;
  for (nr=0;nr<n_r;nr++){ for (nc=1;nc<n_c;nc++){ d_A_rc__[nr+nc*n_r] = d_A_rc__[nr+0*n_r] + nc;}}
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ d_A_cr__[nc+nr*n_c] = d_A_rc__[nr+nc*n_r];}}
  array_printf(d_A_rc__,"double",n_r,n_c," % d_A_rc__: ");
  array_printf(d_A_cr__,"double",n_c,n_r," % d_A_cr__: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_normalize_row(n_r,n_c,f_A_rc__,f_A_cr__,"float",&f_B_rc__,&f_B_cr__);
  array_printf(f_B_rc__,"float",n_r,n_c," % f_B_rc__: ");
  array_printf(f_B_cr__,"float",n_c,n_r," % f_B_cr__: ");
  printf(" %% f_B_rc_ans__ vs f_B_rc__: relative error %0.16f\n",ffnormn(n_r*n_c,f_B_rc_ans__,f_B_rc__));
  printf(" %% f_B_cr_ans__ vs f_B_cr__: relative error %0.16f\n",ffnormn(n_r*n_c,f_B_cr_ans__,f_B_cr__));
  free1(&f_B_rc__);
  free1(&f_B_cr__);
  GLOBAL_toc(0,1," array_extract: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_normalize_row(n_r,n_c,d_A_rc__,d_A_cr__,"double",&d_B_rc__,&d_B_cr__);
  array_printf(d_B_rc__,"double",n_r,n_c," % d_B_rc__: ");
  array_printf(d_B_cr__,"double",n_c,n_r," % d_B_cr__: ");
  printf(" %% d_B_rc_ans__ vs d_B_rc__: relative error %0.16f\n",dfnormn(n_r*n_c,d_B_rc_ans__,d_B_rc__));
  printf(" %% d_B_cr_ans__ vs d_B_cr__: relative error %0.16f\n",dfnormn(n_r*n_c,d_B_cr_ans__,d_B_cr__));
  free1(&d_B_rc__);
  free1(&d_B_cr__);
  GLOBAL_toc(0,1," array_extract: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_normalize_row(n_r,n_c,f_A_rc__,f_A_cr__,"float",&f_A_rc__,&f_A_cr__);
  array_printf(f_A_rc__,"float",n_r,n_c," % f_A_rc__: ");
  array_printf(f_A_cr__,"float",n_c,n_r," % f_A_cr__: ");
  printf(" %% f_B_rc_ans__ vs f_A_rc__: relative error %0.16f\n",ffnormn(n_r*n_c,f_B_rc_ans__,f_A_rc__));
  printf(" %% f_B_cr_ans__ vs f_A_cr__: relative error %0.16f\n",ffnormn(n_r*n_c,f_B_cr_ans__,f_A_cr__));
  GLOBAL_toc(0,1," array_extract: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_normalize_row(n_r,n_c,d_A_rc__,d_A_cr__,"double",&d_A_rc__,&d_A_cr__);
  array_printf(d_A_rc__,"double",n_r,n_c," % d_A_rc__: ");
  array_printf(d_A_cr__,"double",n_c,n_r," % d_A_cr__: ");
  printf(" %% d_B_rc_ans__ vs d_A_rc__: relative error %0.16f\n",dfnormn(n_r*n_c,d_B_rc_ans__,d_A_rc__));
  printf(" %% d_B_cr_ans__ vs d_A_cr__: relative error %0.16f\n",dfnormn(n_r*n_c,d_B_cr_ans__,d_A_cr__));
  GLOBAL_toc(0,1," array_extract: ");
  /* %%%%%%%% */
  free1(&f_A_rc__);
  free1(&f_A_cr__);
  free1(&d_A_rc__);  
  free1(&d_A_cr__);  
}
