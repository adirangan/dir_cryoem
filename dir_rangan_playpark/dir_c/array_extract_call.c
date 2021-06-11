#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

void array_extract_i_from_i(int n_r,int n_c,int *A_rc__,int n_r_rtn,int *r_rtn_,int n_c_rtn,int *c_rtn_,int **B_rc_p_,int **B_cr_p_)
{
  int nr=0,nc=0;
  int nr_0in=0,nc_0in=0,nr_out=0,nc_out=0;
  int *A_r_=NULL;
  int *B_rc__=NULL,*B_cr__=NULL;
  int flag_r_extract=0;
  int flag_c_extract=0;
  int n_r_out=0,n_c_out=0;
  if ( (n_r_rtn>0) && (r_rtn_!=NULL) ){ flag_r_extract=1; n_r_out = n_r_rtn;} else{ flag_r_extract=0; n_r_out = n_r;}
  if ( (n_c_rtn>0) && (c_rtn_!=NULL) ){ flag_c_extract=1; n_c_out = n_c_rtn;} else{ flag_c_extract=0; n_c_out = n_c;}
  B_rc__=NULL;
  if (B_rc_p_!=NULL){
    if ( (*B_rc_p_)==NULL ){ (*B_rc_p_) = (int *) malloc1((unsigned long long int)n_r_out*(unsigned long long int)n_c_out*sizeof(int)); }
    B_rc__ = *B_rc_p_;
    /* if (B_rc_p_!=NULL){ } */}
  B_cr__=NULL;
  if (B_cr_p_!=NULL){
    if ( (*B_cr_p_)==NULL ){ (*B_cr_p_) = (int *) malloc1((unsigned long long int)n_c_out*(unsigned long long int)n_r_out*sizeof(int)); }
    B_cr__ = *B_cr_p_;
    /* if (B_cr_p_!=NULL){ } */}
  for (nc_out=0;nc_out<n_c_out;nc_out++){
    if (flag_c_extract){ nc_0in = c_rtn_[nc_out];} else{ nc_0in = nc_out;}
    A_r_ = A_rc__ + (unsigned long long int)nc_0in * (unsigned long long int)n_r;
    for (nr_out=0;nr_out<n_r_out;nr_out++){
      if (flag_r_extract){ nr_0in = r_rtn_[nr_out];} else{ nr_0in = nr_out;}
      if (B_rc__!=NULL){ B_rc__[nr_out + nc_out*n_r_out] = A_r_[nr_0in];}
      if (B_cr__!=NULL){ B_cr__[nc_out + nr_out*n_c_out] = A_r_[nr_0in];}
      /* for (nr_out=0;nr_out<n_r_out;nr_out++){ } */}
    /* for (nc_out=0;nc_out<n_c_out;nc_out++){ } */}
}

void array_extract_i_from_i_test()
{
  int n_r = 5;
  int n_c = 8;
  int *A_rc__=NULL;
  int *B_rc__=NULL;
  int *B_cr__=NULL;
  int n_r_rtn = 3;
  int n_c_rtn = 5;
  int r_rtn_[3] = { 0 , 2 , 3};
  int c_rtn_[5] = { 1 , 2 , 4 , 5 , 7};
  int nr=0,nc=0;
  A_rc__ = (int *) malloc1(n_r*n_c*sizeof(int));
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ A_rc__[nr+nc*n_r] = (nr+1)*100 + (nc+1);}};
  array_printf(A_rc__,"int",n_r,n_c," % A_rc__: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_i_from_i(n_r,n_c,A_rc__,0,NULL,0,NULL,NULL,NULL);
  GLOBAL_toc(0,1," array_extract_i_from_i: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_i_from_i(n_r,n_c,A_rc__,0,NULL,0,NULL,&B_rc__,NULL);
  array_printf(B_rc__,"int",n_r,n_c," % B_rc__: ");
  free1(&B_rc__);
  GLOBAL_toc(0,1," array_extract_i_from_i: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_i_from_i(n_r,n_c,A_rc__,0,NULL,0,NULL,NULL,&B_cr__);
  array_printf(B_cr__,"int",n_c,n_r," % B_cr__: ");
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_i_from_i: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_i_from_i(n_r,n_c,A_rc__,0,NULL,0,NULL,&B_rc__,&B_cr__);
  array_printf(B_rc__,"int",n_r,n_c," % B_rc__: ");
  array_printf(B_cr__,"int",n_c,n_r," % B_cr__: ");
  free1(&B_rc__);
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_i_from_i: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_i_from_i(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,0,NULL,&B_rc__,NULL);
  array_printf(B_rc__,"int",n_r_rtn,n_c," % B_rc__: ");
  free1(&B_rc__);
  GLOBAL_toc(0,1," array_extract_i_from_i: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_i_from_i(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,0,NULL,&B_rc__,&B_cr__);
  array_printf(B_rc__,"int",n_r_rtn,n_c," % B_rc__: ");
  array_printf(B_cr__,"int",n_c,n_r_rtn," % B_cr__: ");
  free1(&B_rc__);
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_i_from_i: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_i_from_i(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,&B_rc__,NULL);
  array_printf(B_rc__,"int",n_r_rtn,n_c_rtn," % B_rc__: ");
  free1(&B_rc__);
  GLOBAL_toc(0,1," array_extract_i_from_i: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_i_from_i(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,&B_rc__,&B_cr__);
  array_printf(B_rc__,"int",n_r_rtn,n_c_rtn," % B_rc__: ");
  array_printf(B_cr__,"int",n_c_rtn,n_r_rtn," % B_cr__: ");
  free1(&B_rc__);
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_i_from_i: ");
  /* %%%%%%%% */
  free1(&A_rc__);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void array_extract_f_from_d(int n_r,int n_c,double *A_rc__,int n_r_rtn,int *r_rtn_,int n_c_rtn,int *c_rtn_,float **B_rc_p_,float **B_cr_p_)
{
  int nr=0,nc=0;
  int nr_0in=0,nc_0in=0,nr_out=0,nc_out=0;
  double *A_r_=NULL;
  float *B_rc__=NULL,*B_cr__=NULL;
  int flag_r_extract=0;
  int flag_c_extract=0;
  int n_r_out=0,n_c_out=0;
  if ( (n_r_rtn>0) && (r_rtn_!=NULL) ){ flag_r_extract=1; n_r_out = n_r_rtn;} else{ flag_r_extract=0; n_r_out = n_r;}
  if ( (n_c_rtn>0) && (c_rtn_!=NULL) ){ flag_c_extract=1; n_c_out = n_c_rtn;} else{ flag_c_extract=0; n_c_out = n_c;}
  B_rc__=NULL;
  if (B_rc_p_!=NULL){
    if ( (*B_rc_p_)==NULL ){ (*B_rc_p_) = (float *) malloc1((unsigned long long int)n_r_out*(unsigned long long int)n_c_out*sizeof(float)); }
    B_rc__ = *B_rc_p_;
    /* if (B_rc_p_!=NULL){ } */}
  B_cr__=NULL;
  if (B_cr_p_!=NULL){
    if ( (*B_cr_p_)==NULL ){ (*B_cr_p_) = (float *) malloc1((unsigned long long int)n_c_out*(unsigned long long int)n_r_out*sizeof(float)); }
    B_cr__ = *B_cr_p_;
    /* if (B_cr_p_!=NULL){ } */}
  for (nc_out=0;nc_out<n_c_out;nc_out++){
    if (flag_c_extract){ nc_0in = c_rtn_[nc_out];} else{ nc_0in = nc_out;}
    A_r_ = A_rc__ + (unsigned long long int)nc_0in * (unsigned long long int)n_r;
    for (nr_out=0;nr_out<n_r_out;nr_out++){
      if (flag_r_extract){ nr_0in = r_rtn_[nr_out];} else{ nr_0in = nr_out;}
      if (B_rc__!=NULL){ B_rc__[nr_out + nc_out*n_r_out] = A_r_[nr_0in];}
      if (B_cr__!=NULL){ B_cr__[nc_out + nr_out*n_c_out] = A_r_[nr_0in];}
      /* for (nr_out=0;nr_out<n_r_out;nr_out++){ } */}
    /* for (nc_out=0;nc_out<n_c_out;nc_out++){ } */}
}

void array_extract_f_from_d_test()
{
  int n_r = 5;
  int n_c = 8;
  double *A_rc__=NULL;
  float *B_rc__=NULL;
  float *B_cr__=NULL;
  int n_r_rtn = 3;
  int n_c_rtn = 5;
  int r_rtn_[3] = { 0 , 2 , 3};
  int c_rtn_[5] = { 1 , 2 , 4 , 5 , 7};
  int nr=0,nc=0;
  A_rc__ = (double *) malloc1(n_r*n_c*sizeof(double));
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ A_rc__[nr+nc*n_r] = (nr+1) + 0.1*(nc+1);}};
  array_printf(A_rc__,"double",n_r,n_c," % A_rc__: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_d(n_r,n_c,A_rc__,0,NULL,0,NULL,NULL,NULL);
  GLOBAL_toc(0,1," array_extract_f_from_d: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_d(n_r,n_c,A_rc__,0,NULL,0,NULL,&B_rc__,NULL);
  array_printf(B_rc__,"float",n_r,n_c," % B_rc__: ");
  free1(&B_rc__);
  GLOBAL_toc(0,1," array_extract_f_from_d: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_d(n_r,n_c,A_rc__,0,NULL,0,NULL,NULL,&B_cr__);
  array_printf(B_cr__,"float",n_c,n_r," % B_cr__: ");
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_f_from_d: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_d(n_r,n_c,A_rc__,0,NULL,0,NULL,&B_rc__,&B_cr__);
  array_printf(B_rc__,"float",n_r,n_c," % B_rc__: ");
  array_printf(B_cr__,"float",n_c,n_r," % B_cr__: ");
  free1(&B_rc__);
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_f_from_d: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_d(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,0,NULL,&B_rc__,NULL);
  array_printf(B_rc__,"float",n_r_rtn,n_c," % B_rc__: ");
  free1(&B_rc__);
  GLOBAL_toc(0,1," array_extract_f_from_d: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_d(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,0,NULL,&B_rc__,&B_cr__);
  array_printf(B_rc__,"float",n_r_rtn,n_c," % B_rc__: ");
  array_printf(B_cr__,"float",n_c,n_r_rtn," % B_cr__: ");
  free1(&B_rc__);
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_f_from_d: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_d(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,&B_rc__,NULL);
  array_printf(B_rc__,"float",n_r_rtn,n_c_rtn," % B_rc__: ");
  free1(&B_rc__);
  GLOBAL_toc(0,1," array_extract_f_from_d: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_d(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,&B_rc__,&B_cr__);
  array_printf(B_rc__,"float",n_r_rtn,n_c_rtn," % B_rc__: ");
  array_printf(B_cr__,"float",n_c_rtn,n_r_rtn," % B_cr__: ");
  free1(&B_rc__);
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_f_from_d: ");
  /* %%%%%%%% */
  free1(&A_rc__);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void array_extract_d_from_d(int n_r,int n_c,double *A_rc__,int n_r_rtn,int *r_rtn_,int n_c_rtn,int *c_rtn_,double **B_rc_p_,double **B_cr_p_)
{
  int nr=0,nc=0;
  int nr_0in=0,nc_0in=0,nr_out=0,nc_out=0;
  double *A_r_=NULL;
  double *B_rc__=NULL,*B_cr__=NULL;
  int flag_r_extract=0;
  int flag_c_extract=0;
  int n_r_out=0,n_c_out=0;
  if ( (n_r_rtn>0) && (r_rtn_!=NULL) ){ flag_r_extract=1; n_r_out = n_r_rtn;} else{ flag_r_extract=0; n_r_out = n_r;}
  if ( (n_c_rtn>0) && (c_rtn_!=NULL) ){ flag_c_extract=1; n_c_out = n_c_rtn;} else{ flag_c_extract=0; n_c_out = n_c;}
  B_rc__=NULL;
  if (B_rc_p_!=NULL){
    if ( (*B_rc_p_)==NULL ){ (*B_rc_p_) = (double *) malloc1((unsigned long long int)n_r_out*(unsigned long long int)n_c_out*sizeof(double)); }
    B_rc__ = *B_rc_p_;
    /* if (B_rc_p_!=NULL){ } */}
  B_cr__=NULL;
  if (B_cr_p_!=NULL){
    if ( (*B_cr_p_)==NULL ){ (*B_cr_p_) = (double *) malloc1((unsigned long long int)n_c_out*(unsigned long long int)n_r_out*sizeof(double)); }
    B_cr__ = *B_cr_p_;
    /* if (B_cr_p_!=NULL){ } */}
  for (nc_out=0;nc_out<n_c_out;nc_out++){
    if (flag_c_extract){ nc_0in = c_rtn_[nc_out];} else{ nc_0in = nc_out;}
    A_r_ = A_rc__ + (unsigned long long int)nc_0in * (unsigned long long int)n_r;
    for (nr_out=0;nr_out<n_r_out;nr_out++){
      if (flag_r_extract){ nr_0in = r_rtn_[nr_out];} else{ nr_0in = nr_out;}
      if (B_rc__!=NULL){ B_rc__[nr_out + nc_out*n_r_out] = A_r_[nr_0in];}
      if (B_cr__!=NULL){ B_cr__[nc_out + nr_out*n_c_out] = A_r_[nr_0in];}
      /* for (nr_out=0;nr_out<n_r_out;nr_out++){ } */}
    /* for (nc_out=0;nc_out<n_c_out;nc_out++){ } */}
}

void array_extract_d_from_d_test()
{
  int n_r = 5;
  int n_c = 8;
  double *A_rc__=NULL;
  double *B_rc__=NULL;
  double *B_cr__=NULL;
  int n_r_rtn = 3;
  int n_c_rtn = 5;
  int r_rtn_[3] = { 0 , 2 , 3};
  int c_rtn_[5] = { 1 , 2 , 4 , 5 , 7};
  int nr=0,nc=0;
  A_rc__ = (double *) malloc1(n_r*n_c*sizeof(double));
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ A_rc__[nr+nc*n_r] = (nr+1) + 0.1*(nc+1);}};
  array_printf(A_rc__,"double",n_r,n_c," % A_rc__: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_d(n_r,n_c,A_rc__,0,NULL,0,NULL,NULL,NULL);
  GLOBAL_toc(0,1," array_extract_d_from_d: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_d(n_r,n_c,A_rc__,0,NULL,0,NULL,&B_rc__,NULL);
  array_printf(B_rc__,"double",n_r,n_c," % B_rc__: ");
  free1(&B_rc__);
  GLOBAL_toc(0,1," array_extract_d_from_d: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_d(n_r,n_c,A_rc__,0,NULL,0,NULL,NULL,&B_cr__);
  array_printf(B_cr__,"double",n_c,n_r," % B_cr__: ");
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_d_from_d: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_d(n_r,n_c,A_rc__,0,NULL,0,NULL,&B_rc__,&B_cr__);
  array_printf(B_rc__,"double",n_r,n_c," % B_rc__: ");
  array_printf(B_cr__,"double",n_c,n_r," % B_cr__: ");
  free1(&B_rc__);
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_d_from_d: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_d(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,0,NULL,&B_rc__,NULL);
  array_printf(B_rc__,"double",n_r_rtn,n_c," % B_rc__: ");
  free1(&B_rc__);
  GLOBAL_toc(0,1," array_extract_d_from_d: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_d(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,0,NULL,&B_rc__,&B_cr__);
  array_printf(B_rc__,"double",n_r_rtn,n_c," % B_rc__: ");
  array_printf(B_cr__,"double",n_c,n_r_rtn," % B_cr__: ");
  free1(&B_rc__);
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_d_from_d: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_d(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,&B_rc__,NULL);
  array_printf(B_rc__,"double",n_r_rtn,n_c_rtn," % B_rc__: ");
  free1(&B_rc__);
  GLOBAL_toc(0,1," array_extract_d_from_d: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_d(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,&B_rc__,&B_cr__);
  array_printf(B_rc__,"double",n_r_rtn,n_c_rtn," % B_rc__: ");
  array_printf(B_cr__,"double",n_c_rtn,n_r_rtn," % B_cr__: ");
  free1(&B_rc__);
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_d_from_d: ");
  /* %%%%%%%% */
  free1(&A_rc__);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void array_extract_f_from_f(int n_r,int n_c,float *A_rc__,int n_r_rtn,int *r_rtn_,int n_c_rtn,int *c_rtn_,float **B_rc_p_,float **B_cr_p_)
{
  int nr=0,nc=0;
  int nr_0in=0,nc_0in=0,nr_out=0,nc_out=0;
  float *A_r_=NULL;
  float *B_rc__=NULL,*B_cr__=NULL;
  int flag_r_extract=0;
  int flag_c_extract=0;
  int n_r_out=0,n_c_out=0;
  if ( (n_r_rtn>0) && (r_rtn_!=NULL) ){ flag_r_extract=1; n_r_out = n_r_rtn;} else{ flag_r_extract=0; n_r_out = n_r;}
  if ( (n_c_rtn>0) && (c_rtn_!=NULL) ){ flag_c_extract=1; n_c_out = n_c_rtn;} else{ flag_c_extract=0; n_c_out = n_c;}
  B_rc__=NULL;
  if (B_rc_p_!=NULL){
    if ( (*B_rc_p_)==NULL ){ (*B_rc_p_) = (float *) malloc1((unsigned long long int)n_r_out*(unsigned long long int)n_c_out*sizeof(float)); }
    B_rc__ = *B_rc_p_;
    /* if (B_rc_p_!=NULL){ } */}
  B_cr__=NULL;
  if (B_cr_p_!=NULL){
    if ( (*B_cr_p_)==NULL ){ (*B_cr_p_) = (float *) malloc1((unsigned long long int)n_c_out*(unsigned long long int)n_r_out*sizeof(float)); }
    B_cr__ = *B_cr_p_;
    /* if (B_cr_p_!=NULL){ } */}
  for (nc_out=0;nc_out<n_c_out;nc_out++){
    if (flag_c_extract){ nc_0in = c_rtn_[nc_out];} else{ nc_0in = nc_out;}
    A_r_ = A_rc__ + (unsigned long long int)nc_0in * (unsigned long long int)n_r;
    for (nr_out=0;nr_out<n_r_out;nr_out++){
      if (flag_r_extract){ nr_0in = r_rtn_[nr_out];} else{ nr_0in = nr_out;}
      if (B_rc__!=NULL){ B_rc__[nr_out + nc_out*n_r_out] = A_r_[nr_0in];}
      if (B_cr__!=NULL){ B_cr__[nc_out + nr_out*n_c_out] = A_r_[nr_0in];}
      /* for (nr_out=0;nr_out<n_r_out;nr_out++){ } */}
    /* for (nc_out=0;nc_out<n_c_out;nc_out++){ } */}
}

void array_extract_f_from_f_test()
{
  int n_r = 5;
  int n_c = 8;
  float *A_rc__=NULL;
  float *B_rc__=NULL;
  float *B_cr__=NULL;
  int n_r_rtn = 3;
  int n_c_rtn = 5;
  int r_rtn_[3] = { 0 , 2 , 3};
  int c_rtn_[5] = { 1 , 2 , 4 , 5 , 7};
  int nr=0,nc=0;
  A_rc__ = (float *) malloc1(n_r*n_c*sizeof(float));
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ A_rc__[nr+nc*n_r] = (nr+1) + 0.1*(nc+1);}};
  array_printf(A_rc__,"float",n_r,n_c," % A_rc__: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_f(n_r,n_c,A_rc__,0,NULL,0,NULL,NULL,NULL);
  GLOBAL_toc(0,1," array_extract_f_from_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_f(n_r,n_c,A_rc__,0,NULL,0,NULL,&B_rc__,NULL);
  array_printf(B_rc__,"float",n_r,n_c," % B_rc__: ");
  free1(&B_rc__);
  GLOBAL_toc(0,1," array_extract_f_from_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_f(n_r,n_c,A_rc__,0,NULL,0,NULL,NULL,&B_cr__);
  array_printf(B_cr__,"float",n_c,n_r," % B_cr__: ");
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_f_from_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_f(n_r,n_c,A_rc__,0,NULL,0,NULL,&B_rc__,&B_cr__);
  array_printf(B_rc__,"float",n_r,n_c," % B_rc__: ");
  array_printf(B_cr__,"float",n_c,n_r," % B_cr__: ");
  free1(&B_rc__);
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_f_from_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_f(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,0,NULL,&B_rc__,NULL);
  array_printf(B_rc__,"float",n_r_rtn,n_c," % B_rc__: ");
  free1(&B_rc__);
  GLOBAL_toc(0,1," array_extract_f_from_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_f(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,0,NULL,&B_rc__,&B_cr__);
  array_printf(B_rc__,"float",n_r_rtn,n_c," % B_rc__: ");
  array_printf(B_cr__,"float",n_c,n_r_rtn," % B_cr__: ");
  free1(&B_rc__);
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_f_from_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_f(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,&B_rc__,NULL);
  array_printf(B_rc__,"float",n_r_rtn,n_c_rtn," % B_rc__: ");
  free1(&B_rc__);
  GLOBAL_toc(0,1," array_extract_f_from_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_f_from_f(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,&B_rc__,&B_cr__);
  array_printf(B_rc__,"float",n_r_rtn,n_c_rtn," % B_rc__: ");
  array_printf(B_cr__,"float",n_c_rtn,n_r_rtn," % B_cr__: ");
  free1(&B_rc__);
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_f_from_f: ");
  /* %%%%%%%% */
  free1(&A_rc__);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void array_extract_d_from_f(int n_r,int n_c,float *A_rc__,int n_r_rtn,int *r_rtn_,int n_c_rtn,int *c_rtn_,double **B_rc_p_,double **B_cr_p_)
{
  int nr=0,nc=0;
  int nr_0in=0,nc_0in=0,nr_out=0,nc_out=0;
  float *A_r_=NULL;
  double *B_rc__=NULL,*B_cr__=NULL;
  int flag_r_extract=0;
  int flag_c_extract=0;
  int n_r_out=0,n_c_out=0;
  if ( (n_r_rtn>0) && (r_rtn_!=NULL) ){ flag_r_extract=1; n_r_out = n_r_rtn;} else{ flag_r_extract=0; n_r_out = n_r;}
  if ( (n_c_rtn>0) && (c_rtn_!=NULL) ){ flag_c_extract=1; n_c_out = n_c_rtn;} else{ flag_c_extract=0; n_c_out = n_c;}
  B_rc__=NULL;
  if (B_rc_p_!=NULL){
    if ( (*B_rc_p_)==NULL ){ (*B_rc_p_) = (double *) malloc1((unsigned long long int)n_r_out*(unsigned long long int)n_c_out*sizeof(double)); }
    B_rc__ = *B_rc_p_;
    /* if (B_rc_p_!=NULL){ } */}
  B_cr__=NULL;
  if (B_cr_p_!=NULL){
    if ( (*B_cr_p_)==NULL ){ (*B_cr_p_) = (double *) malloc1((unsigned long long int)n_c_out*(unsigned long long int)n_r_out*sizeof(double)); }
    B_cr__ = *B_cr_p_;
    /* if (B_cr_p_!=NULL){ } */}
  for (nc_out=0;nc_out<n_c_out;nc_out++){
    if (flag_c_extract){ nc_0in = c_rtn_[nc_out];} else{ nc_0in = nc_out;}
    A_r_ = A_rc__ + (unsigned long long int)nc_0in * (unsigned long long int)n_r;
    for (nr_out=0;nr_out<n_r_out;nr_out++){
      if (flag_r_extract){ nr_0in = r_rtn_[nr_out];} else{ nr_0in = nr_out;}
      if (B_rc__!=NULL){ B_rc__[nr_out + nc_out*n_r_out] = A_r_[nr_0in];}
      if (B_cr__!=NULL){ B_cr__[nc_out + nr_out*n_c_out] = A_r_[nr_0in];}
      /* for (nr_out=0;nr_out<n_r_out;nr_out++){ } */}
    /* for (nc_out=0;nc_out<n_c_out;nc_out++){ } */}
}

void array_extract_d_from_f_test()
{
  int n_r = 5;
  int n_c = 8;
  float *A_rc__=NULL;
  double *B_rc__=NULL;
  double *B_cr__=NULL;
  int n_r_rtn = 3;
  int n_c_rtn = 5;
  int r_rtn_[3] = { 0 , 2 , 3};
  int c_rtn_[5] = { 1 , 2 , 4 , 5 , 7};
  int nr=0,nc=0;
  A_rc__ = (float *) malloc1(n_r*n_c*sizeof(float));
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ A_rc__[nr+nc*n_r] = (nr+1) + 0.1*(nc+1);}};
  array_printf(A_rc__,"float",n_r,n_c," % A_rc__: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_f(n_r,n_c,A_rc__,0,NULL,0,NULL,NULL,NULL);
  GLOBAL_toc(0,1," array_extract_d_from_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_f(n_r,n_c,A_rc__,0,NULL,0,NULL,&B_rc__,NULL);
  array_printf(B_rc__,"double",n_r,n_c," % B_rc__: ");
  free1(&B_rc__);
  GLOBAL_toc(0,1," array_extract_d_from_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_f(n_r,n_c,A_rc__,0,NULL,0,NULL,NULL,&B_cr__);
  array_printf(B_cr__,"double",n_c,n_r," % B_cr__: ");
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_d_from_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_f(n_r,n_c,A_rc__,0,NULL,0,NULL,&B_rc__,&B_cr__);
  array_printf(B_rc__,"double",n_r,n_c," % B_rc__: ");
  array_printf(B_cr__,"double",n_c,n_r," % B_cr__: ");
  free1(&B_rc__);
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_d_from_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_f(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,0,NULL,&B_rc__,NULL);
  array_printf(B_rc__,"double",n_r_rtn,n_c," % B_rc__: ");
  free1(&B_rc__);
  GLOBAL_toc(0,1," array_extract_d_from_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_f(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,0,NULL,&B_rc__,&B_cr__);
  array_printf(B_rc__,"double",n_r_rtn,n_c," % B_rc__: ");
  array_printf(B_cr__,"double",n_c,n_r_rtn," % B_cr__: ");
  free1(&B_rc__);
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_d_from_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_f(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,&B_rc__,NULL);
  array_printf(B_rc__,"double",n_r_rtn,n_c_rtn," % B_rc__: ");
  free1(&B_rc__);
  GLOBAL_toc(0,1," array_extract_d_from_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract_d_from_f(n_r,n_c,A_rc__,n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,&B_rc__,&B_cr__);
  array_printf(B_rc__,"double",n_r_rtn,n_c_rtn," % B_rc__: ");
  array_printf(B_cr__,"double",n_c_rtn,n_r_rtn," % B_cr__: ");
  free1(&B_rc__);
  free1(&B_cr__);
  GLOBAL_toc(0,1," array_extract_d_from_f: ");
  /* %%%%%%%% */
  free1(&A_rc__);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void array_extract(int n_r,int n_c,void *A_rc__,char *type_A,int n_r_rtn,int *r_rtn_,int n_c_rtn,int *c_rtn_,void *B_rc_p_,void *B_cr_p_,char *type_B)
{
  if ( (strcmp(type_A,"float")==0) && (strcmp(type_B,"double")==0) ){
    array_extract_d_from_f(n_r,n_c,(float *)A_rc__,n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,(double **)B_rc_p_,(double **)B_cr_p_);
    /* float double */}
  if ( (strcmp(type_A,"double")==0) && (strcmp(type_B,"double")==0) ){
    array_extract_d_from_d(n_r,n_c,(double *)A_rc__,n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,(double **)B_rc_p_,(double **)B_cr_p_);
    /* double double */}
  if ( (strcmp(type_A,"float")==0) && (strcmp(type_B,"float")==0) ){
    array_extract_f_from_f(n_r,n_c,(float *)A_rc__,n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,(float **)B_rc_p_,(float **)B_cr_p_);
    /* float float */}
  if ( (strcmp(type_A,"double")==0) && (strcmp(type_B,"float")==0) ){
    array_extract_f_from_d(n_r,n_c,(double *)A_rc__,n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,(float **)B_rc_p_,(float **)B_cr_p_);
    /* double float */}
}  

void array_extract_test()
{
  int n_r = 5;
  int n_c = 8;
  float *f_A_rc__=NULL;
  double *d_A_rc__=NULL;
  float *f_B_rc__=NULL;
  float *f_B_cr__=NULL;
  double *d_B_rc__=NULL;
  double *d_B_cr__=NULL;
  int n_r_rtn = 3;
  int n_c_rtn = 5;
  int r_rtn_[3] = { 0 , 2 , 3};
  int c_rtn_[5] = { 1 , 2 , 4 , 5 , 7};
  float f_B_rc_ans__[15] = { 1.2 , 3.2 , 4.2 , 1.3 , 3.3 , 4.3 , 1.5 , 3.5 , 4.5 , 1.6 , 3.6 , 4.6 , 1.8 , 3.8 , 4.8 };
  float f_B_cr_ans__[15] = { 1.2 , 1.3 , 1.5 , 1.6 , 1.8 , 3.2 , 3.3 , 3.5 , 3.6 , 3.8 , 4.2 , 4.3 , 4.5 , 4.6 , 4.8 };
  double d_B_rc_ans__[15] = { 1.2 , 3.2 , 4.2 , 1.3 , 3.3 , 4.3 , 1.5 , 3.5 , 4.5 , 1.6 , 3.6 , 4.6 , 1.8 , 3.8 , 4.8 };
  double d_B_cr_ans__[15] = { 1.2 , 1.3 , 1.5 , 1.6 , 1.8 , 3.2 , 3.3 , 3.5 , 3.6 , 3.8 , 4.2 , 4.3 , 4.5 , 4.6 , 4.8 };
  int nr=0,nc=0;
  f_A_rc__ = (float *) malloc1(n_r*n_c*sizeof(float));
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ f_A_rc__[nr+nc*n_r] = (nr+1) + 0.1*(nc+1);}};
  array_printf(f_A_rc__,"float",n_r,n_c," % f_A_rc__: ");
  d_A_rc__ = (double *) malloc1(n_r*n_c*sizeof(double));
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ d_A_rc__[nr+nc*n_r] = (nr+1) + 0.1*(nc+1);}};
  array_printf(d_A_rc__,"double",n_r,n_c," % d_A_rc__: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract(n_r,n_c,f_A_rc__,"float",n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,&f_B_rc__,&f_B_cr__,"float");
  array_printf(f_B_rc__,"float",n_r_rtn,n_c_rtn," % f_B_rc__: ");
  array_printf(f_B_cr__,"float",n_c_rtn,n_r_rtn," % f_B_cr__: ");
  printf(" %% f_B_rc_ans__ vs f_B_rc__: relative error %0.16f\n",ffnormn(n_r_rtn*n_c_rtn,f_B_rc_ans__,f_B_rc__));
  printf(" %% f_B_cr_ans__ vs f_B_cr__: relative error %0.16f\n",ffnormn(n_r_rtn*n_c_rtn,f_B_cr_ans__,f_B_cr__));
  free1(&f_B_rc__);
  free1(&f_B_cr__);
  GLOBAL_toc(0,1," array_extract: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract(n_r,n_c,f_A_rc__,"float",n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,&d_B_rc__,&d_B_cr__,"double");
  array_printf(d_B_rc__,"double",n_r_rtn,n_c_rtn," % d_B_rc__: ");
  array_printf(d_B_cr__,"double",n_c_rtn,n_r_rtn," % d_B_cr__: ");
  printf(" %% d_B_rc_ans__ vs d_B_rc__: relative error %0.16f\n",dfnormn(n_r_rtn*n_c_rtn,d_B_rc_ans__,d_B_rc__));
  printf(" %% d_B_cr_ans__ vs d_B_cr__: relative error %0.16f\n",dfnormn(n_r_rtn*n_c_rtn,d_B_cr_ans__,d_B_cr__));
  free1(&d_B_rc__);
  free1(&d_B_cr__);
  GLOBAL_toc(0,1," array_extract: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract(n_r,n_c,d_A_rc__,"double",n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,&f_B_rc__,&f_B_cr__,"float");
  array_printf(f_B_rc__,"float",n_r_rtn,n_c_rtn," % f_B_rc__: ");
  array_printf(f_B_cr__,"float",n_c_rtn,n_r_rtn," % f_B_cr__: ");
  printf(" %% f_B_rc_ans__ vs f_B_rc__: relative error %0.16f\n",ffnormn(n_r_rtn*n_c_rtn,f_B_rc_ans__,f_B_rc__));
  printf(" %% f_B_cr_ans__ vs f_B_cr__: relative error %0.16f\n",ffnormn(n_r_rtn*n_c_rtn,f_B_cr_ans__,f_B_cr__));
  free1(&f_B_rc__);
  free1(&f_B_cr__);
  GLOBAL_toc(0,1," array_extract: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_extract(n_r,n_c,d_A_rc__,"double",n_r_rtn,r_rtn_,n_c_rtn,c_rtn_,&d_B_rc__,&d_B_cr__,"double");
  array_printf(d_B_rc__,"double",n_r_rtn,n_c_rtn," % d_B_rc__: ");
  array_printf(d_B_cr__,"double",n_c_rtn,n_r_rtn," % d_B_cr__: ");
  printf(" %% d_B_rc_ans__ vs d_B_rc__: relative error %0.16f\n",dfnormn(n_r_rtn*n_c_rtn,d_B_rc_ans__,d_B_rc__));
  printf(" %% d_B_cr_ans__ vs d_B_cr__: relative error %0.16f\n",dfnormn(n_r_rtn*n_c_rtn,d_B_cr_ans__,d_B_cr__));
  free1(&d_B_rc__);
  free1(&d_B_cr__);
  GLOBAL_toc(0,1," array_extract: ");
  /* %%%%%%%% */
  free1(&f_A_rc__);
  free1(&d_A_rc__);
}
