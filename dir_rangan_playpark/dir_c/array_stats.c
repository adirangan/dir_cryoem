#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

void icumsum(unsigned long long int ulli_length,int *i_0in_,int **i_out_p_)
{ 
  /* cumulative sum */
  unsigned long long int ulli=0;
  int sum=0;
  if (*i_out_p_==NULL){ *i_out_p_ = (int *) wkspace_all0c(ulli_length*sizeof(int));}
  for (ulli=0;ulli<ulli_length;ulli++){ sum += i_0in_[ulli]; (*i_out_p_)[ulli]=sum;}
}

double ifnormn(unsigned long long int ulli_length,int *i_0_,int *i_1_)
{
  unsigned long long int ulli=0;
  int i_01 = 0;
  int i_0 = 0;
  double d_f01 = 0;
  double d_f0 = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ i_01 = i_0_[ulli]-i_1_[ulli]; d_f01 += i_01*i_01; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  d_f01 = sqrt(d_f01);
  for (ulli=0;ulli<ulli_length;ulli++){ i_0 = i_0_[ulli]; d_f0 += i_0*i_0; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  d_f0 = sqrt(d_f0);
  return d_f01/maximum(1e-12,d_f0);
}

float ffnorm(unsigned long long int ulli_length,float *f_0_,float *f_1_)
{
  float output = 0;
  unsigned long long int ulli=0;
  float f_diff = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ f_diff = f_0_[ulli]-f_1_[ulli]; output += f_diff*f_diff; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  return sqrtf(output);
}

float ffnormn(unsigned long long int ulli_length,float *f_0_,float *f_1_)
{
  unsigned long long int ulli=0;
  float f_01 = 0;
  float f_0 = 0;
  float f_f01 = 0;
  float f_f0 = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ f_01 = f_0_[ulli]-f_1_[ulli]; f_f01 += f_01*f_01; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  f_f01 = sqrt(f_f01);
  for (ulli=0;ulli<ulli_length;ulli++){ f_0 = f_0_[ulli]; f_f0 += f_0*f_0; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  f_f0 = sqrtf(f_f0);
  return f_f01/maximum(1e-12,f_f0);
}

double dfnormn(unsigned long long int ulli_length,double *d_0_,double *d_1_)
{
  unsigned long long int ulli=0;
  double d_01 = 0;
  double d_0 = 0;
  double d_f01 = 0;
  double d_f0 = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ d_01 = d_0_[ulli]-d_1_[ulli]; d_f01 += d_01*d_01; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  d_f01 = sqrt(d_f01);
  for (ulli=0;ulli<ulli_length;ulli++){ d_0 = d_0_[ulli]; d_f0 += d_0*d_0; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  d_f0 = sqrt(d_f0);
  return d_f01/maximum(1e-12,d_f0);
}

double dfnorm(unsigned long long int ulli_length,double *d_0_,double *d_1_)
{
  double output = 0;
  unsigned long long int ulli=0;
  double d_diff = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ d_diff = d_0_[ulli]-d_1_[ulli]; output += d_diff*d_diff; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  return sqrt(output);
}

double cfnorm(unsigned long long int ulli_length,float complex *c_0_,float complex *c_1_)
{
  double output = 0;
  unsigned long long int ulli=0;
  float complex c_diff = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ c_diff = c_0_[ulli]-c_1_[ulli]; output += crealf(c_diff*conjf(c_diff)); /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  return sqrt(output);
}

double zfnorm(unsigned long long int ulli_length,double complex *z_0_,double complex *z_1_)
{
  double output = 0;
  unsigned long long int ulli=0;
  double complex z_diff = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ z_diff = z_0_[ulli]-z_1_[ulli]; output += creal(z_diff*conj(z_diff)); /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  return sqrt(output);
}

void array_stats(void *v_,char *type,unsigned long long int ulli_length,void *max_p,void *min_p,double *mean_p,double *stdev_p)
{
  /* finds the stats of an array */
  unsigned long long int ulli=0;
  double mean=0,stdev=0;
  float f_max=0,f_min=0;
  double d_max=0,d_min=0;
  long long int lli_max=0,lli_min=0;
  int i_max=0,i_min=0;
  if (v_!=NULL){
    if (0){ /* do nothing */}
    else if (strcmp(type,"float")==0){
      for (ulli=0;ulli<ulli_length;ulli++){ mean += ((float *)v_)[ulli];}
      mean /= (double)ulli_length;
      f_max=((float *)v_)[0]; f_min=((float *)v_)[0];
      for (ulli=0;ulli<ulli_length;ulli++){
	if (((float *)v_)[ulli]>f_max){ f_max=((float *)v_)[ulli];}
	if (((float *)v_)[ulli]<f_min){ f_min=((float *)v_)[ulli];}
	stdev += (((float *)v_)[ulli]-mean)*(((float *)v_)[ulli]-mean);}
      if (max_p!=NULL){*(float *)max_p=f_max;}
      if (min_p!=NULL){*(float *)min_p=f_min;}
      stdev /= (double)ulli_length;
      stdev = sqrt(stdev);}
    else if (strcmp(type,"double")==0){
      for (ulli=0;ulli<ulli_length;ulli++){ mean += ((double *)v_)[ulli];}
      mean /= (double)ulli_length;
      d_max=((double *)v_)[0]; d_min=((double *)v_)[0];
      for (ulli=0;ulli<ulli_length;ulli++){
	if (((double *)v_)[ulli]>d_max){ d_max=((double *)v_)[ulli];}
	if (((double *)v_)[ulli]<d_min){ d_min=((double *)v_)[ulli];}
	stdev += (((double *)v_)[ulli]-mean)*(((double *)v_)[ulli]-mean);}
      if (max_p!=NULL){*(double *)max_p=d_max;}
      if (min_p!=NULL){*(double *)min_p=d_min;}
      stdev /= (double)ulli_length;
      stdev = sqrt(stdev);}
    else if (strcmp(type,"long long int")==0){
      for (ulli=0;ulli<ulli_length;ulli++){ mean += ((long long int *)v_)[ulli];}
      mean /= (double)ulli_length;
      lli_max=((long long int *)v_)[0]; lli_min=((long long int *)v_)[0];
      for (ulli=0;ulli<ulli_length;ulli++){
	if (((long long int *)v_)[ulli]>lli_max){ lli_max=((long long int *)v_)[ulli];}
	if (((long long int *)v_)[ulli]<lli_min){ lli_min=((long long int *)v_)[ulli];}
	stdev += (((long long int *)v_)[ulli]-mean)*(((long long int *)v_)[ulli]-mean);}
      if (max_p!=NULL){*(long long int *)max_p=lli_max;}
      if (min_p!=NULL){*(long long int *)min_p=lli_min;}
      stdev /= (double)ulli_length;
      stdev = sqrt(stdev);}
    else if (strcmp(type,"int")==0){
      for (ulli=0;ulli<ulli_length;ulli++){ mean += ((int *)v_)[ulli];}
      mean /= (double)ulli_length;
      i_max=((int *)v_)[0]; i_min=((int *)v_)[0];
      for (ulli=0;ulli<ulli_length;ulli++){
	if (((int *)v_)[ulli]>i_max){ i_max=((int *)v_)[ulli];}
	if (((int *)v_)[ulli]<i_min){ i_min=((int *)v_)[ulli];}
	stdev += (((int *)v_)[ulli]-mean)*(((int *)v_)[ulli]-mean);}
      if (max_p!=NULL){*(int *)max_p=i_max;}
      if (min_p!=NULL){*(int *)min_p=i_min;}
      stdev /= (double)ulli_length;
      stdev = sqrt(stdev);}
    if (mean_p!=NULL){*mean_p=mean;}
    if (stdev_p!=NULL){*stdev_p=stdev;}}
}

void array_mean_center_row_f(int n_r,int n_c,float *f_rc_0in__,float **f_rc_out_p_)
{
  int nr=0,nc=0;
  float *f_rc_out__=NULL;
  float *f_r_=NULL;
  double mean=0,stdev=0;
  f_rc_out__=NULL;
  if (f_rc_out_p_!=NULL){
    if ( (*f_rc_out_p_)==NULL ){ (*f_rc_out_p_) = (float *) malloc1((unsigned long long int)n_r*(unsigned long long int)n_c*sizeof(float));}
    f_rc_out__ = *f_rc_out_p_;
    /* if (f_rc_out_p_!=NULL){ } */}
  if (f_rc_out__!=NULL){
    for (nc=0;nc<n_c;nc++){
      f_r_ = f_rc_0in__ + (unsigned long long int)nc*(unsigned long long int)n_r;
      array_stats(f_r_,"float",n_r,NULL,NULL,&mean,&stdev);
      for (nr=0;nr<n_r;nr++){
	f_rc_out__[nr+nc*n_r] = (f_r_[nr] - mean)/maximum(1e-12,stdev);
	/* for (nr=0;nr<n_r;nr++){ } */}
      /* for (nc=0;nc<n_c;nc++){ } */}
    /* if (f_rc_out__!=NULL){ } */}
}

void array_mean_center_row_d(int n_r,int n_c,double *d_rc_0in__,double **d_rc_out_p_)
{
  int nr=0,nc=0;
  double *d_rc_out__=NULL;
  double *d_r_=NULL;
  double mean=0,stdev=0;
  d_rc_out__=NULL;
  if (d_rc_out_p_!=NULL){
    if ( (*d_rc_out_p_)==NULL ){ (*d_rc_out_p_) = (double *) malloc1((unsigned long long int)n_r*(unsigned long long int)n_c*sizeof(double));}
    d_rc_out__ = *d_rc_out_p_;
    /* if (d_rc_out_p_!=NULL){ } */}
  if (d_rc_out__!=NULL){
    for (nc=0;nc<n_c;nc++){
      d_r_ = d_rc_0in__ + (unsigned long long int)nc*(unsigned long long int)n_r;
      array_stats(d_r_,"double",n_r,NULL,NULL,&mean,&stdev);
      for (nr=0;nr<n_r;nr++){
	d_rc_out__[nr+nc*n_r] = (d_r_[nr] - mean)/maximum(1e-12,stdev);
	/* for (nr=0;nr<n_r;nr++){ } */}
      /* for (nc=0;nc<n_c;nc++){ } */}
    /* if (d_rc_out__!=NULL){ } */}
}

void array_mean_center_row(int n_r,int n_c,void *v_0in__,char *type,void *v_out_p_)
{
  if (strcmp(type,"float")==0){ array_mean_center_row_f(n_r,n_c,(float *)v_0in__,(float **)v_out_p_);}
  if (strcmp(type,"double")==0){ array_mean_center_row_d(n_r,n_c,(double *)v_0in__,(double **)v_out_p_);}
}

void array_mean_center_row_test()
{

  int n_r = 4;
  int n_c = 3;
  float *f_A_rc__=NULL;
  double *d_A_rc__=NULL;
  float *f_B_rc__=NULL;
  double *d_B_rc__=NULL;
  float f_B_rc_ans__[12] = { -1 , -1 , +1 , +1 , -1 , -1 , +1 , +1 , -1 , -1 , +1 , +1 };
  double d_B_rc_ans__[12] = { -1 , -1 , +1 , +1 , -1 , -1 , +1 , +1 , -1 , -1 , +1 , +1 };
  int nr=0,nc=0;
  f_A_rc__ = (float *) malloc1(n_r*n_c*sizeof(float));
  f_A_rc__[0] = -2; f_A_rc__[1] = -2; f_A_rc__[2] = +2; f_A_rc__[3] = +2;
  for (nr=0;nr<n_r;nr++){ for (nc=1;nc<n_c;nc++){ f_A_rc__[nr+nc*n_r] = f_A_rc__[nr+0*n_r] + nc;}}
  array_printf(f_A_rc__,"float",n_r,n_c," % f_A_rc__: ");
  d_A_rc__ = (double *) malloc1(n_r*n_c*sizeof(double));
  d_A_rc__[0] = -2; d_A_rc__[1] = -2; d_A_rc__[2] = +2; d_A_rc__[3] = +2;
  for (nr=0;nr<n_r;nr++){ for (nc=1;nc<n_c;nc++){ d_A_rc__[nr+nc*n_r] = d_A_rc__[nr+0*n_r] + nc;}}
  array_printf(d_A_rc__,"double",n_r,n_c," % d_A_rc__: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_mean_center_row(n_r,n_c,f_A_rc__,"float",&f_B_rc__);
  array_printf(f_B_rc__,"float",n_r,n_c," % f_B_rc__: ");
  printf(" %% f_B_rc_ans__ vs f_B_rc__: relative error %0.16f\n",ffnormn(n_r*n_c,f_B_rc_ans__,f_B_rc__));
  free1(f_B_rc__);f_B_rc__=NULL;
  GLOBAL_toc(0,1," array_extract: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  array_mean_center_row(n_r,n_c,d_A_rc__,"double",&d_B_rc__);
  array_printf(d_B_rc__,"double",n_r,n_c," % d_B_rc__: ");
  printf(" %% d_B_rc_ans__ vs d_B_rc__: relative error %0.16f\n",dfnormn(n_r*n_c,d_B_rc_ans__,d_B_rc__));
  free1(d_B_rc__);d_B_rc__=NULL;
  GLOBAL_toc(0,1," array_extract: ");
  /* %%%%%%%% */
  free1(f_A_rc__);f_A_rc__=NULL;
  free1(d_A_rc__);d_A_rc__=NULL;
  
}
