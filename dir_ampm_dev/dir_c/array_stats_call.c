#ifndef _MONOLITH
#include "ampm_header.h"
#endif /* _MONOLITH */

/* %%%%%%%% */
/* index of max: */
/* %%%%%%%% */
int dmax_index(int n_a,double *a_)
{
  int na=0;
  int index=0;
  double a=0;
  if (n_a>0){
    a = a_[index];
    for (na=0;na<n_a;na++){ if (a_[na]>a){ a=a_[na]; index=na;}}
    /* if (n_a>0){ } */}
  return index;
}

void isum(unsigned long long int ulli_length,int *i_0in_,int *i_out_p)
{ 
  /* sum */
  unsigned long long int ulli=0;
  int sum=0;
  if (i_out_p!=NULL){ for (ulli=0;ulli<ulli_length;ulli++){ sum += i_0in_[ulli];} *i_out_p=sum;}
}

void icumsum(unsigned long long int ulli_length,int *i_0in_,int **i_out_p_)
{ 
  /* cumulative sum */
  unsigned long long int ulli=0;
  int sum=0;
  if (*i_out_p_==NULL){ *i_out_p_ = (int *) malloc1(ulli_length*sizeof(int));}
  if ((*i_out_p_)!=NULL){
    for (ulli=0;ulli<ulli_length;ulli++){ sum += i_0in_[ulli]; (*i_out_p_)[ulli]=sum;}
    /* if ((*i_out_p_)!=NULL){ } */}
}

void i0cumsum(unsigned long long int ulli_length,int *i_0in_,int **i_out_p_)
{ 
  /* cumulative sum, starting at 0 */
  unsigned long long int ulli=0;
  int sum=0;
  if (*i_out_p_==NULL){ *i_out_p_ = (int *) malloc1(((unsigned long long int)1+ulli_length)*sizeof(int));}
  if ((*i_out_p_)!=NULL){
    (*i_out_p_)[0] = 0;
    for (ulli=0;ulli<ulli_length;ulli++){ sum += i_0in_[ulli]; (*i_out_p_)[1+ulli]=sum;}
    /* if ((*i_out_p_)!=NULL){ } */}
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

double ullifnormn(unsigned long long int ulli_length,unsigned long long int *ulli_0_,unsigned long long int *ulli_1_)
{
  unsigned long long int ulli=0;
  unsigned long long int ulli_01 = 0;
  unsigned long long int ulli_0 = 0;
  double d_f01 = 0;
  double d_f0 = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ ulli_01 = ulli_0_[ulli]-ulli_1_[ulli]; d_f01 += ulli_01*ulli_01; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  d_f01 = sqrt(d_f01);
  for (ulli=0;ulli<ulli_length;ulli++){ ulli_0 = ulli_0_[ulli]; d_f0 += ulli_0*ulli_0; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
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

double dfnorm_single(unsigned long long int ulli_length,double *d_0_)
{
  double output = 0;
  unsigned long long int ulli=0;
  double d_curr = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ d_curr = d_0_[ulli]; output += d_curr*d_curr; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  return sqrt(output);
}

double dfnorm(unsigned long long int ulli_length,double *d_0_,double *d_1_)
{
  double output = 0;
  unsigned long long int ulli=0;
  double d_diff = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ d_diff = d_0_[ulli]-d_1_[ulli]; output += d_diff*d_diff; /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  return sqrt(output);
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

#ifdef _COMPLEX
double cfnorm(unsigned long long int ulli_length,float complex *c_0_,float complex *c_1_)
{
  double output = 0;
  unsigned long long int ulli=0;
  float complex c_diff = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ c_diff = c_0_[ulli]-c_1_[ulli]; output += crealf(c_diff*conjf(c_diff)); /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  return sqrt(output);
}
double cfnormn(unsigned long long int ulli_length,float complex *c_0_,float complex *c_1_)
{
  double c_f01 = 0;
  double c_f0 = 0;
  unsigned long long int ulli=0;
  float complex c_01=0,c_0=0;
  for (ulli=0;ulli<ulli_length;ulli++){ c_01 = c_0_[ulli]-c_1_[ulli]; c_f01 += creal(c_01*conj(c_01)); /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  c_f01 = sqrt(c_f01);
  for (ulli=0;ulli<ulli_length;ulli++){ c_0 = c_0_[ulli]; c_f0 += creal(c_0*conj(c_0)); /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  c_f0 = sqrt(c_f0);
  return c_f01/maximum(1e-12,c_f0);
}
#endif /* _COMPLEX */

#ifdef _COMPLEX
double zfnorm(unsigned long long int ulli_length,double complex *z_0_,double complex *z_1_)
{
  double output = 0;
  unsigned long long int ulli=0;
  double complex z_diff = 0;
  for (ulli=0;ulli<ulli_length;ulli++){ z_diff = z_0_[ulli]-z_1_[ulli]; output += creal(z_diff*conj(z_diff)); /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  return sqrt(output);
}
double zfnormn(unsigned long long int ulli_length,double complex *z_0_,double complex *z_1_)
{
  double z_f01 = 0;
  double z_f0 = 0;
  unsigned long long int ulli=0;
  double complex z_01=0,z_0=0;
  for (ulli=0;ulli<ulli_length;ulli++){ z_01 = z_0_[ulli]-z_1_[ulli]; z_f01 += creal(z_01*conj(z_01)); /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  z_f01 = sqrt(z_f01);
  for (ulli=0;ulli<ulli_length;ulli++){ z_0 = z_0_[ulli]; z_f0 += creal(z_0*conj(z_0)); /* for (ulli=0;ulli<ulli_length;ulli++){ } */}
  z_f0 = sqrt(z_f0);
  return z_f01/maximum(1e-12,z_f0);
}
#endif /* _COMPLEX */

void array_maximum_minimum(void *v_,const char *type,unsigned long long int ulli_length,void *max_p,int *index_max_p,void *min_p,int *index_min_p)
{
  /* finds the maximum and minimum of an array */
  unsigned long long int ulli=0;
  double mean=0,stdev=0;
  float f_max=0,f_min=0;
  double d_max=0,d_min=0;
  long long int lli_max=0,lli_min=0;
  int i_max=0,i_min=0;
  int index_max=0,index_min=0;
  if (v_!=NULL){
    if (0){ /* do nothing */}
    else if (strcmp(type,"float")==0){
      for (ulli=0;ulli<ulli_length;ulli++){ mean += ((float *)v_)[ulli];}
      mean /= (double)ulli_length;
      f_max=((float *)v_)[0]; f_min=((float *)v_)[0];
      index_max=0; index_min=0;
      for (ulli=0;ulli<ulli_length;ulli++){
	if (((float *)v_)[ulli]>f_max){ f_max=((float *)v_)[ulli]; index_max = ulli;}
	if (((float *)v_)[ulli]<f_min){ f_min=((float *)v_)[ulli]; index_min = ulli;}
	stdev += (((float *)v_)[ulli]-mean)*(((float *)v_)[ulli]-mean);}
      if (max_p!=NULL){*(float *)max_p=f_max;}
      if (min_p!=NULL){*(float *)min_p=f_min;}
      stdev /= (double)ulli_length;
      stdev = sqrt(stdev);}
    else if (strcmp(type,"double")==0){
      for (ulli=0;ulli<ulli_length;ulli++){ mean += ((double *)v_)[ulli];}
      mean /= (double)ulli_length;
      d_max=((double *)v_)[0]; d_min=((double *)v_)[0];
      index_max=0; index_min=0;
      for (ulli=0;ulli<ulli_length;ulli++){
	if (((double *)v_)[ulli]>d_max){ d_max=((double *)v_)[ulli]; index_max = ulli;}
	if (((double *)v_)[ulli]<d_min){ d_min=((double *)v_)[ulli]; index_min = ulli;}
	stdev += (((double *)v_)[ulli]-mean)*(((double *)v_)[ulli]-mean);}
      if (max_p!=NULL){*(double *)max_p=d_max;}
      if (min_p!=NULL){*(double *)min_p=d_min;}
      stdev /= (double)ulli_length;
      stdev = sqrt(stdev);}
    else if (strcmp(type,"long long int")==0){
      for (ulli=0;ulli<ulli_length;ulli++){ mean += ((long long int *)v_)[ulli];}
      mean /= (double)ulli_length;
      lli_max=((long long int *)v_)[0]; lli_min=((long long int *)v_)[0];
      index_max=0; index_min=0;
      for (ulli=0;ulli<ulli_length;ulli++){
	if (((long long int *)v_)[ulli]>lli_max){ lli_max=((long long int *)v_)[ulli]; index_max = ulli;}
	if (((long long int *)v_)[ulli]<lli_min){ lli_min=((long long int *)v_)[ulli]; index_min = ulli;}
	stdev += (((long long int *)v_)[ulli]-mean)*(((long long int *)v_)[ulli]-mean);}
      if (max_p!=NULL){*(long long int *)max_p=lli_max;}
      if (min_p!=NULL){*(long long int *)min_p=lli_min;}
      stdev /= (double)ulli_length;
      stdev = sqrt(stdev);}
    else if (strcmp(type,"int")==0){
      for (ulli=0;ulli<ulli_length;ulli++){ mean += ((int *)v_)[ulli];}
      mean /= (double)ulli_length;
      i_max=((int *)v_)[0]; i_min=((int *)v_)[0];
      index_max=0; index_min=0;
      for (ulli=0;ulli<ulli_length;ulli++){
	if (((int *)v_)[ulli]>i_max){ i_max=((int *)v_)[ulli]; index_max = ulli;}
	if (((int *)v_)[ulli]<i_min){ i_min=((int *)v_)[ulli]; index_min = ulli;}
	stdev += (((int *)v_)[ulli]-mean)*(((int *)v_)[ulli]-mean);}
      if (max_p!=NULL){*(int *)max_p=i_max;}
      if (min_p!=NULL){*(int *)min_p=i_min;}
      stdev /= (double)ulli_length;
      stdev = sqrt(stdev);}
    /* if (mean_p!=NULL){*mean_p=mean;} */
    /* if (stdev_p!=NULL){*stdev_p=stdev;} */
    /* if (v_!=NULL){ } */}
  if (index_max_p!=NULL){ *index_max_p = index_max;}
  if (index_min_p!=NULL){ *index_min_p = index_min;}
}

void array_stats(void *v_,const char *type,unsigned long long int ulli_length,void *max_p,void *min_p,double *mean_p,double *stdev_p)
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
    if (stdev_p!=NULL){*stdev_p=stdev;}
    /* if (v_!=NULL){ } */}
}
