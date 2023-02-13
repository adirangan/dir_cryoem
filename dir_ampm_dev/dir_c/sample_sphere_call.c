#ifndef _MONOLITH
#include "ampm_header.h"
#endif /* _MONOLITH */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* see sample_shell_5.m with TorL=='L' */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void sample_shell_L(
 double r
,double d
,int *n_all_p
,double **azimu_b_all_p_
,double **polar_a_all_p_
,double **weight_all_p_
,double **k_c_0_all_p_
,double **k_c_1_all_p_
,double **k_c_2_all_p_
,int *n_polar_a_p
,double **polar_a_p_
,int **n_azimu_b_p_
)
{
  int npolar_a=0;
  if (d<=0.0){ d = 1.0/(2.0*PI);}
  int n_equator = 3 + lround(2*PI*r/maximum(1e-16,d));
  int n_polar_a = 3 + lround((double)n_equator/2.0);
  double *lgnd_node_=NULL;
  double *lgnd_weight_=NULL;
  gsl_legpts(n_polar_a,&lgnd_node_,&lgnd_weight_);
  double *polar_a_=NULL;
  polar_a_ = (double *) malloc1(n_polar_a*sizeof(double));
  double *weight_=NULL;
  weight_ = (double *) malloc1(n_polar_a*sizeof(double));
  for (npolar_a=0;npolar_a<n_polar_a;npolar_a++){
    polar_a_[npolar_a] = acos(lgnd_node_[npolar_a]);
    weight_[npolar_a] = lgnd_weight_[npolar_a];
    /* for (npolar_a=0;npolar_a<n_polar_a;npolar_a++){ } */}
  int *n_azimu_b_=NULL;
  n_azimu_b_ = (int *) malloc1(n_polar_a*sizeof(int));
  int n_all=0;
  double polar_a=0,spolar_a=0;
  int n_azimu_b=0;
  for (npolar_a=0;npolar_a<n_polar_a;npolar_a++){
    polar_a = polar_a_[npolar_a];
    spolar_a = sin(polar_a);
    n_azimu_b = 3 + lround(2*PI*spolar_a*r/maximum(1e-16,d));
    n_azimu_b_[npolar_a] = n_azimu_b;
    n_all += n_azimu_b;
    /* for (npolar_a=0;npolar_a<n_polar_a;npolar_a++){ } */}
  double *azimu_b_all_=NULL;
  double *polar_a_all_=NULL;
  double *weight_all_=NULL;
  int ix=0;
  int nazimu_b=0;
  double azimu_b=0,dazimu_b=0;
  double weight=0;
  azimu_b_all_ = (double *) malloc1(n_all*sizeof(double));
  polar_a_all_ = (double *) malloc1(n_all*sizeof(double));
  weight_all_ = (double *) malloc1(n_all*sizeof(double));
  ix=0;
  for (npolar_a=0;npolar_a<n_polar_a;npolar_a++){
    polar_a = polar_a_[npolar_a];
    spolar_a = sin(polar_a);
    n_azimu_b = n_azimu_b_[npolar_a];
    dazimu_b = (2*PI)/maximum(1.0,n_azimu_b);
    weight = pow(r,2) * weight_[npolar_a] * dazimu_b;
    for (nazimu_b=0;nazimu_b<n_azimu_b;nazimu_b++){
      azimu_b = nazimu_b * dazimu_b;
      azimu_b_all_[ix+nazimu_b] = azimu_b;
      polar_a_all_[ix+nazimu_b] = polar_a;
      weight_all_[ix+nazimu_b] = weight;
      /* for (nazimu_b=0;nazimu_b<n_azimu_b;nazimu_b++){ } */}
    ix += n_azimu_b;
    /* for (npolar_a=0;npolar_a<n_polar_a;npolar_a++){ } */}
  if (n_all_p!=NULL){ *n_all_p = n_all;}
  if (azimu_b_all_p_!=NULL){
    if ((*azimu_b_all_p_)==NULL){ (*azimu_b_all_p_) = (double *) malloc1(n_all*sizeof(double)); }
    memcpy((*azimu_b_all_p_),azimu_b_all_,n_all*sizeof(double));
    /* if (azimu_b_all_p_!=NULL){ } */}
  if (polar_a_all_p_!=NULL){
    if ((*polar_a_all_p_)==NULL){ (*polar_a_all_p_) = (double *) malloc1(n_all*sizeof(double)); }
    memcpy((*polar_a_all_p_),polar_a_all_,n_all*sizeof(double));
    /* if (polar_a_all_p_!=NULL){ } */}
  if (weight_all_p_!=NULL){
    if ((*weight_all_p_)==NULL){ (*weight_all_p_) = (double *) malloc1(n_all*sizeof(double)); }
    memcpy((*weight_all_p_),weight_all_,n_all*sizeof(double));
    /* if (weight_all_p_!=NULL){ } */}
  if (k_c_0_all_p_!=NULL){
    if ((*k_c_0_all_p_)==NULL){ (*k_c_0_all_p_) = (double *) malloc1(n_all*sizeof(double)); }
    for (ix=0;ix<n_all;ix++){
      (*k_c_0_all_p_)[ix] = r*cos(azimu_b_all_[ix])*sin(polar_a_all_[ix]);
      /* for (ix=0;ix<n_all;ix++){ } */}
    /* if (k_c_0_all_p_!=NULL){ } */}
  if (k_c_1_all_p_!=NULL){
    if ((*k_c_1_all_p_)==NULL){ (*k_c_1_all_p_) = (double *) malloc1(n_all*sizeof(double)); }
    for (ix=0;ix<n_all;ix++){
      (*k_c_1_all_p_)[ix] = r*sin(azimu_b_all_[ix])*sin(polar_a_all_[ix]);
      /* for (ix=0;ix<n_all;ix++){ } */}
    /* if (k_c_1_all_p_!=NULL){ } */}
  if (k_c_2_all_p_!=NULL){
    if ((*k_c_2_all_p_)==NULL){ (*k_c_2_all_p_) = (double *) malloc1(n_all*sizeof(double)); }
    for (ix=0;ix<n_all;ix++){
      (*k_c_2_all_p_)[ix] = r*cos(polar_a_all_[ix]);
      /* for (ix=0;ix<n_all;ix++){ } */}
    /* if (k_c_2_all_p_!=NULL){ } */}
  if (n_polar_a_p!=NULL){ *n_polar_a_p = n_polar_a;}
  if (polar_a_p_!=NULL){
    if ((*polar_a_p_)==NULL){ (*polar_a_p_) = (double *) malloc1(n_polar_a*sizeof(double)); }
    memcpy((*polar_a_p_),polar_a_,n_polar_a*sizeof(double));
    /* if (polar_a_p_!=NULL){ } */}
  if (n_azimu_b_p_!=NULL){
    if ((*n_azimu_b_p_)==NULL){ (*n_azimu_b_p_) = (int *) malloc1(n_polar_a*sizeof(int)); }
    memcpy((*n_azimu_b_p_),n_azimu_b_,n_polar_a*sizeof(int));
    /* if (n_azimu_b_p_!=NULL){ } */}
  free1(&azimu_b_all_);
  free1(&polar_a_all_);
  free1(&weight_all_);
  free1(&n_azimu_b_);
  free1(&weight_);
  free1(&polar_a_);
  free1(&lgnd_node_);
  free1(&lgnd_weight_);
}

void sample_shell_L_test()
{
  int verbose=1;
  int n_d = 7;
  double ld_[7] = { -0.50 , -0.75 , -1.00 , -1.25 , -1.50 , -1.75 , -2.00 };
  double *d_=NULL;
  d_ = (double *) malloc1(n_d*sizeof(double));
  int nd=0;
  double r = 3.0;
  double d = 0.0;
  int n_equator=0,n_polar_a=0;
  int n_all=0,nall=0;
  double *azimu_b_all_=NULL;
  double *polar_a_all_=NULL;
  double *weight_all_=NULL;
  double *k_c_0_all_=NULL;
  double *k_c_1_all_=NULL;
  double *k_c_2_all_=NULL;
  double *polar_a_=NULL;
  int *n_azimu_b_=NULL;
  double Il=0,Ix=0;
  double f_m = 6;
  // f(polar_a) = exp(-(-f_m*cos(polar_a)).^2)*f_m;
  // F(polar_a) = sqrt(pi)/2*erf(-f_m*cos(polar_a));
  double g_w1=3,g_w2=5;
  // g(azimu_b) = (cos(g_w1*azimu_b) + 1) + (cos(g_w2*azimu_b) + 1);
  // G(azimu_b) = (sin(g_w1*azimu_b)/g_w1 + azimu_b) + (sin(g_w2*azimu_b)/g_w2 + azimu_b);
  double polar_a=0,azimu_b=0,f=0,g=0,h=0,weight=0,F1=0,F0=0,G1=0,G0=0,el=0,lel=0;
  GLOBAL_tic(0);
  for (nd=0;nd<n_d;nd++){
    d_[nd] = pow(2.0,ld_[nd]);
    d = d_[nd];
    n_equator = 3 + lround(2*PI*r/maximum(1e-16,d));
    n_polar_a = 3 + lround((double)n_equator/2.0);
    printf(" %% sampling at equatorial_distance d %0.3f (n_polar_a %d)\n",d,n_polar_a);
    sample_shell_L(r,d,&n_all,&azimu_b_all_,&polar_a_all_,&weight_all_,&k_c_0_all_,&k_c_1_all_,&k_c_2_all_,&n_polar_a,&polar_a_,&n_azimu_b_);
    if (verbose>1){
      printf(" %% r %0.6f d %0.6f n_all %d n_polar_a %d\n",r,d,n_all,n_polar_a);
      array_sub_printf(k_c_0_all_,"double",1,1,n_all,8," %% k_c_0_all_: ");
      array_sub_printf(k_c_1_all_,"double",1,1,n_all,8," %% k_c_1_all_: ");
      array_sub_printf(k_c_2_all_,"double",1,1,n_all,8," %% k_c_2_all_: ");
      array_sub_printf(azimu_b_all_,"double",1,1,n_all,8," %% azimu_b_all_: ");
      array_sub_printf(polar_a_all_,"double",1,1,n_all,8," %% polar_a_all_: ");
      array_sub_printf(weight_all_,"double",1,1,n_all,8," %% weight_all_: ");
      array_sub_printf(polar_a_,"double",1,1,n_polar_a,8," %% polar_a_: ");
      array_sub_printf(n_azimu_b_,"int",1,1,n_polar_a,8," %% n_azimu_b_: ");
      printf(" %%\n");
      printf(" %% exiting\n");
      exit(0);
      /* if (verbose>1){ } */}
    Il=0;
    for (nall=0;nall<n_all;nall++){
      polar_a = polar_a_all_[nall]; azimu_b = azimu_b_all_[nall]; weight = weight_all_[nall];
      f = exp(-pow(-f_m*cos(polar_a),2))*f_m;
      g = (cos(g_w1*azimu_b) + 1) + (cos(g_w2*azimu_b) + 1);
      h = f*g;
      Il += h*weight;
      /* for (nall=0;nall<n_all;nall++){ } */}
    polar_a = 1*PI; F1 = sqrt(PI)/2*erf(-f_m*cos(polar_a));
    polar_a = 0*PI; F0 = sqrt(PI)/2*erf(-f_m*cos(polar_a));
    azimu_b = 2*PI; G1 = (sin(g_w1*azimu_b)/g_w1 + azimu_b) + (sin(g_w2*azimu_b)/g_w2 + azimu_b);
    azimu_b = 0*PI; G0 = (sin(g_w1*azimu_b)/g_w1 + azimu_b) + (sin(g_w2*azimu_b)/g_w2 + azimu_b);
    Ix = pow(r,2) * (F1-F0) * (G1-G0);
    el = fabs(Il-Ix); lel = -log10(el);
    printf(" %% Ix %0.3f Il %0.3f error %0.16f (%0.2f)\n",Ix,Il,el,lel);
    free1(&k_c_0_all_);
    free1(&k_c_1_all_);
    free1(&k_c_2_all_);
    free1(&azimu_b_all_);
    free1(&polar_a_all_);
    free1(&weight_all_);
    free1(&polar_a_);
    free1(&n_azimu_b_);
    /* for (nd=0;nd<n_d;nd++){ } */}
  GLOBAL_toc(0,1," %% sample_sphere_L: ");
  free1(&d_);
}

