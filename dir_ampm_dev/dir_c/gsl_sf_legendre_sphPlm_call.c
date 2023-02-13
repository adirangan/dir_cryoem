#ifndef _MONOLITH
#include "ampm_header.h"
#endif /* _MONOLITH */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* testing gsl_sf_legendre_sphPlm.c */
/* associated with matlab legendre(..,'unnorm') */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void gsl_legendre_unnorm(int l_val,int n_x,double *x_,double **P_p_)
{
  int verbose=0;
  int n_m = 1+l_val;
  unsigned long long int n_P = (unsigned long long int)n_m * (unsigned long long int)n_x;
  double x=0;
  int nm=0,nx=0,na=0,m_val=0;
  double *P_=NULL;
  double P=0;
  double a=1;
  if (P_p_!=NULL){
    if ((*P_p_)==NULL){ (*P_p_) = (double *) malloc1(n_P*sizeof(double));}
    P_ = *P_p_;
    /* if (P_p_!=NULL){ } */}
  if (P_!=NULL){
    for (nx=0;nx<n_x;nx++){
      x = x_[nx];
      for (nm=0;nm<n_m;nm++){
	m_val = nm;
	a = 1.0/sqrt(4.0*PI);
	a *= sqrt(2*l_val + 1);
	a *= exp(0.5*gsl_sf_lnfact((unsigned int)(l_val-m_val)) - 0.5*gsl_sf_lnfact((unsigned int)(l_val+m_val)));
	P = gsl_sf_legendre_sphPlm(l_val,m_val,x);
	if (verbose>1){ printf(" %% l_val %d, m_val %d, x %f: P %f\n",l_val,m_val,x,P); }
	P_[nm + nx*n_m] = P / maximum(1e-16,a) ;
	/* for (nm=0;nm<n_m;nm++){ } */}
      /* for (nx=0;nx<n_x;nx++){ } */}
    /* if (P_!=NULL){ } */}
}

void gsl_legendre_unnorm_test()
{
  int l_val = 4;
  int n_x = 4;
  double x_[4] = { 0.0 , 0.1 , 0.2 , 0.3 };
  double *P_ = NULL;
  //double P_ans_[20] = { +0.7954951288348656 ,  -0.0000000000000000 ,  -0.8385254915624210 ,  +0.0000000000000000 ,  +1.1092649593311779 ,  +0.7168736936016866 ,  -0.3457136165780863 ,  -0.7720304200815212 ,  +0.3090530825157711 ,  +1.0871905866404876 ,  +0.4921463197058368 ,  -0.6320708821010504 ,  -0.5795888197679453 ,  +0.5902243641192727 ,  +1.0222985865196135 ,  +0.1547238025583815 ,  -0.8043064743538375 ,  -0.2823315330090673 ,  +0.8170782140116697 ,  +0.9185823128221485  }; //<-- tmp_ = legendre(4,0:0.1:0.3,'norm'); sprintf(' %+0.16f , ',tmp_);
  double P_ans_[20] = { +0.3749999999999999 ,  +0.0000000000000000 ,  -7.4999999999999991 ,  -0.0000000000000000 ,  +105.0000000000000142 ,  +0.3379375000000000 ,  +0.7288282976805994 ,  -6.9052500000000014 ,  -10.3428944087233159 ,  +102.9105000000000274 ,  +0.2319999999999999 ,  +1.3325224200740491 ,  -5.1839999999999993 ,  -19.7526852858035475 ,  +96.7680000000000007 ,  +0.0729375000000000 ,  +1.6956269305186216 ,  -2.5252500000000015 ,  -27.3446672086167517 ,  +86.9505000000000194 };//<-- tmp_ = legendre(4,0:0.1:0.3,'unnorm'); disp(sprintf(' %+0.16f , ',tmp_));
  /* %%%%%%%%; */
  GLOBAL_tic(0);
  gsl_legendre_unnorm(l_val,n_x,x_,&P_);
  array_printf(P_    ,"double",1+l_val,n_x," %% P_    : ");
  array_printf(P_ans_,"double",1+l_val,n_x," %% P_ans_: ");
  printf(" %% P_ vs P_ans_: %0.16f\n",dfnormn((unsigned long long int)(1+l_val)*n_x,P_,P_ans_));
  free1(&P_);
  GLOBAL_toc(0,1," gsl_legendre_unnorm: ");
  /* %%%%%%%%; */
}

/* 
%%%%%%%%;
% Warning! This has a bug when l_max is greater than about 80. ;
%%%%%%%%;
verbose=1;
n_w = 4; n_viewing_polar_a = 3;
n_x = n_w*n_viewing_polar_a; x_ = transpose(linspace(-1,1,n_x));
l_max = 4; n_m_max = 2*l_max+1;
legendre_evaluate_lmx___ = cell(l_max+1,1);
legendre_normalization_ = cell(l_max+1,1);
for l_val=0:l_max;
tmp_P__ = zeros(1+1*l_val,n_x);
tmp_a1 = ((1+2*l_val)/(4*pi));
m_val_ = -l_val:+l_val;
tmp_a2_ = exp(lfactorial(l_val-abs(m_val_)) - lfactorial(l_val+abs(m_val_)));
tmp_a3_ = sqrt(tmp_a1*tmp_a2_);
tmp_P__ = legendre(l_val,x_,'unnorm');
tmp_t = tic;
legendre_evaluate_lmx___{1+l_val} = reshape(tmp_P__,[1+1*l_val,n_x]);
tmp_t = toc(tmp_t);
if (verbose>1); disp(sprintf(' %% l_val %d/%d legendre_evaluate(%d,%d) %0.2fs',l_val,l_max,n_x,tmp_t)); end;
legendre_normalization_{1+l_val} = tmp_a3_;
end;%for l_val=0:l_max;
%%%%%%%%;
legendre_evaluate_normalized_lmx___ = zeros(l_max,n_m_max,n_x);
if (verbose); disp(sprintf(' %% legendre_evaluate_normalized_lmx___: (%d,%d,%d,%d)=%d (%0.2f GB)',(1+l_max),n_m_max,n_x,(1+l_max)*n_m_max*n_x,(1+l_max)*n_m_max*n_x*8/1e9)); end;
for l_val=0:l_max;
index_m_out_ = l_max + [-l_val:+l_val];
index_m_0in_ = [l_val:-1:1,0:l_val];
legendre_evaluate_normalized_lmx___(1+l_val,1+index_m_out_,:) = bsxfun(@times,legendre_evaluate_lmx___{1+l_val}(1+index_m_0in_,:),reshape(legendre_normalization_{1+l_val},[1+2*l_val,1,1]));
end;%for l_val=0:l_max;
%%%%%%%%;
legendre_evaluate_normalized_lxm___ = reshape(permute(legendre_evaluate_normalized_lmx___,[1,3,2]),[1+l_max,n_x,n_m_max]);
if (verbose); disp(sprintf(' %% legendre_evaluate_normalized_lxm___: (%d,%d,%d)=%d (%0.2f GB)',(1+l_max),n_x,n_m_max,(1+l_max)*n_x*n_m_max,(1+l_max)*n_x*n_m_max*8/1e9)); end;
disp(sprintf(' %+.16f , ',legendre_evaluate_normalized_lxm___(:)));
*/
void gsl_legendre_evaluate_normalized_lxm___(int l_max,int n_x,double *x_,double **P_lxm___p_)
{
  int verbose=0;
  if (verbose){ printf(" %% [entering gsl_legendre_evaluate_normalized_lxm___]\n");}
  int n_m_max = 2*l_max + 1;
  int n_l = 1+l_max;
  double *P_lxm___=NULL;
  unsigned long long int n_P = (unsigned long long int)n_l * (unsigned long long int)n_m_max * (unsigned long long int)n_x;
  double *P_ml__=NULL;
  //unsigned long long int n_lm2 = (unsigned long long int)(n_l*(1+n_l)/2);
  unsigned long long int n_lm2 = (unsigned long long int)(n_l*n_l); //<-- just in case. ;
  unsigned long long int n_lm = (unsigned long long int)(1+l_max) * (unsigned long long int)(1+l_max);
  int nx=0,nl=0,nm=0,l_val=0,m_val=0,m_abs=0,na=0;
  double x=0;
  if (P_lxm___p_!=NULL){
    if ((*P_lxm___p_)==NULL){ (*P_lxm___p_) = (double *) malloc1(n_P*sizeof(double));}
    P_lxm___ = *P_lxm___p_;
    /* if (P_lxm___p_!=NULL){ } */}
  if (P_lxm___!=NULL){
    P_ml__ = (double *) malloc1((unsigned long long int)n_lm2*sizeof(double));
    if (verbose>1){ printf(" n_lm2 %d\n",n_lm2);}
    for (nx=0;nx<n_x;nx++){
      x = x_[nx];
      gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM,l_max,x,-1,P_ml__);
      if (verbose>1){
	array_printf(P_ml__,"double",1,n_lm2," %% P_ml__:");
	/* if (verbose>1){ } */}
      for (nl=0;nl<n_l;nl++){
	l_val = nl;
	nm=0;
	for (m_val=-l_val;m_val<=+l_val;m_val++){
	  m_abs = abs(m_val);
	  na = l_val*(l_val+1)/2 + m_abs;
	  if (verbose>2){ printf(" l_val %d m_val %d m_abs %d: P[%d] -> %f\n",l_val,m_val,m_abs,na,P_ml__[na]);}
	  P_lxm___[nl + nx*n_l + (l_max + m_val)*n_x*n_l] = P_ml__[na];
	  nm++;
	  /* for (m_val=-l_val;m_val<=+l_val;m_val++){ } */}
	/* for (nl=0;nl<n_l;nl++){ } */}
      /* for (nx=0;nx<n_x;nx++){ } */}
    free1(&P_ml__);
    /* if (P_lxm___!=NULL){ } */}
  if (verbose){ printf(" %% [finished gsl_legendre_evaluate_normalized_lxm___]\n");}
}

/* %%%%%%%%; */
/* % Note: some kind of error when l_max is large */
/* %%%%%%%%; */
void gsl_legendre_evaluate_normalized_lxm___test_err()
{
  int verbose=2;
  int l_max = 98;
  int n_l = 1+l_max;
  int n_m_max = 2*l_max + 1,nm=0;
  int n_x = 12;
  unsigned long long int n_P = (unsigned long long int)n_l * (unsigned long long int)n_m_max * (unsigned long long int)n_x;
  double x_[12] = { -1.0000000000000000 , -0.8181818181818182 , -0.6363636363636364 , -0.4545454545454546 , -0.2727272727272727 , -0.0909090909090909 , +0.0909090909090909 , +0.2727272727272727 , +0.4545454545454546 , +0.6363636363636364 , +0.8181818181818182 , +1.0000000000000000 };
  double *P_lxm___=NULL;
  GLOBAL_tic(0);
  if (verbose){ printf(" %% n_P: %d\n",n_P);}
  gsl_legendre_evaluate_normalized_lxm___(l_max,n_x,x_,&P_lxm___);
  if (verbose){ array_printf_margin(P_lxm___,"double",n_l,n_x*n_m_max," %% P_lxm___    : ");}
  if (verbose>1){
    for (nm=0;nm<n_m_max;nm++){
      printf(" %% nm %d m_val %+d\n",nm,nm-l_max);
      array_printf_margin(P_lxm___ + nm*n_l*n_x,"double",n_l,n_x," %% P_lx__    : ");
      /* for (nm=0;nm<n_m_max;nm++){ } */}
    /* if (verbose>1){ } */}
  free1(&P_lxm___);
  GLOBAL_toc(0,1," gsl_legendre_evaluate_normalized_lxm___: ");
}

void gsl_legendre_evaluate_normalized_lxm___test()
{
  int verbose=0;
  int l_max = 4;
  int n_l = 1+l_max;
  int n_m_max = 2*l_max + 1;
  int n_x = 12;
  unsigned long long int n_P = (unsigned long long int)n_l * (unsigned long long int)n_m_max * (unsigned long long int)n_x;
  double x_[12] = { -1.0000000000000000 , -0.8181818181818182 , -0.6363636363636364 , -0.4545454545454546 , -0.2727272727272727 , -0.0909090909090909 , +0.0909090909090909 , +0.2727272727272727 , +0.4545454545454546 , +0.6363636363636364 , +0.8181818181818182 , +1.0000000000000000 };
  double *P_lxm___=NULL;
  double P_lxm_ans___[540] = { +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.0483609253406169 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.1566893981035988 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.2785589299619533 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.3791496546704366 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.4352483280655522 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.4352483280655522 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.3791496546704366 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.2785589299619533 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.1566893981035988 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.0483609253406169 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  -0.0793014321120351 ,  +0.1946489697295407 ,  +0 ,  +0 ,  +0 ,  -0.1915092643488430 ,  +0.3656085955750639 ,  +0 ,  +0 ,  +0 ,  -0.2948478008642978 ,  +0.4020651829967697 ,  +0 ,  +0 ,  +0 ,  -0.3715506944522550 ,  +0.3039960227336631 ,  +0 ,  +0 ,  +0 ,  -0.4120623285930571 ,  +0.1123806350708337 ,  +0 ,  +0 ,  +0 ,  -0.4120623285930571 ,  -0.1123806350708337 ,  +0 ,  +0 ,  +0 ,  -0.3715506944522550 ,  -0.3039960227336631 ,  +0 ,  +0 ,  +0 ,  -0.2948478008642978 ,  -0.4020651829967697 ,  +0 ,  +0 ,  +0 ,  -0.1915092643488430 ,  -0.3656085955750639 ,  +0 ,  +0 ,  +0 ,  -0.0793014321120351 ,  -0.1946489697295407 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.1276939510820461 ,  -0.2764198133102792 ,  +0.4076152700314881 ,  +0 ,  +0 ,  +0.2298491119476831 ,  -0.3869877386343912 ,  +0.3652086589788849 ,  +0 ,  +0 ,  +0.3064654825969107 ,  -0.3685597510803724 ,  +0.1184460515607194 ,  +0 ,  +0 ,  +0.3575430630297292 ,  -0.2579918257562607 ,  -0.1484231386841116 ,  +0 ,  +0 ,  +0.3830818532461384 ,  -0.0921399377700931 ,  -0.3125659693963431 ,  +0 ,  +0 ,  +0.3830818532461384 ,  +0.0921399377700931 ,  -0.3125659693963431 ,  +0 ,  +0 ,  +0.3575430630297292 ,  +0.2579918257562607 ,  -0.1484231386841116 ,  +0 ,  +0 ,  +0.3064654825969107 ,  +0.3685597510803724 ,  +0.1184460515607194 ,  +0 ,  +0 ,  +0.2298491119476831 ,  +0.3869877386343912 ,  +0.3652086589788849 ,  +0 ,  +0 ,  +0.1276939510820461 ,  +0.2764198133102792 ,  +0.4076152700314881 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  -0.1986451691985598 ,  +0.3634233559608495 ,  -0.4361290387784407 ,  +0.3752090158546363 ,  +0 ,  -0.2665104610379564 ,  +0.3792316866517962 ,  -0.2554786437276954 ,  -0.0383852944872703 ,  +0 ,  -0.3077397728442308 ,  +0.3127850233909290 ,  -0.0095161718534443 ,  -0.2976007785548513 ,  +0 ,  -0.3323969450650661 ,  +0.2027078630942048 ,  +0.1952941670415879 ,  -0.3077670442752591 ,  +0 ,  -0.3440635257300219 ,  +0.0699408574645979 ,  +0.3085427036731359 ,  -0.1260117996159369 ,  +0 ,  -0.3440635257300219 ,  -0.0699408574645979 ,  +0.3085427036731359 ,  +0.1260117996159369 ,  +0 ,  -0.3323969450650661 ,  -0.2027078630942048 ,  +0.1952941670415879 ,  +0.3077670442752591 ,  +0 ,  -0.3077397728442308 ,  -0.3127850233909290 ,  -0.0095161718534443 ,  +0.2976007785548513 ,  +0 ,  -0.2665104610379564 ,  -0.3792316866517962 ,  -0.2554786437276954 ,  +0.0383852944872703 ,  +0 ,  -0.1986451691985598 ,  -0.3634233559608495 ,  -0.4361290387784407 ,  -0.3752090158546363 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.2820947917738781 ,  -0.4886025119029199 ,  +0.6307831305050401 ,  -0.7463526651802308 ,  +0.8462843753216345 ,  +0.2820947917738781 ,  -0.3997656915569345 ,  +0.3179981071141112 ,  -0.1059809569639849 ,  -0.1479162431833931 ,  +0.2820947917738781 ,  -0.3109288712109490 ,  +0.0677700884013679 ,  +0.2315880170694482 ,  -0.3606289336542365 ,  +0.2820947917738781 ,  -0.2220920508649636 ,  -0.1199009256331894 ,  +0.3336437534051371 ,  -0.1802855656463478 ,  +0.2820947917738781 ,  -0.1332552305189781 ,  -0.2450149349895610 ,  +0.2674757485281519 ,  +0.1017899586736835 ,  +0.2820947917738781 ,  -0.0444184101729927 ,  -0.3075719396677468 ,  +0.1003734989235623 ,  +0.2913817045281305 ,  +0.2820947917738781 ,  +0.0444184101729927 ,  -0.3075719396677468 ,  -0.1003734989235623 ,  +0.2913817045281305 ,  +0.2820947917738781 ,  +0.1332552305189781 ,  -0.2450149349895610 ,  -0.2674757485281519 ,  +0.1017899586736835 ,  +0.2820947917738781 ,  +0.2220920508649636 ,  -0.1199009256331894 ,  -0.3336437534051371 ,  -0.1802855656463478 ,  +0.2820947917738781 ,  +0.3109288712109490 ,  +0.0677700884013679 ,  -0.2315880170694482 ,  -0.3606289336542365 ,  +0.2820947917738781 ,  +0.3997656915569345 ,  +0.3179981071141112 ,  +0.1059809569639849 ,  -0.1479162431833931 ,  +0.2820947917738781 ,  +0.4886025119029199 ,  +0.6307831305050401 ,  +0.7463526651802308 ,  +0.8462843753216345 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  -0.1986451691985598 ,  +0.3634233559608495 ,  -0.4361290387784407 ,  +0.3752090158546363 ,  +0 ,  -0.2665104610379564 ,  +0.3792316866517962 ,  -0.2554786437276954 ,  -0.0383852944872703 ,  +0 ,  -0.3077397728442308 ,  +0.3127850233909290 ,  -0.0095161718534443 ,  -0.2976007785548513 ,  +0 ,  -0.3323969450650661 ,  +0.2027078630942048 ,  +0.1952941670415879 ,  -0.3077670442752591 ,  +0 ,  -0.3440635257300219 ,  +0.0699408574645979 ,  +0.3085427036731359 ,  -0.1260117996159369 ,  +0 ,  -0.3440635257300219 ,  -0.0699408574645979 ,  +0.3085427036731359 ,  +0.1260117996159369 ,  +0 ,  -0.3323969450650661 ,  -0.2027078630942048 ,  +0.1952941670415879 ,  +0.3077670442752591 ,  +0 ,  -0.3077397728442308 ,  -0.3127850233909290 ,  -0.0095161718534443 ,  +0.2976007785548513 ,  +0 ,  -0.2665104610379564 ,  -0.3792316866517962 ,  -0.2554786437276954 ,  +0.0383852944872703 ,  +0 ,  -0.1986451691985598 ,  -0.3634233559608495 ,  -0.4361290387784407 ,  -0.3752090158546363 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.1276939510820461 ,  -0.2764198133102792 ,  +0.4076152700314881 ,  +0 ,  +0 ,  +0.2298491119476831 ,  -0.3869877386343912 ,  +0.3652086589788849 ,  +0 ,  +0 ,  +0.3064654825969107 ,  -0.3685597510803724 ,  +0.1184460515607194 ,  +0 ,  +0 ,  +0.3575430630297292 ,  -0.2579918257562607 ,  -0.1484231386841116 ,  +0 ,  +0 ,  +0.3830818532461384 ,  -0.0921399377700931 ,  -0.3125659693963431 ,  +0 ,  +0 ,  +0.3830818532461384 ,  +0.0921399377700931 ,  -0.3125659693963431 ,  +0 ,  +0 ,  +0.3575430630297292 ,  +0.2579918257562607 ,  -0.1484231386841116 ,  +0 ,  +0 ,  +0.3064654825969107 ,  +0.3685597510803724 ,  +0.1184460515607194 ,  +0 ,  +0 ,  +0.2298491119476831 ,  +0.3869877386343912 ,  +0.3652086589788849 ,  +0 ,  +0 ,  +0.1276939510820461 ,  +0.2764198133102792 ,  +0.4076152700314881 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  -0.0793014321120351 ,  +0.1946489697295407 ,  +0 ,  +0 ,  +0 ,  -0.1915092643488430 ,  +0.3656085955750639 ,  +0 ,  +0 ,  +0 ,  -0.2948478008642978 ,  +0.4020651829967697 ,  +0 ,  +0 ,  +0 ,  -0.3715506944522550 ,  +0.3039960227336631 ,  +0 ,  +0 ,  +0 ,  -0.4120623285930571 ,  +0.1123806350708337 ,  +0 ,  +0 ,  +0 ,  -0.4120623285930571 ,  -0.1123806350708337 ,  +0 ,  +0 ,  +0 ,  -0.3715506944522550 ,  -0.3039960227336631 ,  +0 ,  +0 ,  +0 ,  -0.2948478008642978 ,  -0.4020651829967697 ,  +0 ,  +0 ,  +0 ,  -0.1915092643488430 ,  -0.3656085955750639 ,  +0 ,  +0 ,  +0 ,  -0.0793014321120351 ,  -0.1946489697295407 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.0483609253406169 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.1566893981035988 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.2785589299619533 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.3791496546704366 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.4352483280655522 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.4352483280655522 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.3791496546704366 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.2785589299619533 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.1566893981035988 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0.0483609253406169 ,  +0 ,  +0 ,  +0 ,  +0 ,  +0 };
  GLOBAL_tic(0);
  if (verbose){ printf(" %% n_P: %d\n",n_P);}
  gsl_legendre_evaluate_normalized_lxm___(l_max,n_x,x_,&P_lxm___);
  if (verbose>1){ array_printf(P_lxm___    ,"double",1,n_P," %% P_lxm___    : ");}
  if (verbose>1){ array_printf(P_lxm_ans___,"double",1,n_P," %% P_lxm_ans___: ");}
  printf(" %% P_lxm___ vs P_lxm_ans___: %0.16f\n",dfnormn(n_P,P_lxm___,P_lxm_ans___));
  free1(&P_lxm___);
  GLOBAL_toc(0,1," gsl_legendre_evaluate_normalized_lxm___: ");
}

void pm_template(
 int verbose
,int l_max
,int n_a
,double *a_k_Y_ya_real__
,double *a_k_Y_ya_imag__
,double viewing_k_eq_d
,double template_k_eq_d
,int n_w_0in
,float **template_waS_real___p_
,float **template_waS_imag___p_
,int *n_w_p
,int *n_viewing_all_p
,double **viewing_azimu_b_all_p_
,double **viewing_polar_a_all_p_
,double **viewing_weight_all_p_
)
{
  /* 
     % uses spherical-harmonic-expansions a_k_Y_ya__ to evaluate templates on a collection of points on spherical shells. ;
     % each spherical-shell has the same resolution, determined by viewing_k_eq_d and template_k_eq_d and/or n_w_max. ;
     % ;
     % inputs: ;
     % ;
     % verbose = integer verbosity_level. ;
     % l_max = spherical harmonic order on each shell. ;
     %         corresponds to n_lm = (l_max+1)^2 coefficients. ;
     % n_a = integer number of shells. ;
     % a_k_Y_ya__ = complex array of size (n_lm,n_a). ;
     %              coefficients are ordered in a row, with m varying quickly, l varying slowly. ;
     % viewing_k_eq_d = real equatorial-distance used for sampling viewing angles and templates. ;
     % template_k_eq_d = real equatorial-distance used for sampling inplane-shifts along each template. ;
     % n_w_0in = integer. used if template_k_eq_d <=0; desired n_w for templates. ;
     % ;
     % outputs: ;
     % ;
     % template_waS___ = complex array of templates for each viewing angle. ;
  */
  if (verbose){ printf(" %% [entering pm_template]\n");}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ printf(" %% l_max %d n_a %d viewing_k_eq_d %0.6f template_k_eq_d %0.6f n_w_0in %d\n",l_max,n_a,viewing_k_eq_d,template_k_eq_d,n_w_0in);}
  unsigned long long int tab=0,tab_0in=0,tab_out=0,na=0,tab_num=0;
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  double elrt_0=0,elrt_1=0,elrt_2=0,elrt_3=0,elrt_4=0,elrt_5=0;
  /* %%%%; */
  int n_lm = (l_max+1)*(l_max+1);
  int n_m_max = 2*l_max+1;
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ printf(" %% n_lm %d n_m_max %d\n",n_lm,n_m_max);}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(a_k_Y_ya_real__,"double",n_lm,n_a," %% a_k_Y_ya_real__: ");}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(a_k_Y_ya_imag__,"double",n_lm,n_a," %% a_k_Y_ya_imag__: ");}
  int *m_max_=NULL;
  m_max_ = (int *) malloc1(n_m_max*sizeof(int));
  int nm=0,nl=0,m_val=0,m_abs=0,l_val=0;
  for (nm=0;nm<n_m_max;nm++){ m_max_[nm] = -l_max + nm;}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(m_max_,"int",1,n_m_max," %% m_max_: ");}
  /* %%%%%%%%; */
  if (verbose>1){ printf(" %% \t \n");}
  if (verbose>1){ printf(" %% First determine the viewing angles.\n");}
  double k_p_r = 1.0;
  int n_viewing_all=0,nviewing_all;
  double *viewing_azimu_b_all_=NULL,viewing_azimu_b=0;
  double *viewing_polar_a_all_=NULL,viewing_polar_a=0;
  double *viewing_weight_all_=NULL,viewing_weight=0;
  int n_viewing_polar_a=0,nviewing_polar_a=0;
  double *viewing_polar_a_=NULL;
  int *n_viewing_azimu_b_=NULL,n_viewing_azimu_b=0,nviewing_azimu_b=0;
  local_tic(0,t_start_,d_start_);
  sample_shell_L(
		 k_p_r
		 ,viewing_k_eq_d
		 ,&n_viewing_all
		 ,&viewing_azimu_b_all_
		 ,&viewing_polar_a_all_
		 ,&viewing_weight_all_
		 ,NULL
		 ,NULL
		 ,NULL
		 ,&n_viewing_polar_a
		 ,&viewing_polar_a_
		 ,&n_viewing_azimu_b_
		 );
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," %% sample_shell_L: ");
  int n_viewing_azimu_b_sum=0;
  isum(n_viewing_polar_a,n_viewing_azimu_b_,&n_viewing_azimu_b_sum);
  int *n_viewing_azimu_b_csum_=NULL;
  i0cumsum(n_viewing_polar_a,n_viewing_azimu_b_,&n_viewing_azimu_b_csum_);
  if (verbose>1){ printf(" %% \t \n");}
  if (verbose>1){ printf(" %% n_viewing_all %d n_viewing_polar_a %d n_viewing_azimu_b_sum %d\n",n_viewing_all,n_viewing_polar_a,n_viewing_azimu_b_sum); /* if (verbose>1){ } */}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(viewing_azimu_b_all_,"double",1,n_viewing_all," %% viewing_azimu_b_all_: ");}
  if (verbose>2){ array_printf_margin(viewing_polar_a_all_,"double",1,n_viewing_all," %% viewing_polar_a_all_: ");}
  if (verbose>2){ array_printf_margin(viewing_weight_all_,"double",1,n_viewing_all," %% viewing_weight_all_: ");}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(viewing_polar_a_,"double",1,n_viewing_polar_a," %% viewing_polar_a_: ");}
  if (verbose>2){ array_printf_margin(n_viewing_azimu_b_,"int",1,n_viewing_polar_a," %% n_viewing_azimu_b_: ");}
  /* %%%%%%%%; */
  if (verbose>1){ printf(" %% \t \n");}
  if (verbose>1){ printf(" %% Now determine the points along each equatorial plane (i.e., the points for each template).\n");}
  int n_w = 0;
  int n_equator = 0;
  int n_polar_a = 0;
  if (template_k_eq_d> 0){
    k_p_r = 1.0;
    n_equator = 3 + lround(2*PI*k_p_r/maximum(1e-16,template_k_eq_d));
    n_polar_a = 3 + lround((double)n_equator/2.0);
    n_w = 2*n_polar_a;
    /* if (template_k_eq_d> 0){ } */}
  if (template_k_eq_d<=0){
    n_w = maximum(6,n_w_0in);
    /* if (template_k_eq_d<=0){ } */}
  if (verbose>1){ printf(" %% \t \n");}
  if (verbose>1){ printf(" %% n_w %d\n",n_w);}
  /* %%%%%%%%; */
  if (verbose>1){ printf(" %% \t \n");}
  if (verbose>1){ printf(" %% Set up inner gamma_z for the templates.\n");}
  double *gamma_z_ = NULL,gamma_z=0;
  double *cc_ = NULL,cc=0;
  double *sc_ = NULL,sc=0;
  gamma_z_ = (double *) malloc1(n_w*sizeof(double));
  cc_ = (double *) malloc1(n_w*sizeof(double));
  sc_ = (double *) malloc1(n_w*sizeof(double));
  int nw=0;
  for (nw=0;nw<n_w;nw++){ 
    gamma_z = (double)nw*2.0*PI/(double)maximum(1,n_w);
    gamma_z_[nw] = gamma_z; cc_[nw] = cos(gamma_z); sc_[nw] = sin(gamma_z);
    /* for (nw=0;nw<n_w;nw++){ } */}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(gamma_z_,"double",1,n_w," %% gamma_z_: ");}
  if (verbose>2){ array_printf_margin(cc_,"double",1,n_w," %% cc_: ");}
  if (verbose>2){ array_printf_margin(sc_,"double",1,n_w," %% sc_: ");}
  /*
    %%%%%%%%;
    % The general formula used here is as follows. ;
    % let sa and ca be sin(polar_a) and cos(polar_a), respectively. ;
    % let sb and cb be sin(azimu_b) and cos(azimu_b), respectively. ;
    % let sc and cc be sin(gamma_z) and cos(gamma_z), respectively. ;
    % And rotation by azimu_b about the +z-axis is represented as: ;
    % Rz(azimu_b) = ;
    % [ +cb -sb 0 ] ;
    % [ +sb +cb 0 ] ;
    % [  0   0  1 ] ;
    % And rotation by polar_a about the +y-axis is represented as: ;
    % Ry(polar_a) = ;
    % [ +ca 0 +sa ] ;
    % [  0  1  0  ] ;
    % [ -sa 0 +ca ] ;
    % And rotation by gamma_z about the +z-axis is represented as: ;
    % Rz(gamma_z) = ;
    % [ +cc -sc 0 ] ;
    % [ +sc +cc 0 ] ;
    % [  0   0  1 ] ;
    % Which, collectively, implies that under the transform: ;
    % Rz(azimu_b) * Ry(polar_a) * Rz(gamma_z), ;
    % Which is the same as: ;
    % [ +cb -sb 0 ] [ +ca*cc -ca*sc +sa ]   [ +cb*ca*cc - sb*sc , -cb*ca*sc -sb*cc , +cb*sa ];
    % [ +sb +cb 0 ] [ +sc    +cc    0   ] = [ +sb*ca*cc + cb*sc , -sb*ca*sc +cb*cc , +sb*sa ];
    % [  0   0  1 ] [ -sa*cc +sa*sc +ca ]   [ -sa*cc            , +sa*sc           , +ca    ];
    % the point [1;0;0] is mapped to: ;
    % [ template_k_c_0 ; template_k_c_1 ; template_k_c_2 ] = [ +cb*ca*cc - sb*sc ; +sb*ca*cc + cb*sc ; -sa*cc ];
    %%%%%%%%;
  */
  if (verbose>1){ printf(" %% \t \n");}
  if (verbose>1){ printf(" %% template: (%d,%d,%d)=%d (%0.2f GB)\n",n_w,n_viewing_all,n_a,n_w*n_viewing_all*n_a,n_w*n_viewing_all*(double)16/(double)1e9);}
  /* %%%%%%%%; */
  if (verbose>1){ printf(" %% \t \n");}
  if (verbose>1){ printf(" %% Now construct array of k_c_?_ values for the templates.\n");}
  double template_k_p_r=1.0;
  double *template_k_c_0__=NULL;
  double *template_k_c_1__=NULL;
  double *template_k_c_2__=NULL;
  double *template_azimu_b__=NULL;
  double *expi_template_azimu_b_real__=NULL;
  double *expi_template_azimu_b_imag__=NULL;
  template_k_c_0__ = (double *) malloc1((unsigned long long int)n_w*(unsigned long long int)n_viewing_all*sizeof(double));
  template_k_c_1__ = (double *) malloc1((unsigned long long int)n_w*(unsigned long long int)n_viewing_all*sizeof(double));
  template_k_c_2__ = (double *) malloc1((unsigned long long int)n_w*(unsigned long long int)n_viewing_all*sizeof(double));
  template_azimu_b__ = (double *) malloc1((unsigned long long int)n_w*(unsigned long long int)n_viewing_all*sizeof(double));
  expi_template_azimu_b_real__ = (double *) malloc1((unsigned long long int)n_w*(unsigned long long int)n_viewing_all*sizeof(double));
  expi_template_azimu_b_imag__ = (double *) malloc1((unsigned long long int)n_w*(unsigned long long int)n_viewing_all*sizeof(double));
  if (verbose>1){ printf(" %% \t \n");}
  if (verbose>1){ printf(" %% template_k_c___: (%d,%d)=%d (%0.2f GB)\n",n_w,n_viewing_all,n_w*n_viewing_all,n_w*n_viewing_all*(double)8/(double)1e9);}
  double ca=0,sa=0,cb=0,sb=0;
  local_tic(0,t_start_,d_start_);
  for (nviewing_all=0;nviewing_all<n_viewing_all;nviewing_all++){
    viewing_polar_a = viewing_polar_a_all_[nviewing_all]; ca = cos(viewing_polar_a); sa = sin(viewing_polar_a);
    viewing_azimu_b = viewing_azimu_b_all_[nviewing_all]; cb = cos(viewing_azimu_b); sb = sin(viewing_azimu_b);
    for (nw=0;nw<n_w;nw++){
      cc=cc_[nw]; sc=sc_[nw];
      tab = (unsigned long long int)nw+(unsigned long long int)nviewing_all*(unsigned long long int)n_w;
      template_k_c_0__[tab] = (+cb*ca*cc - sb*sc)*template_k_p_r;
      template_k_c_1__[tab] = (+sb*ca*cc + cb*sc)*template_k_p_r;
      template_k_c_2__[tab] = (-sa*cc           )*template_k_p_r;
      template_azimu_b__[tab] = atan2(template_k_c_1__[tab],template_k_c_0__[tab]);
      expi_template_azimu_b_real__[tab] = cos(template_azimu_b__[tab]);
      expi_template_azimu_b_imag__[tab] = sin(template_azimu_b__[tab]);
      /* for (nw=0;nw<n_w;nw++){ } */}
    /* for (nviewing_all=0;nviewing_all<n_viewing_all;nviewing_all++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," %% expi_template_azimu_b__: ");
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(template_k_c_0__,"double",n_w,n_viewing_all," %% template_k_c_0__: ");}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(template_k_c_1__,"double",n_w,n_viewing_all," %% template_k_c_1__: ");}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(template_k_c_2__,"double",n_w,n_viewing_all," %% template_k_c_2__: ");}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(template_azimu_b__,"double",n_w,n_viewing_all," %% template_azimu_b__: ");}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(expi_template_azimu_b_real__,"double",n_w,n_viewing_all," %% expi_template_azimu_b_real__: ");}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(expi_template_azimu_b_imag__,"double",n_w,n_viewing_all," %% expi_template_azimu_b_imag__: ");}
  /* %%%%%%%%; */
  if (verbose>1){ printf(" %% \t \n");}
  if (verbose>1){ printf(" %% We also use a condensed array, called condense_k_c_2__, which only depends on the polar_a, and not on azimu_b.\n");}
  double *condense_k_c_2__=NULL;
  condense_k_c_2__ = (double *) malloc1((unsigned long long int)n_w*(unsigned long long int)n_viewing_polar_a*sizeof(double));
  if (verbose>1){ printf(" %% \t \n");}
  if (verbose>1){ printf(" %% condense_k_c_2__: (%d,%d)=%d (%0.2f GB)\n",n_w,n_viewing_polar_a,n_w*n_viewing_polar_a,n_w*n_viewing_polar_a*(double)8/(double)1e9);}
  local_tic(0,t_start_,d_start_);
  for (nviewing_polar_a=0;nviewing_polar_a<n_viewing_polar_a;nviewing_polar_a++){
    viewing_polar_a = viewing_polar_a_[nviewing_polar_a]; ca = cos(viewing_polar_a); sa = sin(viewing_polar_a);
    for (nw=0;nw<n_w;nw++){
      cc=cc_[nw]; sc=sc_[nw];
      tab = (unsigned long long int)nw+(unsigned long long int)nviewing_polar_a*(unsigned long long int)n_w;
      condense_k_c_2__[tab] = -sa*cc ;
      /* for (nw=0;nw<n_w;nw++){ } */}
    /* for (nviewing_polar_a=0;nviewing_polar_a<n_viewing_polar_a;nviewing_polar_a++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," %% condense_k_c_2__: ");
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(condense_k_c_2__,"double",n_w,n_viewing_polar_a," %% condense_k_c_2__: ");}
  /* %%%%%%%%; */
  if (verbose>1){ printf(" %% \t \n");}
  if (verbose>1){ printf(" %% Now evaluate associated legendre polynomials at the various k_c_2 values.\n");}
  /* 
     % Here legendre_evaluate_lmwS____[l_val,l_val+m_val,nw,nviewing_polar_a] contains ;
     % the associated legendre-function of degree l_val and order abs(m_val) (ranging from 0 to +l_val) ;
     % evaluated at the k_c_2 value stored in condense_k_c_2__[nw,nviewing_polar_a]. ;
     % Note that this is associated with viewing_polar_a_[nviewing_polar_a]. ;
     % The legendre_normalization_[l_val,abs(m_val)] contains ;
     % The normalization coefficient for the spherical harmonics associated with l_val and m_val. ;
     % Note that this is somewhat redundant (as it does not depend explicitly on the shell). ;
  */
  unsigned long long int n_lwpm = (unsigned long long int)(1+l_max)*(unsigned long long int)n_w*(unsigned long long int)n_viewing_polar_a*(unsigned long long int)n_m_max;
  double *legendre_evaluate_normalized_lwpm___=NULL;
  local_tic(0,t_start_,d_start_);
  gsl_legendre_evaluate_normalized_lxm___(l_max,n_w*n_viewing_polar_a,condense_k_c_2__,&legendre_evaluate_normalized_lwpm___);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," %% legendre_evaluate_normalized_lwpm___: ");
  if (verbose>1){ printf(" %% \t \n");}
  if (verbose>1){ printf(" %% legendre_evaluate_normalized_lwpm____: (%d,%d,%d,%d)=%d (%0.2f GB)\n",(1+l_max),n_w,n_viewing_polar_a,n_m_max,n_lwpm,n_lwpm*(double)8/(double)1e9);}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(legendre_evaluate_normalized_lwpm___,"double",(1+l_max)*n_w,n_viewing_polar_a*n_m_max," %% legendre_evaluate_normalized_lwpm___: ");}
  /* %%%%%%%%; */
  /* Note: This will store an aligned version of legendre_evaluate_normalized_lwpm___. ; */
  /* Note: If the brute-force alignment (below) is slow we can instead ; */
  /* directly generate an aligned version of legendre_evaluate_normalized_lwpm___ above. ; */
  unsigned long long int n_row_A0 = (unsigned long long int)n_w*(unsigned long long int)n_viewing_polar_a,nrow_A0=0;
  unsigned long long int n_row_A0_pack = (unsigned long long int)n_w*(unsigned long long int)n_viewing_polar_a*(unsigned long long int)n_m_max;
  unsigned long long int n_col_X0 = 1+l_max,ncol_X0=0;
  unsigned long long int n_col_X0_rup = rup(n_col_X0,8ULL);
  unsigned long long int n_col_X0_256 = n_col_X0_rup/8ULL; //%<-- 8 floats per __m256. ;
  __m256 *ps_A0R_trn__ = NULL;
  ps_A0R_trn__ = (__m256 *) _mm_malloc(n_col_X0_256*n_row_A0_pack*sizeof(__m256),32);
  float *f_A0R_trn__ = NULL; 
  f_A0R_trn__ = (float *) ps_A0R_trn__;
  memset(f_A0R_trn__,0,n_col_X0_rup*n_row_A0_pack*sizeof(float));
  unsigned long long int ulli_0in=0,ulli_out=0;
  local_tic(0,t_start_,d_start_);
  for (nrow_A0=0;nrow_A0<n_row_A0_pack;nrow_A0++){ 
    ulli_0in = nrow_A0*n_col_X0; ulli_out = nrow_A0*n_col_X0_rup;
    for (ncol_X0=0;ncol_X0<n_col_X0;ncol_X0++){
      f_A0R_trn__[ulli_out] = (float)(legendre_evaluate_normalized_lwpm___[ulli_0in]);
      ulli_out++; ulli_0in++;
      /* for (ncol_X0=0;ncol_X0<n_col_X0;ncol_X0++){ } */}
    /* for (nrow_A0=0;nrow_A0<n_row_A0_pack;nrow_A0++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," %% f_A0R_trn__: ");
  /* %%%%%%%%; */
  if (verbose>1){ printf(" %% \t \n");}
  if (verbose>1){ printf(" %% unroll a_k_Y_ya__.\n");}
  double *a_k_Y_lam_real___=NULL;
  double *a_k_Y_lam_imag___=NULL;
  unsigned long long int n_lam = (unsigned long long int)(1+l_max)*(unsigned long long int)n_a*(unsigned long long int)n_m_max;
  a_k_Y_lam_real___ = (double *) malloc1(n_lam*sizeof(double));
  a_k_Y_lam_imag___ = (double *) malloc1(n_lam*sizeof(double));
  if (verbose>1){ printf(" %% \t \n");}
  if (verbose>1){ printf(" %% a_k_Y_lam___: (%d,%d,%d)=%lld (%0.2f GB)\n",(1+l_max),n_a,n_m_max,n_lam,n_lam*(double)16/(double)1e9);}
  int index_m_0in=0,index_m_out=0;
  local_tic(0,t_start_,d_start_);
  for (l_val=0;l_val<=l_max;l_val++){
    for (m_val=-l_val;m_val<=+l_val;m_val++){
      index_m_0in = l_val*l_val + l_val + m_val;
      index_m_out = l_max + m_val;
      for (na=0;na<n_a;na++){
	tab_out = (unsigned long long int)l_val + (unsigned long long int)na*(unsigned long long int)(1+l_max) + (unsigned long long int)index_m_out*(unsigned long long int)(1+l_max)*(unsigned long long int)n_a;
	tab_0in = (unsigned long long int)index_m_0in + (unsigned long long int)na*(unsigned long long int)n_lm;
	a_k_Y_lam_real___[tab_out] = a_k_Y_ya_real__[tab_0in];
	a_k_Y_lam_imag___[tab_out] = a_k_Y_ya_imag__[tab_0in];
	/*for (na=0;na<n_a;na++){ } */}
      /* for (m_val=-l_val;m_val<=+l_val;m_val++){ } */}
    /* for (l_val=0;l_val<=l_max;l_val++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," %% a_k_Y_lam___: ");
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(a_k_Y_lam_real___,"double",(1+l_max)*n_a,n_m_max," %% a_k_Y_lam_real___: ");}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(a_k_Y_lam_imag___,"double",(1+l_max)*n_a,n_m_max," %% a_k_Y_lam_imag___: ");}
  /* %%%%%%%%; */
  double tmp_d = (double)n_w*(double)n_viewing_polar_a*(double)n_m_max*(double)8;
  int n_a_per_a_batch = minimum(n_a,maximum(1,ceil((double)0.5e9 / tmp_d ))); //<-- 0.5GB limit. ;
  if (verbose>2){ n_a_per_a_batch = minimum(8,n_a_per_a_batch);}
  int n_a_batch = ceil((double)n_a/(double)n_a_per_a_batch),na_batch=0;
  if (verbose>1){ printf(" %% n_a_per_a_batch %d, n_a_batch %d\n",n_a_per_a_batch,n_a_batch);}
  int *index_a_=NULL;
  index_a_ = (int *) malloc1((unsigned long long int)n_a_per_a_batch*sizeof(int));
  int tmp_n_a=0,tmp_na=0;
  /* %%%%%%%%; */
  /* Note: These will store aligned versions of: ; */
  /* tmp_a_k_Y_lam_real___ and tmp_a_k_Y_lam_imag___, ; */
  /* associated with the batch under consideration. */
  /* These will be formed below. ; */
  unsigned long long int n_row_B0 = (unsigned long long int)n_a_per_a_batch,nrow_B0=0;
  unsigned long long int n_row_B0_pack = (unsigned long long int)n_a_per_a_batch*(unsigned long long int)n_m_max;
  __m256 *ps_B0R_trn__ = NULL;
  __m256 *ps_B0I_trn__ = NULL;
  ps_B0R_trn__ = (__m256 *) _mm_malloc(n_col_X0_256*n_row_B0_pack*sizeof(__m256),32);
  ps_B0I_trn__ = (__m256 *) _mm_malloc(n_col_X0_256*n_row_B0_pack*sizeof(__m256),32);
  float *f_B0R_trn__ = NULL; float *f_B0I_trn__ = NULL;
  f_B0R_trn__ = (float *) ps_B0R_trn__; f_B0I_trn__ = (float *) ps_B0I_trn__;
  /* %%%%%%%%; */
  /* These will store output from f_A0 * f_B0 ; */
  /* i.e., : transpose(legendre_evaluate_normalized_lwpm___(:,:,.)) * tmp_a_k_Y_lam___(:,:,.) for fixed nm. ; */
  float *f_C0R_al__=NULL;
  float *f_C0I_al__=NULL;
  unsigned long long int tmp_n_lam = (unsigned long long int)(1+l_max)*(unsigned long long int)n_a_per_a_batch*(unsigned long long int)n_m_max;
  unsigned long long int n_wpam = (unsigned long long int)n_w*(unsigned long long int)n_viewing_polar_a*(unsigned long long int)n_a_per_a_batch*(unsigned long long int)n_m_max;
  float *spherical_harmonic_unphased_wpam_real____=NULL;
  float *spherical_harmonic_unphased_wpam_imag____=NULL;
  spherical_harmonic_unphased_wpam_real____ = (float *) malloc1(n_wpam*sizeof(float));
  spherical_harmonic_unphased_wpam_imag____ = (float *) malloc1(n_wpam*sizeof(float));
  /* %%%%%%%%; */
  /* Note: This will store aligned versions of: ; */
  /* spherical_harmonic_unphased_mawp_real____ and spherical_harmonic_unphased_mawp_imag____. ; */
  /* We form these arrays (below) directly into the aligned memory. ; */
  unsigned long long int n_row_A1 = (unsigned long long int)n_a_per_a_batch,nrow_A1=0;
  unsigned long long int n_row_A1_pack = (unsigned long long int)n_a_per_a_batch*(unsigned long long int)n_w*(unsigned long long int)n_viewing_polar_a;
  unsigned long long int n_col_X1 = n_m_max,ncol_X1=0;
  unsigned long long int n_col_X1_rup = rup(n_col_X1,8ULL);
  unsigned long long int n_col_X1_256 = n_col_X1_rup/8ULL; //%<-- 8 floats per __m256. ;
  __m256 *ps_A1R_trn__ = NULL;
  __m256 *ps_A1I_trn__ = NULL;
  ps_A1R_trn__ = (__m256 *) _mm_malloc(n_col_X1_256*n_row_A1_pack*sizeof(__m256),32);
  ps_A1I_trn__ = (__m256 *) _mm_malloc(n_col_X1_256*n_row_A1_pack*sizeof(__m256),32);
  float *f_A1R_trn__ = NULL; 
  float *f_A1I_trn__ = NULL; 
  f_A1R_trn__ = (float *) ps_A1R_trn__;
  f_A1I_trn__ = (float *) ps_A1I_trn__;
  memset(f_A1R_trn__,0,n_col_X1_rup*n_row_A1_pack*sizeof(float));
  memset(f_A1I_trn__,0,n_col_X1_rup*n_row_A1_pack*sizeof(float));
  unsigned long long int n_mawp = n_wpam;
  float *spherical_harmonic_unphased_mawp_real____=f_A1R_trn__;
  float *spherical_harmonic_unphased_mawp_imag____=f_A1I_trn__;
  if (verbose>1){ printf(" %% spherical_harmonic_unphased_mawp____: (%d,%d,%d,%d)=%lld (%0.2f GB)\n",n_w,n_viewing_polar_a,n_a_per_a_batch,n_m_max,n_mawp,n_mawp*(double)8/(double)1e9);}
  /* %%%%%%%%; */
  /* Note: This will store aligned versions of: ; */
  /* tmp_expi_mw_real__ and tmp_expi_mw_imag__. ; */
  /* We form these arrays (below) directly into the aligned memory. ; */
  unsigned long long int n_row_B1 = (unsigned long long int)1,nrow_B1=0;
  unsigned long long int n_row_B1_pack = (unsigned long long int)n_w;
  __m256 *ps_B1R_trn__ = NULL;
  __m256 *ps_B1I_trn__ = NULL;
  ps_B1R_trn__ = (__m256 *) _mm_malloc(n_col_X1_256*n_row_B1_pack*sizeof(__m256),32);
  ps_B1I_trn__ = (__m256 *) _mm_malloc(n_col_X1_256*n_row_B1_pack*sizeof(__m256),32);
  float *f_B1R_trn__ = NULL; 
  float *f_B1I_trn__ = NULL; 
  f_B1R_trn__ = (float *) ps_B1R_trn__;
  f_B1I_trn__ = (float *) ps_B1I_trn__;
  memset(f_B1R_trn__,0,n_col_X1_rup*n_row_B1_pack*sizeof(float));
  memset(f_B1I_trn__,0,n_col_X1_rup*n_row_B1_pack*sizeof(float));
  unsigned long long int n_mw = (unsigned long long int)n_m_max*(unsigned long long int)n_w;
  float *tmp_expi_mw_real__=f_B1R_trn__;
  float *tmp_expi_mw_imag__=f_B1I_trn__;
  if (verbose>1){ printf(" %% tmp_expi_mw__: (%d,%d)=%lld (%0.2f GB)\n",n_m_max,n_w,n_mw,n_mw*(double)8/(double)1e9);}
  /* %%%%%%%%; */
  /* These will store output from f_A1 * f_B1 ; */
  /* i.e., : transpose(spherical_harmonic_unphased_mawp____(:,:,.,.)) * tmp_expi_mw__(:,.), for fixed nw, np. ; */
  float *f_C1R_al__=NULL;
  float *f_C1I_al__=NULL;
  /* %%%%%%%%; */
  unsigned long long int n_waS = (unsigned long long int)n_w*(unsigned long long int)n_a*(unsigned long long int)n_viewing_all;
  unsigned long long int n_aw = (unsigned long long int)n_a_per_a_batch*(unsigned long long int)n_w;
  float *spherical_harmonic_evaluate_aw_real__=NULL;
  float *spherical_harmonic_evaluate_aw_imag__=NULL;
  spherical_harmonic_evaluate_aw_real__ = (float *) malloc1(n_aw*sizeof(float));
  spherical_harmonic_evaluate_aw_imag__ = (float *) malloc1(n_aw*sizeof(float));
  if (verbose>1){ printf(" %% spherical_harmonic_evaluate_aw__: (%d,%d)=%lld (%0.2f GB)\n",n_a_per_a_batch,n_w,n_aw,n_aw*(double)8/(double)1e9);}
  double *tmp_expi_sub_real_=(double *) malloc1((unsigned long long int)n_w*sizeof(double));;
  double *tmp_expi_sub_imag_=(double *) malloc1((unsigned long long int)n_w*sizeof(double));;
  double complex z_0in,z_1in,z_out;
  float *template_waS_real___=NULL;
  float *template_waS_imag___=NULL;
  if (template_waS_real___p_!=NULL){
    if ((*template_waS_real___p_)==NULL){ (*template_waS_real___p_) = (float *) malloc1(n_waS*sizeof(float)); }
    template_waS_real___ = (*template_waS_real___p_);
    /* if (template_waS_real___p_!=NULL){ } */}
  if (template_waS_imag___p_!=NULL){
    if ((*template_waS_imag___p_)==NULL){ (*template_waS_imag___p_) = (float *) malloc1(n_waS*sizeof(float)); }
    template_waS_imag___ = (*template_waS_imag___p_);
    /* if (template_waS_imag___p_!=NULL){ } */}
  local_tic(0,t_start_,d_start_);
  elrt_0=0;elrt_1=0;elrt_2=0;elrt_3=0;elrt_4=0;elrt_5=0;  
  for (na_batch=0;na_batch<n_a_batch;na_batch++){
    na=na_batch*n_a_per_a_batch; tmp_na=0;
    while ((tmp_na<n_a_per_a_batch) && (na<n_a)){
      index_a_[tmp_na] = na;
      na++; tmp_na++;
      /* while ((tmp_na<n_a_per_a_batch) && (na<n_a)){ } */}
    if ( (na!=n_a) && (tmp_na!=n_a_per_a_batch) ){ printf(" %% Warning, na %d tmp_na %d in pm_template\n",na,tmp_na);}
    tmp_n_a = tmp_na;
    if (verbose>2){ printf(" %% na_batch %d/%d: tmp_n_a %d\n",na_batch,n_a_batch,tmp_n_a);}
    if (tmp_n_a>0){
      n_row_B0 = tmp_n_a;
      memset(f_B0R_trn__,0,n_col_X0_rup*n_row_B0*sizeof(float));
      memset(f_B0I_trn__,0,n_col_X0_rup*n_row_B0*sizeof(float));
      local_tic(1,t_start_,d_start_);
      for (tmp_na=0;tmp_na<tmp_n_a;tmp_na++){
	na = index_a_[tmp_na];
	if (verbose>3){ printf(" %% na %d\n",na);}
	for (nm=0;nm<n_m_max;nm++){
	  if (verbose>3){ printf(" %% nm %d\n",nm);}
	  ulli_0in = (unsigned long long int)(na + nm*n_a)*(unsigned long long int)(1+l_max);
	  ulli_out = (unsigned long long int)(tmp_na + nm*tmp_n_a)*(unsigned long long int)n_col_X0_rup;
	  if (verbose>3){ printf(" %% ulli_0in %lld ulli_out %lld\n",ulli_0in,ulli_out);}
	  for (ncol_X0=0;ncol_X0<n_col_X0;ncol_X0++){
	    f_B0R_trn__[ulli_out] = (float)(a_k_Y_lam_real___[ulli_0in]);
	    f_B0I_trn__[ulli_out] = (float)(a_k_Y_lam_imag___[ulli_0in]);
	    ulli_out++; ulli_0in++;
	    /* for (ncol_X0=0;ncol_X0<n_col_X0;ncol_X0++){ } */}
	  /* for (nm=0;nm<n_m_max;nm++){ } */}
	/* for (tmp_na=0;tmp_na<tmp_n_a;tmp_na++){ } */}
      local_toc(1,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,0," %% f_B0_trn__: ");
      elrt_0 += elrt_[1];
      /* %%%%%%%%; */
      memset(spherical_harmonic_unphased_wpam_real____,0,n_wpam*sizeof(float));
      memset(spherical_harmonic_unphased_wpam_imag____,0,n_wpam*sizeof(float));
      local_tic(1,t_start_,d_start_);
      for (nm=0;nm<n_m_max;nm++){
	f_C0R_al__ = spherical_harmonic_unphased_wpam_real____ + n_row_A0*n_row_B0*(unsigned long long int)nm;
	f_C0I_al__ = spherical_harmonic_unphased_wpam_imag____ + n_row_A0*n_row_B0*(unsigned long long int)nm;
	rcp_segregated_to_segregated_mult_immintrin_load1_fma
	  (
	   n_row_A0
	   ,n_col_X0
	   ,f_A0R_trn__ + n_row_A0*n_col_X0_rup*(unsigned long long int)nm
	   ,NULL
	   ,n_row_B0
	   ,f_B0R_trn__ + n_row_B0*n_col_X0_rup*(unsigned long long int)nm
	   ,f_B0I_trn__ + n_row_B0*n_col_X0_rup*(unsigned long long int)nm
	   ,&f_C0R_al__
	   ,&f_C0I_al__
	   );
	/* for (nm=0;nm<n_m_max;nm++){ } */}
      local_toc(1,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,0," %% f_B0_trn__: ");
      elrt_1 += elrt_[1];
      /* %%%%%%%%; */
      memset(f_A1R_trn__,0,n_col_X1_rup*n_row_A1_pack*sizeof(float)); //<-- spherical_harmonic_unphased_mawp_real____ ;
      memset(f_A1I_trn__,0,n_col_X1_rup*n_row_A1_pack*sizeof(float)); //<-- spherical_harmonic_unphased_mawp_imag____ ;
      tab_0in=0; tab_out=0;
      local_tic(1,t_start_,d_start_);
      for (nm=0;nm<n_m_max;nm++){
	for (tmp_na=0;tmp_na<tmp_n_a;tmp_na++){
	  for (nviewing_polar_a=0;nviewing_polar_a<n_viewing_polar_a;nviewing_polar_a++){
	    for (nw=0;nw<n_w;nw++){
	      /* tab_0in =  */
	      /* 	(unsigned long long int)nw + */
	      /* 	((unsigned long long int)nviewing_polar_a +  */
	      /* 	 ((unsigned long long int)tmp_na +  */
	      /* 	  ((unsigned long long int)nm */
	      /* 	   )*(unsigned long long int)tmp_n_a */
	      /* 	  )*(unsigned long long int)n_viewing_polar_a */
	      /* 	 )*(unsigned long long int)n_w; */
	      tab_out = 
		(unsigned long long int)nm +
		((unsigned long long int)tmp_na + 
		 ((unsigned long long int)nw +
		  ((unsigned long long int)nviewing_polar_a
		   )*(unsigned long long int)n_w
		  )*(unsigned long long int)tmp_n_a
		 )*(unsigned long long int)n_col_X1_rup;
	      spherical_harmonic_unphased_mawp_real____[tab_out] = spherical_harmonic_unphased_wpam_real____[tab_0in];
	      spherical_harmonic_unphased_mawp_imag____[tab_out] = spherical_harmonic_unphased_wpam_imag____[tab_0in];
	      tab_0in++;
	      /* for (nw=0;nw<n_w;nw++){ } */}
	    /* for (nviewing_polar_a=0;nviewing_polar_a<n_viewing_polar_a;nviewing_polar_a++){ } */}
	  /* for (tmp_na=0;tmp_na<tmp_n_a;tmp_na++){ } */}
	/* for (nm=0;nm<n_m_max;nm++){ } */}
      local_toc(1,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,0," %% spherical_harmonic_unphased_mawp_imag____: ");
      elrt_2 += elrt_[1];
      if (verbose>2){ printf(" %% \t \n");}
      if (verbose>2){ farray_aligned_printf_margin(spherical_harmonic_unphased_mawp_real____,n_m_max,tmp_n_a*n_w*n_viewing_polar_a," %% spherical_harmonic_unphased_mawp_real____ (aligned): ");}
      if (verbose>2){ printf(" %% \t \n");}
      if (verbose>2){ farray_aligned_printf_margin(spherical_harmonic_unphased_mawp_imag____,n_m_max,tmp_n_a*n_w*n_viewing_polar_a," %% spherical_harmonic_unphased_mawp_imag____ (aligned): ");}
      /* %%%%%%%%; */
      memset(spherical_harmonic_evaluate_aw_real__,0,n_aw*sizeof(float));
      memset(spherical_harmonic_evaluate_aw_imag__,0,n_aw*sizeof(float));
      memset(f_B1R_trn__,0,n_col_X1_rup*n_row_B1_pack*sizeof(float)); //<-- tmp_expi_mw_real__. ;
      memset(f_B1I_trn__,0,n_col_X1_rup*n_row_B1_pack*sizeof(float)); //<-- tmp_expi_mw_imag__. ;
      nviewing_all=0;
      for (nviewing_polar_a=0;nviewing_polar_a<n_viewing_polar_a;nviewing_polar_a++){
	n_viewing_azimu_b = n_viewing_azimu_b_[nviewing_polar_a];
	for (nviewing_azimu_b=0;nviewing_azimu_b<n_viewing_azimu_b;nviewing_azimu_b++){
	  if (verbose>3){ printf(" %% nviewing_polar_a %d/%d <-- n_viewing_azimu_b %d/%d\n",nviewing_polar_a,n_viewing_polar_a,nviewing_azimu_b,n_viewing_azimu_b);}
	  memcpy(tmp_expi_sub_real_,expi_template_azimu_b_real__ + nviewing_all*n_w,n_w*sizeof(double));
	  memcpy(tmp_expi_sub_imag_,expi_template_azimu_b_imag__ + nviewing_all*n_w,n_w*sizeof(double));
	  local_tic(1,t_start_,d_start_);
	  nm=0;
	  for (nw=0;nw<n_w;nw++){
	    z_0in = CMPLX( tmp_expi_sub_real_[nw] , tmp_expi_sub_imag_[nw] );
	    z_1in = CMPLX( (double) -l_max , 0.0 );
	    z_out = cpow(z_0in,z_1in);
	    ulli_out = nm + nw*n_col_X1_rup;
	    tmp_expi_mw_real__[ulli_out] = creal(z_out);
	    tmp_expi_mw_imag__[ulli_out] = cimag(z_out);
	    /* for (nw=0;nw<n_w;nw++){ } */}
	  for (nm=1;nm<n_m_max;nm++){
	    for (nw=0;nw<n_w;nw++){
	      ulli_0in = (nm-1) + nw*n_col_X1_rup;
	      z_0in = CMPLX( tmp_expi_mw_real__[ulli_0in] , tmp_expi_mw_imag__[ulli_0in] );
	      z_1in = CMPLX( tmp_expi_sub_real_[nw] , tmp_expi_sub_imag_[nw] );
	      z_out = z_0in * z_1in;
	      ulli_out = nm + nw*n_col_X1_rup;
	      tmp_expi_mw_real__[ulli_out] = creal(z_out);
	      tmp_expi_mw_imag__[ulli_out] = cimag(z_out);
	      /* for (nw=0;nw<n_w;nw++){ } */}
	    /* for (nm=1;nm<n_m_max;nm++){ } */}
	  local_toc(1,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,0," %% tmp_expi_mw__: ");
	  elrt_3 += elrt_[1];
	  local_tic(1,t_start_,d_start_);
	  for (nw=0;nw<n_w;nw++){
	    f_C1R_al__ = spherical_harmonic_evaluate_aw_real__ + (unsigned long long int)nw*(unsigned long long int)tmp_n_a;
	    f_C1I_al__ = spherical_harmonic_evaluate_aw_imag__ + (unsigned long long int)nw*(unsigned long long int)tmp_n_a;
	    nhp_segregated_to_segregated_mult_immintrin_load1_fma
	      (
	        tmp_n_a
	       ,n_col_X1
	       ,f_A1R_trn__ + (unsigned long long int)tmp_n_a*n_col_X1_rup*((unsigned long long int)nw + (unsigned long long int)nviewing_polar_a*(unsigned long long int)n_w)
	       ,f_A1I_trn__ + (unsigned long long int)tmp_n_a*n_col_X1_rup*((unsigned long long int)nw + (unsigned long long int)nviewing_polar_a*(unsigned long long int)n_w)
	       ,n_row_B1
	       ,f_B1R_trn__ + n_row_B1*n_col_X1_rup*(unsigned long long int)nw
	       ,f_B1I_trn__ + n_row_B1*n_col_X1_rup*(unsigned long long int)nw
	       ,&f_C1R_al__
	       ,&f_C1I_al__
	       );
	    /* for (nw=0;nw<n_w;nw++){ } */}
	  local_toc(1,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,0," %% spherical_harmonic_evaluate_aw__: ");
	  elrt_4 += elrt_[1];
	  if ( (verbose>2) && (nviewing_all==0) ){
	    printf(" %% \t \n");
	    farray_aligned_printf_margin(tmp_expi_mw_real__,n_m_max,n_w," %% tmp_expi_mw_real__: ");
	    printf(" %% \t \n");
	    farray_aligned_printf_margin(tmp_expi_mw_imag__,n_m_max,n_w," %% tmp_expi_mw_imag__: ");
	    printf(" %% \t \n");
	    farray_printf_margin(spherical_harmonic_evaluate_aw_real__,tmp_n_a,n_w," %% spherical_harmonic_evaluate_aw_real__: ");
	    printf(" %% \t \n");
	    farray_printf_margin(spherical_harmonic_evaluate_aw_imag__,tmp_n_a,n_w," %% spherical_harmonic_evaluate_aw_imag__: ");
	    /* if ( (verbose>2) && (nviewing_all==0) ){ } */}
	  if ( (template_waS_real___!=NULL) && (template_waS_imag___!=NULL) ){
	    local_tic(1,t_start_,d_start_);
	    for (tmp_na=0;tmp_na<tmp_n_a;tmp_na++){
	      na = index_a_[tmp_na];
	      for (nw=0;nw<n_w;nw++){
	  	tab_0in = (unsigned long long int)tmp_na + (unsigned long long int)nw*(unsigned long long int)tmp_n_a;
	  	tab_out =
	  	  (unsigned long long int)nw +
	  	  ((unsigned long long int)na +
	  	   ((unsigned long long int)nviewing_all
	  	    )*(unsigned long long int)n_a
	  	   )*(unsigned long long int)n_w;
	  	template_waS_real___[tab_out] = spherical_harmonic_evaluate_aw_real__[tab_0in];
	  	template_waS_imag___[tab_out] = spherical_harmonic_evaluate_aw_imag__[tab_0in];
	  	/* for (nw=0;nw<n_w;nw++){ } */}
	      /* for (tmp_na=0;tmp_na<tmp_n_a;tmp_na++){ } */}
	    local_toc(1,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,0," %% template_waS__: ");
	    elrt_5 += elrt_[1];
	    /* if ( (template_waS_real___!=NULL) && (template_waS_imag___!=NULL) ){ } */}
	  nviewing_all++;
	  /* for (nviewing_azimu_b=0;nviewing_azimu_b<n_viewing_azimu_b;nviewing_azimu_b++){ } */}
	/* for (nviewing_polar_a=0;nviewing_polar_a<n_viewing_polar_a;nviewing_polar_a++){ } */}
      /* if (tmp_n_a>0){ } */}
    /* for (na_batch=0;na_batch<n_a_batch;na_batch++){ } */}
  if (verbose>0){
    printf(" %% elrt_0 %0.6fs\n",elrt_0);
    printf(" %% elrt_1 %0.6fs\n",elrt_1);
    printf(" %% elrt_2 %0.6fs\n",elrt_2);
    printf(" %% elrt_3 %0.6fs\n",elrt_3);
    printf(" %% elrt_4 %0.6fs\n",elrt_4);
    printf(" %% elrt_5 %0.6fs\n",elrt_5);
    /* if (verbose>0){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," %% all batches: ");
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(template_waS_real___,"float",n_w*n_a,n_viewing_all," %% template_waS_real___: ");}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(template_waS_imag___,"float",n_w*n_a,n_viewing_all," %% template_waS_imag___: ");}
  /* %%%%%%%%; */
  local_tic(0,t_start_,d_start_);
  if (n_w_p!=NULL){ (*n_w_p) = n_w;}
  if (n_viewing_all_p!=NULL){ (*n_viewing_all_p) = n_viewing_all;}
  if (viewing_azimu_b_all_p_!=NULL){
    tab = (unsigned long long int)n_viewing_all;
    if ((*viewing_azimu_b_all_p_)==NULL){ (*viewing_azimu_b_all_p_) = (double *) malloc1(tab*sizeof(double));}
    memcpy((*viewing_azimu_b_all_p_),viewing_azimu_b_all_,tab*sizeof(double));
    /* if (viewing_azimu_b_all_p_!=NULL){ } */}
  if (viewing_polar_a_all_p_!=NULL){
    tab = (unsigned long long int)n_viewing_all;
    if ((*viewing_polar_a_all_p_)==NULL){ (*viewing_polar_a_all_p_) = (double *) malloc1(tab*sizeof(double));}
    memcpy((*viewing_polar_a_all_p_),viewing_polar_a_all_,tab*sizeof(double));
    /* if (viewing_polar_a_all_p_!=NULL){ } */}
  if (viewing_weight_all_p_!=NULL){
    tab = (unsigned long long int)n_viewing_all;
    if ((*viewing_weight_all_p_)==NULL){ (*viewing_weight_all_p_) = (double *) malloc1(tab*sizeof(double));}
    memcpy((*viewing_weight_all_p_),viewing_weight_all_,tab*sizeof(double));
    /* if (viewing_weight_all_p_!=NULL){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," %% copy output: ");
  _mm_free(ps_A1R_trn__); ps_A1R_trn__ = NULL;
  _mm_free(ps_A1R_trn__); ps_A1R_trn__ = NULL;
  _mm_free(ps_B1R_trn__); ps_B1R_trn__ = NULL;
  _mm_free(ps_B1I_trn__); ps_B1I_trn__ = NULL;
  free1(&tmp_expi_sub_real_);
  free1(&tmp_expi_sub_imag_);
  free1(&spherical_harmonic_evaluate_aw_real__);
  free1(&spherical_harmonic_evaluate_aw_imag__);
  _mm_free(ps_A0R_trn__); ps_A0R_trn__ = NULL;
  _mm_free(ps_B0R_trn__); ps_B0R_trn__ = NULL;
  _mm_free(ps_B0I_trn__); ps_B0I_trn__ = NULL;
  free1(&spherical_harmonic_unphased_wpam_real____);
  free1(&spherical_harmonic_unphased_wpam_imag____);
  free1(&index_a_);
  free1(&a_k_Y_lam_real___);
  free1(&a_k_Y_lam_imag___);
  free1(&legendre_evaluate_normalized_lwpm___);
  free1(&condense_k_c_2__);
  free1(&template_azimu_b__);
  free1(&expi_template_azimu_b_real__);
  free1(&expi_template_azimu_b_imag__);
  free1(&template_k_c_0__);
  free1(&template_k_c_1__);
  free1(&template_k_c_2__);
  free1(&gamma_z_);
  free1(&cc_);
  free1(&sc_);
  free1(&n_viewing_azimu_b_csum_); //<-- not stored. ;
  free1(&viewing_azimu_b_all_); //<-- stored. ;
  free1(&viewing_polar_a_all_); //<-- stored. ;
  free1(&viewing_weight_all_); //<-- stored. ;
  free1(&viewing_polar_a_); //<-- not stored. ;
  free1(&n_viewing_azimu_b_); //<-- not stored. ;
  free1(&m_max_);
  if (verbose>1){ printf(" %% \t \n");}
  if (verbose){ printf(" %% [finished pm_template]\n");}
}

void pm_template_test()
{
  int verbose=1,verbose_local=2;
  if (verbose){ printf(" %% [entering pm_template_test]\n");}
  int n_k_p_r = 49,nk_p_r=0;
  double k_p_r_max = 1.0;
  double *k_p_r_=NULL;
  k_p_r_ = (double *) malloc1((unsigned long long int)n_k_p_r*sizeof(double));
  for (nk_p_r=0;nk_p_r<n_k_p_r;nk_p_r++){ k_p_r_[nk_p_r]=1.0;}
  double *weight_k_p_r_=NULL;
  weight_k_p_r_ = (double *) malloc1((unsigned long long int)n_k_p_r*sizeof(double));
  for (nk_p_r=0;nk_p_r<n_k_p_r;nk_p_r++){ weight_k_p_r_[nk_p_r]=1.0;}
  int l_max = 80;
  int n_lm = (1+l_max)*(1+l_max);
  int *l_max_=NULL;
  l_max_ = (int *) malloc1((unsigned long long int)n_k_p_r*sizeof(int));
  for (nk_p_r=0;nk_p_r<n_k_p_r;nk_p_r++){ l_max_[nk_p_r]=l_max;}
  int n_lm_sum=n_lm*n_k_p_r,ny=0;
  double *a_k_Y_ya_real__=NULL;
  double *a_k_Y_ya_imag__=NULL;
  a_k_Y_ya_real__ = (double *) malloc1((unsigned long long int)n_lm_sum*sizeof(double));
  a_k_Y_ya_imag__ = (double *) malloc1((unsigned long long int)n_lm_sum*sizeof(double));
  for (ny=0;ny<n_lm_sum;ny++){
    a_k_Y_ya_real__[ny] = (double) ( (ny%89) - 44 ) / (double)89;
    a_k_Y_ya_imag__[ny] = (double) ( (ny%97) - 48 ) / (double)97;
    /* for (ny=0;ny<n_lm_sum;ny++){ } */}
  double viewing_k_eq_d = 1.0/(4.0*PI);
  double template_k_eq_d = -1.0;
  int n_w_max = 98;
  float *template_waS_real___=NULL;
  float *template_waS_imag___=NULL;
  int n_w;
  int n_viewing_all;
  double *viewing_azimu_b_all_=NULL;
  double *viewing_polar_a_all_=NULL;
  double *viewing_weight_all_=NULL;
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ printf(" %% l_max %d n_k_p_r %d viewing_k_eq_d %0.6f template_k_eq_d %0.6f n_w_max %d\n",l_max,n_k_p_r,viewing_k_eq_d,template_k_eq_d,n_w_max);}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(a_k_Y_ya_real__,"double",n_lm,n_k_p_r," %% a_k_Y_ya_real__: ");}
  if (verbose>2){ printf(" %% \t \n");}
  if (verbose>2){ array_printf_margin(a_k_Y_ya_imag__,"double",n_lm,n_k_p_r," %% a_k_Y_ya_imag__: ");}
  GLOBAL_tic(0);
  pm_template(
	      verbose_local
	      ,l_max
	      ,n_k_p_r
	      ,a_k_Y_ya_real__
	      ,a_k_Y_ya_imag__
	      ,viewing_k_eq_d
	      ,template_k_eq_d
	      ,n_w_max
	      ,&template_waS_real___
	      ,&template_waS_imag___
	      ,&n_w
	      ,&n_viewing_all
	      ,&viewing_azimu_b_all_
	      ,&viewing_polar_a_all_
	      ,&viewing_weight_all_
	      );
  if (verbose>0){ printf(" %% \t \n"); printf(" %% n_viewing_all %d\n",n_viewing_all);}
  if (verbose>0){ printf(" %% \t \n");}
  if (verbose>0){ array_printf_margin(template_waS_real___,"float",n_w*n_k_p_r,n_viewing_all," %% template_waS_real___: ");}
  if (verbose>0){ printf(" %% \t \n");}
  if (verbose>0){ array_printf_margin(template_waS_imag___,"float",n_w*n_k_p_r,n_viewing_all," %% template_waS_imag___: ");}
  free1(&template_waS_real___);
  free1(&template_waS_imag___);
  free1(&viewing_azimu_b_all_);
  free1(&viewing_polar_a_all_);
  free1(&viewing_weight_all_);
  GLOBAL_toc(0,1," %% pm_template_test: ");
  free1(&k_p_r_);
  free1(&weight_k_p_r_);
  free1(&l_max_);
  free1(&a_k_Y_ya_real__);
  free1(&a_k_Y_ya_imag__);
  if (verbose){ printf(" %% [finished pm_template_test]\n");}
}
