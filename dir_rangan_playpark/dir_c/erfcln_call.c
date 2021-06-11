#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

float erfcln_single_f(float f_0in)
{
 /* approximates log(erfc(input)). ; */
 /* This is quite accurate (errors < 1e-6) ; */
 /* when the input is > 3 or so. ; */
  float f_out=0;
  float x1=0,x2=0,pi=PI;
  if (f_0in<=10){ f_out = logf(erfcf(f_0in));}
  if (f_0in> 10){
    x1 = f_0in;
    x2 = x1*x1;
    f_out = -x2 + logf(x1) - 0.5*logf(pi)
      - logf( x2 + 0.5
	     / (1.0 + 1.0
		 / (x2 + 1.5
		     / ( 1.0 + 2.0
			  / ( x2 + 2.5
			      ) ) ) ) ) ;
    /* if (f_0in> 10){ } */}
  return f_out;
}

void erfcln_f(int n_r,float *f_0in_,float **f_out_p_)
{
  int nr=0;
  float *f_out_=NULL;
  f_out_=NULL;
  if (f_out_p_!=NULL){
    if ( (*f_out_p_)==NULL ){ (*f_out_p_) = (float *) malloc1(n_r*sizeof(float));}
    f_out_ = *f_out_p_;
    /* if (f_out_p_!=NULL){ } */}
  if (f_out_!=NULL){
    for (nr=0;nr<n_r;nr++){ f_out_[nr] = erfcln_single_f(f_0in_[nr]);}
    /* if (f_out_!=NULL){ } */}
}

void erfcln_f_test()
{
  int n_r=7;
  float f_0in_[7] = { -3 , -1 , 3 , 9 , 10 , 11 , 50 };
  float *f_out_=NULL;
  float f_out_ans_[7] = { 0.6931361352504468 , 0.6112323176780705 , -10.7203630419811127 , -83.7756698795290902 , -102.8798890248540658 , -123.9743506042298691 , -2504.4845878484511559 };
  GLOBAL_tic(0);
  erfcln_f(n_r,f_0in_,&f_out_);
  array_printf(f_out_,"float",1,n_r," % f_out_: ");
  printf(" %% f_out_ans_ vs f_out_: relative error %0.16f\n",ffnormn(n_r,f_out_ans_,f_out_));
  free1(&f_out_);
  GLOBAL_toc(0,1," % erfcln_f: ");
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

double erfcln_single_d(double d_0in)
{
 /* approximates log(erfc(input)). ; */
 /* This is quite accurate (errors < 1e-6) ; */
 /* when the input is > 3 or so. ; */
  double d_out=0;
  double x1=0,x2=0,pi=PI;
  if (d_0in<=10){ d_out = log(erfc(d_0in));}
  if (d_0in> 10){
    x1 = d_0in;
    x2 = x1*x1;
    d_out = -x2 + log(x1) - 0.5*log(pi)
      - log( x2 + 0.5
	     / (1.0 + 1.0
		 / (x2 + 1.5
		     / ( 1.0 + 2.0
			  / ( x2 + 2.5
			      ) ) ) ) ) ;
    /* if (d_0in> 10){ } */}
  return d_out;
}

void erfcln_d(int n_r,double *d_0in_,double **d_out_p_)
{
  int nr=0;
  double *d_out_=NULL;
  d_out_=NULL;
  if (d_out_p_!=NULL){
    if ( (*d_out_p_)==NULL ){ (*d_out_p_) = (double *) malloc1(n_r*sizeof(double));}
    d_out_ = *d_out_p_;
    /* if (d_out_p_!=NULL){ } */}
  if (d_out_!=NULL){
    for (nr=0;nr<n_r;nr++){ d_out_[nr] = erfcln_single_d(d_0in_[nr]);}
    /* if (d_out_!=NULL){ } */}
}

void erfcln_d_test()
{
  int n_r=7;
  double d_0in_[7] = { -3 , -1 , 3 , 9 , 10 , 11 , 50 };
  double *d_out_=NULL;
  double d_out_ans_[7] = { 0.6931361352504468 , 0.6112323176780705 , -10.7203630419811127 , -83.7756698795290902 , -102.8798890248540658 , -123.9743506042298691 , -2504.4845878484511559 };
  GLOBAL_tic(0);
  erfcln_d(n_r,d_0in_,&d_out_);
  array_printf(d_out_,"double",1,n_r," % d_out_: ");
  printf(" %% d_out_ans_ vs d_out_: relative error %0.16f\n",dfnormn(n_r,d_out_ans_,d_out_));
  free1(&d_out_);
  GLOBAL_toc(0,1," % erfcln_d: ");
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

double z_to_lp_single_d(double d_0in)
{
  double d_out=0;
  d_out = log((double)0.5) + erfcln_single_d(d_0in/sqrt(2));
  return d_out;
}

void z_to_lp_d(int n_r,double *d_0in_,double **d_out_p_)
{
  int nr=0;
  double *d_out_=NULL;
  d_out_=NULL;
  if (d_out_p_!=NULL){
    if ( (*d_out_p_)==NULL ){ (*d_out_p_) = (double *) malloc1(n_r*sizeof(double));}
    d_out_ = *d_out_p_;
    /* if (d_out_p_!=NULL){ } */}
  if (d_out_!=NULL){
    for (nr=0;nr<n_r;nr++){ d_out_[nr] = z_to_lp_single_d(d_0in_[nr]);}
    /* if (d_out_!=NULL){ } */}
}

void z_to_lp_d_test()
{
  int n_r=7;
  double d_0in_[7] = { -3 , -1 , 3 , 9 , 10 , 11 , 50 };
  double *d_out_=NULL;
  double d_out_ans_[7] = { -0.0013508099647482 , -0.1727537790234499 , -6.6077262215103492 , -43.6281491133321140 , -53.2312851505124627 , -63.8249340944237105 , -1254.8313611394194140 };
  GLOBAL_tic(0);
  z_to_lp_d(n_r,d_0in_,&d_out_);
  array_printf(d_out_,"double",1,n_r," % d_out_: ");
  printf(" %% d_out_ans_ vs d_out_: relative error %0.16f\n",dfnormn(n_r,d_out_ans_,d_out_));
  free1(&d_out_);
  GLOBAL_toc(0,1," % z_to_lp_d: ");
}
