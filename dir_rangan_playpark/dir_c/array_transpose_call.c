#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

void dtranspose_bruteforce(int n_row_A,int n_col_A,double *d_0in__,double *d_out__)
{
  /* transposes double array */
  int nrow_A=0,ncol_A=0;
#pragma omp parallel for
  for (ncol_A=0;ncol_A<n_col_A;ncol_A++){
    for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
      d_out__[ncol_A+nrow_A*n_col_A] = d_0in__[nrow_A+ncol_A*n_row_A];
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
    /* for (ncol_A=0;ncol_A<n_col_A;ncol_A++{ } */}
}

inline void dtranspose_block_AtoB(const int n_row_A,const int n_col_A,const double* A_,double* B_,const int block_size) {
  int nrow_A=0,ncol_A=0;
  int nrow_A_set=0,ncol_A_set=0;
  int n_row_B = n_col_A;
  int n_col_B = n_row_A;
  int nrow_B=0,ncol_B=0;
  int nrow_B_set=0,ncol_B_set=0;
  int n_row_A_sub=0,n_col_A_sub=0;
  int nrow_A_sub=0,ncol_A_sub=0;
  int n_row_B_sub=0,n_col_B_sub=0;
  int nrow_B_sub=0,ncol_B_sub=0;
#pragma omp parallel default(shared) private(ncol_A_set,nrow_A_set,n_row_A_sub,n_col_A_sub,ncol_A_sub,nrow_A_sub,nrow_A,ncol_A,nrow_B,ncol_B)
  { /* begin omp parallel */
    ncol_A_set = 0; nrow_A_set=0; n_row_A_sub=0; n_col_A_sub=0; ncol_A_sub=0; nrow_A_sub=0; nrow_A=0; ncol_A=0; nrow_B=0; ncol_B=0;
#pragma omp for schedule(dynamic)
    for (ncol_A_set=0;ncol_A_set<n_col_A;ncol_A_set+=block_size){
      n_col_A_sub = minimum(block_size,n_col_A-ncol_A_set);
      for (nrow_A_set=0;nrow_A_set<n_row_A;nrow_A_set+=block_size){
	n_row_A_sub = minimum(block_size,n_row_A-nrow_A_set);
	for (ncol_A_sub=0;ncol_A_sub<n_col_A_sub;ncol_A_sub++){
	  ncol_A = ncol_A_set+ncol_A_sub;
	  for (nrow_A_sub=0;nrow_A_sub<n_row_A_sub;nrow_A_sub++){
	    nrow_A = nrow_A_set+nrow_A_sub;
	    nrow_B = ncol_A;
	    ncol_B = nrow_A;
	    B_[nrow_B + ncol_B*n_row_B] = A_[nrow_A + ncol_A*n_row_A];
	    /* for (nrow_A_sub=0;nrow_A_sub<n_row_A_sub;nrow_A_sub++){ } */}
	  /* for (ncol_A_sub=0;ncol_A_sub<n_col_A_sub;ncol_A_sub++){ } */}
	/* for (nrow_A_set=0;nrow_A_set<n_row_A;nrow_A_set+=block_size){ } */}
      /* for (ncol_A_set=0;ncol_A_set<n_col_A;ncol_A_set+=block_size){ } */}
    /* end omp parallel */}
}

inline void dtranspose_block_BtoA(const int n_row_A,const int n_col_A,const double* A_,double* B_,const int block_size) {
  int nrow_A=0,ncol_A=0;
  int nrow_A_set=0,ncol_A_set=0;
  int n_row_B = n_col_A;
  int n_col_B = n_row_A;
  int nrow_B=0,ncol_B=0;
  int nrow_B_set=0,ncol_B_set=0;
  int n_row_A_sub=0,n_col_A_sub=0;
  int nrow_A_sub=0,ncol_A_sub=0;
  int n_row_B_sub=0,n_col_B_sub=0;
  int nrow_B_sub=0,ncol_B_sub=0;
#pragma omp parallel default(shared) private(ncol_A_set,nrow_A_set,n_row_A_sub,n_col_A_sub,ncol_A_sub,nrow_A_sub,nrow_A,ncol_A,nrow_B,ncol_B)
  { /* begin omp parallel */
    ncol_A_set = 0; nrow_A_set=0; n_row_A_sub=0; n_col_A_sub=0; ncol_A_sub=0; nrow_A_sub=0; nrow_A=0; ncol_A=0; nrow_B=0; ncol_B=0;
#pragma omp for schedule(dynamic)
    for (nrow_A_set=0;nrow_A_set<n_row_A;nrow_A_set+=block_size){
      n_row_A_sub = minimum(block_size,n_row_A-nrow_A_set);
      for (ncol_A_set=0;ncol_A_set<n_col_A;ncol_A_set+=block_size){
	n_col_A_sub = minimum(block_size,n_col_A-ncol_A_set);
	for (nrow_A_sub=0;nrow_A_sub<n_row_A_sub;nrow_A_sub++){
	  nrow_A = nrow_A_set+nrow_A_sub;
	  for (ncol_A_sub=0;ncol_A_sub<n_col_A_sub;ncol_A_sub++){
	    ncol_A = ncol_A_set+ncol_A_sub;
	    nrow_B = ncol_A;
	    ncol_B = nrow_A;
	    B_[nrow_B + ncol_B*n_row_B] = A_[nrow_A + ncol_A*n_row_A];
	    /* for (ncol_A_sub=0;ncol_A_sub<n_col_A_sub;ncol_A_sub++){ } */}
	  /* for (nrow_A_sub=0;nrow_A_sub<n_row_A_sub;nrow_A_sub++){ } */}
	/* for (ncol_A_set=0;ncol_A_set<n_col_A;ncol_A_set+=block_size){ } */}
      /* for (nrow_A_set=0;nrow_A_set<n_row_A;nrow_A_set+=block_size){ } */}
    /* end omp parallel */}
}

inline void dtranspose(const int n_row_A,const int n_col_A,const double* A_,const double* B_)
{
  dtranspose_block_AtoB(n_row_A,n_col_A,A_,B_,32);
}

void dtranspose_test()
{
  int n_row_A = 1024*10 + 723;
  int n_col_A = 1024*6 + 817;
  /* int n_row_A = 4; */
  /* int n_col_A = 4; */
  int n_row_A_sub = 4;
  int n_col_A_sub = 4;
  int nrow_A=0,ncol_A=0;
  unsigned long long int ulli_total = (unsigned long long int) n_row_A * (unsigned long long int) n_col_A ;
  unsigned long long int ulli=0;
  double *d_0in_ = NULL;
  double *d_tru_ = NULL;
  double *d_out_ = NULL;
  double *d_0in_sub_ = NULL;
  double *d_tru_sub_ = NULL;
  double *d_out_sub_ = NULL;
  double derror=0;
  int block_size=0;
  d_0in_ = (double *) malloc1(n_row_A*n_col_A*sizeof(double));
  d_tru_ = (double *) malloc1(n_row_A*n_col_A*sizeof(double));
  d_out_ = (double *) malloc1(n_row_A*n_col_A*sizeof(double));
  d_0in_sub_ = (double *) malloc1(n_row_A_sub*n_col_A_sub*sizeof(double));
  d_tru_sub_ = (double *) malloc1(n_row_A_sub*n_col_A_sub*sizeof(double));
  d_out_sub_ = (double *) malloc1(n_row_A_sub*n_col_A_sub*sizeof(double));
  GLOBAL_tic(0);
  //for (ulli=0;ulli<ulli_total;ulli++){ d_0in_[ulli] = rand01;}
  for (ulli=0;ulli<ulli_total;ulli++){ d_0in_[ulli] = (double)ulli;}
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_A=0;ncol_A<n_col_A_sub;ncol_A++){ d_0in_sub_[nrow_A+ncol_A*n_row_A_sub] = d_0in_[nrow_A+ncol_A*n_row_A];}}
  printf(" %% upper corner of d_0in_: \n");
  array_printf(d_0in_sub_,"double",n_row_A_sub,n_col_A_sub," % d_0in_sub_: ");
  GLOBAL_toc(0,1," initialize: ");
  GLOBAL_tic(0);
  dtranspose_bruteforce(n_row_A,n_col_A,d_0in_,d_tru_);
  GLOBAL_toc(0,1," dtranspose_bruteforce: ");
  /* %%%%%%%% */
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_A=0;ncol_A<n_col_A_sub;ncol_A++){ d_tru_sub_[nrow_A+ncol_A*n_row_A_sub] = d_tru_[nrow_A+ncol_A*n_col_A];}}
  printf(" %% upper corner of d_tru_: \n");
  array_printf(d_tru_sub_,"double",n_row_A_sub,n_col_A_sub," % d_tru_sub_: ");
  derror = dfnorm(ulli_total,d_tru_,d_tru_);
  printf(" %% derror %0.16f\n",derror);
  /* %%%%%%%% */
  for (block_size=1;block_size<=128;block_size*=2){
    GLOBAL_tic(0);
    dtranspose_block_AtoB(n_row_A,n_col_A,d_0in_,d_out_,block_size);
    GLOBAL_toc(0,0," dtranspose_block_AtoB: ");
    derror = dfnorm(ulli_total,d_tru_,d_out_);
    printf(" %% dtranspose_block_AtoB: block_size %d --> elct/elrt %0.3f/%0.3f --> error %0.16f\n",block_size,GLOBAL_elct[0],GLOBAL_elrt[0],derror);
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_A=0;ncol_A<n_col_A_sub;ncol_A++){ d_out_sub_[nrow_A+ncol_A*n_row_A_sub] = d_out_[nrow_A+ncol_A*n_col_A];}}
    //printf(" %% upper corner of d_out_: \n");
    //array_printf(d_out_sub_,"double",n_row_A_sub,n_col_A_sub," % d_out_sub_: ");
    /* for (block_size=1;block_size<128;block_size*=2){ } */}
  /* %%%%%%%% */
  for (block_size=1;block_size<=128;block_size*=2){
    GLOBAL_tic(0);
    dtranspose_block_BtoA(n_row_A,n_col_A,d_0in_,d_out_,block_size);
    GLOBAL_toc(0,0," dtranspose_block_BtoA: ");
    derror = dfnorm(ulli_total,d_tru_,d_out_);
    printf(" %% dtranspose_block_BtoA: block_size %d --> elct/elrt %0.3f/%0.3f --> error %0.16f\n",block_size,GLOBAL_elct[0],GLOBAL_elrt[0],derror);
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_A=0;ncol_A<n_col_A_sub;ncol_A++){ d_out_sub_[nrow_A+ncol_A*n_row_A_sub] = d_out_[nrow_A+ncol_A*n_col_A];}}
    //printf(" %% upper corner of d_out_: \n");
    //array_printf(d_out_sub_,"double",n_row_A_sub,n_col_A_sub," % d_out_sub_: ");
    /* for (block_size=1;block_size<128;block_size*=2){ } */}
  /* %%%%%%%% */
  free1(&d_0in_);
  free1(&d_tru_);
  free1(&d_out_);
  free1(&d_0in_sub_);
  free1(&d_tru_sub_);
  free1(&d_out_sub_);
  //wkspace_printf();
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void cntranspose_bruteforce(int n_row_A,int n_col_A,float complex *c_0in__,float complex *c_out__)
{
  /* transposes float complex array */
  int nrow_A=0,ncol_A=0;
#pragma omp parallel for
  for (ncol_A=0;ncol_A<n_col_A;ncol_A++){
    for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
      c_out__[ncol_A+nrow_A*n_col_A] = c_0in__[nrow_A+ncol_A*n_row_A];
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
    /* for (ncol_A=0;ncol_A<n_col_A;ncol_A++{ } */}
}

inline void cntranspose_block_AtoB(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_,const int block_size) {
  int nrow_A=0,ncol_A=0;
  int nrow_A_set=0,ncol_A_set=0;
  int n_row_B = n_col_A;
  int n_col_B = n_row_A;
  int nrow_B=0,ncol_B=0;
  int nrow_B_set=0,ncol_B_set=0;
  int n_row_A_sub=0,n_col_A_sub=0;
  int nrow_A_sub=0,ncol_A_sub=0;
  int n_row_B_sub=0,n_col_B_sub=0;
  int nrow_B_sub=0,ncol_B_sub=0;
#pragma omp parallel default(shared) private(ncol_A_set,nrow_A_set,n_row_A_sub,n_col_A_sub,ncol_A_sub,nrow_A_sub,nrow_A,ncol_A,nrow_B,ncol_B)
  { /* begin omp parallel */
    ncol_A_set = 0; nrow_A_set=0; n_row_A_sub=0; n_col_A_sub=0; ncol_A_sub=0; nrow_A_sub=0; nrow_A=0; ncol_A=0; nrow_B=0; ncol_B=0;
#pragma omp for schedule(dynamic)
    for (ncol_A_set=0;ncol_A_set<n_col_A;ncol_A_set+=block_size){
      n_col_A_sub = minimum(block_size,n_col_A-ncol_A_set);
      for (nrow_A_set=0;nrow_A_set<n_row_A;nrow_A_set+=block_size){
	n_row_A_sub = minimum(block_size,n_row_A-nrow_A_set);
	for (ncol_A_sub=0;ncol_A_sub<n_col_A_sub;ncol_A_sub++){
	  ncol_A = ncol_A_set+ncol_A_sub;
	  for (nrow_A_sub=0;nrow_A_sub<n_row_A_sub;nrow_A_sub++){
	    nrow_A = nrow_A_set+nrow_A_sub;
	    nrow_B = ncol_A;
	    ncol_B = nrow_A;
	    B_[nrow_B + ncol_B*n_row_B] = A_[nrow_A + ncol_A*n_row_A];
	    /* for (nrow_A_sub=0;nrow_A_sub<n_row_A_sub;nrow_A_sub++){ } */}
	  /* for (ncol_A_sub=0;ncol_A_sub<n_col_A_sub;ncol_A_sub++){ } */}
	/* for (nrow_A_set=0;nrow_A_set<n_row_A;nrow_A_set+=block_size){ } */}
      /* for (ncol_A_set=0;ncol_A_set<n_col_A;ncol_A_set+=block_size){ } */}
    /* end omp parallel */}
}

inline void cntranspose_block_BtoA(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_,const int block_size) {
  int nrow_A=0,ncol_A=0;
  int nrow_A_set=0,ncol_A_set=0;
  int n_row_B = n_col_A;
  int n_col_B = n_row_A;
  int nrow_B=0,ncol_B=0;
  int nrow_B_set=0,ncol_B_set=0;
  int n_row_A_sub=0,n_col_A_sub=0;
  int nrow_A_sub=0,ncol_A_sub=0;
  int n_row_B_sub=0,n_col_B_sub=0;
  int nrow_B_sub=0,ncol_B_sub=0;
#pragma omp parallel default(shared) private(ncol_A_set,nrow_A_set,n_row_A_sub,n_col_A_sub,ncol_A_sub,nrow_A_sub,nrow_A,ncol_A,nrow_B,ncol_B)
  { /* begin omp parallel */
    ncol_A_set = 0; nrow_A_set=0; n_row_A_sub=0; n_col_A_sub=0; ncol_A_sub=0; nrow_A_sub=0; nrow_A=0; ncol_A=0; nrow_B=0; ncol_B=0;
#pragma omp for schedule(dynamic)
    for (nrow_A_set=0;nrow_A_set<n_row_A;nrow_A_set+=block_size){
      n_row_A_sub = minimum(block_size,n_row_A-nrow_A_set);
      for (ncol_A_set=0;ncol_A_set<n_col_A;ncol_A_set+=block_size){
	n_col_A_sub = minimum(block_size,n_col_A-ncol_A_set);
	for (nrow_A_sub=0;nrow_A_sub<n_row_A_sub;nrow_A_sub++){
	  nrow_A = nrow_A_set+nrow_A_sub;
	  for (ncol_A_sub=0;ncol_A_sub<n_col_A_sub;ncol_A_sub++){
	    ncol_A = ncol_A_set+ncol_A_sub;
	    nrow_B = ncol_A;
	    ncol_B = nrow_A;
	    B_[nrow_B + ncol_B*n_row_B] = A_[nrow_A + ncol_A*n_row_A];
	    /* for (ncol_A_sub=0;ncol_A_sub<n_col_A_sub;ncol_A_sub++){ } */}
	  /* for (nrow_A_sub=0;nrow_A_sub<n_row_A_sub;nrow_A_sub++){ } */}
	/* for (ncol_A_set=0;ncol_A_set<n_col_A;ncol_A_set+=block_size){ } */}
      /* for (nrow_A_set=0;nrow_A_set<n_row_A;nrow_A_set+=block_size){ } */}
    /* end omp parallel */}
}

inline void cntranspose(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_)
{
  cntranspose_block_AtoB(n_row_A,n_col_A,A_,B_,32);
}

void cntranspose_test()
{
  int n_row_A = 1024*10 + 723;
  int n_col_A = 1024*6 + 817;
  /* int n_row_A = 4; */
  /* int n_col_A = 4; */
  int n_row_A_sub = 4;
  int n_col_A_sub = 4;
  int nrow_A=0,ncol_A=0;
  unsigned long long int ulli_total = (unsigned long long int) n_row_A * (unsigned long long int) n_col_A ;
  unsigned long long int ulli=0;
  float complex *c_0in_ = NULL;
  float complex *c_tru_ = NULL;
  float complex *c_out_ = NULL;
  float complex *c_0in_sub_ = NULL;
  float complex *c_tru_sub_ = NULL;
  float complex *c_out_sub_ = NULL;
  double derror=0;
  int block_size=0;
  c_0in_ = (float complex *) malloc1(n_row_A*n_col_A*sizeof(float complex));
  c_tru_ = (float complex *) malloc1(n_row_A*n_col_A*sizeof(float complex));
  c_out_ = (float complex *) malloc1(n_row_A*n_col_A*sizeof(float complex));
  c_0in_sub_ = (float complex *) malloc1(n_row_A_sub*n_col_A_sub*sizeof(float complex));
  c_tru_sub_ = (float complex *) malloc1(n_row_A_sub*n_col_A_sub*sizeof(float complex));
  c_out_sub_ = (float complex *) malloc1(n_row_A_sub*n_col_A_sub*sizeof(float complex));
  GLOBAL_tic(0);
  //for (ulli=0;ulli<ulli_total;ulli++){ c_0in_[ulli] = rand01;}
  for (ulli=0;ulli<ulli_total;ulli++){ c_0in_[ulli] = (float complex)ulli + _Complex_I * (float complex)(ulli%7);}
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_A=0;ncol_A<n_col_A_sub;ncol_A++){ c_0in_sub_[nrow_A+ncol_A*n_row_A_sub] = c_0in_[nrow_A+ncol_A*n_row_A];}}
  printf(" %% upper corner of c_0in_: \n");
  array_printf(c_0in_sub_,"float complex",n_row_A_sub,n_col_A_sub," % c_0in_sub_: ");
  GLOBAL_toc(0,1," initialize: ");
  GLOBAL_tic(0);
  cntranspose_bruteforce(n_row_A,n_col_A,c_0in_,c_tru_);
  GLOBAL_toc(0,1," cntranspose_bruteforce: ");
  /* %%%%%%%% */
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_A=0;ncol_A<n_col_A_sub;ncol_A++){ c_tru_sub_[nrow_A+ncol_A*n_row_A_sub] = c_tru_[nrow_A+ncol_A*n_col_A];}}
  printf(" %% upper corner of c_tru_: \n");
  array_printf(c_tru_sub_,"float complex",n_row_A_sub,n_col_A_sub," % c_tru_sub_: ");
  derror = cfnorm(ulli_total,c_tru_,c_tru_);
  printf(" %% derror %0.16f\n",derror);
  /* %%%%%%%% */
  for (block_size=1;block_size<=128;block_size*=2){
    GLOBAL_tic(0);
    cntranspose_block_AtoB(n_row_A,n_col_A,c_0in_,c_out_,block_size);
    GLOBAL_toc(0,0," cntranspose_block_AtoB: ");
    derror = cfnorm(ulli_total,c_tru_,c_out_);
    printf(" %% cntranspose_block_AtoB: block_size %d --> elct/elrt %0.3f/%0.3f --> error %0.16f\n",block_size,GLOBAL_elct[0],GLOBAL_elrt[0],derror);
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_A=0;ncol_A<n_col_A_sub;ncol_A++){ c_out_sub_[nrow_A+ncol_A*n_row_A_sub] = c_out_[nrow_A+ncol_A*n_col_A];}}
    //printf(" %% upper corner of c_out_: \n");
    //array_printf(c_out_sub_,"float complex",n_row_A_sub,n_col_A_sub," % c_out_sub_: ");
    /* for (block_size=1;block_size<128;block_size*=2){ } */}
  /* %%%%%%%% */
  for (block_size=1;block_size<=128;block_size*=2){
    GLOBAL_tic(0);
    cntranspose_block_BtoA(n_row_A,n_col_A,c_0in_,c_out_,block_size);
    GLOBAL_toc(0,0," cntranspose_block_BtoA: ");
    derror = cfnorm(ulli_total,c_tru_,c_out_);
    printf(" %% cntranspose_block_BtoA: block_size %d --> elct/elrt %0.3f/%0.3f --> error %0.16f\n",block_size,GLOBAL_elct[0],GLOBAL_elrt[0],derror);
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_A=0;ncol_A<n_col_A_sub;ncol_A++){ c_out_sub_[nrow_A+ncol_A*n_row_A_sub] = c_out_[nrow_A+ncol_A*n_col_A];}}
    //printf(" %% upper corner of c_out_: \n");
    //array_printf(c_out_sub_,"float complex",n_row_A_sub,n_col_A_sub," % c_out_sub_: ");
    /* for (block_size=1;block_size<128;block_size*=2){ } */}
  /* %%%%%%%% */
  free1(&c_0in_);
  free1(&c_tru_);
  free1(&c_out_);
  free1(&c_0in_sub_);
  free1(&c_tru_sub_);
  free1(&c_out_sub_);
  //wkspace_printf();
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void cctranspose_bruteforce(int n_row_A,int n_col_A,float complex *c_0in__,float complex *c_out__)
{
  /* transposes float complex array */
  int nrow_A=0,ncol_A=0;
#pragma omp parallel for
  for (ncol_A=0;ncol_A<n_col_A;ncol_A++){
    for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
      c_out__[ncol_A+nrow_A*n_col_A] = conjf(c_0in__[nrow_A+ncol_A*n_row_A]);
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
    /* for (ncol_A=0;ncol_A<n_col_A;ncol_A++{ } */}
}

inline void cctranspose_block_AtoB(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_,const int block_size) {
  int nrow_A=0,ncol_A=0;
  int nrow_A_set=0,ncol_A_set=0;
  int n_row_B = n_col_A;
  int n_col_B = n_row_A;
  int nrow_B=0,ncol_B=0;
  int nrow_B_set=0,ncol_B_set=0;
  int n_row_A_sub=0,n_col_A_sub=0;
  int nrow_A_sub=0,ncol_A_sub=0;
  int n_row_B_sub=0,n_col_B_sub=0;
  int nrow_B_sub=0,ncol_B_sub=0;
#pragma omp parallel default(shared) private(ncol_A_set,nrow_A_set,n_row_A_sub,n_col_A_sub,ncol_A_sub,nrow_A_sub,nrow_A,ncol_A,nrow_B,ncol_B)
  { /* begin omp parallel */
    ncol_A_set = 0; nrow_A_set=0; n_row_A_sub=0; n_col_A_sub=0; ncol_A_sub=0; nrow_A_sub=0; nrow_A=0; ncol_A=0; nrow_B=0; ncol_B=0;
#pragma omp for schedule(dynamic)
    for (ncol_A_set=0;ncol_A_set<n_col_A;ncol_A_set+=block_size){
      n_col_A_sub = minimum(block_size,n_col_A-ncol_A_set);
      for (nrow_A_set=0;nrow_A_set<n_row_A;nrow_A_set+=block_size){
	n_row_A_sub = minimum(block_size,n_row_A-nrow_A_set);
	for (ncol_A_sub=0;ncol_A_sub<n_col_A_sub;ncol_A_sub++){
	  ncol_A = ncol_A_set+ncol_A_sub;
	  for (nrow_A_sub=0;nrow_A_sub<n_row_A_sub;nrow_A_sub++){
	    nrow_A = nrow_A_set+nrow_A_sub;
	    nrow_B = ncol_A;
	    ncol_B = nrow_A;
	    B_[nrow_B + ncol_B*n_row_B] = conjf(A_[nrow_A + ncol_A*n_row_A]);
	    /* for (nrow_A_sub=0;nrow_A_sub<n_row_A_sub;nrow_A_sub++){ } */}
	  /* for (ncol_A_sub=0;ncol_A_sub<n_col_A_sub;ncol_A_sub++){ } */}
	/* for (nrow_A_set=0;nrow_A_set<n_row_A;nrow_A_set+=block_size){ } */}
      /* for (ncol_A_set=0;ncol_A_set<n_col_A;ncol_A_set+=block_size){ } */}
    /* end omp parallel */}
}

inline void cctranspose_block_BtoA(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_,const int block_size) {
  int nrow_A=0,ncol_A=0;
  int nrow_A_set=0,ncol_A_set=0;
  int n_row_B = n_col_A;
  int n_col_B = n_row_A;
  int nrow_B=0,ncol_B=0;
  int nrow_B_set=0,ncol_B_set=0;
  int n_row_A_sub=0,n_col_A_sub=0;
  int nrow_A_sub=0,ncol_A_sub=0;
  int n_row_B_sub=0,n_col_B_sub=0;
  int nrow_B_sub=0,ncol_B_sub=0;
#pragma omp parallel default(shared) private(ncol_A_set,nrow_A_set,n_row_A_sub,n_col_A_sub,ncol_A_sub,nrow_A_sub,nrow_A,ncol_A,nrow_B,ncol_B)
  { /* begin omp parallel */
    ncol_A_set = 0; nrow_A_set=0; n_row_A_sub=0; n_col_A_sub=0; ncol_A_sub=0; nrow_A_sub=0; nrow_A=0; ncol_A=0; nrow_B=0; ncol_B=0;
#pragma omp for schedule(dynamic)
    for (nrow_A_set=0;nrow_A_set<n_row_A;nrow_A_set+=block_size){
      n_row_A_sub = minimum(block_size,n_row_A-nrow_A_set);
      for (ncol_A_set=0;ncol_A_set<n_col_A;ncol_A_set+=block_size){
	n_col_A_sub = minimum(block_size,n_col_A-ncol_A_set);
	for (nrow_A_sub=0;nrow_A_sub<n_row_A_sub;nrow_A_sub++){
	  nrow_A = nrow_A_set+nrow_A_sub;
	  for (ncol_A_sub=0;ncol_A_sub<n_col_A_sub;ncol_A_sub++){
	    ncol_A = ncol_A_set+ncol_A_sub;
	    nrow_B = ncol_A;
	    ncol_B = nrow_A;
	    B_[nrow_B + ncol_B*n_row_B] = conjf(A_[nrow_A + ncol_A*n_row_A]);
	    /* for (ncol_A_sub=0;ncol_A_sub<n_col_A_sub;ncol_A_sub++){ } */}
	  /* for (nrow_A_sub=0;nrow_A_sub<n_row_A_sub;nrow_A_sub++){ } */}
	/* for (ncol_A_set=0;ncol_A_set<n_col_A;ncol_A_set+=block_size){ } */}
      /* for (nrow_A_set=0;nrow_A_set<n_row_A;nrow_A_set+=block_size){ } */}
    /* end omp parallel */}
}

inline void cctranspose(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_)
{
  cctranspose_block_AtoB(n_row_A,n_col_A,A_,B_,32);
}
  

void cctranspose_test()
{
  int n_row_A = 1024*10 + 723;
  int n_col_A = 1024*6 + 817;
  /* int n_row_A = 4; */
  /* int n_col_A = 4; */
  int n_row_A_sub = 4;
  int n_col_A_sub = 4;
  int nrow_A=0,ncol_A=0;
  unsigned long long int ulli_total = (unsigned long long int) n_row_A * (unsigned long long int) n_col_A ;
  unsigned long long int ulli=0;
  float complex *c_0in_ = NULL;
  float complex *c_tru_ = NULL;
  float complex *c_out_ = NULL;
  float complex *c_0in_sub_ = NULL;
  float complex *c_tru_sub_ = NULL;
  float complex *c_out_sub_ = NULL;
  double derror=0;
  int block_size=0;
  c_0in_ = (float complex *) malloc1(n_row_A*n_col_A*sizeof(float complex));
  c_tru_ = (float complex *) malloc1(n_row_A*n_col_A*sizeof(float complex));
  c_out_ = (float complex *) malloc1(n_row_A*n_col_A*sizeof(float complex));
  c_0in_sub_ = (float complex *) malloc1(n_row_A_sub*n_col_A_sub*sizeof(float complex));
  c_tru_sub_ = (float complex *) malloc1(n_row_A_sub*n_col_A_sub*sizeof(float complex));
  c_out_sub_ = (float complex *) malloc1(n_row_A_sub*n_col_A_sub*sizeof(float complex));
  GLOBAL_tic(0);
  //for (ulli=0;ulli<ulli_total;ulli++){ c_0in_[ulli] = rand01;}
  for (ulli=0;ulli<ulli_total;ulli++){ c_0in_[ulli] = (float complex)ulli + _Complex_I * (float complex)(ulli%7);}
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_A=0;ncol_A<n_col_A_sub;ncol_A++){ c_0in_sub_[nrow_A+ncol_A*n_row_A_sub] = c_0in_[nrow_A+ncol_A*n_row_A];}}
  printf(" %% upper corner of c_0in_: \n");
  array_printf(c_0in_sub_,"float complex",n_row_A_sub,n_col_A_sub," % c_0in_sub_: ");
  GLOBAL_toc(0,1," initialize: ");
  GLOBAL_tic(0);
  cctranspose_bruteforce(n_row_A,n_col_A,c_0in_,c_tru_);
  GLOBAL_toc(0,1," cctranspose_bruteforce: ");
  /* %%%%%%%% */
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_A=0;ncol_A<n_col_A_sub;ncol_A++){ c_tru_sub_[nrow_A+ncol_A*n_row_A_sub] = c_tru_[nrow_A+ncol_A*n_col_A];}}
  printf(" %% upper corner of c_tru_: \n");
  array_printf(c_tru_sub_,"float complex",n_row_A_sub,n_col_A_sub," % c_tru_sub_: ");
  derror = cfnorm(ulli_total,c_tru_,c_tru_);
  printf(" %% derror %0.16f\n",derror);
  /* %%%%%%%% */
  for (block_size=1;block_size<=128;block_size*=2){
    GLOBAL_tic(0);
    cctranspose_block_AtoB(n_row_A,n_col_A,c_0in_,c_out_,block_size);
    GLOBAL_toc(0,0," cctranspose_block_AtoB: ");
    derror = cfnorm(ulli_total,c_tru_,c_out_);
    printf(" %% cctranspose_block_AtoB: block_size %d --> elct/elrt %0.3f/%0.3f --> error %0.16f\n",block_size,GLOBAL_elct[0],GLOBAL_elrt[0],derror);
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_A=0;ncol_A<n_col_A_sub;ncol_A++){ c_out_sub_[nrow_A+ncol_A*n_row_A_sub] = c_out_[nrow_A+ncol_A*n_col_A];}}
    //printf(" %% upper corner of c_out_: \n");
    //array_printf(c_out_sub_,"float complex",n_row_A_sub,n_col_A_sub," % c_out_sub_: ");
    /* for (block_size=1;block_size<128;block_size*=2){ } */}
  /* %%%%%%%% */
  for (block_size=1;block_size<=128;block_size*=2){
    GLOBAL_tic(0);
    cctranspose_block_BtoA(n_row_A,n_col_A,c_0in_,c_out_,block_size);
    GLOBAL_toc(0,0," cctranspose_block_BtoA: ");
    derror = cfnorm(ulli_total,c_tru_,c_out_);
    printf(" %% cctranspose_block_BtoA: block_size %d --> elct/elrt %0.3f/%0.3f --> error %0.16f\n",block_size,GLOBAL_elct[0],GLOBAL_elrt[0],derror);
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_A=0;ncol_A<n_col_A_sub;ncol_A++){ c_out_sub_[nrow_A+ncol_A*n_row_A_sub] = c_out_[nrow_A+ncol_A*n_col_A];}}
    //printf(" %% upper corner of c_out_: \n");
    //array_printf(c_out_sub_,"float complex",n_row_A_sub,n_col_A_sub," % c_out_sub_: ");
    /* for (block_size=1;block_size<128;block_size*=2){ } */}
  /* %%%%%%%% */
  free1(&c_0in_);
  free1(&c_tru_);
  free1(&c_out_);
  free1(&c_0in_sub_);
  free1(&c_tru_sub_);
  free1(&c_out_sub_);
  //wkspace_printf();
}
