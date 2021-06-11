#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

void fftf_1d_bruteforce(int n_row,float complex *c_0in_,float complex *c_out_,int sgn)
{
  int n_col = n_row;
  int nrow=0,ncol=0;
  float twopiovern = 0;
  float complex w=0;
  float complex wj=0;
  float complex wjk=0;
  twopiovern = 2.0*PI/(double)n_row;
  w = +cos(twopiovern) + _Complex_I*sgn*sin(twopiovern);
  wj = +1.0;
  for (nrow=0;nrow<n_row;nrow++){
    wjk = +1.0;
    c_out_[nrow] = 0.0;
    for (ncol=0;ncol<n_col;ncol++){
      c_out_[nrow] += wjk * c_0in_[ncol];
      wjk *= wj;
      /* for (ncol=0;ncol<n_col;ncol++){ } */}
    wj *= w;
    /* for (nrow=0;nrow<n_row;nrow++){ } */}
  if (sgn>0){ for (nrow=0;nrow<n_row;nrow++){ c_out_[nrow] /= n_row; }}
}

void fftf__1d_bruteforce(int n_row,int n_col,float complex *c_0in_,float complex *c_out_,int sgn)
{
  int ncol=0;
  for (ncol=0;ncol<n_col;ncol++){ fftf_1d_bruteforce(n_row,&(c_0in_[ncol*n_row]),&(c_out_[ncol*n_row]),sgn);}
}

void fftw3f__1d_andplan(int n_row,int n_col,float complex *c_0in_,float complex *c_out_,int sgn)
{
  int rank = 1;
  int n = n_row;
  int howmany = n_col;
  int idist = n_row;
  int odist = n_row;
  int istride = 1;
  int ostride = 1;
  int *inembed = NULL, *onembed = NULL;
  fftwf_complex *fftwf_0in_ = NULL;
  fftwf_complex *fftwf_out_ = NULL;
  fftw_plan tmp_fftw_plan;
  unsigned flags=0;
  unsigned long long int ulli=0,ulli_total=0;
  fftwf_0in_ = (fftwf_complex *) fftw_malloc(n_row*n_col*sizeof(fftwf_complex));
  fftwf_out_ = (fftwf_complex *) fftw_malloc(n_row*n_col*sizeof(fftwf_complex));
  tmp_fftw_plan = fftwf_plan_many_dft(rank,&n,howmany,fftwf_0in_,inembed,istride,idist,fftwf_out_,onembed,ostride,odist,sgn,flags);
  ulli_total = (unsigned long long int) n_row * (unsigned long long int) n_col;
  for (ulli=0;ulli<ulli_total;ulli++){ fftwf_0in_[ulli] = c_0in_[ulli];}
  fftw_execute(tmp_fftw_plan);
  for (ulli=0;ulli<ulli_total;ulli++){ c_out_[ulli] = fftwf_out_[ulli];}
  fftw_destroy_plan(tmp_fftw_plan);
  fftw_free(fftwf_0in_);
  fftw_free(fftwf_out_);
}

void fftw3f__1d_test()
{
  int n_row_A = 5;
  int n_col_A = 4;
  int n_row_A_sub = 5;
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
  c_0in_ = (float complex *) malloc1(n_row_A*n_col_A*sizeof(float complex));
  c_tru_ = (float complex *) malloc1(n_row_A*n_col_A*sizeof(float complex));
  c_out_ = (float complex *) malloc1(n_row_A*n_col_A*sizeof(float complex));
  c_0in_sub_ = (float complex *) malloc1(n_row_A_sub*n_col_A_sub*sizeof(float complex));
  c_tru_sub_ = (float complex *) malloc1(n_row_A_sub*n_col_A_sub*sizeof(float complex));
  c_out_sub_ = (float complex *) malloc1(n_row_A_sub*n_col_A_sub*sizeof(float complex));
  for (ulli=0;ulli<ulli_total;ulli++){ c_0in_[ulli] = 0.25*(float complex)ulli + 0.25*_Complex_I * (float complex)(ulli%7);}
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_A=0;ncol_A<n_col_A_sub;ncol_A++){ c_0in_sub_[nrow_A+ncol_A*n_row_A_sub] = c_0in_[nrow_A+ncol_A*n_row_A];}}
  printf(" %% upper corner of c_0in_: \n");
  array_printf(c_0in_sub_,"float complex",n_row_A_sub,n_col_A_sub," % c_0in_sub_: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  fftf__1d_bruteforce(n_row_A,n_col_A,c_0in_,c_tru_,-1);
  GLOBAL_toc(0,1," fftf__1d_bruteforce: ");
   for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_A=0;ncol_A<n_col_A_sub;ncol_A++){ c_tru_sub_[nrow_A+ncol_A*n_row_A_sub] = c_tru_[nrow_A+ncol_A*n_row_A];}}
  printf(" %% upper corner of c_tru_: \n");
  array_printf(c_tru_sub_,"float complex",n_row_A_sub,n_col_A_sub," % c_tru_sub_: ");
  derror = cfnorm(ulli_total,c_tru_,c_tru_);
  printf(" %% derror %0.16f\n",derror);
  /* %%%%%%%% */
  /* %%%%%%%% */
  GLOBAL_tic(0);
  fftw3f__1d_andplan(n_row_A,n_col_A,c_0in_,c_out_,-1);
  GLOBAL_toc(0,1," fftw3f__1d_andplan: ");
   for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_A=0;ncol_A<n_col_A_sub;ncol_A++){ c_out_sub_[nrow_A+ncol_A*n_row_A_sub] = c_out_[nrow_A+ncol_A*n_row_A];}}
  printf(" %% upper corner of c_out_: \n");
  array_printf(c_out_sub_,"float complex",n_row_A_sub,n_col_A_sub," % c_out_sub_: ");
  derror = cfnorm(ulli_total,c_tru_,c_out_);
  printf(" %% derror %0.16f\n",derror);
  free1(&c_0in_);
  free1(&c_tru_);
  free1(&c_out_);
  free1(&c_0in_sub_);
  free1(&c_tru_sub_);
  free1(&c_out_sub_);
}
