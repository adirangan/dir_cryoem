/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* The computational routine */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

#include "common.h"
#include "ampmh_omp.h"

int main(int argc, char *argv[])
{
  int verbose = 0;
  int flag_test = 0;
  /* 0in */
  int n_M_per_Mbatch = 1 + ceil(1024 / 35);
  int n_S_per_Sbatch = 24;
  int flag_optimize_over_gamma_z = 0;
  int flag_compute_I_value = 1;
  double tolerance_master = 0.01;
  int FTK_n_svd_l = 10;
  int FTK_n_delta_v = 29;
  double *FTK_svd_U_d_expiw_s_real__ = NULL, *FTK_svd_U_d_expiw_s_imag__ = NULL;
  double *FTK_delta_x_ = NULL;
  double *FTK_delta_y_ = NULL;
  int n_w_max = 98;
  int pm_n_UX_rank = 18;
  int n_S = 3 * n_S_per_Sbatch + 5; // int n_S=993;
  n_S = 993;
  double *CTF_UX_S_k_q_wnS_real__ = NULL, *CTF_UX_S_k_q_wnS_imag__ = NULL;
  double *CTF_UX_S_l2_ = NULL;
  int n_M = 2 * n_M_per_Mbatch + 7; // int n_M=1024;
  n_M = 1024;
  double *svd_VUXM_lwnM_real____ = NULL, *svd_VUXM_lwnM_imag____ = NULL;
  double *UX_M_l2_dM__ = NULL;
  /* out */
  double *X_wSM___ = NULL;
  double *delta_x_wSM___ = NULL;
  double *delta_y_wSM___ = NULL;
  double *gamma_z_wSM___ = NULL;
  double *I_value_wSM___ = NULL;
  /* %%%% */
  int ntick = 0, n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick], l_ssec_[n_tick], l_usec_[n_tick];
  double elct_[n_tick], elrt_[n_tick];
  /* %%%% */
  unsigned long long int na = 0, n_a = 0;
  /* %%%% */
  if (verbose)
  {
    printf(" %% [entering main]\n");
  }
  /* %%%%; */
  n_a = (unsigned long long int)FTK_n_delta_v * (unsigned long long int)FTK_n_svd_l;
  FTK_svd_U_d_expiw_s_real__ = (double *)calloc(n_a, sizeof(double));
  for (na = 0; na < n_a; na++)
  {
    FTK_svd_U_d_expiw_s_real__[na] = (double)(na % 123) - 61;
  }
  FTK_svd_U_d_expiw_s_imag__ = (double *)calloc(n_a, sizeof(double));
  for (na = 0; na < n_a; na++)
  {
    FTK_svd_U_d_expiw_s_imag__[na] = (double)(na % 125) - 62;
  }
  n_a = FTK_n_delta_v;
  FTK_delta_x_ = (double *)calloc(n_a, sizeof(double));
  for (na = 0; na < n_a; na++)
  {
    FTK_delta_x_[na] = (double)(na % 127) - 63;
  }
  FTK_delta_y_ = (double *)calloc(n_a, sizeof(double));
  for (na = 0; na < n_a; na++)
  {
    FTK_delta_y_[na] = (double)(na % 129) - 65;
  }
  n_a = (unsigned long long int)n_w_max * (unsigned long long int)pm_n_UX_rank * (unsigned long long int)n_S;
  CTF_UX_S_k_q_wnS_real__ = (double *)calloc(n_a, sizeof(double));
  for (na = 0; na < n_a; na++)
  {
    CTF_UX_S_k_q_wnS_real__[na] = (double)(na % 123) - 61;
  }
  CTF_UX_S_k_q_wnS_imag__ = (double *)calloc(n_a, sizeof(double));
  for (na = 0; na < n_a; na++)
  {
    CTF_UX_S_k_q_wnS_imag__[na] = (double)(na % 125) - 62;
  }
  n_a = n_S;
  CTF_UX_S_l2_ = (double *)calloc(n_a, sizeof(double));
  for (na = 0; na < n_a; na++)
  {
    CTF_UX_S_l2_[na] = (double)(1 + na % 123);
  }
  n_a = (unsigned long long int)FTK_n_svd_l * (unsigned long long int)n_w_max * (unsigned long long int)pm_n_UX_rank * (unsigned long long int)n_M;
  svd_VUXM_lwnM_real____ = (double *)calloc(n_a, sizeof(double));
  for (na = 0; na < n_a; na++)
  {
    svd_VUXM_lwnM_real____[na] = (double)(na % 123) - 61;
  }
  svd_VUXM_lwnM_imag____ = (double *)calloc(n_a, sizeof(double));
  for (na = 0; na < n_a; na++)
  {
    svd_VUXM_lwnM_imag____[na] = (double)(na % 125) - 62;
  }
  n_a = FTK_n_delta_v * n_M;
  UX_M_l2_dM__ = (double *)calloc(n_a, sizeof(double));
  for (na = 0; na < n_a; na++)
  {
    UX_M_l2_dM__[na] = (double)(1 + na % 123);
  }
  /* %%%%; */
  n_a = (unsigned long long int)n_w_max * (unsigned long long int)n_S * (unsigned long long int)n_M;
  if (flag_optimize_over_gamma_z)
  {
    n_a = (unsigned long long int)n_S * (unsigned long long int)n_M;
  }
  X_wSM___ = (double *)calloc(n_a, sizeof(double));
  delta_x_wSM___ = (double *)calloc(n_a, sizeof(double));
  delta_y_wSM___ = (double *)calloc(n_a, sizeof(double));
  gamma_z_wSM___ = (double *)calloc(n_a, sizeof(double));
  I_value_wSM___ = (double *)calloc(n_a, sizeof(double));
  /* %%%%; */
  if (flag_test)
  {
    // printf(" %% calling MDA_io_test()\n"); MDA_io_test();
    //  printf(" %% calling test_cblas()\n"); test_cblas();
    // printf(" %% calling test_transpose_float()\n"); test_transpose_float();
    // printf(" %% calling test_transpose_float_complex()\n"); test_transpose_float_complex();
    // printf(" %% calling mex_ampmh_X_wSM___fftw_slow_test()\n"); mex_ampmh_X_wSM___fftw_slow_test();
    // printf(" %% calling mex_ampmh_X_wSM___fftwf_fast_test()\n"); mex_ampmh_X_wSM___fftwf_fast_test();
    // printf(" %% calling mex_ampmh_X_wSM___fft_mkl_prealloc_test()\n"); mex_ampmh_X_wSM___fft_mkl_prealloc_test();
    // printf(" %% calling mex_ampmh_X_wSM___fft_mkl_test()\n"); mex_ampmh_X_wSM___fft_mkl_test();
    // printf(" %% calling dp_ps_mult_immintrin_test()\n"); dp_ps_mult_immintrin_test();
    // printf(" %% calling hp_ps_mult_immintrin_test()\n"); hp_ps_mult_immintrin_test();
    // printf(" %% calling nhp_ps_mult_immintrin_test()\n"); nhp_ps_mult_immintrin_test();
    // printf(" %% calling nhpr_ps_mult_immintrin_test()\n"); nhpr_ps_mult_immintrin_test();
    // printf(" %% finished testing, exiting before calling mex_ampmh_X_wSM___16_omp \n");
    // exit(0);
    /* if (flag_test){ } */}
    /* %%%%; */
    // local_tic(0, t_start_, d_start_);
    time_t t_start, t_end;
    time(&t_start);
    mex_ampmh_X_wSM___16_omp(n_M_per_Mbatch, n_S_per_Sbatch, flag_optimize_over_gamma_z, flag_compute_I_value, tolerance_master, FTK_n_svd_l, FTK_n_delta_v, FTK_svd_U_d_expiw_s_real__, FTK_svd_U_d_expiw_s_imag__, FTK_delta_x_, FTK_delta_y_, n_w_max, pm_n_UX_rank, n_S, CTF_UX_S_k_q_wnS_real__, CTF_UX_S_k_q_wnS_imag__, CTF_UX_S_l2_, n_M, svd_VUXM_lwnM_real____, svd_VUXM_lwnM_imag____, UX_M_l2_dM__, X_wSM___, delta_x_wSM___, delta_y_wSM___, gamma_z_wSM___, I_value_wSM___);
    time(&t_end);
    printf(" %% mex_ampmh_X_wSM___16_omp: took %d seconds.\n", (int)difftime(t_end, t_start));
    if (verbose)
    {
      printf(" %% [finished main]\n");
    }
    return 0;
}
