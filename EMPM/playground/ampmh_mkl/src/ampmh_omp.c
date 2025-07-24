#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <errno.h>
// #include <cblas.h>
#include <string.h>
#include <complex.h> /* <-- must include before fftw3.h to ensure that fftw_complex is compatible */
#include "fftw3.h"	 /* <-- must include after complex.h to ensure that fftw_complex is compatible */
#include <omp.h>
#include <emmintrin.h>
#include <immintrin.h>
#define MKL_Complex8 float complex
#include "mkl.h"

#define PI (3.141592653589793)
#define rup(A, B) ((A) + !!((A) % (B)) * ((B) - ((A) % (B))))
#define maximum(A, B) ((A) > (B) ? (A) : (B))
#define minimum(A, B) ((A) < (B) ? (A) : (B))
#define periodize(A, B, C) ((A) < (B) ? (A) + (C) - (B) : ((A) >= (C) ? (A) - (C) + (B) : (A)))
#define crop(A, B, C) ((A) < (B) ? (B) : ((A) > (C) ? (C) : (A)))

unsigned long long dmax_index(unsigned long long n_a, double *a_)
{
  if (n_a <= 0)
  {
    return 0;
  }
  unsigned long long index = 0;
  double amax = a_[0];
  for (unsigned long long i = 1; i < n_a; i++)
  {
    if (a_[i] > amax)
    {
      amax = a_[i];
      index = i;
    }
  }
  return index;
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* The helper function */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
MKL_LONG mex_ampmh_X_wSM___16_omp_helper(
	unsigned long long n_M_per_Mbatch,
	unsigned long long n_S_per_Sbatch,
	int flag_optimize_over_gamma_z,
	int flag_compute_I_value,
	double tolerance_master,
	unsigned long long FTK_n_svd_l,
	unsigned long long FTK_n_delta_v,
	double *FTK_delta_x_,
	double *FTK_delta_y_,
	unsigned long long n_w_max,
	unsigned long long pm_n_UX_rank,
	unsigned long long n_S,
	double *d_CTF_UX_S_k_q_wnS_real__,
	double *d_CTF_UX_S_k_q_wnS_imag__,
	double *d_CTF_UX_S_l2_,
	unsigned long long n_M,
	double *d_svd_VUXM_lwnM_real____,
	double *d_svd_VUXM_lwnM_imag____,
	double *d_UX_M_l2_dM__,
	double *X_wSM___,
	double *delta_x_wSM___,
	double *delta_y_wSM___,
	double *gamma_z_wSM___,
	double *I_value_wSM___,
	int verbose,
	float complex *c_FTK_svd_U_d_expiw_s_dl__,
	float complex *c_FTK_svd_U_d_expiw_s_ld__,
	float complex *c_CTF_UX_S_k_q_wnS___,
	float complex *c_CTF_UX_S_k_q_nSw___,
	float complex *c_svd_VUXM_lwnM____,
	double *gamma_z_,
	unsigned long long nMbatch,
	unsigned long long n_Mbatch,
	fftw_plan fftw_plan_guru_split_dft_plan,
	fftwf_plan fftwf_plan_guru_split_dft_plan,
	DFTI_DESCRIPTOR *mkl_dfti_fft1d_p_)
{
	int flag_check = 0;
	double d_fnormn = 0.0;
	const float complex c_cblas_alpha = (float complex)1.0;
	const float complex c_cblas_beta = (float complex)0.0;
	double d_l2 = 0;
	MKL_LONG mkl_err = 0;
	/* %%%% */
	if (verbose)
	{
		printf(" %% [entering mex_ampmh_X_wSM___16_omp_helper]\n");
	}
	unsigned long long n_M_sub = (nMbatch == n_Mbatch - 1) ? n_M_per_Mbatch : n_M - n_M_per_Mbatch * nMbatch;
	if (n_M_sub < 1)
	{
		printf(" %% Warning! n_M_sub %llu < 1 in mex_ampmh_X_wSM___16_omp_helper. Exiting...\n", n_M_sub);
		return;
	}
	unsigned long long *index_M_in_Mbatch_ = (unsigned long long *)malloc(n_M_sub * sizeof(unsigned long long));
	for (unsigned long long nM_sub = 0; nM_sub < n_M_sub; nM_sub++)
	{
		index_M_in_Mbatch_[nM_sub] = n_M_per_Mbatch * nMbatch + nM_sub;
	}
	if (verbose)
	{
		printf(" %% nMbatch %llu/%llu: index_M_in_Mbatch_ %llu --> %llu\n", nMbatch, n_Mbatch, index_M_in_Mbatch_[0], index_M_in_Mbatch_[n_M_sub - 1]);
	}

	unsigned long long n_Sbatch = n_S % n_S_per_Sbatch ? n_S / n_S_per_Sbatch + 1 : n_S / n_S_per_Sbatch;

	unsigned long long n_nS = pm_n_UX_rank * n_S;
	unsigned long long n_nM = pm_n_UX_rank * n_M;
	unsigned long long n_Sm = n_S * n_M_sub;
	unsigned long long n_lw = FTK_n_svd_l * n_w_max;
	unsigned long long n_nm = pm_n_UX_rank * n_M_sub;
	unsigned long long n_dw = FTK_n_delta_v * n_w_max;

	unsigned long long n_dSm = FTK_n_delta_v * n_Sm;
	unsigned long long n_lSm = FTK_n_svd_l * n_Sm;
	unsigned long long n_wSm = n_w_max * n_Sm;
	unsigned long long n_lwn = n_lw * pm_n_UX_rank;
	unsigned long long n_lwm = n_lw * n_M_sub;

	unsigned long long n_lwnm = n_lw * n_nm;
	unsigned long long n_lwSm = n_lw * n_Sm;
	unsigned long long n_dwSm = n_w_max * n_dSm;

	float complex *c_svd_VUXM_nmlw____ = (float complex *)mkl_malloc(n_lwnm * sizeof(float complex), 64);
	mkl_comatcopy('C', 'T', n_lw, n_nm, (float)1.0, c_svd_VUXM_lwnM____ + n_lwn * nMbatch * n_M_per_Mbatch,
				  n_lw, c_svd_VUXM_nmlw____, n_nm);
	float complex *c_svd_SVUXM_Smlw____ = (float complex *)mkl_malloc(n_lwSm * sizeof(float complex), 64);
	for (unsigned long long nl = 0; nl < FTK_n_svd_l; nl++)
	{
		for (unsigned long long nw = 0; nw < n_w_max; nw++)
		{
			unsigned long long tabA = nw * n_nS;
			unsigned long long tabB = (nl + nw * FTK_n_svd_l) * n_nm;
			unsigned long long tabC = (nl + nw * FTK_n_svd_l) * n_Sm;
			cblas_cgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, n_S, n_M_sub, pm_n_UX_rank, &c_cblas_alpha, c_CTF_UX_S_k_q_nSw___ + tabA, pm_n_UX_rank, c_svd_VUXM_nmlw____ + tabB, pm_n_UX_rank, &c_cblas_beta, c_svd_SVUXM_Smlw____ + tabC, n_S);
		} /* for (nw=0;nw<n_w_max;nw++){ } */
	}	  /* for (nl=0;nl<FTK_n_svd_l;nl++){ } */
	mkl_free(c_svd_VUXM_nmlw____);
	float complex *c_svd_SVUXM_0inout_wSml____ = (float complex *)mkl_malloc(n_lwSm * sizeof(float complex), 64);
	mkl_comatcopy('C', 'T', n_lSm, n_w_max, (MKL_Complex8)1.0, c_svd_SVUXM_Smlw____, n_lSm, c_svd_SVUXM_0inout_wSml____, n_w_max);
	mkl_free(c_svd_SVUXM_Smlw____);
	mkl_err = DftiComputeBackward(mkl_dfti_fft1d_p_, c_svd_SVUXM_0inout_wSml____);
	if (mkl_err != 0)
	{
		printf(" %% Warning, mkl_err %llu in mex_ampmh_X_wSM____16_omp_helper\n", (int)mkl_err);
		return mkl_err;
	}
	float complex *c_svd_SVUXM_lwSm____ = (float complex *)mkl_malloc(n_lwSm * sizeof(float complex), 64);
	mkl_comatcopy('C', 'T', n_wSm, FTK_n_svd_l, (MKL_Complex8)1.0, c_svd_SVUXM_0inout_wSml____, n_wSm, c_svd_SVUXM_lwSm____, FTK_n_svd_l);
	mkl_free(c_svd_SVUXM_0inout_wSml____);

	unsigned long long n_sm = n_S_per_Sbatch * n_M_sub;
	unsigned long long n_lwsm = n_lw * n_sm;
	unsigned long long n_dwsm = n_dw * n_sm;
	unsigned long long n_dsm = FTK_n_delta_v * n_sm;
	unsigned long long n_wsm = n_w_max * n_sm;

	unsigned long long *index_S_in_Sbatch_ = (unsigned long long *)malloc(n_S_per_Sbatch * sizeof(unsigned long long));
	float complex *c_svd_SVUXM_lwsm____ = (float complex *)mkl_malloc(n_lwsm * sizeof(float complex), 64);
	float complex *c_svd_USESVUXM_dwsm____ = (float complex *)mkl_malloc(n_dwsm * sizeof(float complex), 64);

	double *d_l2_dsm___ = (double *)malloc(n_dsm * sizeof(double));
	double *d_n2_dsm___ = (double *)malloc(n_dsm * sizeof(double));
	double *d_f2_dsm___ = (double *)malloc(n_dsm * sizeof(double));

	double *d_X_sub_dwsm____ = (double *)malloc(n_dwsm * sizeof(double));
	double *I_value_use_dwsm____ = flag_compute_I_value ? (double *)malloc(n_dwsm * sizeof(double)) : NULL;

	for (unsigned long long nSbatch = 0; nSbatch < n_Sbatch; nSbatch++)
	{
		unsigned long long n_S_sub = (nSbatch == n_Sbatch - 1) ? n_S - n_S_per_Sbatch * nSbatch : n_S_per_Sbatch;
		if (n_S_sub < 1)
		{
			continue;
		}
		if (nSbatch == n_Sbatch - 1)
		{
			n_sm = n_S_sub * n_M_sub;
			n_lwsm = n_lw * n_sm;
			n_dwsm = n_dw * n_sm;
			n_dsm = FTK_n_delta_v * n_sm;
			n_wsm = n_w_max * n_sm;
		}
		for (unsigned long long nS_sub = 0; nS_sub < n_S_sub; nS_sub++)
		{
			index_S_in_Sbatch_[nS_sub] = n_S_per_Sbatch * nSbatch + nS_sub;
		}
		// /* %%%% */
		// /* svd_SVUXM_lwsM____ = svd_SVUXM_lwSM____(:,:,1+index_S_in_Sbatch_,:); */
		// /* %%%% */
		for (unsigned long long nM_sub = 0; nM_sub < n_M_sub; nM_sub++)
		{
			for (unsigned long long nS_sub = 0; nS_sub < n_S_sub; nS_sub++)
			{
				unsigned long long nS = index_S_in_Sbatch_[nS_sub];
				unsigned long long tabA = (nS_sub + nM_sub * n_S_sub) * n_lw;
				unsigned long long tabB = (nS + nM_sub * n_S) * n_lw;
				mkl_comatcopy('C', 'N', n_lw, 1, (MKL_Complex8)1.0, c_svd_SVUXM_lwSm____ + tabB, n_lw, c_svd_SVUXM_lwsm____ + tabA, n_lw);
			} /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */
		}	  /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */

		cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, FTK_n_delta_v, n_wsm, FTK_n_svd_l, &c_cblas_alpha, c_FTK_svd_U_d_expiw_s_dl__, FTK_n_delta_v, c_svd_SVUXM_lwsm____, FTK_n_svd_l, &c_cblas_beta, c_svd_USESVUXM_dwsm____, FTK_n_delta_v);

		for (unsigned long long nS_sub = 0; nS_sub < n_S_sub; nS_sub++)
		{
			unsigned long long nS = index_S_in_Sbatch_[nS_sub];
			for (unsigned long long nM_sub = 0; nM_sub < n_M_sub; nM_sub++)
			{
				unsigned long long nM = index_M_in_Mbatch_[nM_sub];
				for (unsigned long long ndelta_v = 0; ndelta_v < FTK_n_delta_v; ndelta_v++)
				{
					unsigned long long tabA = ndelta_v + nM * FTK_n_delta_v;
					unsigned long long tabB = ndelta_v + (nS_sub + nM_sub * n_S_sub) * FTK_n_delta_v;
					d_l2_dsm___[tabB] = sqrt(d_CTF_UX_S_l2_[nS]) * sqrt(d_UX_M_l2_dM__[tabA]);
					d_n2_dsm___[tabB] = 1.0 / maximum(1.0e-14, d_l2_dsm___[tabB]);
					d_f2_dsm___[tabB] = sqrt(d_CTF_UX_S_l2_[nS]) / maximum(1.0e-14, sqrt(d_UX_M_l2_dM__[tabA]));
				} /* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */
			}	  /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */
		}		  /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */

		unsigned long long na = 0;
		for (unsigned long long nM_sub = 0; nM_sub < n_M_sub; nM_sub++)
		{
			for (unsigned long long nS_sub = 0; nS_sub < n_S_sub; nS_sub++)
			{
				for (unsigned long long nw = 0; nw < n_w_max; nw++)
				{
					for (unsigned long long ndelta_v = 0; ndelta_v < FTK_n_delta_v; ndelta_v++)
					{
						unsigned long long tabA = ndelta_v + (nS_sub + nM_sub * n_S_sub) * FTK_n_delta_v;
						d_X_sub_dwsm____[na] = d_n2_dsm___[tabA] * creal(c_svd_USESVUXM_dwsm____[na]);
						na++;
					} /* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */
				}	  /* for (nw=0;nw<n_w_max;nw++){ } */
			}		  /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */
		}			  /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */

		/* %%%% */

		if (flag_compute_I_value)
		{
			unsigned long long na = 0;
			for (unsigned long long nM_sub = 0; nM_sub < n_M_sub; nM_sub++)
			{
				for (unsigned long long nS_sub = 0; nS_sub < n_S_sub; nS_sub++)
				{
					for (unsigned long long nw = 0; nw < n_w_max; nw++)
					{
						for (unsigned long long ndelta_v = 0; ndelta_v < FTK_n_delta_v; ndelta_v++)
						{
							unsigned long long tabA = ndelta_v + (nS_sub + nM_sub * n_S_sub) * FTK_n_delta_v;
							I_value_use_dwsm____[na] = maximum(0.0, d_f2_dsm___[tabA] * d_X_sub_dwsm____[na]);
							na++;
						} /* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */
					}	  /* for (nw=0;nw<n_w_max;nw++){ } */
				}		  /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */
			}			  /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */
		}				  /* if (flag_compute_I_value){ } */

		if (flag_optimize_over_gamma_z)
		{
			na = 0;
			for (unsigned long long nM_sub = 0; nM_sub < n_M_sub; nM_sub++)
			{
				unsigned long long nM = index_M_in_Mbatch_[nM_sub];
				for (unsigned long long nS_sub = 0; nS_sub < n_S_sub; nS_sub++)
				{
					unsigned long long nS = index_S_in_Sbatch_[nS_sub];
					unsigned long long tabA = nS_sub + n_sm;
					unsigned long long ndw = dmax_index(n_dw, d_X_sub_dwsm____ + tabA * n_dw);
					unsigned long long nw = ndw / maximum(1, FTK_n_delta_v);
					unsigned long long ndelta_v = ndw % FTK_n_delta_v;
					unsigned long long tabB = nS + nM * n_S;
					unsigned long long tabC = ndw + tabA * n_dw;
					X_wSM___[tabB] = d_X_sub_dwsm____[tabC];
					delta_x_wSM___[tabB] = FTK_delta_x_[ndelta_v];
					delta_y_wSM___[tabB] = FTK_delta_y_[ndelta_v];
					gamma_z_wSM___[tabB] = gamma_z_[nw];
					if (flag_compute_I_value)
					{
						I_value_wSM___[tabB] = I_value_use_dwsm____[tabC];
					}
					na += n_dw;
				} /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */
			}	  /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */
		}
		else /* if (flag_optimize_over_gamma_z){ } */
		{
			na = 0;
			for (unsigned long long nM_sub = 0; nM_sub < n_M_sub; nM_sub++)
			{
				unsigned long long nM = index_M_in_Mbatch_[nM_sub];
				for (unsigned long long nS_sub = 0; nS_sub < n_S_sub; nS_sub++)
				{
					unsigned long long nS = index_S_in_Sbatch_[nS_sub];
					for (unsigned long long nw = 0; nw < n_w_max; nw++)
					{
						unsigned long long tabA = nw + (nS_sub + nM_sub * n_S_sub) * n_w_max;
						unsigned long long ndelta_v = dmax_index(FTK_n_delta_v, d_X_sub_dwsm____ + tabA * FTK_n_delta_v);
						unsigned long long tabB = nw + (nS + nM * n_S) * n_w_max;
						unsigned long long tabC = ndelta_v + tabA * FTK_n_delta_v;
						X_wSM___[tabB] = d_X_sub_dwsm____[tabC];
						delta_x_wSM___[tabB] = FTK_delta_x_[ndelta_v];
						delta_y_wSM___[tabB] = FTK_delta_y_[ndelta_v];
						gamma_z_wSM___[tabB] = gamma_z_[nw];
						if (flag_compute_I_value)
						{
							I_value_wSM___[tabB] = I_value_use_dwsm____[tabC];
						}
						na += FTK_n_delta_v;
					} /* for (nw=0;nw<n_w_max;nw++){ } */
				}	  /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */
			}		  /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */
		}			  /* if (flag_optimize_over_gamma_z==0){ } */
	}				  /* for (nSbatch=0;nSbatch<n_Sbatch;nSbatch++){ } */

	free(index_S_in_Sbatch_);
	mkl_free(c_svd_SVUXM_lwsm____);
	mkl_free(c_svd_USESVUXM_dwsm____);
	free(d_l2_dsm___);
	free(d_n2_dsm___);
	free(d_f2_dsm___);
	free(d_X_sub_dwsm____);
	if (flag_compute_I_value)
	{
		free(I_value_use_dwsm____);
	}
	mkl_free(c_svd_SVUXM_lwSm____);
	free(index_M_in_Mbatch_);

	/* %%%%%%%%%%%%%%%% */
	if (verbose)
	{
		printf(" %% [finished mex_ampmh_X_wSM___16_omp_helper for batch %llu]\n", nMbatch);
	}
	return mkl_err;
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* The computational routine with mkl openmp support. */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

MKL_LONG mex_ampmh_X_wSM___16_omp(
	const unsigned long long n_M_per_Mbatch,
	const unsigned long long n_S_per_Sbatch,
	int flag_optimize_over_gamma_z,
	int flag_compute_I_value,
	double tolerance_master,
	const unsigned long long FTK_n_svd_l,
	const unsigned long long FTK_n_delta_v,
	double *d_FTK_svd_U_d_expiw_s_real__,
	double *d_FTK_svd_U_d_expiw_s_imag__,
	double *FTK_delta_x_,
	double *FTK_delta_y_,
	const unsigned long long n_w_max,
	const unsigned long long pm_n_UX_rank,
	const unsigned long long n_S,
	double *d_CTF_UX_S_k_q_wnS_real__,
	double *d_CTF_UX_S_k_q_wnS_imag__,
	double *d_CTF_UX_S_l2_,
	const unsigned long long n_M,
	double *d_svd_VUXM_lwnM_real____,
	double *d_svd_VUXM_lwnM_imag____,
	double *d_UX_M_l2_dM__,
	double *X_wSM___,
	double *delta_x_wSM___,
	double *delta_y_wSM___,
	double *gamma_z_wSM___,
	double *I_value_wSM___)
{
	int flag_check = 0;
	double d_fnormn = 0.0;
	int flag_omp = 1;
	int verbose = 1;
	fftw_plan fftw_plan_guru_split_dft_plan;
	fftw_iodim fftw_iodim_use;
	fftwf_plan fftwf_plan_guru_split_dft_plan;
	fftwf_iodim fftwf_iodim_use;
	MKL_LONG mkll_n_dfti = 0; //%<-- DFTI_NUMBER_OF_TRANSFORMS ;
	MKL_LONG mkll_i_dfti = 0; //%<-- DFTI_INPUT_DISTANCE ;
	MKL_LONG mkl_err = 0;
	/* %%%% */
	const unsigned long long n_wnS = n_w_max * pm_n_UX_rank * n_S;
	const unsigned long long n_dl = FTK_n_delta_v * FTK_n_svd_l;
	const unsigned long long n_vM = FTK_n_delta_v * n_M;
	const unsigned long long n_lwnM = FTK_n_svd_l * n_w_max * pm_n_UX_rank * n_M;
	const unsigned long long n_wSM = flag_optimize_over_gamma_z ? n_S * n_M : n_w_max * n_S * n_M;
	printf("n_wSM = %llu\n", n_wSM);

	time_t t_start, t_end;

	if (verbose)
	{
		printf(" %% [entering mex_ampmh_X_wSM___16_omp]\n");
	}
	omp_set_max_active_levels(1);
	mkl_set_dynamic(1); //%<-- directives to control threading. Do not seem necessary. ;
	
	double *gamma_z_ = (double *)mkl_malloc(n_w_max * sizeof(double), 64);
	double w_step_size = 2 * PI / (double)(int)n_w_max;
	for (int nw = 0; nw < (int)n_w_max; nw++)
	{
		gamma_z_[nw] = w_step_size * (double)nw;
	}
	printf("gamma_z_ assigned.\n");

	float complex *c_FTK_svd_U_d_expiw_s_dl__ = (float complex *)mkl_malloc(n_dl * sizeof(float complex), 64);
	for (unsigned long long na = 0; na < n_dl; na++)
	{
		c_FTK_svd_U_d_expiw_s_dl__[na] = (float complex)(d_FTK_svd_U_d_expiw_s_real__[na] + _Complex_I * d_FTK_svd_U_d_expiw_s_imag__[na]);
	}
	printf("c_FTK_svd_U_d_expiw_s_dl__ assigned.\n");

	float complex *c_CTF_UX_S_k_q_wnS___ = (float complex *)mkl_malloc(n_wnS * sizeof(float complex), 64);
	for (unsigned long long na = 0; na < n_wnS; na++)
	{
		c_CTF_UX_S_k_q_wnS___[na] = (float complex)(d_CTF_UX_S_k_q_wnS_real__[na] + _Complex_I * d_CTF_UX_S_k_q_wnS_imag__[na]);
	}
	printf("c_CTF_UX_S_k_q_wnS___ assigned.\n");

	float complex *c_svd_VUXM_lwnM____ = (float complex *)mkl_malloc(n_lwnM * sizeof(float complex), 64);
	for (unsigned long long na = 0; na < n_lwnM; na++)
	{
		c_svd_VUXM_lwnM____[na] = (float complex)(d_svd_VUXM_lwnM_real____[na] + _Complex_I * d_svd_VUXM_lwnM_imag____[na]);
	}
	printf("c_svd_VUXM_lwnM____ assigned.\n");

	/* %%%%%%%%%%%%%%%% */

	float complex *c_FTK_svd_U_d_expiw_s_ld__ = (float complex *)mkl_malloc(n_dl * sizeof(MKL_Complex8), 64);
	mkl_comatcopy('C', 'T', (size_t)FTK_n_svd_l, (size_t)FTK_n_delta_v, (MKL_Complex8)1.0, c_FTK_svd_U_d_expiw_s_dl__, FTK_n_svd_l, c_FTK_svd_U_d_expiw_s_ld__, FTK_n_delta_v);
	printf("c_FTK_svd_U_d_expiw_s_ld__ assigned.\n");

	if (verbose)
	{
		printf(" %% n_M_per_Mbatch %llu\n", n_M_per_Mbatch);
		printf(" %% n_S_per_Sbatch %llu\n", n_S_per_Sbatch);
		printf(" %% flag_optimize_over_gamma_z %d\n", flag_optimize_over_gamma_z);
		printf(" %% flag_compute_I_value %d\n", flag_compute_I_value);
		printf(" %% tolerance_master %0.16f\n", tolerance_master);
		printf(" %% FTK_n_svd_l %llu\n", FTK_n_svd_l);
		printf(" %% FTK_n_delta_v %llu\n", FTK_n_delta_v);
		printf(" %% FTK_delta_x_ %0.16f --> %0.16f\n", FTK_delta_x_[0], FTK_delta_x_[FTK_n_delta_v - 1]);
		printf(" %% FTK_delta_y_ %0.16f --> %0.16f\n", FTK_delta_y_[0], FTK_delta_y_[FTK_n_delta_v - 1]);
		printf(" %% n_w_max %llu\n", n_w_max);
		printf(" %% pm_n_UX_rank %llu\n", pm_n_UX_rank);
		printf(" %% n_S %llu\n", n_S);
		printf(" %% CTF_UX_S_l2_ %0.16f --> %0.16f\n", d_CTF_UX_S_l2_[0], d_CTF_UX_S_l2_[n_S - 1]);
		printf(" %% n_M %llu\n", n_M);
		printf(" %% UX_M_l2_dM__ %0.16f --> %0.16f\n", d_UX_M_l2_dM__[0], d_UX_M_l2_dM__[FTK_n_delta_v * n_M - 1]);
	}

	if (flag_compute_I_value)
	{
		for (unsigned long long na = 0; na < n_wSM; na++)
		{
			I_value_wSM___[na] = 1.0;
		}
		printf("I_value_wSM___ of size %llu assigned. \n", n_wSM);
	} /* if (flag_compute_I_value){ } */

	/* %%%% */
	/* CTF_UX_S_k_q_nSw___ = permute(CTF_UX_S_k_q_wnS___,[2,3,1]); */
	/* %%%% */

	MKL_Complex8 *c_CTF_UX_S_k_q_nSw___ = (MKL_Complex8 *)mkl_malloc(n_wnS * sizeof(MKL_Complex8), 64);
	unsigned long long n_nS = pm_n_UX_rank * n_S;
	printf("pm_n_UX_rank %llu\n", pm_n_UX_rank);
	printf("n_S %llu\n", n_S);
	printf("n_w_max %llu\n", n_w_max);

	time(&t_start);
	mkl_comatcopy('C', 'T', n_w_max, n_nS, (MKL_Complex8)1.0, c_CTF_UX_S_k_q_wnS___, n_w_max, c_CTF_UX_S_k_q_nSw___, n_nS);
	time(&t_end);
	printf("time to transpose c_CTF_UX_S_k_q_nSw___: %0.2f seconds\n", difftime(t_end, t_start));

	/* %%%%%%%%%%%%%%%% */
	unsigned long long n_Mbatch = ceil((double)n_M / (double)n_M_per_Mbatch);
	if (verbose)
	{
		printf("n_Mbatch %llu\n", n_Mbatch);
	}

	time(&t_start);
	DFTI_DESCRIPTOR **mkl_dfti_fft1d_p__ = (DFTI_DESCRIPTOR **)mkl_malloc(n_Mbatch * sizeof(DFTI_DESCRIPTOR *), 64);
	for (unsigned long long nMbatch = 0; nMbatch < n_Mbatch; nMbatch++)
	{
		unsigned long long n_M_sub = n_M_per_Mbatch;
		if (nMbatch == n_Mbatch - 1)
		{
			n_M_sub = n_M - n_M_per_Mbatch * nMbatch;
		}
		mkl_err = 0;
		mkll_i_dfti = n_w_max;
		mkll_n_dfti = FTK_n_svd_l * n_S * n_M_sub; //%<-- do all at once. ;
		mkl_err += DftiCreateDescriptor(&(mkl_dfti_fft1d_p__[nMbatch]), DFTI_SINGLE, DFTI_COMPLEX, 1, mkll_i_dfti);
		mkl_err += DftiSetValue((mkl_dfti_fft1d_p__[nMbatch]), DFTI_NUMBER_OF_TRANSFORMS, mkll_n_dfti);
		mkl_err += DftiSetValue((mkl_dfti_fft1d_p__[nMbatch]), DFTI_INPUT_DISTANCE, mkll_i_dfti);
		mkl_err += DftiSetValue((mkl_dfti_fft1d_p__[nMbatch]), DFTI_OUTPUT_DISTANCE, mkll_i_dfti);
		mkl_err += DftiSetValue((mkl_dfti_fft1d_p__[nMbatch]), DFTI_PLACEMENT, DFTI_INPLACE);
		mkl_err += DftiCommitDescriptor((mkl_dfti_fft1d_p__[nMbatch]));
		if (mkl_err != 0)
		{
			printf(" %% Warning, nMbatch %llu, mkl_err %llu in mex_ampmh_X_wSM___16_omp DftiCommitDescriptor\n", nMbatch, (int)mkl_err);
		}
	} /* for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){ } */
	time(&t_end);
	printf(" %% time to create mkl_dfti_fft1d_p__: %0.2f seconds\n", difftime(t_end, t_start));

	time(&t_start);
	if (flag_omp)
	{
		unsigned long long nMbatch = 0;
#pragma omp parallel private(nMbatch)
		{ /* begin omp parallel */
#pragma omp for schedule(dynamic)
			for (nMbatch = 0; nMbatch < n_Mbatch; nMbatch++)
			{
				MKL_LONG mkl_err_nM = mex_ampmh_X_wSM___16_omp_helper(
					n_M_per_Mbatch,
					n_S_per_Sbatch,
					flag_optimize_over_gamma_z,
					flag_compute_I_value,
					tolerance_master,
					FTK_n_svd_l,
					FTK_n_delta_v,
					FTK_delta_x_,
					FTK_delta_y_,
					n_w_max,
					pm_n_UX_rank,
					n_S,
					d_CTF_UX_S_k_q_wnS_real__,
					d_CTF_UX_S_k_q_wnS_imag__,
					d_CTF_UX_S_l2_,
					n_M,
					d_svd_VUXM_lwnM_real____,
					d_svd_VUXM_lwnM_imag____,
					d_UX_M_l2_dM__,
					X_wSM___,
					delta_x_wSM___,
					delta_y_wSM___,
					gamma_z_wSM___,
					I_value_wSM___,
					verbose,
					c_FTK_svd_U_d_expiw_s_dl__,
					c_FTK_svd_U_d_expiw_s_ld__,
					c_CTF_UX_S_k_q_wnS___,
					c_CTF_UX_S_k_q_nSw___,
					c_svd_VUXM_lwnM____,
					gamma_z_,
					nMbatch,
					n_Mbatch,
					fftw_plan_guru_split_dft_plan,
					fftwf_plan_guru_split_dft_plan,
					mkl_dfti_fft1d_p__[nMbatch]);
				if (mkl_err_nM != 0)
				{
					printf(" %% Warning, nMbatch %llu, mkl_error %llu in mex_ampmh_X_wSM___16_omp_helper\n", nMbatch, (int)mkl_err_nM);
					mkl_err = mkl_err_nM;
				}
			} /* for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){ } */
		}	  /* end omp parallel */
	}
	else /* if (flag_omp){ } */
	{
		for (unsigned long long nMbatch = 0; nMbatch < n_Mbatch; nMbatch++)
		{
			MKL_LONG mkl_err_nM = mex_ampmh_X_wSM___16_omp_helper(
				n_M_per_Mbatch,
				n_S_per_Sbatch,
				flag_optimize_over_gamma_z,
				flag_compute_I_value,
				tolerance_master,
				FTK_n_svd_l,
				FTK_n_delta_v,
				FTK_delta_x_,
				FTK_delta_y_,
				n_w_max,
				pm_n_UX_rank,
				n_S,
				d_CTF_UX_S_k_q_wnS_real__,
				d_CTF_UX_S_k_q_wnS_imag__,
				d_CTF_UX_S_l2_,
				n_M,
				d_svd_VUXM_lwnM_real____,
				d_svd_VUXM_lwnM_imag____,
				d_UX_M_l2_dM__,
				X_wSM___,
				delta_x_wSM___,
				delta_y_wSM___,
				gamma_z_wSM___,
				I_value_wSM___,
				verbose,
				c_FTK_svd_U_d_expiw_s_dl__,
				c_FTK_svd_U_d_expiw_s_ld__,
				c_CTF_UX_S_k_q_wnS___,
				c_CTF_UX_S_k_q_nSw___,
				c_svd_VUXM_lwnM____,
				gamma_z_,
				nMbatch,
				n_Mbatch,
				fftw_plan_guru_split_dft_plan,
				fftwf_plan_guru_split_dft_plan,
				mkl_dfti_fft1d_p__[nMbatch]);
			if (mkl_err_nM != 0)
			{
				printf(" %% Warning, nMbatch %llu, mkl_error %llu in mex_ampmh_X_wSM___16_omp_helper\n", nMbatch, (int)mkl_err_nM);
				mkl_err = mkl_err_nM;
			}
		} /* for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){ } */
	}	  /* if (!flag_omp){ } */
	time(&t_end);
	printf(" %% time to run mex_ampmh_X_wSM___16_omp_helper: %0.2f seconds\n", difftime(t_end, t_start));
	/* %%%%%%%%%%%%%%%% */

	mkl_free(gamma_z_);
	mkl_free(c_FTK_svd_U_d_expiw_s_dl__);
	mkl_free(c_CTF_UX_S_k_q_wnS___);
	mkl_free(c_svd_VUXM_lwnM____);
	mkl_free(c_FTK_svd_U_d_expiw_s_ld__);
	mkl_free(c_CTF_UX_S_k_q_nSw___);

	for (unsigned long long nMbatch = 0; nMbatch < n_Mbatch; nMbatch++)
	{
		mkl_err = DftiFreeDescriptor(&(mkl_dfti_fft1d_p__[nMbatch]));
		if (mkl_err != 0)
		{
			printf(" %% Warning, nMbatch %llu, mkl_err %llu in mex_ampmh_X_wSM___16_omp DftiFreeDescriptor\n", nMbatch, (int)mkl_err);
		}
	} /* for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){ } */
	mkl_free(mkl_dfti_fft1d_p__);

	if (verbose)
	{
		printf(" %% [finished mex_ampmh_X_wSM___16_omp]\n");
	}
	return mkl_err;
}

int main(int argc, char *argv[])
{
	int verbose = 1;
	int flag_test = 0;

	const unsigned long long n_S_per_Sbatch = 24;
	const unsigned long long flag_optimize_over_gamma_z = 0;
	const unsigned long long flag_compute_I_value = 1;
	const double tolerance_master = 0.01;
	const unsigned long long FTK_n_svd_l = 10;
	const unsigned long long FTK_n_delta_v = 29;
	const unsigned long long n_w_max = 98;
	const unsigned long long pm_n_UX_rank = 18;
	const unsigned long long n_S = 993;
	const unsigned long long n_M = 1024;
	const unsigned long long NUM_BATCH_M = 8;// getenv("OMP_NUM_THREADS") ? atoi(getenv("OMP_NUM_THREADS")) : 24;
	const unsigned long long n_M_per_Mbatch = n_M % NUM_BATCH_M ? n_M / NUM_BATCH_M + 1 : n_M / NUM_BATCH_M;

	const unsigned long long n_dl = FTK_n_delta_v * FTK_n_svd_l;
	double *FTK_svd_U_d_expiw_s_real__ = (double *)mkl_malloc(n_dl * sizeof(double), 64);
	double *FTK_svd_U_d_expiw_s_imag__ = (double *)mkl_malloc(n_dl * sizeof(double), 64);
	for (unsigned long long na = 0; na < n_dl; na++)
	{
		FTK_svd_U_d_expiw_s_real__[na] = (double)(na % 123) - 61;
		FTK_svd_U_d_expiw_s_imag__[na] = (double)(na % 125) - 62;
	}
	printf("FTK_svd_U_d_expiw_s_real__ and FTK_svd_U_d_expiw_s_imag__ assigned.\n");

	double *FTK_delta_x_ = (double *)mkl_malloc(FTK_n_delta_v * sizeof(double), 64);
	double *FTK_delta_y_ = (double *)mkl_malloc(FTK_n_delta_v * sizeof(double), 64);
	for (unsigned long long na = 0; na < FTK_n_delta_v; na++)
	{
		FTK_delta_x_[na] = (double)(na % 127) - 63;
		FTK_delta_y_[na] = (double)(na % 129) - 65;
	}
	printf("FTK_delta_x_ and FTK_delta_y_ assigned.\n");

	const unsigned long long n_wnS = n_w_max * pm_n_UX_rank * n_S;
	printf("n_wnS = %llu\n", n_wnS);

	double *CTF_UX_S_k_q_wnS_real__ = (double *)mkl_malloc(n_wnS * sizeof(double), 64);
	double *CTF_UX_S_k_q_wnS_imag__ = (double *)mkl_malloc(n_wnS * sizeof(double), 64);
	for (unsigned long long na = 0; na < n_wnS; na++)
	{
		CTF_UX_S_k_q_wnS_real__[na] = (double)(na % 123) - 61;
		CTF_UX_S_k_q_wnS_imag__[na] = (double)(na % 125) - 62;
	}
	printf("CTF_UX_S_k_q_wnS_real__ and CTF_UX_S_k_q_wnS_imag__ assigned.\n");

	double *CTF_UX_S_l2_ = (double *)mkl_malloc(n_S * sizeof(double), 64);
	for (unsigned long long na = 0; na < n_S; na++)
	{
		CTF_UX_S_l2_[na] = (double)(1 + na % 123);
	}
	printf("CTF_UX_S_l2_ assigned.\n");

	const unsigned long long n_lwnM = FTK_n_svd_l * n_w_max * pm_n_UX_rank * n_M;
	double *svd_VUXM_lwnM_real____ = (double *)mkl_malloc(n_lwnM * sizeof(double), 64);
	double *svd_VUXM_lwnM_imag____ = (double *)mkl_malloc(n_lwnM * sizeof(double), 64);
	for (unsigned long long na = 0; na < n_lwnM; na++)
	{
		svd_VUXM_lwnM_real____[na] = (double)(na % 123) - 61;
		svd_VUXM_lwnM_imag____[na] = (double)(na % 125 - 62);
	}
	printf("svd_VUXM_lwnM_real____ and svd_VUXM_lwnM_imag____ assigned.\n");

	const unsigned long long n_vM = FTK_n_delta_v * n_M;
	double *UX_M_l2_dM__ = (double *)mkl_malloc(n_vM * sizeof(double), 64);
	for (unsigned long long na = 0; na < n_vM; na++)
	{
		UX_M_l2_dM__[na] = (double)(1 + na % 123);
	}
	printf("UX_M_l2_dM__ assigned.\n");

	/* out */
	const unsigned long long n_wSM = flag_optimize_over_gamma_z ? n_S * n_M : n_w_max * n_S * n_M;
	printf("n_wSM = %llu\n", n_wSM);
	double *X_wSM___ = (double *)mkl_malloc(n_wSM * sizeof(double), 64);
	double *delta_x_wSM___ = (double *)mkl_malloc(n_wSM * sizeof(double), 64);
	double *delta_y_wSM___ = (double *)mkl_malloc(n_wSM * sizeof(double), 64);
	double *gamma_z_wSM___ = (double *)mkl_malloc(n_wSM * sizeof(double), 64);
	double *I_value_wSM___ = (double *)mkl_malloc(n_wSM * sizeof(double), 64);
	printf("X_wSM___, delta_x_wSM___, delta_y_wSM___, gamma_z_wSM___, I_value_wSM___ assigned.\n");

	time_t t_start, t_end;
	time(&t_start);
	mex_ampmh_X_wSM___16_omp(n_M_per_Mbatch, n_S_per_Sbatch, flag_optimize_over_gamma_z, flag_compute_I_value, tolerance_master, FTK_n_svd_l, FTK_n_delta_v, FTK_svd_U_d_expiw_s_real__, FTK_svd_U_d_expiw_s_imag__, FTK_delta_x_, FTK_delta_y_, n_w_max, pm_n_UX_rank, n_S, CTF_UX_S_k_q_wnS_real__, CTF_UX_S_k_q_wnS_imag__, CTF_UX_S_l2_, n_M, svd_VUXM_lwnM_real____, svd_VUXM_lwnM_imag____, UX_M_l2_dM__, X_wSM___, delta_x_wSM___, delta_y_wSM___, gamma_z_wSM___, I_value_wSM___);
	time(&t_end);
	printf(" %% mex_ampmh_X_wSM___16_omp: took %d seconds.\n", (int)difftime(t_end, t_start));

	// mkl_free(FTK_svd_U_d_expiw_s_real__);
	// mkl_free(FTK_svd_U_d_expiw_s_imag__);
	// mkl_free(FTK_delta_x_);
	// mkl_free(FTK_delta_y_);
	// mkl_free(CTF_UX_S_k_q_wnS_real__);
	// mkl_free(CTF_UX_S_k_q_wnS_imag__);
	// mkl_free(CTF_UX_S_l2_);
	// mkl_free(svd_VUXM_lwnM_real____);
	// mkl_free(svd_VUXM_lwnM_imag____);
	// mkl_free(UX_M_l2_dM__);
	// mkl_free(X_wSM___);
	// mkl_free(delta_x_wSM___);
	// mkl_free(delta_y_wSM___);
	// mkl_free(gamma_z_wSM___);
	// mkl_free(I_value_wSM___);

	return 0;
}