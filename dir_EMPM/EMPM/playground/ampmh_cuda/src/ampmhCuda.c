#include <stdio.h>
#include <complex.h>
#include <time.h>
#include <math.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cublas_v2.h>
#include <cufft.h>

#include "ampmhKernel.h"

#define PI 3.14159265358979323846

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* The helper function */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void ampmh_cuda(

	cublasHandle_t handle,
	cufftHandle handlecufftPlan,

	int flag_optimize_over_gamma_z,
	int flag_compute_I_value,
	double tolerance_master,

	const unsigned long long FTK_n_svd_l,
	const unsigned long long FTK_n_delta_v,
	const unsigned long long n_w_max,
	const unsigned long long pm_n_UX_rank,
	const unsigned long long n_S,
	const unsigned long long n_M_batch,

	float *CTF_UX_S_l2_cuda_,

	cuFloatComplex *svd_VUXM_lwnm_cuda____,
	float *UX_M_l2_dm_cuda__,
	float *FTK_delta_x_cuda_,
	float *FTK_delta_y_cuda_,
	float *gamma_z_cuda_,

	float *X_wSm_cuda___,
	float *delta_x_wSm_cuda___,
	float *delta_y_wSm_cuda___,
	float *gamma_z_wSm_cuda___,
	float *I_value_wSm_cuda___,

	cuFloatComplex *FTK_svd_U_d_expiw_s_dl_cuda__,
	cuFloatComplex *CTF_UX_S_k_q_nSw_cuda___,

	cuFloatComplex *svd_VUXM_nmlw_cuda____,			// memory allocated for the transposed array
	cuFloatComplex *svd_SVUXM_Smlw_cuda____,		// memory allocated for the result array
	cuFloatComplex *svd_SVUXM_0inout_wSml_cuda____, // memory allocated for the result array
	cuFloatComplex *svd_SVUXM_lwSm_cuda____,		// memory allocated for the result array
	cuFloatComplex *svd_USESVUXM_dwSm_cuda____,		// memory allocated for the result array

	float *l2_dSm_cuda___,
	float *n2_dSm_cuda___,
	float *f2_dSm_cuda___,
	float *X_sub_dwSm_cuda___,
	float *I_value_use_dwSm_cuda___,

	int verbose)
{

	cublasStatus_t status;
	cudaError_t cudaErrorOut = cudaGetLastError();

	const cuFloatComplex alpha = make_cuFloatComplex(1.0f, 0.0f);
	const cuFloatComplex beta = make_cuFloatComplex(0.0f, 0.0f);

	unsigned long long n_nS = pm_n_UX_rank * n_S;
	unsigned long long n_Sm = n_S * n_M_batch;
	unsigned long long n_lw = FTK_n_svd_l * n_w_max;
	unsigned long long n_nm = pm_n_UX_rank * n_M_batch;
	unsigned long long n_dw = FTK_n_delta_v * n_w_max;

	unsigned long long n_dSm = FTK_n_delta_v * n_Sm;
	unsigned long long n_lSm = FTK_n_svd_l * n_Sm;
	unsigned long long n_wSm = n_w_max * n_Sm;
	unsigned long long n_lwn = n_lw * pm_n_UX_rank;
	unsigned long long n_lwm = n_lw * n_M_batch;

	unsigned long long n_lwnm = n_lw * n_nm;
	unsigned long long n_lwSm = n_lw * n_Sm;
	unsigned long long n_dwSm = n_w_max * n_dSm;

	printf("Transposing the array svd_VUXM_lwnm_cuda____ to svd_VUXM_nml...\n");
	/* Transpose the array svd_VUXM_lwnm_cuda____ to svd_VUXM_nmlw_cuda____ */
	status = cublasCgeam(handle, CUBLAS_OP_T, CUBLAS_OP_N, n_lw, n_nm, &alpha, svd_VUXM_lwnm_cuda____, n_nm, &beta, svd_VUXM_nmlw_cuda____, n_lw, svd_VUXM_nmlw_cuda____, n_lw);
	if (status != CUBLAS_STATUS_SUCCESS)
	{
		printf("Failed to transpose the array.\n");
	}
	else
	{
		printf("Succeeded to transpose the array svd_VUXM_lwnm_cuda____ to svd_VUXM_nmlw_cuda____\n");
	}

	calculate_svd_SVUXM_Smlw_wrapper(FTK_n_svd_l, FTK_n_delta_v, n_w_max, n_S, n_M_batch, pm_n_UX_rank, n_nS, n_nm, n_Sm, CTF_UX_S_k_q_nSw_cuda___, svd_VUXM_nmlw_cuda____, svd_SVUXM_Smlw_cuda____);
	if (cudaErrorOut != cudaSuccess)
	{
		printf("CUDA error: %s\n", cudaGetErrorString(cudaErrorOut));
	}
	else
	{
		printf("Succeeded to calculate dSm.\n");
	}

	/* transpose the array svd_SVUXM_Smlw_cuda____ to svd_SVUXM_0inout_wSml_cuda____ */
	status = cublasCgeam(handle, CUBLAS_OP_T, CUBLAS_OP_N, n_lSm, n_w_max, &alpha, svd_SVUXM_Smlw_cuda____, n_w_max, &beta, svd_SVUXM_0inout_wSml_cuda____, n_lSm, svd_SVUXM_0inout_wSml_cuda____, n_lSm);
	if (status != CUBLAS_STATUS_SUCCESS)
	{
		printf("Failed to transpose the array.\n");
	}
	else
	{
		printf("Succeeded to transpose the array svd_SVUXM_Smlw_cuda____ to svd_SVUXM_0inout_wSml_cuda____\n");
	}

	/* perform inversed FFT */
	cufftResult resultcufftPlan;
	resultcufftPlan = cufftExecC2C(handlecufftPlan, svd_SVUXM_0inout_wSml_cuda____, svd_SVUXM_0inout_wSml_cuda____, CUFFT_INVERSE);
	if (resultcufftPlan != CUFFT_SUCCESS)
	{
		printf("cufftExecC2C failed!\n");
	}
	else
	{
		printf("Succeeded to perform inversed FFT\n");
	}

	// transpose the array svd_SVUXM_0inout_wSml_cuda____ to svd_SVUXM_lwSm_cuda____
	status = cublasCgeam(handle, CUBLAS_OP_T, CUBLAS_OP_N, n_wSm, FTK_n_svd_l, &alpha, svd_SVUXM_0inout_wSml_cuda____, FTK_n_svd_l, &beta, svd_SVUXM_lwSm_cuda____, n_wSm, svd_SVUXM_lwSm_cuda____, n_wSm);
	if (status != CUBLAS_STATUS_SUCCESS)
	{
		printf("Failed to transpose the array svd_SVUXM_0inout_wSml_cuda____ to svd_SVUXM_lwSm_cuda____\n");
		printf("status: %d\n", status);
	}
	else
	{
		printf("Succeeded to transpose the array svd_SVUXM_0inout_wSml_cuda____ to svd_SVUXM_lwSm_cuda____\n");
	}

	status = cublasCgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, FTK_n_delta_v, n_wSm, FTK_n_svd_l, &alpha, FTK_svd_U_d_expiw_s_dl_cuda__, FTK_n_delta_v, svd_SVUXM_lwSm_cuda____, FTK_n_svd_l, &beta, svd_USESVUXM_dwSm_cuda____, FTK_n_delta_v);
	if (status != CUBLAS_STATUS_SUCCESS)
	{
		printf("Failed to perform cublasCgemm on svd_SVUXM_lwSm_cuda____ and FTK_svd_U_d_expiw_s_dl_cuda__.\n");
	}
	else
	{
		printf("Succeeded to perform cublasCgemm on svd_SVUXM_lwSm_cuda____ and FTK_svd_U_d_expiw_s_dl_cuda__.\n");
	}

	// calculate l2_dSm
	calculate_dSm_wrapper(n_S, n_M_batch, FTK_n_delta_v, n_w_max, CTF_UX_S_l2_cuda_, UX_M_l2_dm_cuda__, svd_USESVUXM_dwSm_cuda____, l2_dSm_cuda___, n2_dSm_cuda___, f2_dSm_cuda___, X_sub_dwSm_cuda___);
	cudaErrorOut = cudaGetLastError();
	if (cudaErrorOut != cudaSuccess)
	{
		printf("CUDA error: %s\n", cudaGetErrorString(cudaErrorOut));
	}
	else
	{
		printf("Succeeded to calculate dSm.\n");
	}

	if (flag_compute_I_value)
	{
		compute_I_value_wrapper(n_M_batch, n_S, n_w_max, FTK_n_delta_v, f2_dSm_cuda___, X_sub_dwSm_cuda___, I_value_use_dwSm_cuda___);
		cudaErrorOut = cudaGetLastError();
		if (cudaErrorOut != cudaSuccess)
		{
			printf("CUDA error: %s\n", cudaGetErrorString(cudaErrorOut));
		}
		else
		{
			printf("Succeeded to compute I_value.\n");
		}
	}

	optimized_output_wrapper(n_M_batch, n_S, n_Sm, n_w_max, FTK_n_delta_v, n_dw, X_sub_dwSm_cuda___, FTK_delta_x_cuda_, FTK_delta_y_cuda_, gamma_z_cuda_, I_value_use_dwSm_cuda___, X_wSm_cuda___, delta_x_wSm_cuda___, delta_y_wSm_cuda___, gamma_z_wSm_cuda___, I_value_wSm_cuda___, flag_optimize_over_gamma_z);
	cudaErrorOut = cudaGetLastError();
	if (cudaErrorOut != cudaSuccess)
	{
		printf("CUDA error: %s\n", cudaGetErrorString(cudaErrorOut));
	}
	else
	{
		printf("Succeeded to compute maxlikelihood_output.\n");
	}

	return;
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* The computational routine with mkl openmp support. */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
int ampmh_cuda_wrapper(
	int flag_optimize_over_gamma_z,
	int flag_compute_I_value,
	float tolerance_master,
	const unsigned long long FTK_n_svd_l,
	const unsigned long long FTK_n_delta_v,
	cuFloatComplex *FTK_svd_U_d_expiw_s_dl__,
	float *FTK_delta_x_,
	float *FTK_delta_y_,
	const unsigned long long n_w_max,
	const unsigned long long pm_n_UX_rank,
	const unsigned long long n_S,
	cuFloatComplex *CTF_UX_S_k_q_wnS___,
	float *CTF_UX_S_l2_,
	const unsigned long long n_M,
	cuFloatComplex *svd_VUXM_lwnM____,
	float *UX_M_l2_dM__,
	float *X_wSM___,
	float *delta_x_wSM___,
	float *delta_y_wSM___,
	float *gamma_z_wSM___,
	float *I_value_wSM___)
{

	int flag_check = 0;
	double d_fnormn = 0.0;
	int flag_omp = 1;
	int verbose = 1;
	cudaError_t cuda_err;

	const cuFloatComplex alpha = make_cuFloatComplex(1.0f, 0.0f);
	const cuFloatComplex beta = make_cuFloatComplex(0.0f, 0.0f);

	cudaError_t cudaErrorOut;

	const unsigned long long n_wnS = n_w_max * pm_n_UX_rank * n_S;
	const unsigned long long n_dl = FTK_n_delta_v * FTK_n_svd_l;
	const unsigned long long n_dM = FTK_n_delta_v * n_M;
	const unsigned long long n_lwnM = FTK_n_svd_l * n_w_max * pm_n_UX_rank * n_M;
	const unsigned long long n_wSM = flag_optimize_over_gamma_z ? n_S * n_M : n_w_max * n_S * n_M;
	const unsigned long long n_nS = pm_n_UX_rank * n_S;
	const unsigned long long n_nM = pm_n_UX_rank * n_M;

	time_t t_start, t_end;

	if (verbose)
	{
		printf(" %% [entering mex_ampmh_X_wSM___16_omp]\n");
	}

	float *gamma_z_ = (float *)malloc(n_w_max * sizeof(float));
	float w_step_size = 2 * PI / (float)n_w_max;
	for (int nw = 0; nw < (int)n_w_max; nw++)
	{
		gamma_z_[nw] = w_step_size * (float)nw;
	}
	// copy to gpu
	float *gamma_z_cuda_ = NULL;
	cudaMalloc((void **)&gamma_z_cuda_, n_w_max * sizeof(float));
	cudaMemcpy(gamma_z_cuda_, gamma_z_, n_w_max * sizeof(float), cudaMemcpyHostToDevice);
	free(gamma_z_);

	float *CTF_UX_S_l2_cuda_ = NULL;
	cudaMalloc((void **)&CTF_UX_S_l2_cuda_, n_S * sizeof(float));
	cudaMemcpy(CTF_UX_S_l2_cuda_, CTF_UX_S_l2_, n_S * sizeof(float), cudaMemcpyHostToDevice);

	cuFloatComplex *FTK_svd_U_d_expiw_s_dl_cuda__ = NULL;
	cudaMalloc((void **)&FTK_svd_U_d_expiw_s_dl_cuda__, n_dl * sizeof(cuFloatComplex));
	cudaMemcpy(FTK_svd_U_d_expiw_s_dl_cuda__, FTK_svd_U_d_expiw_s_dl__, n_dl * sizeof(cuFloatComplex), cudaMemcpyHostToDevice);

	if (verbose)
	{
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
		printf(" %% CTF_UX_S_l2_ %0.16f --> %0.16f\n", CTF_UX_S_l2_[0], CTF_UX_S_l2_[n_S - 1]);
		printf(" %% n_M %llu\n", n_M);
		printf(" %% UX_M_l2_dM__ %0.16f --> %0.16f\n", UX_M_l2_dM__[0], UX_M_l2_dM__[FTK_n_delta_v * n_M - 1]);
	}

	if (flag_compute_I_value)
	{
		for (unsigned long long i = 0; i < n_wSM; i++)
		{
			I_value_wSM___[i] = 1.0;
		}
		printf("I_value_wSM___ of size %llu assigned. \n", n_wSM);
	} /* if (flag_compute_I_value){ } */

	cublasHandle_t handle;
	cublasCreate(&handle);
	cublasStatus_t status;

	time(&t_start);
	cuFloatComplex *CTF_UX_S_k_q_nSw_cuda___ = NULL;
	cuFloatComplex *CTF_UX_S_k_q_wnS_cuda___ = NULL;
	cudaMalloc((void **)&CTF_UX_S_k_q_nSw_cuda___, n_wnS * sizeof(cuFloatComplex));
	cudaMalloc((void **)&CTF_UX_S_k_q_wnS_cuda___, n_wnS * sizeof(cuFloatComplex));
	cudaMemcpy(CTF_UX_S_k_q_wnS_cuda___, CTF_UX_S_k_q_wnS___, n_wnS * sizeof(cuFloatComplex), cudaMemcpyHostToDevice);
	// Transpose the array c_CTF_UX_S_k_q_wnS___ to c_CTF_UX_S_k_q_nSw___
	status = cublasCgeam(handle, CUBLAS_OP_T, CUBLAS_OP_N, n_nS, n_w_max, &alpha, CTF_UX_S_k_q_wnS_cuda___, n_w_max, &beta, CTF_UX_S_k_q_nSw_cuda___, n_nS, CTF_UX_S_k_q_nSw_cuda___, n_nS);
	if (status != CUBLAS_STATUS_SUCCESS)
	{
		printf("Failed to transpose the array.\n");
	}
	cudaFree(CTF_UX_S_k_q_wnS_cuda___);
	time(&t_end);
	printf("time to transpose the array: %0.2f seconds\n", difftime(t_end, t_start));

	float *FTK_delta_x_cuda_ = NULL;
	float *FTK_delta_y_cuda_ = NULL;
	cudaMalloc((void **)&FTK_delta_x_cuda_, FTK_n_delta_v * sizeof(float));
	cudaMalloc((void **)&FTK_delta_y_cuda_, FTK_n_delta_v * sizeof(float));
	cudaMemcpy(FTK_delta_x_cuda_, FTK_delta_x_, FTK_n_delta_v * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(FTK_delta_y_cuda_, FTK_delta_y_, FTK_n_delta_v * sizeof(float), cudaMemcpyHostToDevice);

	cuFloatComplex *CTF_UX_S_k_q_nSw_cuda__ = NULL;
	cudaMalloc((void **)&CTF_UX_S_k_q_nSw_cuda__, n_wnS * sizeof(cuFloatComplex));
	cudaMemcpy(CTF_UX_S_k_q_nSw_cuda__, CTF_UX_S_k_q_wnS___, n_wnS * sizeof(cuFloatComplex), cudaMemcpyHostToDevice);

	cuFloatComplex *svd_VUXM_lwnm_cuda____ = NULL;
	float *UX_M_l2_dm_cuda__ = NULL;
	float *X_wSm_cuda___ = NULL;
	float *delta_x_wSm_cuda___ = NULL;
	float *delta_y_wSm_cuda___ = NULL;
	float *gamma_z_wSm_cuda___ = NULL;
	float *I_value_wSm_cuda___ = NULL;
	cuFloatComplex *svd_VUXM_nmlw_cuda____ = NULL;
	cuFloatComplex *svd_SVUXM_Smlw_cuda____ = NULL;
	cuFloatComplex *svd_SVUXM_0inout_wSml_cuda____ = NULL;
	cuFloatComplex *svd_SVUXM_lwSm_cuda____ = NULL;
	cuFloatComplex *svd_USESVUXM_dwSm_cuda____ = NULL;
	float *l2_dSm_cuda___ = NULL;
	float *n2_dSm_cuda___ = NULL;
	float *f2_dSm_cuda___ = NULL;
	float *X_sub_dwSm_cuda___ = NULL;
	float *I_value_use_dwSm_cuda___ = NULL;

	// check available memory left on the GPU
	size_t free_byte;
	size_t total_byte;
	cudaMemGetInfo(&free_byte, &total_byte);

	// batch data according to the available memory
	unsigned long long n_M_per_Mbatch;

	unsigned long long n_lwn = FTK_n_svd_l * n_w_max * pm_n_UX_rank;
	unsigned long long n_wS = n_w_max * n_S;
	unsigned long long n_lwS = FTK_n_svd_l * n_wS;
	unsigned long long n_dwS = FTK_n_delta_v * n_wS;
	unsigned long long n_dS = FTK_n_delta_v * n_S;

	// memory allocation for the batched data
	size_t memory_size = (n_lwS * 3 + n_lwn * 2 + n_dwS) * n_M * sizeof(cuFloatComplex) + (FTK_n_delta_v + n_wS * 4 + n_dS * 3 + n_dwS * 2) * n_M * sizeof(float);
	if (flag_compute_I_value)
	{
		memory_size += n_wSM * sizeof(float);
	}
	int n_Mbatch = 1;
	n_M_per_Mbatch = n_M;
	int flag_memory_allocated = 0;
	while ((n_Mbatch <= n_M) && (!flag_memory_allocated))
	{
		n_M_per_Mbatch = n_M % n_Mbatch ? n_M / n_Mbatch + 1 : n_M / n_Mbatch;
		memory_size = (n_lwS * 3 + n_lwn * 2) * n_M_per_Mbatch * sizeof(cuFloatComplex) + (FTK_n_delta_v + n_wS * 4 + n_dS * 3 + n_dwS * 2) * n_M_per_Mbatch * sizeof(float);
		if (flag_compute_I_value)
		{
			memory_size += n_wS * n_M_per_Mbatch * sizeof(float);
		}
		if (memory_size > free_byte / 4)
		{
			n_Mbatch++;
			continue;
		}

		// try to allocate memory on the GPU
		printf("Trying to allocate memory on the GPU for n_Mbatch %d n_M_per_Mbatch %d \n", n_Mbatch, n_M_per_Mbatch);
		cudaErrorOut = cudaMalloc((void **)&svd_VUXM_lwnm_cuda____, n_lwn * n_M_per_Mbatch * sizeof(cuFloatComplex));
		if (cudaErrorOut != cudaSuccess)
		{
			printf("Failed to allocate memory on the GPU for svd_VUXM_lwnm_cuda____.\n");
			flag_memory_allocated = 0;
			n_Mbatch++;
			continue;
		}
		else
		{
			printf("Memory allocated on the GPU for svd_VUXM_lwnm_cuda____.\n");
		}
		cudaErrorOut = cudaMalloc((void **)&UX_M_l2_dm_cuda__, FTK_n_delta_v * n_M_per_Mbatch * sizeof(float));
		if (cudaErrorOut != cudaSuccess)
		{
			printf("Failed to allocate memory on the GPU for UX_M_l2_dm_cuda__.\n");
			flag_memory_allocated = 0;
			cudaFree(svd_VUXM_lwnm_cuda____);
			n_Mbatch++;
			continue;
		}
		else
		{
			printf("Memory allocated on the GPU for UX_M_l2_dm_cuda__.\n");
		}
		cudaErrorOut = cudaMalloc((void **)&X_wSm_cuda___, n_wS * n_M_per_Mbatch * sizeof(float));
		if (cudaErrorOut != cudaSuccess)
		{
			printf("Failed to allocate memory on the GPU for X_wSm_cuda___.\n");
			flag_memory_allocated = 0;
			cudaFree(svd_VUXM_lwnm_cuda____);
			cudaFree(UX_M_l2_dm_cuda__);
			n_Mbatch++;
			continue;
		}
		else
		{
			printf("Memory allocated on the GPU for X_wSm_cuda___.\n");
		}
		cudaErrorOut = cudaMalloc((void **)&delta_x_wSm_cuda___, n_wS * n_M_per_Mbatch * sizeof(float));
		if (cudaErrorOut != cudaSuccess)
		{
			printf("Failed to allocate memory on the GPU for delta_x_wSm_cuda___.\n");
			flag_memory_allocated = 0;
			cudaFree(svd_VUXM_lwnm_cuda____);
			cudaFree(UX_M_l2_dm_cuda__);
			cudaFree(X_wSm_cuda___);
			n_Mbatch++;
			continue;
		}
		else
		{
			printf("Memory allocated on the GPU for delta_x_wSm_cuda___.\n");
		}
		cudaErrorOut = cudaMalloc((void **)&delta_y_wSm_cuda___, n_wS * n_M_per_Mbatch * sizeof(float));
		if (cudaErrorOut != cudaSuccess)
		{
			printf("Failed to allocate memory on the GPU for delta_y_wSm_cuda___.\n");
			flag_memory_allocated = 0;
			cudaFree(svd_VUXM_lwnm_cuda____);
			cudaFree(UX_M_l2_dm_cuda__);
			cudaFree(X_wSm_cuda___);
			cudaFree(delta_x_wSm_cuda___);
			n_Mbatch++;
			continue;
		}
		else
		{
			printf("Memory allocated on the GPU for delta_y_wSm_cuda___.\n");
		}
		cudaErrorOut = cudaMalloc((void **)&gamma_z_wSm_cuda___, n_wS * n_M_per_Mbatch * sizeof(float));
		if (cudaErrorOut != cudaSuccess)
		{
			printf("Failed to allocate memory on the GPU for gamma_z_wSm_cuda___.\n");
			flag_memory_allocated = 0;
			cudaFree(svd_VUXM_lwnm_cuda____);
			cudaFree(UX_M_l2_dm_cuda__);
			cudaFree(X_wSm_cuda___);
			cudaFree(delta_x_wSm_cuda___);
			cudaFree(delta_y_wSm_cuda___);
			n_Mbatch++;
			continue;
		}
		else
		{
			printf("Memory allocated on the GPU for gamma_z_wSm_cuda___.\n");
		}
		if (flag_compute_I_value)
		{
			cudaErrorOut = cudaMalloc((void **)&I_value_wSm_cuda___, n_wS * n_M_per_Mbatch * sizeof(float));
			if (cudaErrorOut != cudaSuccess)
			{
				printf("Failed to allocate memory on the GPU for I_value_wSm_cuda___.\n");
				flag_memory_allocated = 0;
				cudaFree(svd_VUXM_lwnm_cuda____);
				cudaFree(UX_M_l2_dm_cuda__);
				cudaFree(X_wSm_cuda___);
				cudaFree(delta_x_wSm_cuda___);
				cudaFree(delta_y_wSm_cuda___);
				cudaFree(gamma_z_wSm_cuda___);
				n_Mbatch++;
				continue;
			}
		}
		else
		{
			printf("Memory allocated on the GPU for I_value_wSm_cuda___.\n");
		}
		// place holder for transposed array
		cudaErrorOut = cudaMalloc((void **)&svd_VUXM_nmlw_cuda____, n_lwn * n_M_per_Mbatch * sizeof(cuFloatComplex));
		if (cudaErrorOut != cudaSuccess)
		{
			printf("Failed to allocate memory on the GPU for svd_VUXM_nmlw_cuda____.\n");
			flag_memory_allocated = 0;
			cudaFree(svd_VUXM_lwnm_cuda____);
			cudaFree(UX_M_l2_dm_cuda__);
			cudaFree(X_wSm_cuda___);
			cudaFree(delta_x_wSm_cuda___);
			cudaFree(delta_y_wSm_cuda___);
			cudaFree(gamma_z_wSm_cuda___);
			if (flag_compute_I_value)
			{
				cudaFree(I_value_wSm_cuda___);
			}
			n_Mbatch++;
			continue;
		}
		else
		{
			printf("Memory allocated on the GPU for svd_VUXM_nmlw_cuda____.\n");
		}
		cudaErrorOut = cudaMalloc((void **)&svd_SVUXM_Smlw_cuda____, n_lwS * n_M_per_Mbatch * sizeof(cuFloatComplex));
		if (cudaErrorOut != cudaSuccess)
		{
			printf("Failed to allocate memory on the GPU for svd_SVUXM_Smlw_cuda____.\n");
			flag_memory_allocated = 0;
			cudaFree(svd_VUXM_lwnm_cuda____);
			cudaFree(UX_M_l2_dm_cuda__);
			cudaFree(X_wSm_cuda___);
			cudaFree(delta_x_wSm_cuda___);
			cudaFree(delta_y_wSm_cuda___);
			cudaFree(gamma_z_wSm_cuda___);
			if (flag_compute_I_value)
			{
				cudaFree(I_value_wSm_cuda___);
			}
			cudaFree(svd_VUXM_nmlw_cuda____);
			n_Mbatch++;
			continue;
		}
		else
		{
			printf("Memory allocated on the GPU for svd_SVUXM_Smlw_cuda____.\n");
		}
		cudaErrorOut = cudaMalloc((void **)&svd_SVUXM_0inout_wSml_cuda____, n_lwS * n_M_per_Mbatch * sizeof(cuFloatComplex));
		if (cudaErrorOut != cudaSuccess)
		{
			printf("Failed to allocate memory on the GPU for svd_SVUXM_0inout_wSml_cuda____.\n");
			flag_memory_allocated = 0;
			cudaFree(svd_VUXM_lwnm_cuda____);
			cudaFree(UX_M_l2_dm_cuda__);
			cudaFree(X_wSm_cuda___);
			cudaFree(delta_x_wSm_cuda___);
			cudaFree(delta_y_wSm_cuda___);
			cudaFree(gamma_z_wSm_cuda___);
			if (flag_compute_I_value)
			{
				cudaFree(I_value_wSm_cuda___);
			}
			cudaFree(svd_VUXM_nmlw_cuda____);
			cudaFree(svd_SVUXM_Smlw_cuda____);
			n_Mbatch++;
			continue;
		}
		else
		{
			printf("Memory allocated on the GPU for svd_SVUXM_0inout_wSml_cuda____.\n");
		}
		cudaErrorOut = cudaMalloc((void **)&svd_SVUXM_lwSm_cuda____, n_lwS * n_M_per_Mbatch * sizeof(cuFloatComplex));
		if (cudaErrorOut != cudaSuccess)
		{
			printf("Failed to allocate memory on the GPU for svd_SVUXM_lwSm_cuda____.\n");
			flag_memory_allocated = 0;
			cudaFree(svd_VUXM_lwnm_cuda____);
			cudaFree(UX_M_l2_dm_cuda__);
			cudaFree(X_wSm_cuda___);
			cudaFree(delta_x_wSm_cuda___);
			cudaFree(delta_y_wSm_cuda___);
			cudaFree(gamma_z_wSm_cuda___);
			if (flag_compute_I_value)
			{
				cudaFree(I_value_wSm_cuda___);
			}
			cudaFree(svd_VUXM_nmlw_cuda____);
			cudaFree(svd_SVUXM_Smlw_cuda____);
			cudaFree(svd_SVUXM_0inout_wSml_cuda____);
			n_Mbatch++;
			continue;
		}
		else
		{
			printf("Memory allocated on the GPU for svd_SVUXM_lwSm_cuda____.\n");
		}
		cudaErrorOut = cudaMalloc((void **)&svd_USESVUXM_dwSm_cuda____, n_dwS * n_M_per_Mbatch * sizeof(cuFloatComplex));
		if (cudaErrorOut != cudaSuccess)
		{
			printf("Failed to allocate memory on the GPU for svd_USESVUXM_dwSm_cuda____.\n");
			flag_memory_allocated = 0;
			cudaFree(svd_VUXM_lwnm_cuda____);
			cudaFree(UX_M_l2_dm_cuda__);
			cudaFree(X_wSm_cuda___);
			cudaFree(delta_x_wSm_cuda___);
			cudaFree(delta_y_wSm_cuda___);
			cudaFree(gamma_z_wSm_cuda___);
			if (flag_compute_I_value)
			{
				cudaFree(I_value_wSm_cuda___);
			}
			cudaFree(svd_VUXM_nmlw_cuda____);
			cudaFree(svd_SVUXM_Smlw_cuda____);
			cudaFree(svd_SVUXM_0inout_wSml_cuda____);
			cudaFree(svd_SVUXM_lwSm_cuda____);
			n_Mbatch++;
			continue;
		}
		else
		{
			printf("Memory allocated on the GPU for svd_USESVUXM_dwSm_cuda____.\n");
		}
		cudaErrorOut = cudaMalloc((void **)&l2_dSm_cuda___, n_dS * n_M_per_Mbatch * sizeof(float));
		if (cudaErrorOut != cudaSuccess)
		{
			printf("Failed to allocate memory on the GPU for l2_dSm_cuda___.\n");
			flag_memory_allocated = 0;
			cudaFree(svd_VUXM_lwnm_cuda____);
			cudaFree(UX_M_l2_dm_cuda__);
			cudaFree(X_wSm_cuda___);
			cudaFree(delta_x_wSm_cuda___);
			cudaFree(delta_y_wSm_cuda___);
			cudaFree(gamma_z_wSm_cuda___);
			if (flag_compute_I_value)
			{
				cudaFree(I_value_wSm_cuda___);
			}
			cudaFree(svd_VUXM_nmlw_cuda____);
			cudaFree(svd_SVUXM_Smlw_cuda____);
			cudaFree(svd_SVUXM_0inout_wSml_cuda____);
			cudaFree(svd_SVUXM_lwSm_cuda____);
			cudaFree(svd_USESVUXM_dwSm_cuda____);
			n_Mbatch++;
			continue;
		}
		else
		{
			printf("Memory allocated on the GPU for l2_dSm_cuda___.\n");
		}
		cudaErrorOut = cudaMalloc((void **)&n2_dSm_cuda___, n_dS * n_M_per_Mbatch * sizeof(float));
		if (cudaErrorOut != cudaSuccess)
		{
			printf("Failed to allocate memory on the GPU for n2_dSm_cuda___.\n");
			flag_memory_allocated = 0;
			cudaFree(svd_VUXM_lwnm_cuda____);
			cudaFree(UX_M_l2_dm_cuda__);
			cudaFree(X_wSm_cuda___);
			cudaFree(delta_x_wSm_cuda___);
			cudaFree(delta_y_wSm_cuda___);
			cudaFree(gamma_z_wSm_cuda___);
			if (flag_compute_I_value)
			{
				cudaFree(I_value_wSm_cuda___);
			}
			cudaFree(svd_VUXM_nmlw_cuda____);
			cudaFree(svd_SVUXM_Smlw_cuda____);
			cudaFree(svd_SVUXM_0inout_wSml_cuda____);
			cudaFree(svd_SVUXM_lwSm_cuda____);
			cudaFree(svd_USESVUXM_dwSm_cuda____);
			cudaFree(l2_dSm_cuda___);
			n_Mbatch++;
			continue;
		}
		else
		{
			printf("Memory allocated on the GPU for n2_dSm_cuda___.\n");
		}
		cudaErrorOut = cudaMalloc((void **)&f2_dSm_cuda___, n_dS * n_M_per_Mbatch * sizeof(float));
		if (cudaErrorOut != cudaSuccess)
		{
			printf("Failed to allocate memory on the GPU for f2_dSm_cuda___.\n");
			flag_memory_allocated = 0;
			cudaFree(svd_VUXM_lwnm_cuda____);
			cudaFree(UX_M_l2_dm_cuda__);
			cudaFree(X_wSm_cuda___);
			cudaFree(delta_x_wSm_cuda___);
			cudaFree(delta_y_wSm_cuda___);
			cudaFree(gamma_z_wSm_cuda___);
			if (flag_compute_I_value)
			{
				cudaFree(I_value_wSm_cuda___);
			}
			cudaFree(svd_VUXM_nmlw_cuda____);
			cudaFree(svd_SVUXM_Smlw_cuda____);
			cudaFree(svd_SVUXM_0inout_wSml_cuda____);
			cudaFree(svd_SVUXM_lwSm_cuda____);
			cudaFree(svd_USESVUXM_dwSm_cuda____);
			cudaFree(l2_dSm_cuda___);
			cudaFree(n2_dSm_cuda___);
			n_Mbatch++;
			continue;
		}
		else
		{
			printf("Memory allocated on the GPU for f2_dSm_cuda___.\n");
		}
		cudaErrorOut = cudaMalloc((void **)&X_sub_dwSm_cuda___, n_dwS * n_M_per_Mbatch * sizeof(float));
		if (cudaErrorOut != cudaSuccess)
		{
			printf("Failed to allocate memory on the GPU for X_sub_dwSm_cuda___.\n");
			flag_memory_allocated = 0;
			cudaFree(svd_VUXM_lwnm_cuda____);
			cudaFree(UX_M_l2_dm_cuda__);
			cudaFree(X_wSm_cuda___);
			cudaFree(delta_x_wSm_cuda___);
			cudaFree(delta_y_wSm_cuda___);
			cudaFree(gamma_z_wSm_cuda___);
			if (flag_compute_I_value)
			{
				cudaFree(I_value_wSm_cuda___);
			}
			cudaFree(svd_VUXM_nmlw_cuda____);
			cudaFree(svd_SVUXM_Smlw_cuda____);
			cudaFree(svd_SVUXM_0inout_wSml_cuda____);
			cudaFree(svd_SVUXM_lwSm_cuda____);
			cudaFree(svd_USESVUXM_dwSm_cuda____);
			cudaFree(l2_dSm_cuda___);
			cudaFree(n2_dSm_cuda___);
			cudaFree(f2_dSm_cuda___);
			n_Mbatch++;
			continue;
		}
		else
		{
			printf("Memory allocated on the GPU for X_sub_dwSm_cuda___.\n");
		}
		if (flag_compute_I_value)
		{
			cudaErrorOut = cudaMalloc((void **)&I_value_use_dwSm_cuda___, n_dwS * sizeof(float));
			if (cudaErrorOut != cudaSuccess)
			{
				printf("Failed to allocate memory on the GPU for I_value_use_dwSm_cuda___.\n");
				flag_memory_allocated = 0;
				cudaFree(svd_VUXM_lwnm_cuda____);
				cudaFree(UX_M_l2_dm_cuda__);
				cudaFree(X_wSm_cuda___);
				cudaFree(delta_x_wSm_cuda___);
				cudaFree(delta_y_wSm_cuda___);
				cudaFree(gamma_z_wSm_cuda___);
				cudaFree(I_value_wSm_cuda___);
				cudaFree(svd_VUXM_nmlw_cuda____);
				cudaFree(svd_SVUXM_Smlw_cuda____);
				cudaFree(svd_SVUXM_0inout_wSml_cuda____);
				cudaFree(svd_SVUXM_lwSm_cuda____);
				cudaFree(svd_USESVUXM_dwSm_cuda____);
				cudaFree(l2_dSm_cuda___);
				cudaFree(n2_dSm_cuda___);
				cudaFree(f2_dSm_cuda___);
				cudaFree(X_sub_dwSm_cuda___);
				n_Mbatch++;
				continue;
			}
			else
			{
				printf("Memory allocated on the GPU for I_value_use_dwSm_cuda___.\n");
			}
		}
		flag_memory_allocated = 1;
		break;
	}
	if ((n_Mbatch > n_M) || (!flag_memory_allocated))
	{
		printf("Not enough memory on the GPU to allocate the arrays.\n");
		printf("n_Mbatch %llu\n", n_Mbatch);
		return 1;
	}

	// initialize the arrays
	cudaMemset(svd_VUXM_lwnm_cuda____, 0, n_lwn * n_M_per_Mbatch * sizeof(cuFloatComplex));
	cudaMemset(UX_M_l2_dm_cuda__, 0, FTK_n_delta_v * n_M_per_Mbatch * sizeof(float));
	cudaMemset(X_wSm_cuda___, 0, n_wS * n_M_per_Mbatch * sizeof(float));
	cudaMemset(delta_x_wSm_cuda___, 0, n_wS * n_M_per_Mbatch * sizeof(float));
	cudaMemset(delta_y_wSm_cuda___, 0, n_wS * n_M_per_Mbatch * sizeof(float));
	cudaMemset(gamma_z_wSm_cuda___, 0, n_wS * n_M_per_Mbatch * sizeof(float));
	if (flag_compute_I_value)
	{
		cudaMemset(I_value_wSm_cuda___, 0, n_wS * n_M_per_Mbatch * sizeof(float));
	}
	cudaMemset(svd_VUXM_nmlw_cuda____, 0, n_lwn * n_M_per_Mbatch * sizeof(cuFloatComplex));
	cudaMemset(svd_SVUXM_Smlw_cuda____, 0, n_lwS * n_M_per_Mbatch * sizeof(cuFloatComplex));
	cudaMemset(svd_SVUXM_0inout_wSml_cuda____, 0, n_lwS * n_M_per_Mbatch * sizeof(cuFloatComplex));
	cudaMemset(svd_SVUXM_lwSm_cuda____, 0, n_lwS * n_M_per_Mbatch * sizeof(cuFloatComplex));
	cudaMemset(svd_USESVUXM_dwSm_cuda____, 0, n_dwS * n_M_per_Mbatch * sizeof(cuFloatComplex));
	cudaMemset(l2_dSm_cuda___, 0, n_dS * n_M_per_Mbatch * sizeof(float));
	cudaMemset(n2_dSm_cuda___, 0, n_dS * n_M_per_Mbatch * sizeof(float));
	cudaMemset(f2_dSm_cuda___, 0, n_dS * n_M_per_Mbatch * sizeof(float));
	cudaMemset(X_sub_dwSm_cuda___, 0, n_dwS * n_M_per_Mbatch * sizeof(float));
	if (flag_compute_I_value)
	{
		cudaMemset(I_value_use_dwSm_cuda___, 0, n_dwS * sizeof(float));
	}

	/* create a plan for the batched FFT */
	printf("Creating a plan for the batched FFT\n");
	unsigned long long n_lSm = FTK_n_svd_l * n_S * n_M_per_Mbatch;
	cufftHandle handlecufftPlan;
	cufftResult resultcufftPlan;
	resultcufftPlan = cufftPlan1d(&handlecufftPlan, n_w_max, CUFFT_C2C, n_lSm);

	if (resultcufftPlan != CUFFT_SUCCESS)
	{
		printf("cufftMakePlan1d failed!\n");
	}

	if (verbose)
	{
		printf("n_Mbatch %llu\n", n_Mbatch);
		printf("n_M_per_Mbatch %llu\n", n_M_per_Mbatch);
		printf("memory_size %llu out of %llu\n", memory_size, free_byte);
	}

	unsigned long long n_processed_M = 0;
	// for (unsigned long long iMbatch = 0; iMbatch < n_Mbatch; iMbatch++)
	unsigned long long iMbatch = 0;
	while (n_processed_M < n_M)
	{

		/* batch M into ampmh program */
		// unsigned long long n_M_batch = (iMbatch == n_Mbatch - 1) ? n_M_per_Mbatch : n_M - n_M_per_Mbatch * iMbatch;
		unsigned long long n_M_batch = n_M_per_Mbatch;
		if (n_M_batch > n_M - n_processed_M)
		{
			n_M_batch = n_M - n_processed_M;
		}
		if (n_M_batch < 1)
			continue;
		printf("Processing batch %llu with size %llu starting %llu/%llu \n", iMbatch, n_M_batch, n_processed_M, n_M);

		/* copy to gpu */
		// cudaMemset(svd_VUXM_lwnm_cuda____, 0, n_lwn * n_M_per_Mbatch * sizeof(cuFloatComplex));
		// cudaMemset(UX_M_l2_dm_cuda__, 0, FTK_n_delta_v * n_M_per_Mbatch * sizeof(float));

		// printf("svd_VUXM_lwnM____ + n_lwn * n_processed_M = %.16f %.16f\n", svd_VUXM_lwnM____[n_lwn * n_processed_M].x, svd_VUXM_lwnM____[n_lwn * n_processed_M].y);
		cuda_err = cudaMemcpy(svd_VUXM_lwnm_cuda____, svd_VUXM_lwnM____ + n_lwn * n_processed_M, n_lwn * n_M_batch * sizeof(cuFloatComplex), cudaMemcpyHostToDevice);
		if (cuda_err != cudaSuccess)
		{
			printf("Failed to copy svd_VUXM_lwnm_cuda____ to the GPU, %s.\n", cudaGetErrorString(cuda_err));
			return 1;
		}
		else
		{
			printf("svd_VUXM_lwnm_cuda____ of size %llu copied. \n", n_lwn * n_M_batch);
		}
		cuda_err = cudaMemcpy(UX_M_l2_dm_cuda__, UX_M_l2_dM__ + FTK_n_delta_v * n_processed_M, FTK_n_delta_v * n_M_batch * sizeof(float), cudaMemcpyHostToDevice);
		if (cuda_err != cudaSuccess)
		{
			printf("Failed to copy UX_M_l2_dm_cuda__ to the GPU.\n");
			return 1;
		}
		else
		{
			/* check if first elements equal */
			// printf("UX_M_l2_dM__ + FTK_n_delta_v * n_processed_M = %.16f\n", UX_M_l2_dM__[0]);
			// float test = 0.0;
			// cudaMemcpy(&test, UX_M_l2_dm_cuda__, sizeof(float), cudaMemcpyDeviceToHost);
			// printf("UX_M_l2_dm_cuda__ + FTK_n_delta_v * n_processed_M = %.16f\n", test);
			printf("UX_M_l2_dm_cuda__ of size %llu copied. \n", FTK_n_delta_v * n_M_batch);
		}

		ampmh_cuda(
			handle, handlecufftPlan,
			flag_optimize_over_gamma_z, flag_compute_I_value, tolerance_master,
			FTK_n_svd_l, FTK_n_delta_v, n_w_max, pm_n_UX_rank, n_S, n_M_batch,
			CTF_UX_S_l2_cuda_, svd_VUXM_lwnm_cuda____, UX_M_l2_dm_cuda__,
			FTK_delta_x_cuda_, FTK_delta_y_cuda_, gamma_z_cuda_,
			X_wSm_cuda___, delta_x_wSm_cuda___, delta_y_wSm_cuda___,
			gamma_z_wSm_cuda___, I_value_wSm_cuda___,
			FTK_svd_U_d_expiw_s_dl_cuda__, CTF_UX_S_k_q_nSw_cuda___,
			svd_VUXM_nmlw_cuda____, svd_SVUXM_Smlw_cuda____, svd_SVUXM_0inout_wSml_cuda____,
			svd_SVUXM_lwSm_cuda____, svd_USESVUXM_dwSm_cuda____,
			l2_dSm_cuda___, n2_dSm_cuda___, f2_dSm_cuda___, X_sub_dwSm_cuda___, I_value_use_dwSm_cuda___,
			verbose);

		/* copy back to host */
		cuda_err = cudaMemcpy(X_wSM___ + n_wS * n_processed_M, X_wSm_cuda___, n_wS * n_M_batch * sizeof(float), cudaMemcpyDeviceToHost);
		if (cuda_err != cudaSuccess)
		{
			printf("Failed to copy X_wSm_cuda___ back to the Host. %s.\n", cudaGetErrorString(cuda_err));
			return 1;
		}
		else
		{
			printf("X_wSM___ of size %llu copied. \n", n_wS * n_M_batch);
		}
		cuda_err = cudaMemcpy(delta_x_wSM___ + n_wS * n_processed_M, delta_x_wSm_cuda___, n_wS * n_M_batch * sizeof(float), cudaMemcpyDeviceToHost);
		if (cuda_err != cudaSuccess)
		{
			printf("Failed to copy delta_x_wSm_cuda___ back to the Host. %s.\n", cudaGetErrorString(cuda_err));
			return 1;
		}
		else
		{
			printf("delta_x_wSM___ of size %llu copied. \n", n_wS * n_M_batch);
		}
		cuda_err = cudaMemcpy(delta_y_wSM___ + n_wS * n_processed_M, delta_y_wSm_cuda___, n_wS * n_M_batch * sizeof(float), cudaMemcpyDeviceToHost);
		if (cuda_err != cudaSuccess)
		{
			printf("Failed to copy delta_y_wSm_cuda___ back to the Host. %s.\n", cudaGetErrorString(cuda_err));
			return 1;
		}
		else
		{
			printf("delta_y_wSM___ of size %llu copied. \n", n_wS * n_M_batch);
		}
		cuda_err = cudaMemcpy(gamma_z_wSM___ + n_wS * n_processed_M, gamma_z_wSm_cuda___, n_wS * n_M_batch * sizeof(float), cudaMemcpyDeviceToHost);
		if (cuda_err != cudaSuccess)
		{
			printf("Failed to copy gamma_z_wSm_cuda___ back to the Host. %s.\n", cudaGetErrorString(cuda_err));
			return 1;
		}
		else
		{
			printf("gamma_z_wSM___ of size %llu copied. \n", n_wS * n_M_batch);
		}
		if (flag_compute_I_value)
		{
			cuda_err = cudaMemcpy(I_value_wSM___ + n_wS * n_processed_M, I_value_wSm_cuda___, n_wS * n_M_batch * sizeof(float), cudaMemcpyDeviceToHost);
			if (cuda_err != cudaSuccess)
			{
				printf("Failed to copy I_value_wSm_cuda___ back to the Host. %s.\n", cudaGetErrorString(cuda_err));
				return 1;
			}
			else
			{
				printf("I_value_wSM___ of size %llu copied. \n", n_wS * n_M_batch);
			}
		}
		n_processed_M += n_M_batch;
		iMbatch++;
	} /* for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){ } */
	time(&t_end);
	printf(" %% time to run mex_ampmh_X_wSM___16_omp_helper: %0.2f seconds\n", difftime(t_end, t_start));
	/* %%%%%%%%%%%%%%%% */

	/* free */

	status = cublasDestroy(handle);
	if (status != CUBLAS_STATUS_SUCCESS)
	{
		printf("cublasDestroy failed\n");
	}
	resultcufftPlan = cufftDestroy(handlecufftPlan);
	if (resultcufftPlan != CUFFT_SUCCESS)
	{
		printf("cufftDestroy failed\n");
	}

	cudaFree(gamma_z_cuda_);
	cudaFree(FTK_svd_U_d_expiw_s_dl_cuda__);
	cudaFree(CTF_UX_S_k_q_nSw_cuda__);
	cudaFree(FTK_delta_x_cuda_);
	cudaFree(FTK_delta_y_cuda_);
	cudaFree(CTF_UX_S_k_q_nSw_cuda__);
	cudaFree(svd_VUXM_lwnm_cuda____);
	cudaFree(UX_M_l2_dm_cuda__);
	cudaFree(X_wSm_cuda___);
	cudaFree(delta_x_wSm_cuda___);
	cudaFree(delta_y_wSm_cuda___);
	cudaFree(gamma_z_wSm_cuda___);
	if (flag_compute_I_value)
	{
		cudaFree(I_value_wSm_cuda___);
	}
	cudaFree(svd_VUXM_nmlw_cuda____);
	cudaFree(svd_SVUXM_Smlw_cuda____);
	cudaFree(svd_SVUXM_0inout_wSml_cuda____);
	cudaFree(svd_SVUXM_lwSm_cuda____);
	cudaFree(svd_USESVUXM_dwSm_cuda____);

	if (verbose)
	{
		printf(" %% [finished mex_ampmh_X_wSM___16_omp]\n");
	}
	return 0;
}

int main(int argc, char **argv)
{
	int verbose = 1;
	int flag_test = 0;

	// const unsigned long long n_S_per_Sbatch = 24;
	const unsigned long long flag_optimize_over_gamma_z = 0;
	const unsigned long long flag_compute_I_value = 1;
	const double tolerance_master = 0.01;
	const unsigned long long FTK_n_svd_l = 10;
	const unsigned long long FTK_n_delta_v = 64;//29;
	const unsigned long long n_w_max = 98;
	const unsigned long long pm_n_UX_rank = 48;//18;
	const unsigned long long n_S = 993;
	const unsigned long long n_M = 1024;
	// const unsigned long long NUM_BATCH_M = 24;
	// const unsigned long long n_M_per_Mbatch = n_M % NUM_BATCH_M ? n_M / NUM_BATCH_M + 1 : n_M / NUM_BATCH_M;

	const unsigned long long n_dl = FTK_n_delta_v * FTK_n_svd_l;
	const unsigned long long n_wnS = n_w_max * pm_n_UX_rank * n_S;
	const unsigned long long n_lwnM = FTK_n_svd_l * n_w_max * pm_n_UX_rank * n_M;
	const unsigned long long n_dM = FTK_n_delta_v * n_M;
	const unsigned long long n_wSM = flag_optimize_over_gamma_z ? n_S * n_M : n_w_max * n_S * n_M;

	cuFloatComplex *FTK_svd_U_d_expiw_s_dl__ = (cuFloatComplex *)malloc(sizeof(cuFloatComplex) * n_dl);
	float *FTK_delta_x_ = (float *)malloc(sizeof(float) * FTK_n_delta_v);
	float *FTK_delta_y_ = (float *)malloc(sizeof(float) * FTK_n_delta_v);
	cuFloatComplex *CTF_UX_S_k_q_wnS___ = (cuFloatComplex *)malloc(sizeof(cuFloatComplex) * n_wnS);
	float *CTF_UX_S_l2_ = (float *)malloc(sizeof(float) * n_S);
	cuFloatComplex *svd_VUXM_lwnM____ = (cuFloatComplex *)malloc(sizeof(cuFloatComplex) * n_lwnM);
	float *UX_M_l2_dM__ = (float *)malloc(sizeof(float) * n_dM);

	for (unsigned long long i = 0; i < n_dl; i++)
	{
		float f_real = (float)(i % 123) - 61;
		float f_imag = (float)(i % 125) - 62;
		// FTK_svd_U_d_expiw_s_dl__[i] = (float complex)(f_real, f_imag);
		FTK_svd_U_d_expiw_s_dl__[i] = make_cuFloatComplex(f_real, f_imag);
	}
	for (unsigned long long i = 0; i < FTK_n_delta_v; i++)
	{
		FTK_delta_x_[i] = (float)(i % 127) - 63;
		FTK_delta_y_[i] = (float)(i % 129) - 65;
	}
	for (unsigned long long i = 0; i < n_wnS; i++)
	{
		float f_real = (float)(i % 123) - 61;
		float f_imag = (float)(i % 125) - 62;
		// CTF_UX_S_k_q_wnS___[i] = (float complex)(f_real, f_imag);
		CTF_UX_S_k_q_wnS___[i] = make_cuFloatComplex(f_real, f_imag);
	}
	for (unsigned long long i = 0; i < n_S; i++)
	{
		CTF_UX_S_l2_[i] = (float)(1 + i % 123);
	}
	for (unsigned long long i = 0; i < n_lwnM; i++)
	{
		float f_real = (float)(i % 123) - 61;
		float f_imag = (float)(i % 125) - 62;
		// svd_VUXM_lwnM____[i] = (float complex)(f_real, f_imag);
		svd_VUXM_lwnM____[i] = make_cuFloatComplex(f_real, f_imag);
	}
	for (unsigned long long i = 0; i < n_dM; i++)
	{
		UX_M_l2_dM__[i] = (float)(1 + i % 123);
	}

	/* output */
	float *X_wSM___ = (float *)malloc(sizeof(float) * n_wSM);
	float *delta_x_wSM___ = (float *)malloc(sizeof(float) * n_wSM);
	float *delta_y_wSM___ = (float *)malloc(sizeof(float) * n_wSM);
	float *gamma_z_wSM___ = (float *)malloc(sizeof(float) * n_wSM);
	float *I_value_wSM___ = (float *)malloc(sizeof(float) * n_wSM);

	time_t t_start, t_end;
	time(&t_start);
	ampmh_cuda_wrapper(flag_optimize_over_gamma_z, flag_compute_I_value, tolerance_master, FTK_n_svd_l, FTK_n_delta_v, FTK_svd_U_d_expiw_s_dl__, FTK_delta_x_, FTK_delta_y_, n_w_max, pm_n_UX_rank, n_S, CTF_UX_S_k_q_wnS___, CTF_UX_S_l2_, n_M, svd_VUXM_lwnM____, UX_M_l2_dM__, X_wSM___, delta_x_wSM___, delta_y_wSM___, gamma_z_wSM___, I_value_wSM___);
	time(&t_end);
	printf(" %% ampmh_cuda_wrapper: took %d seconds.\n", (int)difftime(t_end, t_start));

	/* free */
	free(FTK_svd_U_d_expiw_s_dl__);
	free(FTK_delta_x_);
	free(FTK_delta_y_);
	free(CTF_UX_S_k_q_wnS___);
	free(CTF_UX_S_l2_);
	free(svd_VUXM_lwnM____);
	free(UX_M_l2_dM__);

	free(X_wSM___);
	free(delta_x_wSM___);
	free(delta_y_wSM___);
	free(gamma_z_wSM___);
	free(I_value_wSM___);

	return 0;
}