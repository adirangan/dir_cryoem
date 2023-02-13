// Written by Aaditya Rangan. ;
// This program is free software: ;
// you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
// See the GNU General Public License for more details.
// You should have received a copy of the GNU General Public License along with this program.  
// If not, see <http://www.gnu.org/licenses/>.

#ifndef _MONOLITH
#include "ampm_header.h"
#endif /* _MONOLITH */
#ifdef _MONOLITH
#include "ampm_header.h"
#endif /* _MONOLITH */

int GLOBAL_verbose=0;
char GLOBAL_mode[FNAMESIZE] = "\0";
double GLOBAL_tolerance=0.000000000001; //%<-- 1e-12;
char GLOBAL_dir_trunk[PNAMESIZE] = "\0";
char GLOBAL_prefix_base[FNAMESIZE] = "\0";
int GLOBAL_flag_force_create=0;
int GLOBAL_flag_omp_use=1;

/* ---------------------------------------------------------------- */
/* global variables used for timing */
clock_t GLOBAL_t_start[GLOBAL_NTICKS],GLOBAL_t_final[GLOBAL_NTICKS];
struct timeval GLOBAL_d_start[GLOBAL_NTICKS],GLOBAL_d_final[GLOBAL_NTICKS];
long GLOBAL_l_msec[GLOBAL_NTICKS],GLOBAL_l_ssec[GLOBAL_NTICKS],GLOBAL_l_usec[GLOBAL_NTICKS]; 
double GLOBAL_elct[GLOBAL_NTICKS],GLOBAL_elrt[GLOBAL_NTICKS];
/* ---------------------------------------------------------------- */
/* global variables used for malloc1 */
int GLOBAL_malloc1_notupdate=0;
unsigned long long int GLOBAL_n_malloc1=0;
int GLOBAL_n_malloc1_[GLOBAL_NTICKS];
/* ---------------------------------------------------------------- */
/* RAND functions */
unsigned long int POW2RPOWPLUSRADD=35L;
unsigned long int POW22RPOWMINUSONE=2147483647LL;
/*---------------------------------------------------------------- */

#ifdef _MONOLITH
#include "ampm_function.c"
#endif /* _MONOLITH */

#ifdef _WITHOUT_MEX
int main(int argc, char** argv) {
  int cur_arg = 1,ii=0;
  char* argptr;
  int retval=0;
  int cmdline_param=0;
  int flag_test=0;
  while (cur_arg < argc) {
    argptr = argv[cur_arg];
    if (0){ /* do nothing */}
    else if (!strcmp(argptr, "--test")) {
      if (cur_arg == argc - 1) { printf("Error: Missing --test parameter."); exit(EXIT_FAILURE);}
      ii = atoi(argv[cur_arg + 1]);
      if (ii < 1) { printf("Error: Invalid --test parameter.\n"); exit(EXIT_FAILURE);}
      flag_test = ii; printf(" %% flag_test %d\n",flag_test);
      cur_arg += 2; /* else if flag_test */}
    else { printf("Error: Invalid argument (%s).", argv[cur_arg]); exit(EXIT_FAILURE);}
    /* while (cur_arg < argc) { } */}
  read_input();
  if (strstr(GLOBAL_mode,"MDA_io_test")!=NULL){ MDA_io_test();}
  if (strstr(GLOBAL_mode,"R01GET_test")!=NULL){ R01GET_test();}
  if (strstr(GLOBAL_mode,"dtranspose_test")!=NULL){ dtranspose_test();}
  if (strstr(GLOBAL_mode,"cblas_test")!=NULL){ cblas_test();}
  if (strstr(GLOBAL_mode,"transpose_float_test")!=NULL){ transpose_float_test();}
  if (strstr(GLOBAL_mode,"dp_ps_mult_immintrin_test")!=NULL){ dp_ps_mult_immintrin_test();}
  if (strstr(GLOBAL_mode,"hp_ps_mult_immintrin_test")!=NULL){ hp_ps_mult_immintrin_test();}
  if (strstr(GLOBAL_mode,"rcp_ps_mult_immintrin_test")!=NULL){ rcp_ps_mult_immintrin_test();}
  if (strstr(GLOBAL_mode,"array_extract_i_from_i_test")!=NULL){ array_extract_i_from_i_test();}
  if (strstr(GLOBAL_mode,"array_extract_f_from_d_test")!=NULL){ array_extract_f_from_d_test();}
  if (strstr(GLOBAL_mode,"array_extract_d_from_d_test")!=NULL){ array_extract_d_from_d_test();}
  if (strstr(GLOBAL_mode,"array_extract_f_from_f_test")!=NULL){ array_extract_f_from_f_test();}
  if (strstr(GLOBAL_mode,"array_extract_d_from_f_test")!=NULL){ array_extract_d_from_f_test();}
  if (strstr(GLOBAL_mode,"array_extract_test")!=NULL){ array_extract_test();}
  if (strstr(GLOBAL_mode,"gsl_legpts_test")!=NULL){ gsl_legpts_test();}
  if (strstr(GLOBAL_mode,"gsl_legendre_unnorm_test")!=NULL){ gsl_legendre_unnorm_test();}
  if (strstr(GLOBAL_mode,"gsl_legendre_evaluate_normalized_lxm___test")!=NULL){ gsl_legendre_evaluate_normalized_lxm___test();}
  if (strstr(GLOBAL_mode,"sample_shell_L_test")!=NULL){ sample_shell_L_test();}
  if (strstr(GLOBAL_mode,"pm_template_test")!=NULL){ pm_template_test();}
  if (strstr(GLOBAL_mode,"mex_ampmh_X_wSM_test")!=NULL){ mex_ampmh_X_wSM_test();}
  if (GLOBAL_verbose>0){ printf("exiting successfully\n");}
  return 0;
}
#endif /* _WITHOUT_MEX */
