// Written by Aaditya Rangan, with thanks to Chris Chang (see https://www.cog-genomics.org/plink2/). ;
// This program is free software: ;
// you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
// See the GNU General Public License for more details.
// You should have received a copy of the GNU General Public License along with this program.  
// If not, see <http://www.gnu.org/licenses/>.

#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */
#ifdef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

/* global variables used for timing */
clock_t GLOBAL_t_start[GLOBAL_NTICKS],GLOBAL_t_final[GLOBAL_NTICKS];
struct timeval GLOBAL_d_start[GLOBAL_NTICKS],GLOBAL_d_final[GLOBAL_NTICKS];
long GLOBAL_l_msec[GLOBAL_NTICKS],GLOBAL_l_ssec[GLOBAL_NTICKS],GLOBAL_l_usec[GLOBAL_NTICKS]; 
double GLOBAL_elct[GLOBAL_NTICKS],GLOBAL_elrt[GLOBAL_NTICKS];

char GLOBAL_CWD[FNAMESIZE];
char GLOBAL_dir_base[FNAMESIZE]="\0";
char GLOBAL_dir_name[FNAMESIZE]="\0";
char GLOBAL_dir_xpre[FNAMESIZE]="\0";
char GLOBAL_out_name[FNAMESIZE]="\0";
int GLOBAL_verbose=0; // set to 1 to see sysconf(_SC_NPROCESSORS_ONLN) ;
int GLOBAL_thread_count=1; // Set >1 to use pthreads to parallelize across nbins ;
int GLOBAL_omp_type=GLOBAL_omp_off; // manual parallelization using pthreads ;
int GLOBAL_1_icache_linesize=64;
int GLOBAL_1_icache_size=1024;
int GLOBAL_1_icache_assoc=16;
int GLOBAL_1_dcache_linesize=64;
int GLOBAL_1_dcache_size=1024;
int GLOBAL_1_dcache_assoc=16;
int GLOBAL_2_cache_linesize=64;
int GLOBAL_2_cache_size=1024;
int GLOBAL_2_cache_assoc=16;
int GLOBAL_3_cache_linesize=64;
int GLOBAL_3_cache_size=1024;
int GLOBAL_3_cache_assoc=16;
int GLOBAL_4_cache_linesize=64;
int GLOBAL_4_cache_size=1024;
int GLOBAL_4_cache_assoc=16;
double GLOBAL_tolerance=0.000000000001;
unsigned int GLOBAL_recursion_limit=1024*32;
int addressable_1=1;
int addressable_0=0;
int addressable_int_length=128;
int addressable_int[128] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127};

/* thread management */
int GLOBAL_nf=0; int GLOBAL_nf_cur=0; int GLOBAL_nf_ind=0; int GLOBAL_nf_opn=0;
int GLOBAL_tint[MAX_THREADS];
void *GLOBAL_tvp[MAX_THREADS][128];
pthread_t GLOBAL_threads[MAX_THREADS];
unsigned long long int GLOBAL_ops_f_[MAX_THREADS];
unsigned long long int GLOBAL_ops_f_sum=0;
unsigned long long int GLOBAL_ops_b_[MAX_THREADS];
unsigned long long int GLOBAL_ops_b_sum=0;

/* ---------------------------------------------------------------- */
unsigned long long int GLOBAL_n_malloc1=0;
int GLOBAL_n_malloc1_[GLOBAL_NTICKS];
/* manually managed memory stack */
unsigned char* wkspace=NULL;
unsigned char* wkspace_base=NULL;
int GLOBAL_wkspace_point=1;
struct wkspace_point *wkspace_point_0=NULL,*wkspace_point_t=NULL;
unsigned long long int GLOBAL_memory_gb=GLOBAL_MEMORY_GB_DEFAULT;
unsigned long long int GLOBAL_memory_mb=GLOBAL_MEMORY_GB_DEFAULT*(unsigned long long int)1024;
unsigned long long int GLOBAL_memory_kb=GLOBAL_MEMORY_GB_DEFAULT*(unsigned long long int)1048576;
long long int wkspace_left=0;
long long int wkspace_used=0;

/* ---------------------------------------------------------------- */

/* RAND functions */
unsigned long int POW2RPOWPLUSRADD=35L;
unsigned long int POW22RPOWMINUSONE=2147483647LL;
int RCYCLENUM=7;

/*---------------------------------------------------------------- */

#ifdef _MONOLITH
#include "playpark_function.c"
#endif /* _MONOLITH */

inline void ping(){ printf(" %% ping\n");}
inline void pong(){ printf(" %% pong\n");}

void set_globals()
{
  int verbose=0;
  int np=0,nb=0;
  struct stat st_tmp = {0};
  if (verbose){ printf(" %% setting GLOBAL_dir_name\n");}
  if (!strcmp(GLOBAL_out_name,"\0")){ sprintf(GLOBAL_out_name,"test"); /* if (!strcmp(GLOBAL_out_name,"\0")){ } */}
  if (strcmp(GLOBAL_dir_name,"\0")){ /* do nothing */ }
  else /* if GLOBAL_dir_name not defined */{
    if (strcmp(GLOBAL_out_name,"\0")){
      if (!strcmp(GLOBAL_dir_base,"\0")){ sprintf(GLOBAL_dir_name,"%s/dir_%s",GLOBAL_CWD,GLOBAL_out_name);}
      else /* if GLOBAL_dir_base */{ sprintf(GLOBAL_dir_name,"%s/dir_%s",GLOBAL_dir_base,GLOBAL_out_name);}
      /* if (strcmp(GLOBAL_out_name,"\0"))){ } */}
    /* if not defined */}
  if (verbose>-1){ printf(" %% setting GLOBAL_dir_name: %s\n",GLOBAL_dir_name);}
  if (strcmp(GLOBAL_dir_name,"\0") && stat(GLOBAL_dir_name, &st_tmp) == -1) {  printf(" %% mkdir %s;\n",GLOBAL_dir_name); mkdir(GLOBAL_dir_name, 0755);}
  if (verbose){ printf(" %% initializing GLOBAL_tint\n");}
  for (GLOBAL_nf=0;GLOBAL_nf<MAX_THREADS;GLOBAL_nf++){ GLOBAL_tint[GLOBAL_nf]=GLOBAL_nf;}
}

int main(int argc, char** argv) {
  unsigned char* wkspace_ua;
  int cur_arg = 1,ii=0;
  char* argptr;
  char* bubble;
  int retval=0;
  int cmdline_param=0;
  if (omp_get_max_threads()>MAX_THREADS){ printf(" %% Warning! global variable MAX_THREADS set to %d when it should be %d\n",MAX_THREADS,omp_get_max_threads()); exit(0);}
  if (!getcwd(GLOBAL_CWD,FNAMESIZE)){ printf(" %% Warning! GLOBAL_CWD not read in main\n");} else{ if (GLOBAL_verbose>0){ printf(" %% GLOBAL_CWD: %s\n",GLOBAL_CWD);};} 
  GLOBAL_thread_count = sysconf(_SC_NPROCESSORS_ONLN); if (GLOBAL_verbose>-1){ printf("sysconf(_SC_NPROCESSORS_ONLN) = %d\n",GLOBAL_thread_count);}
  GLOBAL_1_icache_linesize = sysconf(_SC_LEVEL1_ICACHE_LINESIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL1_ICACHE_LINESIZE) = %d\n",GLOBAL_1_icache_linesize);}
  GLOBAL_1_icache_size = sysconf(_SC_LEVEL1_ICACHE_SIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL1_ICACHE_SIZE) = %d\n",GLOBAL_1_icache_size);}
  GLOBAL_1_icache_assoc = sysconf(_SC_LEVEL1_ICACHE_ASSOC); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL1_ICACHE_ASSOC) = %d\n",GLOBAL_1_icache_assoc);}
  GLOBAL_1_dcache_linesize = sysconf(_SC_LEVEL1_DCACHE_LINESIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL1_DCACHE_LINESIZE) = %d\n",GLOBAL_1_dcache_linesize);}
  GLOBAL_1_dcache_size = sysconf(_SC_LEVEL1_DCACHE_SIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL1_DCACHE_SIZE) = %d\n",GLOBAL_1_dcache_size);}
  GLOBAL_1_dcache_assoc = sysconf(_SC_LEVEL1_DCACHE_ASSOC); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL1_DCACHE_ASSOC) = %d\n",GLOBAL_1_dcache_assoc);}
  GLOBAL_2_cache_linesize = sysconf(_SC_LEVEL2_CACHE_LINESIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL2_CACHE_LINESIZE) = %d\n",GLOBAL_2_cache_linesize);}
  GLOBAL_2_cache_size = sysconf(_SC_LEVEL2_CACHE_SIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL2_CACHE_SIZE) = %d\n",GLOBAL_2_cache_size);}
  GLOBAL_2_cache_assoc = sysconf(_SC_LEVEL2_CACHE_ASSOC); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL2_CACHE_ASSOC) = %d\n",GLOBAL_2_cache_assoc);}
  GLOBAL_3_cache_linesize = sysconf(_SC_LEVEL3_CACHE_LINESIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL3_CACHE_LINESIZE) = %d\n",GLOBAL_3_cache_linesize);}
  GLOBAL_3_cache_size = sysconf(_SC_LEVEL3_CACHE_SIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL3_CACHE_SIZE) = %d\n",GLOBAL_3_cache_size);}
  GLOBAL_3_cache_assoc = sysconf(_SC_LEVEL3_CACHE_ASSOC); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL3_CACHE_ASSOC) = %d\n",GLOBAL_3_cache_assoc);}
  GLOBAL_4_cache_linesize = sysconf(_SC_LEVEL4_CACHE_LINESIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL4_CACHE_LINESIZE) = %d\n",GLOBAL_4_cache_linesize);}
  GLOBAL_4_cache_size = sysconf(_SC_LEVEL4_CACHE_SIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL4_CACHE_SIZE) = %d\n",GLOBAL_4_cache_size);}
  GLOBAL_4_cache_assoc = sysconf(_SC_LEVEL4_CACHE_ASSOC); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL4_CACHE_ASSOC) = %d\n",GLOBAL_4_cache_assoc);}
  if (GLOBAL_thread_count == -1) { GLOBAL_thread_count = 1; } else if (GLOBAL_thread_count > MAX_THREADS) { GLOBAL_thread_count = MAX_THREADS; }
  while (cur_arg < argc) {
    argptr = argv[cur_arg];
    if (0){ /* do nothing */}
    else if (!strcmp(argptr, "--GLOBAL_memory_gb")) {
      if (cur_arg == argc - 1) { printf("Error: Missing --GLOBAL_memory_gb parameter."); exit(RET_INVALID_CMDLINE);}
      ii = atoi(argv[cur_arg + 1]);
      if (ii < 1) { printf("Error: Invalid --GLOBAL_memory_gb parameter.\n"); exit(RET_INVALID_CMDLINE);}
      GLOBAL_memory_gb = ii; printf(" %% GLOBAL_memory_gb %d\n",GLOBAL_memory_gb);
      cur_arg += 2; /* else if GLOBAL_memory_gb */}
    else if (!strcmp(argptr, "--nthreads")) {
      if (cur_arg == argc - 1) { printf("Error: Missing --nthreads parameter."); exit(RET_INVALID_CMDLINE);}
      ii = atoi(argv[cur_arg + 1]);
      if (ii < 1) { printf("Error: Invalid --nthreads parameter.\n"); exit(RET_INVALID_CMDLINE);}
      GLOBAL_thread_count = ii; printf(" %% nthreads %d\n",GLOBAL_thread_count);
      cur_arg += 2; /* else if nthreads */}
    else if (!strcmp(argptr, "--cmdline_param")) {
      if (cur_arg == argc - 1) { printf("Error: Missing --cmdline_param parameter."); exit(RET_INVALID_CMDLINE);}
      ii = atoi(argv[cur_arg + 1]);
      if (ii < 1) { printf("Error: Invalid parameter.\n"); exit(RET_INVALID_CMDLINE);}
      cmdline_param = ii; printf(" %% cmdline_param %d\n",cmdline_param);
      cur_arg += 2; /* else if other */}
    else { printf("Error: Invalid argument (%s).", argv[cur_arg]); exit(RET_INVALID_CMDLINE);}
    /* while (cur_arg < argc) { } */}
  bubble = (char*)malloc1((unsigned long long int)67108864 * sizeof(char));
  if (!bubble) { printf("Error: not enough memory!\n");}
  GLOBAL_memory_mb = GLOBAL_memory_gb * (unsigned long long int)1024; GLOBAL_memory_kb = GLOBAL_memory_mb * (unsigned long long int)1024;
  wkspace_ua = (unsigned char*)malloc1(GLOBAL_memory_kb * (unsigned long long int)1024 * sizeof(char));
  while (!wkspace_ua) { 
    if (GLOBAL_verbose>-2){ if (!wkspace_ua){ printf("Could not allocate %ld GB; trying %ld GB instead...\n",GLOBAL_memory_gb,GLOBAL_memory_gb - 1);}}
    GLOBAL_memory_gb -= 1; GLOBAL_memory_mb = GLOBAL_memory_gb * (unsigned long long int)1024; GLOBAL_memory_kb = GLOBAL_memory_mb * (unsigned long long int)1024;
    wkspace_ua = (unsigned char*)malloc1(GLOBAL_memory_kb * (unsigned long long int)1024 * sizeof(char));
    if (wkspace_ua) { printf("Allocated %ld KB = %ld MB = %ld GB successfully.\n", GLOBAL_memory_kb,GLOBAL_memory_mb,GLOBAL_memory_gb);}
    /* while (!wkspace_ua) { } */}
  // force 64-byte align on OS X to make cache line sensitivity work
  wkspace = (unsigned char*)CACHEALIGN((unsigned long)wkspace_ua);
  wkspace_base = wkspace;
  wkspace_left = GLOBAL_memory_kb * (unsigned long long int)1024 - (unsigned long long int)(wkspace - wkspace_ua);
  free1(bubble);
  if (GLOBAL_verbose>0){ wkspace_printf();}
  wkspace_point_0 = wkspace_make_point();
  wkspace_point_0->parent = NULL; wkspace_point_0->child = NULL;
  wkspace_point_0->check = 0; *(wkspace_point_0->point) = wkspace_point_0->check;
  wkspace_point_t = wkspace_point_0;
  GLOBAL_verbose=0;
  read_input();
  set_globals();
  //MDA_io_test();
  //R01GET_test();
  //fftw3f__1d_test();
  //finufft_1d_test();
  //dtranspose_test();
  //dp_ps_single_test();
  //dp_ps_mult_immintrin_test();
  //dp_pd_mult_immintrin_test();
  //get_xdrop_logscale_array_test();
  //iquicksort_index_driver_test();
  //fquicksort_index_driver_test();
  //dquicksort_index_driver_test();
  //iquicksort_index_index_driver_test();
  //fquicksort_index_index_driver_test();
  //dquicksort_index_index_driver_test();
  //irandperm_test();
  //dexcluster_nonbinary_test_error();
  //dexcluster_nonbinary_rdrop_test_error();
  //gumbel_nll_test();
  //nelder_mead_test();
  //gumbel_fit_test();
  //array_extract_i_from_i_test();
  //array_extract_f_from_d_test();
  //array_extract_d_from_d_test();
  //array_extract_f_from_f_test();
  //array_extract_d_from_f_test();
  //array_extract_test();
  //array_mean_center_row_test();
  //array_normalize_row_test();
  //array_orth_f_test_error();
  //array_orth_f_test_speed();
  //array_orth_d_test_error();
  //array_orth_d_test_speed();
  //erfcln_f_test();
  //erfcln_d_test();
  //z_to_lp_d_test();
  //find_internal_maximum_test();
  dexcluster_nonbinary_f_recursive_test();
  if (GLOBAL_verbose>-1){ printf("exiting successfully\n");}
  free1(&wkspace_ua);
  exit(RET_SUCCESS);
  return 0;
}
