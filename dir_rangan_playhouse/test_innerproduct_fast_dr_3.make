# Barnett adjusted for global include flag options. 6/24/16

# get the global options (barnett 6/23/16)
include make.inc.dist
-include make.inc

PROJECT=int2         # for historical reasons we always like to 
                     # call the executable int2, but you could set it to 
                     # something more descriptive
OBJSUF=o
MODSUF=mod
FLINK=$(FC) -o $(PROJECT)

.f.$(OBJSUF):
	$(FC) -c $(FFLAGS) $<

.f.$(MODSUF):
	$(FC) -c $(FFLAGS) $<

.SUFFIXES: $(MODSUF) .$(OBJSUF) .f .c

vpath %.f .:../spharmproj

FMODS = 

FSRCS = test_innerproduct_fast_dr_3.f  \
	adi_rand_f.f \
	get_template_size.f \
	get_F2_x_c.f \
	get_F2_x_c_.f \
	get_ctf_star_k_p_.f \
	get_ctf_ones_k_p_.f \
	stdlim_c16.f \
	pearson_c16.f \
	cl1_i4.f \
	cl1_c16.f \
	cp1_r8.f \
	cp1_c16.f \
	cps_c16.f \
	af1_c16.f \
	afs_c16.f \
	xx1_c16.f \
	xc1_c16.f \
	trn0_r8.f \
	trn0_c16.f \
	write_xxx_xx.f \
	randperm.f \
	adi_randperm.f \
	randinclude.f \
	linspace.f \
	interp1_c16.f \
	interp2_c16.f \
	interp_c_to_p.f \
	interp_p_to_q.f \
	interp_p_to_q_fftw.f \
	periodize_r8.f \
	periodize_i.f \
	transl_c_to_c.f \
	transf_c_to_c.f \
	transf_p_to_p.f \
	transf_svd_q_to_q.f \
	rotate_c_to_c.f \
	rotate_p_to_p.f \
	rotate_q_to_q.f \
	innerproduct_c.f \
	innerproduct_p.f \
	innerproduct_q_stretch.f \
	innerproduct_q__k_only.f \
	innerproduct_q__k_svdr.f \
	innerproduct_q__k_svdd.f \
	get_interchange_delta.f \
	get_C_S_use_.f \
	get_CTF_R_S_use_.f \
	get_delta_0.f \
	get_gamma_0.f \
	get_svd_0.f \
	dfftw_block_many_0.f \
	dfftw_block_many_1.f \
	test_innerproduct_fast_wrapper_0.f \
	test_innerproduct_fast_wrapper_1.f \
	test_innerproduct_fast_wrapper_2.f \
	test_innerproduct_fast_stage_0.f \
	test_innerproduct_fast_omp_stage_0.f \
	test_innerproduct_fast_stage_1.f \
	test_innerproduct_fast_stage_2.f \
	test_innerproduct_fast_MS_0.f \
	test_innerproduct_fast_MS_sort_0.f \
	test_innerproduct_fast_MS_update_0.f \
	test_innerproduct_fast_MS_update_1.f \
	test_innerproduct_batch_SM_stage_1a.f \
	test_innerproduct_batch_SM_stage_2a.f \
	test_innerproduct_fast_CTF_R_S_0.f \
	test_innerproduct_fast_CTF_R_S_1.f \
	test_innerproduct_fast_T_S_q_0.f \
	test_innerproduct_fast_T_S_q_1.f \
	test_innerproduct_fast_T_S_q_2.f \
	test_innerproduct_fast_Z_S_q_0.f \
	test_innerproduct_fast_Z_S_q_1.f \
	test_innerproduct_fast_Z_S_q_2.f \
	test_innerproduct_fast_T_R_CTF_M_q_0.f \
	test_innerproduct_fast_T_R_CTF_M_q_1.f \
	test_innerproduct_fast_T_R_CTF_M_q_2.f \
	test_innerproduct_fast_STxTRM_0.f \
	test_innerproduct_fast_STxTRM_1.f \
	test_innerproduct_fast_SZxTRM_0.f \
	test_innerproduct_fast_SZxTRM_1.f \
	test_innerproduct_fast_ZZ_0.f \
	test_innerproduct_fast_ZZ_1.f \
	test_innerproduct_fast_ZZ_2.f \
	test_innerproduct_fast_ZZ_scan_0.f \
	test_innerproduct_batch_excerpt_1.f \
	test_innerproduct_batch_excerpt_2.f \
	test_innerproduct_batch_excerpt_3.f \
	test_innerproduct_batch_excerpt_4.f \
	test_innerproduct_batch_excerpt_5.f \
	test_innerproduct_batch_excerpt_6.f \
	test_innerproduct_batch_excerpt_7.f \
	test_innerproduct_batch_excerpt_8.f \
	test_innerproduct_batch_sort_a.f \
	max_i4_f.f \
	max_r8_f.f \
	avg_r8_f.f \
	al2_r8_f.f \
	sum_i4_f.f \
	cp1_i4.f \
	unique_i4.f \
	cross_r8.f \
	dot_r8.f \
	dot_c16.f \
	distance_r8.f \
	normalize_r8.f \
	normalize_c16.f \
	gramschmidt_c16.f \
	get_angle_to_vp_.f \
	tesselation_size.f \
	tesselation_index.f \
	tesselation_neighborhood.f \
	tesselation_define_octants.f \
	tesselation_get_nl_nm_ll.f \
	tesselation_index_wrapper_0.f \
	tesselation_neighborhood_wrapper_0.f \
	hsv2rgb.f \
	colorscale.f \
	Fig_header.f \
	Fig_text.f \
	recenter_c16.f \
	Fig_c16_carte.f \
	Fig_c16_polar.f \
	Fig_c16_bessl.f \
	Fig_c16_bess2.f \
	Fig_gen_ver0.f \
	Fig_gen_ver2.f \
	Fig_gen_ver3.f \
	Fig_gen_ver4.f \
	Fig_gen_ctf_ver1.f \
	Fig_gen_ctf_ver2.f \
	adi_fft1.f \
	adi_fft2.f \
	quicksort_c16.f \
	quicksort_r8.f \
	quicksort_i4.f \
	source_end.f 

#
# object files list
MODS    =  $(FMODS:.f=.$(MODSUF)) 
OBJS    =  $(FMODS:.f=.$(OBJSUF)) $(FSRCS:.f=.$(OBJSUF)) 
#
#
#  needs finufft.a library in working directory
#  eventually should be fixed with proper intallation/makefile
#
#
$(PROJECT):   $(MODS)   $(OBJS)
	rm -f $(PROJECT)
	$(FLINK) $(OBJS) $(LFLAGS) -L/home/rangan/finufft-master/lib/ -lfinufft -lfftw3 -lfftw3_threads -fopenmp -lopenblas -lstdc++
	./$(PROJECT)
#
clean: 
	rm -f $(OBJS)
# 
list: $(FSRCS)
	echo $^
#
pack: $(FSRCS)
	cat $^ > _tmp_.pack.f
#
