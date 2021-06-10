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

FSRCS = test_innerproduct_7_dr.f  \
        prini.f \
        chebexps.f \
        gridrouts3d.f \
	ping.f \
	adi_rand_f.f \
	get_template_size.f \
	get_F2_x_c.f \
	get_F2_x_c_.f \
	get_ctf_star_k_p_.f \
	get_ctf_ones_k_p_.f \
	stdlim_c16.f \
	pearson_c16.f \
	cl1_l2.f \
	cl1_i4.f \
	cl1_r8.f \
	cl1_c16.f \
	cs1_l2.f \
	cs1_i4.f \
	cs1_r8.f \
	cs1_c16.f \
	cxs_l2.f \
	cxs_i4.f \
	cxs_r8.f \
	cxs_c16.f \
	cp1_l2.f \
	cp1_i4.f \
	cp1_r8.f \
	cp1_c16.f \
	cps_c16.f \
	af1_c16.f \
	afs_c16.f \
	xx1_c16.f \
	xc1_c16.f \
	trn0_r8.f \
	trn0_c16.f \
	trn0_1_c16.f \
	write_xxx_xx.f \
	MDA_write_i4.f \
	MDA_write_r8.f \
	MDA_write_c16.f \
	randperm.f \
	adi_randperm.f \
	randinclude.f \
	adi_randinclude.f \
	linspace.f \
	block_0.f \
	block_1.f \
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
	transf_svd_q_to_q_3.f \
	transf_svd_q_to_q_4.f \
	transf_svd_q_to_q_5.f \
	transf_svd_q_to_q_6.f \
	rotate_c_to_c.f \
	rotate_p_to_p.f \
	rotate_p_to_p_lagrange.f \
	rotate_p_to_p_spectral.f \
	rotate_q_to_q.f \
	rotate_q_to_q_0.f \
	innerproduct_c.f \
	innerproduct_p.f \
	innerproduct_q_stretch.f \
	innerproduct_q__k_only.f \
	innerproduct_q__k_svdr.f \
	innerproduct_q__k_svdr_3.f \
	innerproduct_q__k_svdr_4.f \
	innerproduct_q__k_svdr_5.f \
	innerproduct_q__k_svdr_6.f \
	innerproduct_q__k_svdd.f \
	innerproduct_q__k_svdd_3.f \
	innerproduct_q__k_svdd_4.f \
	get_interchange_delta.f \
	get_interchange_delta_RTRT_vs_RTTR.f \
	get_C_S_use_.f \
	get_CTF_R_S_use_.f \
	get_delta_0.f \
	get_delta_1.f \
	get_gamma_0.f \
	get_svd_2.f \
	get_svd_polyval_U_d_0.f \
	get_svd_polyval_V_r_0.f \
	get_CTF_R_S_periodize_0.f \
	dfftw_block_many_0.f \
	dfftw_block_many_1.f \
	alpha_SM_update_1.f \
	alpha_SM_sort_0.f \
	alpha_SM_write_0.f \
	test_innerproduct_bruteforce_0.f \
	test_innerproduct_bruteforce_1.f \
	test_innerproduct_fast_CTF_R_S_1.f \
	test_innerproduct_7.f \
	test_innerproduct_5_T_S_q_4.f \
	test_innerproduct_5_CTF_R_S_2.f \
	test_innerproduct_5_O_S_q_1.f \
	test_innerproduct_5_O_T_R_CTF_M_q_3.f \
	test_innerproduct_5_T_T_R_CTF_M_q_2.f \
	ti7_build_CTF_R_S_3.f \
	ti7_build_O_S_q_1.f \
	ti7_build_T_S_q_4.f \
	ti7_build_Z_S_q_5.f \
	test_innerproduct_7_global_0.f \
	test_innerproduct_7_local_0.f \
	ti7_build_O_T_R_CTF_M_q_3.f \
	ti7_build_T_T_R_CTF_M_q_2.f \
	ti7_build_Z_T_R_CTF_M_q_2.f \
	ti7_check_STxTRM_3.f \
	ti7_check_SZxTRM_4.f \
	ti7_check_SxTTRM_3.f \
	ti7_check_SxZTRM_4.f \
	test_innerproduct_5_Zstore_1.f \
	test_innerproduct_5_Zstore_3a.f \
	test_innerproduct_5_Zstore_3b.f \
	test_innerproduct_5_Zstore_3c.f \
	test_innerproduct_5_Zstore_3x.f \
	test_innerproduct_5_Zstore_3d.f \
	test_innerproduct_5_Zstore_3e.f \
	test_innerproduct_5_Zstore_3.f \
	test_innerproduct_5_Zstore_4.f \
	test_innerproduct_5_ZZ_scan_2.f \
	test_innerproduct_5_ZZ_scan_4x.f \
	test_innerproduct_5_ZZ_scan_4.f \
	test_alpha_update_SM_3.f \
	test_alpha_update_MS_3.f \
	test_max_0.f \
	min_i4_f.f \
	max_i4_f.f \
	max_r8_f.f \
	maxr_c16_f.f \
	avg_r8_f.f \
	al2_r8_f.f \
	al2_c16_f.f \
	polyval_r8_reverse_0.f \
	costheta_c16_f.f \
	sum_i4_f.f \
	sum_l2_f.f \
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
	Fig_gen_ver8.f \
	Fig_gen_ctf_ver0.f \
	Fig_gen_ctf_ver4.f \
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
