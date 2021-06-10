# Barnett adjusted for global include flag options. 6/24/16

# get the global options (barnett 6/23/16)
include make.inc.dist
-include make.inc

PROJECT=test_ver15.out
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

FSRCS = test_ver15.f  \
	ac1_c16.f \
	additive_noise_k_p.f \
	adi_fft1.f \
	adi_fft2.f \
	adi_rand_f.f \
	adi_randinclude.f \
	adi_randperm.f \
	af1_c16.f \
	afs_c16.f \
	al1_c16_f.f \
	al2_c16_f.f \
	al2_r8_f.f \
	alpha_SM_sort_0.f \
	alpha_SM_update_1.f \
	alpha_SM_write_0.f \
	avg_r8_f.f \
	block_0.f \
	block_1.f \
        ccongrlsq.f \
        chebexps.f \
	cl1_c16.f \
	cl1_i4.f \
	cl1_l2.f \
	cl1_r8.f \
	colorscale.f \
	costheta_c16_f.f \
	cp1_c16.f \
	cp1_i4.f \
	cp1_l2.f \
	cp1_r8.f \
	cps_c16.f \
	cross_r8.f \
	cs1_c16.f \
	cs1_i4.f \
	cs1_l2.f \
	cs1_r8.f \
	cxs_c16.f \
	cxs_i4.f \
	cxs_l2.f \
	cxs_r8.f \
        dfftpack.f \
	dfftw_block_many_0.f \
	dfftw_block_many_1.f \
	distance_r8.f \
	dot_c16.f \
	dot_c16_f.f \
	dot_r8.f \
	Fig_c16_bess2.f \
	Fig_c16_bessl.f \
	Fig_c16_carte.f \
	Fig_c16_carte_no_recenter.f \
	Fig_c16_polar.f \
	Fig_gen_ctf_ver0.f \
	Fig_gen_ctf_ver3.f \
	Fig_gen_ctf_ver4.f \
	Fig_gen_ver5.f \
	Fig_gen_ver8.f \
	Fig_header.f \
	Fig_R_ver1.f \
	Fig_text.f \
        gaussq.f \
	get_angle_to_vp_.f \
	get_C_S_use_.f \
	get_ctf_k_p_.f \
	get_ctf_ones_k_p_.f \
	get_CTF_R_S_periodize_0.f \
	get_CTF_R_S_use_.f \
	get_ctf_star_k_p_.f \
	get_ctf_x_c_.f \
	get_ctf_x_c.f \
	get_delta_0.f \
	get_delta_1.f \
	get_F2_x_c_.f \
	get_F2_x_c.f \
	get_gamma_0.f \
	get_interchange_delta.f \
	get_interchange_delta_RTRT_vs_RTTR.f \
        get_lsqdata_1.f \
        get_lsqdata_2.f \
	get_svd_0.f \
	get_svd_2.f \
	get_svd_polyval_U_d_0.f \
	get_svd_polyval_V_r_0.f \
	gramschmidt_c16.f \
        gridrouts3d.f \
	hsv2rgb.f \
	innerproduct_c.f \
	innerproduct_p.f \
	innerproduct_q__k_only.f \
	innerproduct_q__k_svdd_3.f \
	innerproduct_q__k_svdd_4.f \
	innerproduct_q__k_svdd.f \
	innerproduct_q__k_svdr_3.f \
	innerproduct_q__k_svdr_4.f \
	innerproduct_q__k_svdr_5.f \
	innerproduct_q__k_svdr.f \
	innerproduct_q_stretch.f \
	interp1_c16.f \
	interp2_c16.f \
	interp_c_to_p.f \
	interp_p_to_q.f \
	interp_p_to_q_fftw.f \
        kspace_model_norm.f \
	linspace.f \
        lsqctfsolvespharm.f \
	max_i4_f.f \
	max_r8.f \
	max_r8_f.f \
	maxr_c16_f.f \
	MDA_write_c16.f \
	MDA_write_i4.f \
	MDA_write_r8.f \
	min_i4_f.f \
        mk_simulated_slices_alpha2d.f \
        next235.f \
	normalize_c16.f \
	normalize_r8.f \
	pearson_c16.f \
	periodize_i.f \
	periodize_r8.f \
	ping.f \
	polyval_r8_reverse_0.f \
        prini.f \
	quicksort_c16.f \
	quicksort_i4.f \
	quicksort_r8.f \
	randinclude.f \
	randperm.f \
	rebuild_model_1.f \
	rebuild_model_2.f \
	recenter_c16.f \
	rotate_c_to_c.f \
	rotate_p_to_p_1.f \
	rotate_p_to_p.f \
	rotate_p_to_p_lagrange.f \
	rotate_p_to_p_single.f \
	rotate_p_to_p_single_spectral.f \
	rotate_p_to_p_spectral.f \
	rotate_q_to_q_0.f \
	rotate_q_to_q.f \
        spharmroutsall.f \
        spharmroutsnewp.f \
	stdlim_c16.f \
	sum_i4_f.f \
	sum_l2_f.f \
        template_gen.f \
	tesselation_define_octants.f \
	tesselation_get_nl_nm_ll.f \
	tesselation_index.f \
	tesselation_index_wrapper_0.f \
	tesselation_neighborhood.f \
	tesselation_neighborhood_wrapper_0.f \
	tesselation_size.f \
	test_alpha_update_5.f \
	test_alpha_update_MS_3.f \
	test_alpha_update_SM_3.f \
        testfourier.f \
        testgreatcircles.f \
	test_innerproduct_5_CTF_R_S_2.f \
	test_innerproduct_5_O_S_q_1.f \
	test_innerproduct_5_O_T_R_CTF_M_q_3.f \
	test_innerproduct_5_T_S_q_4.f \
	test_innerproduct_5_T_T_R_CTF_M_q_2.f \
	test_innerproduct_5_Zstore_1.f \
	test_innerproduct_5_Zstore_3a.f \
	test_innerproduct_5_Zstore_3b.f \
	test_innerproduct_5_Zstore_3c.f \
	test_innerproduct_5_Zstore_3d.f \
	test_innerproduct_5_Zstore_3e.f \
	test_innerproduct_5_Zstore_3.f \
	test_innerproduct_5_Zstore_3x.f \
	test_innerproduct_5_Zstore_4.f \
	test_innerproduct_5_ZZ_scan_2.f \
	test_innerproduct_5_ZZ_scan_4.f \
	test_innerproduct_5_ZZ_scan_4x.f \
	test_innerproduct_7.f \
	test_innerproduct_7_global_0.f \
	test_innerproduct_7_local_0.f \
	test_innerproduct_bruteforce_0.f \
	test_innerproduct_bruteforce_1.f \
	test_innerproduct_fast_CTF_R_S_1.f \
	test_max_0.f \
        testprojeval.f \
        test_residual_multa.f \
        test_residual_multah.f \
	test_residual_4.f \
	ti7_build_CTF_R_S_3.f \
	ti7_build_O_S_q_1.f \
	ti7_build_O_T_R_CTF_M_q_3.f \
	ti7_build_T_S_q_4.f \
	ti7_build_T_T_R_CTF_M_q_2.f \
	ti7_build_Z_S_q_5.f \
	ti7_build_Z_T_R_CTF_M_q_2.f \
	ti7_check_STxTRM_3.f \
	ti7_check_SxTTRM_3.f \
	ti7_check_SxZTRM_4.f \
	ti7_check_SZxTRM_4.f \
	transf_c_to_c.f \
	transf_p_to_p_1.f \
	transf_p_to_p.f \
	transf_p_to_p_single.f \
	transf_svd_q_to_q_3.f \
	transf_svd_q_to_q_5.f \
	transf_svd_q_to_q.f \
	transl_c_to_c.f \
	trn0_1_c16.f \
	trn0_c16.f \
	trn0_r8.f \
	unique_i4.f \
        utils.f \
	write_xxx_xx.f \
	xc1_c16.f \
	xx1_c16.f \
        yrecursion.f \
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
#	./$(PROJECT)
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
