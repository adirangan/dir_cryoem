clear; test_pm_trpv1_1;

a_k_Y_true_ = a_k_Y_quad_;
euler_polar_a_true_ = +euler_polar_a_star_(1+(0:n_M-1));
euler_azimu_b_true_ = +euler_azimu_b_star_(1+(0:n_M-1));
euler_gamma_z_true_ = -euler_gamma_z_star_(1+(0:n_M-1)); %<-- note sign change. ;
image_delta_x_true_ = +image_delta_x_star_plus_M_abs_x_c_0_avg_(1+(0:n_M-1));
image_delta_y_true_ = +image_delta_y_star_plus_M_abs_x_c_1_avg_(1+(0:n_M-1));

delta_sigma_use = delta_sigma;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
dat_n_UX_rank_ = [1,2,4,8,16];  n_dat_n_UX_rank = numel(dat_n_UX_rank_);
delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
ndat_rseed = 0;%for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
ndelta_r_max_factor = 2;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2)); %<-- allow large accumulated translations. ;
ndat_n_UX_rank=n_dat_n_UX_rank-1;%for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
dat_n_UX_rank = dat_n_UX_rank_(1+ndat_n_UX_rank);
%%%%%%%%;

%%%%%%%%;
% original run (no exclusion). ;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.rseed = dat_rseed;
parameter.delta_r_max = delta_r_max_use;
parameter.delta_r_upb = delta_r_max_upb;
parameter.dir_pm = dir_pm;
ampmut_wrap_wrap_3( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,dat_n_UX_rank ...
,UX__ ...
,X_weight_r_ ...
,n_M ...
,M_k_p__ ...
,n_ctf ...
,n_CTF_rank ...
,CTF_index_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_true_ ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
);

%%%%%%%%;
% modified run (local exclusion). ;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.rseed = dat_rseed;
parameter.delta_r_max = delta_r_max_use;
parameter.delta_r_upb = delta_r_max_upb;
parameter.dir_pm = dir_pm;
parameter.flag_local_exclusion = 1;
ampmut_wrap_wrap_3( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,dat_n_UX_rank ...
,UX__ ...
,X_weight_r_ ...
,n_M ...
,M_k_p__ ...
,n_ctf ...
,n_CTF_rank ...
,CTF_index_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_true_ ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
);

ut_fname_pre = sprintf('%s_mat/ut%.4dn%.2dr%d',dir_pm,floor(1000*delta_r_max_use),dat_n_UX_rank,dat_rseed);
fname_ut_align_a_CTF_avg_UX_Y_mat = sprintf('%s_align_a_CTF_avg_UX_Y_.mat',ut_fname_pre);
fname_ut_align_a_k_Y_mat = sprintf('%s_align_a_k_Y_.mat',ut_fname_pre);
vt_fname_pre = sprintf('%s_mat/vt%.4dn%.2dr%d',dir_pm,floor(1000*delta_r_max_use),dat_n_UX_rank,dat_rseed);
fname_vt_align_a_CTF_avg_UX_Y_mat = sprintf('%s_align_a_CTF_avg_UX_Y_.mat',vt_fname_pre);
fname_vt_align_a_k_Y_mat = sprintf('%s_align_a_k_Y_.mat',vt_fname_pre);
tmp_ut_align_a_CTF_avg_UX_Y_ = load(fname_ut_align_a_CTF_avg_UX_Y_mat);
tmp_ut_align_a_k_Y_ = load(fname_ut_align_a_k_Y_mat);
tmp_vt_align_a_CTF_avg_UX_Y_ = load(fname_vt_align_a_CTF_avg_UX_Y_mat);
tmp_vt_align_a_k_Y_ = load(fname_vt_align_a_k_Y_mat);
xt_align_a_CTF_avg_UX_Y_ = [ tmp_ut_align_a_CTF_avg_UX_Y_.X_best_ ; tmp_vt_align_a_CTF_avg_UX_Y_.X_best_ ] ;
xt_align_a_k_Y_ = [ tmp_ut_align_a_k_Y_.X_best_ ; tmp_vt_align_a_k_Y_.X_best_ ] ;

ue1t_fname_pre = sprintf('%s_mat/ue1t%.4dn%.2dr%d',dir_pm,floor(1000*delta_r_max_use),dat_n_UX_rank,dat_rseed);
fname_ue1t_align_a_CTF_avg_UX_Y_mat = sprintf('%s_align_a_CTF_avg_UX_Y_.mat',ue1t_fname_pre);
fname_ue1t_align_a_k_Y_mat = sprintf('%s_align_a_k_Y_.mat',ue1t_fname_pre);
ve1t_fname_pre = sprintf('%s_mat/ve1t%.4dn%.2dr%d',dir_pm,floor(1000*delta_r_max_use),dat_n_UX_rank,dat_rseed);
fname_ve1t_align_a_CTF_avg_UX_Y_mat = sprintf('%s_align_a_CTF_avg_UX_Y_.mat',ve1t_fname_pre);
fname_ve1t_align_a_k_Y_mat = sprintf('%s_align_a_k_Y_.mat',ve1t_fname_pre);
tmp_ue1t_align_a_CTF_avg_UX_Y_ = load(fname_ue1t_align_a_CTF_avg_UX_Y_mat);
tmp_ue1t_align_a_k_Y_ = load(fname_ue1t_align_a_k_Y_mat);
tmp_ve1t_align_a_CTF_avg_UX_Y_ = load(fname_ve1t_align_a_CTF_avg_UX_Y_mat);
tmp_ve1t_align_a_k_Y_ = load(fname_ve1t_align_a_k_Y_mat);
xe1t_align_a_CTF_avg_UX_Y_ = [ tmp_ue1t_align_a_CTF_avg_UX_Y_.X_best_ ; tmp_ve1t_align_a_CTF_avg_UX_Y_.X_best_ ] ;
xe1t_align_a_k_Y_ = [ tmp_ue1t_align_a_k_Y_.X_best_ ; tmp_ve1t_align_a_k_Y_.X_best_ ] ;
