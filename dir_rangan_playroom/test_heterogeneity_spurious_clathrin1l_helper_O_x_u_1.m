%%%%%%%%;
% Now form images. ;
% Even though the montage is 10-by-10 (in terms of cells), ;
% we will take subsets of the montage so that we end up with 1024 images. ;
%%%%%%%%;
n_M_mon = prod(size(P_x_u_pack_)/n_x_u_pack); %<-- initial montage images. ;
% assume a 10-x-10 grid. ;
M_mon_x_c_xxM___ = permute(reshape(permute(reshape(P_x_u_pack_,[n_x_u_pack*10,n_x_u_pack,10]),[2,1,3]),[n_x_u_pack,n_x_u_pack,10*10]),[2,1,3]);
tmp_f=0;
for nx0=0:10-1; for nx1=0:10-1;
tmp_f = tmp_f + fnorm(M_mon_x_c_xxM___(:,:,1+nx0+10*nx1) - P_x_u_pack_(n_x_u_pack*nx0+[1:64],n_x_u_pack*nx1+[1:64]));
end;end;%for nx0=0:10-1; for nx1=0:10-1;
if (verbose>-1); disp(sprintf(' %% error: %0.16f',tmp_f)); end;
%%%%;
% redefine above to include 1024=32-x-32 images. ;
%%%%;
n_M_mon = 1024; %<-- standard number. ;
M_mon_x_c_xxM___ = zeros(n_x_u_pack,n_x_u_pack,n_M_mon);
tmp_stride0 = (size(P_x_u_pack_,1+0)-n_x_u_pack)/32;
tmp_stride1 = (size(P_x_u_pack_,1+1)-n_x_u_pack)/32;
for tmp_nx0=0:32-1;for tmp_nx1=0:32-1;
nM_mon = tmp_nx0 + 32*tmp_nx1;
tmp_nx0_ = tmp_nx0*tmp_stride0 + [0:n_x_u_pack-1];
tmp_nx1_ = tmp_nx1*tmp_stride1 + [0:n_x_u_pack-1];
M_mon_x_c_xxM___(:,:,1+nM_mon) = P_x_u_pack_(1+tmp_nx0_,1+tmp_nx1_);
end;end;%for tmp_nx0=0:32-1;for tmp_nx1=0:32-1;
%%%%%%%%;
n_x_M_u = n_x_u_pack;
dx = diameter_x_c/max(1,n_x_M_u);
M_mon_k_p_wkM__ = zeros(n_w_sum,n_M_mon);
M_mon_k_q_wkM__ = zeros(n_w_sum,n_M_mon);
for nM_mon=0:n_M_mon-1;
M_x_c_ = M_mon_x_c_xxM___(:,:,1+nM_mon);
M_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2 ;
M_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_);
M_mon_k_p_wkM__(:,1+nM_mon) = M_k_p_;
M_mon_k_q_wkM__(:,1+nM_mon) = M_k_q_;
clear M_x_c_ M_k_p_ M_k_q_ ;
end;%for nM_mon=0:n_M_mon-1;
%%%%;

%%%%%%%%;
% correct image-indices, assuming only a single CTF-cluster. ;
%%%%%%%%;
n_M = n_M_mon;
n_index_nM_from_ncluster = n_M;
n_index_nM_from_ncluster_ = n_M;
index_nCTF_from_nM_ = zeros(n_M,1);
index_nM_from_ncluster_ = transpose(0:n_M-1);
index_nM_from_ncluster__ = {transpose(0:n_M-1)};
index_ncluster_from_nCTF_ = 0;
index_ncluster_from_nM_ = zeros(n_M,1);
%%%%%%%%;
% align to a_p123_k_Y_norm_yk_. ;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.fname_align_a_k_Y_pre = sprintf('%s_mat/test_heterogeneity_spurious_clathrin1l_p123_from_mon',dir_pm);
tmp_fname_mat = sprintf('%s.mat',parameter.fname_align_a_k_Y_pre);
if flag_realign &  exist(tmp_fname_mat,'file'); delete(tmp_fname_mat); end;
parameter.n_iteration = 4;
tmp_n_cluster = 1; %<-- just use the first cluster. ;
tmp_index_ncluster_from_nCTF_ = 0;
tmp_n_CTF = 1;
tmp_index_nCTF_from_nM_ = zeros(n_M,1);
tmp_CTF_k_p_r_kC__ = CTF_k_p_r_kC__(:,1+0); %<-- just use the first cluster. ;
tmp_t = tic();
[ ...
 parameter ...
,p123_from_mon_a_k_Y_reco_yki__ ...
,p123_from_mon_corr_a_k_Y_i_ ...
,p123_from_mon_euler_polar_a_Mi__ ...
,p123_from_mon_euler_azimu_b_Mi__ ...
,p123_from_mon_euler_gamma_z_Mi__ ...
,p123_from_mon_image_delta_x_acc_Mi__ ...
,p123_from_mon_image_delta_y_acc_Mi__ ...
,p123_from_mon_image_delta_x_upd_Mi__ ...
,p123_from_mon_image_delta_y_upd_Mi__ ...
,p123_from_mon_flag_image_delta_upd_Mi__ ...
,p123_from_mon_image_I_value_Mi__ ...
,p123_from_mon_image_X_value_Mi__ ...
,p123_from_mon_image_S_index_Mi__ ...
,p123_from_mon_n_S ...
,p123_from_mon_template_viewing_azimu_b_all_ ...
,p123_from_mon_template_viewing_polar_a_all_ ...
,p123_from_mon_X_SMi___ ...
,p123_from_mon_delta_x_SMi___ ...
,p123_from_mon_delta_y_SMi___ ...
,p123_from_mon_gamma_z_SMi___ ...
,p123_from_mon_I_value_SMi___ ...
] = ...
pm_align_M_k_p_to_a_k_Y_3( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_mon_k_p_wkM__ ...
,tmp_n_cluster ...
,tmp_index_ncluster_from_nCTF_ ...
,tmp_n_CTF ...
,tmp_index_nCTF_from_nM_(1:n_M) ...
,tmp_CTF_k_p_r_kC__ ...
,l_max_ ...
,a_p123_k_Y_norm_yk_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% pm_align_M_k_p_to_a_k_Y_3: %0.3fs',tmp_t)); end;
%%%%%%%%;
% extract best alignment. ;
%%%%%%%%;
[~,ij_use] = max(real(p123_from_mon_corr_a_k_Y_i_)); niteration_use = ij_use-1;
p123_from_mon_a_k_Y_reco_yk__ = zeros(n_lm_max,n_k_p_r);
p123_from_mon_a_k_Y_reco_yk_ = p123_from_mon_a_k_Y_reco_yki__(:,1+niteration_use);
p123_from_mon_corr_a_k_Y = p123_from_mon_corr_a_k_Y_i_(1+niteration_use);
p123_from_mon_euler_polar_a_M_ = p123_from_mon_euler_polar_a_Mi__(:,1+niteration_use);
p123_from_mon_euler_azimu_b_M_ = p123_from_mon_euler_azimu_b_Mi__(:,1+niteration_use);
p123_from_mon_euler_gamma_z_M_ = p123_from_mon_euler_gamma_z_Mi__(:,1+niteration_use);
p123_from_mon_image_delta_x_acc_M_ = p123_from_mon_image_delta_x_acc_Mi__(:,1+niteration_use);
p123_from_mon_image_delta_y_acc_M_ = p123_from_mon_image_delta_y_acc_Mi__(:,1+niteration_use);
p123_from_mon_image_delta_x_upd_M_ = p123_from_mon_image_delta_x_upd_Mi__(:,1+niteration_use);
p123_from_mon_image_delta_y_upd_M_ = p123_from_mon_image_delta_y_upd_Mi__(:,1+niteration_use);
p123_from_mon_image_delta_x_M_ = p123_from_mon_image_delta_x_acc_M_ + p123_from_mon_image_delta_x_acc_M_;
p123_from_mon_image_delta_y_M_ = p123_from_mon_image_delta_y_acc_M_ + p123_from_mon_image_delta_y_acc_M_;
p123_from_mon_flag_image_delta_upd_M_ = p123_from_mon_flag_image_delta_upd_Mi__(:,1+niteration_use);
p123_from_mon_image_I_value_M_ = p123_from_mon_image_I_value_Mi__(:,1+niteration_use);
p123_from_mon_image_X_value_M_ = p123_from_mon_image_X_value_Mi__(:,1+niteration_use);
p123_from_mon_image_S_index_M_ = p123_from_mon_image_S_index_Mi__(:,1+niteration_use);
p123_from_mon_X_SM__ = p123_from_mon_X_SMi___(:,:,1+niteration_use);
p123_from_mon_delta_x_SM__ = p123_from_mon_delta_x_SMi___(:,:,1+niteration_use);
p123_from_mon_delta_y_SM__ = p123_from_mon_delta_y_SMi___(:,:,1+niteration_use);
p123_from_mon_gamma_z_SM__ = p123_from_mon_gamma_z_SMi___(:,:,1+niteration_use);
p123_from_mon_I_value_SM__ = p123_from_mon_I_value_SMi___(:,:,1+niteration_use);
%%%%%%%%;

%%%%%%%%;
% set snr: ;
%%%%%%%%;
sim_sigma = 0.50;

%%%%%%%%;
% generate simulated images. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/M_sim_sigma_%.2d_k_p_wkM__.mat',dir_pm,floor(100*sim_sigma));
if flag_recalc | ~exist(fname_mat,'file');
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
% Set up residual R_mon_k_p__. ;
%%%%%%%%;
R_mon_k_p_wkM__ = zeros(n_w_sum,n_M_mon);
M_sim_k_p_wkM__ = zeros(n_w_sum,n_M_mon);
T_mon_x_c_l2_M_ = zeros(n_M_mon,1);
R_mon_x_c_l2_M_ = zeros(n_M_mon,1);
N_mon_x_c_l2_M_ = zeros(n_M_mon,1);
M_mon_x_c_l2_M_ = zeros(n_M_mon,1);
M_sim_x_c_l2_M_ = zeros(n_M_mon,1);
viewing_k_c_0_all_ = sin(viewing_polar_a_all_).*cos(viewing_azimu_b_all_);
viewing_k_c_1_all_ = sin(viewing_polar_a_all_).*sin(viewing_azimu_b_all_);
viewing_k_c_2_all_ = cos(viewing_polar_a_all_);
for nM_mon=0:n_M_mon-1;
if (mod(nM_mon,32)==0); disp(sprintf(' %% nM_mon %d/%d',nM_mon,n_M_mon)); end;
tmp_euler_polar_a = +p123_from_mon_euler_polar_a_M_(1+nM_mon);
tmp_euler_azimu_b = +p123_from_mon_euler_azimu_b_M_(1+nM_mon);
tmp_euler_gamma_z = +p123_from_mon_euler_gamma_z_M_(1+nM_mon);
tmp_image_delta_x = +1.0*p123_from_mon_image_delta_x_M_(1+nM_mon);
tmp_image_delta_y = +1.0*p123_from_mon_image_delta_y_M_(1+nM_mon);
M_k_p_ = M_mon_k_p_wkM__(:,1+nM_mon);
tmp_k_c_0 = sin(tmp_euler_polar_a)*cos(tmp_euler_azimu_b);
tmp_k_c_1 = sin(tmp_euler_polar_a)*sin(tmp_euler_azimu_b);
tmp_k_c_2 = cos(tmp_euler_polar_a);
nS = knnsearch([viewing_k_c_0_all_,viewing_k_c_1_all_,viewing_k_c_2_all_],[tmp_k_c_0,tmp_k_c_1,tmp_k_c_2]) - 1;
S_k_p_ = S_p123_k_p__(:,1+nS);
S_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u_pack^2) * n_w_sum;
T_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_,+tmp_euler_gamma_z);
T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_,-tmp_image_delta_x,-tmp_image_delta_y);
T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+nM_mon)),n_w_sum,n_w_sum)*T_k_p_;
tmp_TT_k_p = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,T_k_p_)/(2*pi);
tmp_TM_k_p = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,M_k_p_)/(2*pi);
tmp_MM_k_p = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,M_k_p_,M_k_p_)/(2*pi);
T_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u_pack^2) * n_w_sum;
T_mon_x_c_l2 = sum(abs(T_x_c_).^2,'all')*dx_u_pack.^2;
M_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u_pack^2) * n_w_sum;
tmp_TM_x_c = sum(conj(T_x_c_).*M_x_c_,'all')*dx_u_pack.^2;
tmp_TT_x_c = sum(conj(T_x_c_).*T_x_c_,'all')*dx_u_pack.^2;
tmp_MM_x_c = sum(conj(M_x_c_).*M_x_c_,'all')*dx_u_pack.^2;
N_x_c_ = tmp_TM_x_c/max(1e-12,tmp_TT_x_c)*T_x_c_; %<-- projection of M_x_c_ onto T_x_c_. ;
R_x_c_ = M_x_c_ - N_x_c_;
R_mon_x_c_l2 = sum(abs(R_x_c_).^2,'all')*dx_u_pack.^2;
N_mon_x_c_l2 = sum(abs(N_x_c_).^2,'all')*dx_u_pack.^2;
M_mon_x_c_l2 = sum(abs(M_x_c_).^2,'all')*dx_u_pack.^2;
tmp_f = sqrt(R_mon_x_c_l2)/max(1e-12,sqrt(T_mon_x_c_l2));
M_sim_x_c_ = R_x_c_ + T_x_c_*sim_sigma*tmp_f; %<-- snr roughly sim_sigma. ;
M_sim_x_c_l2 = sum(abs(M_sim_x_c_).^2,'all')*dx_u_pack.^2;
M_sim_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,M_sim_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_u_pack^2)*dx_u_pack^2;
M_sim_k_p_wkM__(:,1+nM_mon) = M_sim_k_p_;
T_mon_x_c_l2_M_(1+nM_mon) = T_mon_x_c_l2;
R_mon_x_c_l2_M_(1+nM_mon) = R_mon_x_c_l2;
N_mon_x_c_l2_M_(1+nM_mon) = N_mon_x_c_l2;
M_mon_x_c_l2_M_(1+nM_mon) = M_mon_x_c_l2;
M_sim_x_c_l2_M_(1+nM_mon) = M_sim_x_c_l2;
R_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,R_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_u_pack^2)*dx_u_pack^2;
R_mon_k_p_wkM__(:,1+nM_mon) = R_k_p_;
clear tmp_euler_polar_a tmp_euler_azimu_b tmp_euler_gamma_z tmp_image_delta_x tmp_image_delta_y M_k_p_ tmp_k_c_0 tmp_k_c_1 tmp_k_c_2 nS S_k_p_ T_k_p_ T_k_p_ T_k_p_ tmp_TT_k_p tmp_TM_k_p tmp_MM_k_p T_x_c_ M_x_c_ tmp_TM_x_c tmp_TT_x_c tmp_MM_x_c N_x_c_ R_x_c_ R_k_p_ M_sim_x_c_ M_xim_k_p_ ;
end;%for nM_mon=0:n_M_mon-1;
%%%%%%%%;
% R_mon_k_q_wkM__. ;
%%%%%%%%;
M_sim_k_q_wkM__ = zeros(n_w_sum,n_M_mon);
for nM_mon=0:n_M_mon-1;
M_sim_k_q_wkM__(:,1+nM_mon) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_sim_k_p_wkM__(:,1+nM_mon));
end;%for nM_mon=0:n_M_mon-1;
R_mon_k_q_wkM__ = zeros(n_w_sum,n_M_mon);
for nM_mon=0:n_M_mon-1;
R_mon_k_q_wkM__(:,1+nM_mon) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,R_mon_k_p_wkM__(:,1+nM_mon));
end;%for nM_mon=0:n_M_mon-1;
%%%%%%%%;
save(fname_mat ...
,'sim_sigma' ...
,'n_M_mon' ...
,'R_mon_k_p_wkM__' ...
,'M_sim_k_p_wkM__' ...
,'T_mon_x_c_l2_M_' ...
,'R_mon_x_c_l2_M_' ...
,'N_mon_x_c_l2_M_' ...
,'M_mon_x_c_l2_M_' ...
,'M_sim_x_c_l2_M_' ...
,'M_sim_k_q_wkM__' ...
,'R_mon_k_q_wkM__' ...
);
%%%%%%%%;
end;%if flag_recalc | ~exist(fname_mat,'file');
if  exist(fname_mat,'file');
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if  exist(fname_mat,'file');

fname_fig_pre = sprintf('%s_jpg/test_heterogeneity_spurious_clathrin1l_sim_sigma_%.2d_FIGA',dir_pm,floor(100*sim_sigma));
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_jpg_stripped = sprintf('%s_stripped.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 5; p_col = 7; np=0;
for pcol=0:p_col-1;
%%%%%%%%;
% plot some of the residuals. ;
%%%%%%%%;
viewing_k_c_0_all_ = sin(viewing_polar_a_all_).*cos(viewing_azimu_b_all_);
viewing_k_c_1_all_ = sin(viewing_polar_a_all_).*sin(viewing_azimu_b_all_);
viewing_k_c_2_all_ = cos(viewing_polar_a_all_);
nM_mon = max(0,min(n_M_mon-1,floor(n_M_mon*pcol/p_col)));
tmp_euler_polar_a = +p123_from_mon_euler_polar_a_M_(1+nM_mon);
tmp_euler_azimu_b = +p123_from_mon_euler_azimu_b_M_(1+nM_mon);
tmp_euler_gamma_z = +p123_from_mon_euler_gamma_z_M_(1+nM_mon);
tmp_image_delta_x = +1.0*p123_from_mon_image_delta_x_M_(1+nM_mon);
tmp_image_delta_y = +1.0*p123_from_mon_image_delta_y_M_(1+nM_mon);
M_k_p_ = M_mon_k_p_wkM__(:,1+nM_mon);
tmp_k_c_0 = sin(tmp_euler_polar_a)*cos(tmp_euler_azimu_b);
tmp_k_c_1 = sin(tmp_euler_polar_a)*sin(tmp_euler_azimu_b);
tmp_k_c_2 = cos(tmp_euler_polar_a);
nS = knnsearch([viewing_k_c_0_all_,viewing_k_c_1_all_,viewing_k_c_2_all_],[tmp_k_c_0,tmp_k_c_1,tmp_k_c_2]) - 1;
S_k_p_ = S_p123_k_p__(:,1+nS);
S_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u_pack^2) * n_w_sum;
T_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_,+tmp_euler_gamma_z);
T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_,-tmp_image_delta_x,-tmp_image_delta_y);
S_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u_pack^2) * n_w_sum;
T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+nM_mon)),n_w_sum,n_w_sum)*T_k_p_;
tmp_TT_k_p = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,T_k_p_)/(2*pi);
tmp_TM_k_p = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,M_k_p_)/(2*pi);
tmp_MM_k_p = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,M_k_p_,M_k_p_)/(2*pi);
T_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u_pack^2) * n_w_sum;
T_mon_x_c_l2 = sum(abs(T_x_c_).^2,'all')*dx_u_pack.^2;
M_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u_pack^2) * n_w_sum;
tmp_TM_x_c = sum(conj(T_x_c_).*M_x_c_,'all')*dx_u_pack.^2;
tmp_TT_x_c = sum(conj(T_x_c_).*T_x_c_,'all')*dx_u_pack.^2;
tmp_MM_x_c = sum(conj(M_x_c_).*M_x_c_,'all')*dx_u_pack.^2;
N_x_c_ = tmp_TM_x_c/max(1e-12,tmp_TT_x_c)*T_x_c_; %<-- projection of M_x_c_ onto T_x_c_. ;
R_x_c_ = M_x_c_ - N_x_c_;
R_mon_x_c_l2 = sum(abs(R_x_c_).^2,'all')*dx_u_pack.^2;
N_mon_x_c_l2 = sum(abs(N_x_c_).^2,'all')*dx_u_pack.^2;
M_mon_x_c_l2 = sum(abs(M_x_c_).^2,'all')*dx_u_pack.^2;
tmp_f = sqrt(R_mon_x_c_l2)/max(1e-12,sqrt(T_mon_x_c_l2));
M_sim_x_c_ = R_x_c_ + T_x_c_*sim_sigma*tmp_f; %<-- snr roughly sim_sigma. ;
M_sim_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,M_sim_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_u_pack^2)*dx_u_pack^2;
R_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,R_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_u_pack^2)*dx_u_pack^2;
R_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,R_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u_pack^2) * n_w_sum;
%%%%;
subplot(p_row,p_col,1+pcol+0*p_col);
imagesc_c(n_x_u_pack,x_u_pack_0_,n_x_u_pack,x_u_pack_1_,real(S_x_c_),[],[]);
axis image; axisnotick; title('real(S_mon_x_c_)','Interpreter','none');
subplot(p_row,p_col,1+pcol+1*p_col);
imagesc_c(n_x_u_pack,x_u_pack_0_,n_x_u_pack,x_u_pack_1_,real(T_x_c_),[],[]);
axis image; axisnotick; title('real(T_mon_x_c_)','Interpreter','none');
subplot(p_row,p_col,1+pcol+2*p_col);
imagesc_c(n_x_u_pack,x_u_pack_0_,n_x_u_pack,x_u_pack_1_,real(M_x_c_),[],[]);
axis image; axisnotick; title('real(M_mon_x_c_)','Interpreter','none');
subplot(p_row,p_col,1+pcol+3*p_col);
imagesc_c(n_x_u_pack,x_u_pack_0_,n_x_u_pack,x_u_pack_1_,real(R_x_c_),[],[]);
axis image; axisnotick; title('real(R_mon_x_c_)','Interpreter','none');
subplot(p_row,p_col,1+pcol+4*p_col);
imagesc_c(n_x_u_pack,x_u_pack_0_,n_x_u_pack,x_u_pack_1_,real(M_sim_x_c_),[],[]);
axis image; axisnotick; title('real(M_sim_x_c_)','Interpreter','none');
%%%%;
end;%for pcol=0:p_col-1;
%%%%%%%%;
print('-djpeg',fname_fig_jpg_stripped);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
%%%%%%%%;
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_jpg));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;
