%%%%%%%%;
% tests heterogeneity_spurious on trpv1c. ;
%%%%%%%%;

clear;
%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;
test_pm_trpv1c_9b;

CTF_k_p_r_kC__ = CTF_k_p_r__;
a_k_Y_true_yk_ = a_k_Y_true_;
delta_r_max_factor = 0.25;
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max = delta_r_max_use;
tolerance_pm = tolerance_master;
svd_eps = min(1e-2,tolerance_master);
n_delta_v_requested = 24;

%%%%%%%%;
% test innerproduct. ;
%%%%%%%%;
dx_u = diameter_x_c/n_x_u;
x_c_r__ = sqrt(x_c_0__.^2 + x_c_1__.^2);
S_x_c_ = 1/(2*pi)/0.1.^2 * exp(-x_c_r__.^2/(2*0.1.^2));
S_x_c_l2 = sum(abs(S_x_c_).^2,'all')*dx_u^2;
S_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,S_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_u^2)*dx_u^2;
S_k_p_l2 = sum(abs(S_k_p_).^2 .* weight_2d_k_all_) * (2*pi)^2;
S_k_p_l2 = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,S_k_p_,S_k_p_)/(2*pi);
R_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u^2) * n_w_sum;
R_x_c_l2 = sum(abs(R_x_c_).^2,'all')*dx_u^2;
if (verbose); disp(sprintf(' %% S_x_c_l2: %0.6f',S_x_c_l2)); end;
if (verbose); disp(sprintf(' %% S_k_p_l2: %0.6f',S_k_p_l2)); end;
if (verbose); disp(sprintf(' %% R_x_c_l2: %0.6f',R_x_c_l2)); end;

%%%%%%%%;
% Recalculate templates using unit norm volume. ;
%%%%%%%%;
[ ...
 X0_quad ...
,C0_quad ...
] = ...
register_spharm_to_spharm_3( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_quad_ ...
,a_k_Y_quad_ ...
);
%%%%%%%%;
a_k_Y_norm_yk_ = a_k_Y_quad_/max(1e-12,sqrt(X0_quad));
a_k_Y_norm_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_norm_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_norm_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
[ ...
 X0_norm ...
,C0_norm ...
] = ...
register_spharm_to_spharm_3( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_norm_yk_ ...
,a_k_Y_norm_yk_ ...
);
%%%%%%%%;
template_k_eq_d = -1;
viewing_k_eq_d = 1.0; %<-- subsample just a few templates for visualization. ;
tmp_t = tic();
[ ...
 S_k_p_wkS__ ...
,~ ...
,n_S ...
,template_viewing_azimu_b_all_ ...
,template_viewing_polar_a_all_ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max_max ...
,n_k_p_r ...
,reshape(a_k_Y_norm_yk__,[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% pm_template_2: %0.6fs',tmp_t)); end;
S_k_p_wkS__ = reshape(S_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
%%%%%%%%;
S_k_q_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_q_wkS__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS));
end;%for nS=0:n_S-1;
%%%%%%%%;

flag_disp=1;
%%%%%%%%;
% examine S_k_p_. ;
%%%%%%%%;
if flag_disp;
tmp_n_S = min(n_S,24);
p_row=6; p_col=4; np=0;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024,1024*2]);
for nS=0:tmp_n_S-1;
flag_drawnow=0;
subplot(p_row,p_col,1+np); np=np+1; if np==p_row*p_col;np=0; flag_drawnow=0; end;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(abs(S_k_p_wkS__(:,1+nS))));
axis image; axisnotick; title(sprintf('nS %d',nS)); 
if flag_drawnow; drawnow(); end;
end;%for nS=0:tmp_n_S-1;
end;%if flag_disp;

flag_disp=1;
%%%%%%%%;
% examine M_k_p_. ;
%%%%%%%%;
if flag_disp;
tmp_n_M = min(n_M,6);
p_row=6; p_col=4; np=0;
n_S = n_viewing_all;
viewing_k_c_0_all_ = sin(viewing_polar_a_all_).*cos(viewing_azimu_b_all_);
viewing_k_c_1_all_ = sin(viewing_polar_a_all_).*sin(viewing_azimu_b_all_);
viewing_k_c_2_all_ = cos(viewing_polar_a_all_);
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024,1024*2]);
for nM=0:tmp_n_M-1;
tmp_euler_polar_a = +euler_polar_a_true_(1+nM);
tmp_euler_azimu_b = +euler_azimu_b_true_(1+nM);
tmp_euler_gamma_z = +euler_gamma_z_true_(1+nM);
tmp_image_delta_x = +1.0*image_delta_x_true_(1+nM);
tmp_image_delta_y = +1.0*image_delta_y_true_(1+nM);
M_k_p_ = M_k_p__(:,1+nM);
tmp_k_c_0 = sin(tmp_euler_polar_a)*cos(tmp_euler_azimu_b);
tmp_k_c_1 = sin(tmp_euler_polar_a)*sin(tmp_euler_azimu_b);
tmp_k_c_2 = cos(tmp_euler_polar_a);
nS = knnsearch([viewing_k_c_0_all_,viewing_k_c_1_all_,viewing_k_c_2_all_],[tmp_k_c_0,tmp_k_c_1,tmp_k_c_2]) - 1;
S_k_p_ = S_k_p_wkS__(:,1+nS);
T_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_,+tmp_euler_gamma_z);
T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_,-tmp_image_delta_x,-tmp_image_delta_y);
T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,CTF_k_p__(:,1+index_nCTF_from_nM_(1+nM)),n_w_sum,n_w_sum)*T_k_p_;
tmp_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,T_k_p_)/(2*pi);
tmp_MM = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,M_k_p_,M_k_p_)/(2*pi);
tmp_TM = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,M_k_p_)/(2*pi);
N_k_p_ = tmp_TM/max(1e-12,tmp_TT)*T_k_p_; %<-- projection of M_k_p_ onto T_k_p_. ;
R_k_p_ = M_k_p_ - T_k_p_;
Rlim_ = max([max(abs(T_k_p_)),max(abs(M_k_p_)),max(abs(R_k_p_)),max(abs(N_k_p_))])*[-1.0,+1.0];
flag_drawnow = 0;
%%%%;
subplot(p_row,p_col,1+np); np=np+1; if np==p_row*p_col;np=0; flag_drawnow = 1; end;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_),Rlim_);
axis image; axisnotick; title(sprintf('nM %d (T)',nM)); 
%%%%;
subplot(p_row,p_col,1+np); np=np+1; if np==p_row*p_col;np=0; flag_drawnow = 1; end;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_),Rlim_);
axis image; axisnotick; title(sprintf('nM %d (M)',nM)); 
%%%%;
subplot(p_row,p_col,1+np); np=np+1; if np==p_row*p_col;np=0; flag_drawnow = 1; end;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(N_k_p_),Rlim_);
axis image; axisnotick; title(sprintf('nM %d (N)',nM)); 
%%%%;
subplot(p_row,p_col,1+np); np=np+1; if np==p_row*p_col;np=0; flag_drawnow = 1; end;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(R_k_p_),Rlim_);
axis image; axisnotick; title(sprintf('nM %d (R)',nM)); 
%%%%;
if flag_drawnow; drawnow(); end;
end;%for nM=0:tmp_n_M-1;
end;%if flag_disp;
 
flag_disp=1;
%%%%%%%%;
% examine M_x_c_. ;
%%%%%%%%;
if flag_disp;
tmp_n_M = min(n_M,6);
p_row=6; p_col=4; np=0;
n_S = n_viewing_all;
viewing_k_c_0_all_ = sin(viewing_polar_a_all_).*cos(viewing_azimu_b_all_);
viewing_k_c_1_all_ = sin(viewing_polar_a_all_).*sin(viewing_azimu_b_all_);
viewing_k_c_2_all_ = cos(viewing_polar_a_all_);
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024,1024*2]);
%R_x_c_cat_ = [];
for nM=0:tmp_n_M-1;
tmp_euler_polar_a = +euler_polar_a_true_(1+nM);
tmp_euler_azimu_b = +euler_azimu_b_true_(1+nM);
tmp_euler_gamma_z = +euler_gamma_z_true_(1+nM);
tmp_image_delta_x = +1.0*image_delta_x_true_(1+nM);
tmp_image_delta_y = +1.0*image_delta_y_true_(1+nM);
M_k_p_ = M_k_p__(:,1+nM);
tmp_k_c_0 = sin(tmp_euler_polar_a)*cos(tmp_euler_azimu_b);
tmp_k_c_1 = sin(tmp_euler_polar_a)*sin(tmp_euler_azimu_b);
tmp_k_c_2 = cos(tmp_euler_polar_a);
nS = knnsearch([viewing_k_c_0_all_,viewing_k_c_1_all_,viewing_k_c_2_all_],[tmp_k_c_0,tmp_k_c_1,tmp_k_c_2]) - 1;
S_k_p_ = S_k_p_wkS__(:,1+nS);
T_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_,+tmp_euler_gamma_z);
T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_,-tmp_image_delta_x,-tmp_image_delta_y);
T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,CTF_k_p__(:,1+index_nCTF_from_nM_(1+nM)),n_w_sum,n_w_sum)*T_k_p_;
tmp_TT_k_p = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,T_k_p_)/(2*pi);
tmp_TM_k_p = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,M_k_p_)/(2*pi);
tmp_MM_k_p = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,M_k_p_,M_k_p_)/(2*pi);
T_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u^2) * n_w_sum;
M_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u^2) * n_w_sum;
tmp_TM_x_c = sum(conj(T_x_c_).*M_x_c_,'all')*dx_u.^2;
tmp_TT_x_c = sum(conj(T_x_c_).*T_x_c_,'all')*dx_u.^2;
tmp_MM_x_c = sum(conj(M_x_c_).*M_x_c_,'all')*dx_u.^2;
N_x_c_ = tmp_TM_x_c/max(1e-12,tmp_TT_x_c)*T_x_c_; %<-- projection of M_x_c_ onto T_x_c_. ;
R_x_c_ = M_x_c_ - N_x_c_;
Rlim_ = max([max(abs(T_x_c_)),max(abs(M_x_c_)),max(abs(R_x_c_)),max(abs(N_x_c_))])*[-1.0,+1.0];
%R_x_c_cat_ = [R_x_c_cat_;real(R_x_c_(:))];
flag_drawnow = 0;
%%%%;
subplot(p_row,p_col,1+np); np=np+1; if np==p_row*p_col;np=0; flag_drawnow = 1; end;
imagesc_c(n_x_u,x_c_0_,n_x_u,x_c_1_,real(T_x_c_),Rlim_);
axis image; axisnotick; title(sprintf('nM %d (T)',nM)); 
%%%%;
subplot(p_row,p_col,1+np); np=np+1; if np==p_row*p_col;np=0; flag_drawnow = 1; end;
imagesc_c(n_x_u,x_c_0_,n_x_u,x_c_1_,real(M_x_c_),Rlim_);
axis image; axisnotick; title(sprintf('nM %d (M)',nM)); 
%%%%;
subplot(p_row,p_col,1+np); np=np+1; if np==p_row*p_col;np=0; flag_drawnow = 1; end;
imagesc_c(n_x_u,x_c_0_,n_x_u,x_c_1_,real(N_x_c_),Rlim_);
axis image; axisnotick; title(sprintf('nM %d (N)',nM)); 
%%%%;
subplot(p_row,p_col,1+np); np=np+1; if np==p_row*p_col;np=0; flag_drawnow = 1; end;
imagesc_c(n_x_u,x_c_0_,n_x_u,x_c_1_,real(R_x_c_),Rlim_);
axis image; axisnotick; title(sprintf('nM %d (R)',nM)); 
%%%%;
if flag_drawnow; drawnow(); end;
end;%for nM=0:tmp_n_M-1;
end;%if flag_disp;

%%%%%%%%;
% Now estimate noise variance for M_k_p__. ;
% Verdict: per-pixel noise std is about 0.10, ;
% while per-pixel signal (when it exists) is ~0.05 or so. ;
%%%%%%%%;
tmp_n_M = min(n_M,1024);
N_factor_M_ = zeros(n_M,1);
R_sum1_k_p_ = zeros(n_w_sum,1);
R_sum2_k_p_ = zeros(n_w_sum,1);
N_sum1_k_p_ = zeros(n_w_sum,1);
N_sum2_k_p_ = zeros(n_w_sum,1);
R_sum1_x_c_ = zeros(n_x_u,n_x_u);
R_sum2_x_c_ = zeros(n_x_u,n_x_u);
N_sum1_x_c_ = zeros(n_x_u,n_x_u);
N_sum2_x_c_ = zeros(n_x_u,n_x_u);
n_S = n_viewing_all;
viewing_k_c_0_all_ = sin(viewing_polar_a_all_).*cos(viewing_azimu_b_all_);
viewing_k_c_1_all_ = sin(viewing_polar_a_all_).*sin(viewing_azimu_b_all_);
viewing_k_c_2_all_ = cos(viewing_polar_a_all_);
for nM=0:tmp_n_M-1;
if verbose; if (mod(nM,32)==0); disp(sprintf(' %% nM %d/%d',nM,tmp_n_M)); end; end;
tmp_euler_polar_a = +euler_polar_a_true_(1+nM);
tmp_euler_azimu_b = +euler_azimu_b_true_(1+nM);
tmp_euler_gamma_z = +euler_gamma_z_true_(1+nM);
tmp_image_delta_x = +1.0*image_delta_x_true_(1+nM);
tmp_image_delta_y = +1.0*image_delta_y_true_(1+nM);
M_k_p_ = M_k_p__(:,1+nM);
tmp_k_c_0 = sin(tmp_euler_polar_a)*cos(tmp_euler_azimu_b);
tmp_k_c_1 = sin(tmp_euler_polar_a)*sin(tmp_euler_azimu_b);
tmp_k_c_2 = cos(tmp_euler_polar_a);
nS = knnsearch([viewing_k_c_0_all_,viewing_k_c_1_all_,viewing_k_c_2_all_],[tmp_k_c_0,tmp_k_c_1,tmp_k_c_2]) - 1;
S_k_p_ = S_k_p_wkS__(:,1+nS);
T_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_,+tmp_euler_gamma_z);
T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_,-tmp_image_delta_x,-tmp_image_delta_y);
T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,CTF_k_p__(:,1+index_nCTF_from_nM_(1+nM)),n_w_sum,n_w_sum)*T_k_p_;
tmp_TT_k_p = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,T_k_p_)/(2*pi);
tmp_TM_k_p = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,M_k_p_)/(2*pi);
tmp_MM_k_p = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,M_k_p_,M_k_p_)/(2*pi);
N_factor_k_p = real(tmp_TM_k_p/max(1e-12,tmp_TT_k_p));
N_k_p_ = N_factor_k_p*T_k_p_; %<-- projection of M_k_p_ onto T_k_p_. ;
R_k_p_ = M_k_p_ - N_k_p_;
R_sum1_k_p_ = R_sum1_k_p_ + real(R_k_p_).^1;
R_sum2_k_p_ = R_sum2_k_p_ + real(R_k_p_).^2;
N_sum1_k_p_ = N_sum1_k_p_ + real(N_k_p_).^1;
N_sum2_k_p_ = N_sum2_k_p_ + real(N_k_p_).^2;
T_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u^2) * n_w_sum;
M_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u^2) * n_w_sum;
tmp_TM_x_c = sum(conj(T_x_c_).*M_x_c_,'all')*dx_u.^2;
tmp_TT_x_c = sum(conj(T_x_c_).*T_x_c_,'all')*dx_u.^2;
tmp_MM_x_c = sum(conj(M_x_c_).*M_x_c_,'all')*dx_u.^2;
N_factor_x_c = real(tmp_TM_x_c/max(1e-12,tmp_TT_x_c));
N_factor_M_(1+nM) = N_factor_x_c;
N_x_c_ = N_factor_x_c*T_x_c_; %<-- projection of M_x_c_ onto T_x_c_. ;
R_x_c_ = M_x_c_ - N_x_c_;
R_sum1_x_c_ = R_sum1_x_c_ + real(R_x_c_).^1;
R_sum2_x_c_ = R_sum2_x_c_ + real(R_x_c_).^2;
N_sum1_x_c_ = N_sum1_x_c_ + real(N_x_c_).^1;
N_sum2_x_c_ = N_sum2_x_c_ + real(N_x_c_).^2;
end;%for nM=0:tmp_n_M-1;
R_avg_k_p_ = R_sum1_k_p_/max(1,tmp_n_M);
R_var_k_p_=  R_sum2_k_p_/max(1,tmp_n_M) - R_avg_k_p_.^2;
R_std_k_p_ = sqrt(R_var_k_p_);
N_avg_k_p_ = N_sum1_k_p_/max(1,tmp_n_M);
N_var_k_p_=  N_sum2_k_p_/max(1,tmp_n_M) - N_avg_k_p_.^2;
N_std_k_p_ = sqrt(N_var_k_p_);
R_avg_x_c_ = R_sum1_x_c_/max(1,tmp_n_M);
R_var_x_c_=  R_sum2_x_c_/max(1,tmp_n_M) - R_avg_x_c_.^2;
R_std_x_c_ = sqrt(R_var_x_c_);
N_avg_x_c_ = N_sum1_x_c_/max(1,tmp_n_M);
N_var_x_c_=  N_sum2_x_c_/max(1,tmp_n_M) - N_avg_x_c_.^2;
N_std_x_c_ = sqrt(N_var_x_c_);
RR_x_c = sum(R_var_x_c_,'all')*dx_u.^2; %<-- covers an area of 4. ;
M_var = RR_x_c/diameter_x_c^2; %<-- mean(R_var_x_c_). ;
M_std = sqrt(M_var);
%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1);imagesc_c(n_x_u,x_c_0_,n_x_u,x_c_1_,R_avg_x_c_); title('R_avg_x_c_','Interpreter','none'); axis image; axisnotick;
subplot(2,2,2);imagesc_c(n_x_u,x_c_0_,n_x_u,x_c_1_,R_std_x_c_); title('R_std_x_c_','Interpreter','none'); axis image; axisnotick;
subplot(2,2,3);imagesc_c(n_x_u,x_c_0_,n_x_u,x_c_1_,N_avg_x_c_); title('N_avg_x_c_','Interpreter','none'); axis image; axisnotick;
subplot(2,2,4);imagesc_c(n_x_u,x_c_0_,n_x_u,x_c_1_,N_std_x_c_); title('N_std_x_c_','Interpreter','none'); axis image; axisnotick;
%%%%%%%%;

flag_check=0;
if flag_check;
%%%%%%%%;
% Now let us see what the noise looks like when converting from fourier to real space. ;
% Verdict: close to iid with variance ~ M_sigma^2*n_w_max/(pi*k_p_r_max^3)/k_eq_d/mean(RR_k_p_vs_x_c_) ;
% Note that this is a complicated relationship because the integration weights are designed for smooth functions. ;
%%%%%%%%;
M_sigma = 1.0;
n_shuffle = 1024;
R_sum1_k_p_ = zeros(n_w_sum,1);
R_sum2_k_p_ = zeros(n_w_sum,1);
R_sum1_x_c_ = zeros(n_x_u,n_x_u);
R_sum2_x_c_ = zeros(n_x_u,n_x_u);
RR_k_p_vs_x_c_ = zeros(n_shuffle,1);
for nshuffle=0:n_shuffle-1;
if verbose; if (mod(nshuffle,32)==0); disp(sprintf(' %% nshuffle %d/%d',nshuffle,n_shuffle)); end; end;
rng(nshuffle);
R_k_p_ = M_sigma*randn(n_w_sum,1);
tmp_RR_k_p = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,R_k_p_,R_k_p_)/(2*pi);
R_sum1_k_p_ = R_sum1_k_p_ + R_k_p_;
R_sum2_k_p_ = R_sum2_k_p_ + abs(R_k_p_).^2;
R_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,R_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u^2) * n_w_sum;
tmp_RR_x_c = sum(conj(R_x_c_).*R_x_c_,'all')*dx_u.^2;
RR_k_p_vs_x_c_(1+nshuffle) = tmp_RR_k_p/max(1e-12,tmp_RR_x_c);
R_sum1_x_c_ = R_sum1_x_c_ + R_x_c_;
R_sum2_x_c_ = R_sum2_x_c_ + abs(R_x_c_).^2;
end;%for nshuffle=0:n_shuffle-1;
R_avg_k_p_ = R_sum1_k_p_/max(1,n_shuffle);
R_var_k_p_=  R_sum2_k_p_/max(1,n_shuffle) - abs(R_avg_k_p_).^2;
R_std_k_p_ = sqrt(R_var_k_p_);
R_avg_x_c_ = R_sum1_x_c_/max(1,n_shuffle);
R_var_x_c_=  R_sum2_x_c_/max(1,n_shuffle) - abs(R_avg_x_c_).^2;
R_std_x_c_ = sqrt(R_var_x_c_);
%%%%;
if (verbose); disp(sprintf(' %% R_var_k_p %0.6f, R_var_x_c*n_w_max/(pi*k_p_r_max^3)/k_eq_d/mean(RR_k_p_vs_x_c_) %0.6f',mean(R_var_k_p_,'all'),mean(R_var_x_c_,'all')*n_w_max/(pi*k_p_r_max^3)/k_eq_d/mean(RR_k_p_vs_x_c_))); end;
%%%%%%%%;
end;%if flag_check;

flag_check=0;
if flag_check;
%%%%%%%%;
% Now let us see what the noise looks like when converting from fourier to real space. ;
%%%%%%%%;
M_sigma = 1.0;
n_shuffle = 1024;
%%%%%%%%;
k_factor = 1.2;
k_2_p_r_max = k_factor*48/(2*pi); k_2_eq_d = 0.6/(2*pi);
[ ...
 ~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,n_k_2_p_r ...
,k_2_p_r_ ...
,~ ...
] = ...
sample_sphere_7( ...
 verbose ...
,k_2_p_r_max ...
,k_2_eq_d ...
,'L' ...
) ;
%%%%%%%%;
n_w_factor = 2; n_w_2_max = n_w_factor*n_w_max; n_w_2_ = n_w_2_max*ones(n_k_2_p_r,1); n_w_2_sum = sum(n_w_2_);
[ ...
 n_w_2_ ...
,weight_2d_k_2_p_r_ ...
,weight_2d_k_2_all_ ...
] = ...
get_weight_2d_1( ...
 0*verbose ...
,n_k_2_p_r ...
,k_2_p_r_ ...
,k_2_p_r_max ...
,-1 ...
,n_w_2_ ...
);
assert(abs(sum(weight_2d_k_2_all_)*(4*pi^2) - (pi*k_2_p_r_max^2))<1e-6);
%%%%%%%%;
n_x_factor = 1.0; n_x_2_u  = n_x_factor*n_x_u; dx_2_u = dx_u/max(1,n_x_factor);
R2_sum1_k_p_ = zeros(n_w_2_sum,1);
R2_sum2_k_p_ = zeros(n_w_2_sum,1);
R2_sum1_x_c_ = zeros(n_x_2_u,n_x_2_u);
R2_sum2_x_c_ = zeros(n_x_2_u,n_x_2_u);
R2R2_k_p_vs_x_c_ = zeros(n_shuffle,1);
for nshuffle=0:n_shuffle-1;
if verbose; if (mod(nshuffle,32)==0); disp(sprintf(' %% nshuffle %d/%d',nshuffle,n_shuffle)); end; end;
rng(nshuffle);
R2_k_p_ = M_sigma*randn(n_w_2_sum,1);
R2_sum1_k_p_ = R2_sum1_k_p_ + R2_k_p_;
R2_sum2_k_p_ = R2_sum2_k_p_ + abs(R2_k_p_).^2;
R2_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_2_u,diameter_x_c,n_x_2_u,diameter_x_c,n_k_2_p_r,k_2_p_r_,n_w_2_,R2_k_p_.*weight_2d_k_2_all_*(2*pi)^2)*sqrt(n_x_2_u^2) * n_w_2_sum;
tmp_RR_k_p = innerproduct_p_quad(n_k_2_p_r,k_2_p_r_,weight_2d_k_2_p_r_,n_w_2_,n_w_2_sum,R2_k_p_,R2_k_p_)/(2*pi);
tmp_RR_x_c = sum(conj(R2_x_c_).*R2_x_c_,'all')*dx_2_u.^2;
R2R2_k_p_vs_x_c_(1+nshuffle) = tmp_RR_k_p/max(1e-12,tmp_RR_x_c);
R2_sum1_x_c_ = R2_sum1_x_c_ + R2_x_c_;
R2_sum2_x_c_ = R2_sum2_x_c_ + abs(R2_x_c_).^2;
end;%for nshuffle=0:n_shuffle-1;
R2_avg_k_p_ = R2_sum1_k_p_/max(1,n_shuffle);
R2_var_k_p_=  R2_sum2_k_p_/max(1,n_shuffle) - abs(R2_avg_k_p_).^2;
R2_std_k_p_ = sqrt(R2_var_k_p_);
R2_avg_x_c_ = R2_sum1_x_c_/max(1,n_shuffle);
R2_var_x_c_=  R2_sum2_x_c_/max(1,n_shuffle) - abs(R2_avg_x_c_).^2;
R2_std_x_c_ = sqrt(R2_var_x_c_);
%%%%;
%R_var_k_p_from_x_c = n_w_max/(pi*k_p_r_max^3)/k_eq_d;
R_var_k_p_from_x_c = mean(RR_k_p_vs_x_c_)/k_p_r_max^2 * 4/pi;
%R2_var_k_p_from_x_c = n_w_2_max/(pi*k_2_p_r_max^3)/k_2_eq_d;
R2_var_k_p_from_x_c = mean(R2R2_k_p_vs_x_c_)/k_2_p_r_max^2 * 4/pi;
if (verbose); disp(sprintf(' %% R_var_x_c*R_var_k_p_from_x_c %0.6f <-- R2_var_x_c*R2_var_k_p_from_x_c %0.6f ',mean(R_var_x_c_,'all')*R_var_k_p_from_x_c,mean(R2_var_x_c_,'all')*R2_var_k_p_from_x_c)); end;
%%%%%%%%;
end;%if flag_check;

flag_check=0;
if flag_check;
%%%%%%%%;
% Now let us see what the noise looks like when converting from real to fourier space. ;
% Verdict: iid with variance M_sigma^2/(diameter_x_c.^2*dx_u^2);
% This is quite a bit simpler than the inverse operation because there are no integration weights involved. ;
%%%%%%%%;
M_sigma = 1.0;
n_shuffle = 1024;
R_sum1_x_c_ = zeros(n_x_u,n_x_u);
R_sum2_x_c_ = zeros(n_x_u,n_x_u);
R_sum1_k_p_ = zeros(n_w_sum,1);
R_sum2_k_p_ = zeros(n_w_sum,1);
for nshuffle=0:n_shuffle-1;
if verbose; if (mod(nshuffle,32)==0); disp(sprintf(' %% nshuffle %d/%d',nshuffle,n_shuffle)); end; end;
rng(nshuffle);
R_x_c_ = M_sigma*randn(n_x_u,n_x_u);
R_sum1_x_c_ = R_sum1_x_c_ + R_x_c_;
R_sum2_x_c_ = R_sum2_x_c_ + abs(R_x_c_).^2;
R_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,R_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_u^2)*dx_u^2;
R_sum1_k_p_ = R_sum1_k_p_ + R_k_p_;
R_sum2_k_p_ = R_sum2_k_p_ + abs(R_k_p_).^2;
end;%for nshuffle=0:n_shuffle-1;
R_avg_x_c_ = R_sum1_x_c_/max(1,n_shuffle);
R_var_x_c_=  R_sum2_x_c_/max(1,n_shuffle) - abs(R_avg_x_c_).^2;
R_std_x_c_ = sqrt(R_var_x_c_);
R_avg_k_p_ = R_sum1_k_p_/max(1,n_shuffle);
R_var_k_p_=  R_sum2_k_p_/max(1,n_shuffle) - abs(R_avg_k_p_).^2;
R_std_k_p_ = sqrt(R_var_k_p_);
%%%%;
if (verbose); disp(sprintf(' %% R_var_x_c %0.6f <-- R_var_k_p/(diameter_x_c.^2*dx_u.^2) %0.6f ',mean(R_var_x_c_,'all'),mean(R_var_k_p_,'all')/(diameter_x_c.^2*dx_u.^2))); end;
%%%%%%%%;
end;%if flag_check;

%%%%%%%%;
% Now cluster the CTF based on tolerance_cluster. ;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.tolerance_master = tolerance_master;
[ ...
 parameter ...
,index_ncluster_from_nCTF_ ...
] = ...
knn_cluster_CTF_k_p_r_kC__0( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_CTF ...
,CTF_k_p_r_kC__ ...
);
%%%%%%%%;
n_cluster = 1+max(index_ncluster_from_nCTF_);
index_ncluster_from_nM_ = index_ncluster_from_nCTF_(1+index_nCTF_from_nM_);
index_nM_from_ncluster__ = cell(n_cluster,1);
n_index_nM_from_ncluster_ = zeros(n_cluster,1);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster__{1+ncluster} = efind(index_ncluster_from_nM_==ncluster);
n_index_nM_from_ncluster_(1+ncluster) = numel(index_nM_from_ncluster__{1+ncluster});
end;%for ncluster=0:n_cluster-1;

%%%%%%%%;
% Now calculate average CTFs for each cluster. ;
%%%%%%%%;
CTF_k_p_r_xavg_kc__ = zeros(n_k_p_r,n_cluster);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),2);
CTF_k_p_r_xavg_kc__(:,1+ncluster) = CTF_k_p_r_xavg_k_;
end;%for ncluster=0:n_cluster-1;

%%%%%%%%;
% Now calculate idealized principal-modes for each nCTF. ;
%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
a_k_Y_norm_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_norm_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_norm_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_t = tic();
X_2d_xavg_dx_kkc___ = zeros(n_k_p_r,n_k_p_r,n_cluster);
X_2d_xavg_dx_weight_rc__ = zeros(n_k_p_r,n_cluster);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
CTF_k_p_r_xavg_kk__ = CTF_k_p_r_xavg_k_*transpose(CTF_k_p_r_xavg_k_);
tmp_delta_sigma = delta_r_max;
if (tmp_delta_sigma> 0);
[X_2d_xavg_dx_kk__,X_2d_xavg_dx_weight_r_] = principled_marching_cost_matrix_5(n_k_p_r,k_p_r_,weight_2d_k_p_r_,l_max_,[],[],a_k_Y_norm_yk_,CTF_k_p_r_xavg_kk__,tmp_delta_sigma);
end;%if (tmp_delta_sigma> 0);
if (tmp_delta_sigma==0);
[X_2d_xavg_dx_kk__,X_2d_xavg_dx_weight_r_] = principled_marching_cost_matrix_3(n_k_p_r,weight_2d_k_p_r_,l_max_max,a_k_Y_norm_yk__,CTF_k_p_r_xavg_kk__);
end;%if (tmp_delta_sigma==0);
X_2d_xavg_dx_kkc___(:,:,1+ncluster) = X_2d_xavg_dx_kk__;
X_2d_xavg_dx_weight_rc__(:,1+ncluster) = X_2d_xavg_dx_weight_r_;
clear X_2d_xavg_dx_kk__ X_2d_xavg_dx_weight_r_ ;
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_2d_xavg_dx_kkc___: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
UX_knc___ = zeros(n_k_p_r,n_UX_rank,n_cluster);
SX_kc__ = zeros(n_UX_rank,n_cluster);
pm_n_UX_rank_c_ = zeros(n_cluster,1);
for ncluster=0:n_cluster-1;
tmp_X__ = X_2d_xavg_dx_kkc___(:,:,1+ncluster);
[tmp_UX__,tmp_SX__,tmp_VX__] = svds(tmp_X__,n_UX_rank); tmp_SX_ = diag(tmp_SX__);
pm_n_UX_rank = max(find(tmp_SX_/max(tmp_SX_)> tolerance_pm));
UX_knc___(:,:,1+ncluster) = tmp_UX__(:,1+[0:n_UX_rank-1]);
SX_kc__(:,1+ncluster) = tmp_SX_(1+[0:n_UX_rank-1]);
pm_n_UX_rank_c_(1+ncluster) = pm_n_UX_rank;
if (verbose>1); disp(sprintf(' %% ncluster %.2d/%.2d: pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,n_UX_rank)); end;
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_knc___: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'X__',tmp_t);
%%%%%%%%;

pm_n_UX_rank_max = max(pm_n_UX_rank_c_);
pm_n_lm_max_max = (1+l_max_max)^2;
if (verbose);
for ncluster=0:n_cluster-1;
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
disp(sprintf(' %% ncluster %.2d/%.2d --> pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,pm_n_UX_rank_max));
end;%for ncluster=0:n_cluster-1;
end;%if (verbose);

FTK = [];
if isempty(FTK);
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_FTK_1',tmp_t);
end;%if isempty(FTK);
assert(FTK.svd_d_max>=delta_r_max);
assert(FTK.n_delta_v>=n_delta_v_requested);

%%%%%%%%;
% Now estimate distances between templates. ;
% Use one of the CTF clusters for now. ;
%%%%%%%%;
ncluster = 0; %<-- pick one of the CTF clusters. ;
UX_kn__ = UX_knc___(:,:,1+ncluster);
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
T_k_p_wkS__ = reshape(bsxfun(@times,reshape(S_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_k_q_wkS__ = reshape(bsxfun(@times,reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
tmp_t = tic();
svd_VUXT_lwnS____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_S,T_k_q_wkS__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXT_lwnS____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the templates. ;
%%%%%%%%;
tmp_t = tic();
UX_T_l2_dS__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_S,n_UX_rank,svd_VUXT_lwnS____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_l2_dS__: %0.3fs',tmp_t)); end;
if (verbose); disp(sprintf(' %% average l2-norm of templates: %0.16f',mean(UX_T_l2_dS__(:))/(pi*k_p_r_max^2))); end;
tmp_TT_S_ = zeros(n_S,1);
for nS=0:n_S-1;
tmp_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_wkS__(:,1+nS),T_k_p_wkS__(:,1+nS))/(2*pi);
tmp_TT_S_(1+nS) = tmp_TT;
end;%for nS=0:n_S-1;
tmp_ij = find((FTK.delta_x_.^2 + FTK.delta_y_.^2)<1e-6);
UX_T_l2_S_ = transpose(UX_T_l2_dS__(1+tmp_ij,:));
if (verbose); disp(sprintf(' %% tmp_TT_S_ vs UX_T_l2_S_: %0.16f',fnorm(tmp_TT_S_ - UX_T_l2_S_)/fnorm(tmp_TT_S_))); end;
flag_plot=0;
if flag_plot;
tmp_ij = find((FTK.delta_x_.^2 + FTK.delta_y_.^2)<1e-6);
subplot(1,2,1); hold on; 
plot(0:n_S-1,UX_T_l2_S_/(pi*k_p_r_max^2),'rx'); xlabel('nS'); ylabel('l2');
plot(0:n_S-1,tmp_TT_S_/(pi*k_p_r_max^2),'bo'); xlabel('nS'); ylabel('l2');
hold off;
subplot(1,2,2); plot(UX_T_l2_S_/(pi*k_p_r_max^2),tmp_TT_S_/(pi*k_p_r_max^2),'g.'); xlabel('l2_A'); ylabel('l2_B');
end;%if flag_plot;
%%%%%%%%;
tmp_t = tic();
[UX_T_k_q_wnS___,UX_T_k_p_wnS___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,n_UX_rank,n_S,svd_VUXT_lwnS____,zeros(n_S,1),zeros(n_S,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_k_q_wnS___: %0.6fs',tmp_t)); end;
UX_T_k_p_wnS__ = reshape(UX_T_k_p_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
UX_T_k_q_wnS__ = reshape(UX_T_k_q_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
%%%%%%%%;
% Visualize: ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
figure(3);clf;figbig;fig80s;
subplot(2,2,1);imagesc(reshape(permute(log10(abs(svd_VUXT_lwnS____)),[1,3,4,2]),[FTK.n_svd_l*n_UX_rank*n_S,n_w_max]));axisnotick; colorbar;
subplot(2,2,2);plot(0:n_w_max-1,log10(sqrt(mean(reshape(permute(abs(svd_VUXT_lwnS____).^2,[1,3,4,2]),[FTK.n_svd_l*n_UX_rank*n_S,n_w_max]),1))));
subplot(2,2,3);imagesc(reshape(permute(reshape(log10(abs(UX_T_k_q_wnS__)),[n_w_max,n_UX_rank,n_S]),[2,3,1]),[n_UX_rank*n_S,n_w_max]));axisnotick;colorbar;
subplot(2,2,4);plot(0:n_w_max-1,log10(sqrt(mean(reshape(permute(reshape(abs(UX_T_k_q_wnS__).^2,[n_w_max,n_UX_rank,n_S]),[2,3,1]),[n_UX_rank*n_S,n_w_max]),1))));
end;%if flag_plot;
%%%%%%%%;
% Calculate ampmh_X_wSM___8. ;
%%%%%%%%;
tmp_t = tic();
TT_k_q_ = svd(UX_T_k_q_wnS__);
n_T_rank = min(find(TT_k_q_/TT_k_q_(1)<1e-3)); if isempty(n_T_rank); n_T_rank = min(size(UX_T_k_q_wnS__)); end;
[UT_k_q__,ST_k_q__,VT_k_q__] = svds(UX_T_k_q_wnS__,n_T_rank);
if (verbose>1); disp(sprintf(' %% n_S %d --> n_T_rank %d',n_S,n_T_rank)); end;
%%%%%%%%;
tmp_ij = find((FTK.delta_x_.^2 + FTK.delta_y_.^2)<1e-6);
parameter = struct('type','parameter');
parameter.flag_optimize_over_gamma_z = 1;
parameter.flag_compute_I_value = 0;
parameter.tolerance_master = tolerance_master;
parameter.pm_n_UX_rank_use = n_UX_rank;
tmp_t = tic();
[ ...
 parameter ...
,X_TT__ ...
,delta_x_TT__ ...
,delta_y_TT__ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK ...
,n_w_max ...
,n_UX_rank ...
,n_S ...
,UX_T_k_q_wnS__ ...
,UX_T_l2_S_ ...
,n_S ...
,svd_VUXT_lwnS____ ...
,UX_T_l2_dS__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wTT___: %0.3fs',tmp_t)); end;

%%%%%%%%;
% Now let us assume that we are given an image M_{B} drawn from CTF\times template T_{B}. ;
% We can assume that M_{B} = T_{B} + iid noise, ;
% where the noise has variance \sigma^{2} per unit-area in real-space. ;
% This means that over an area of dx^{2} the noise has a variance of (\sigma^{2}/dx^{2}). ;
% Le us also assume that we have another CTF\times template T_{A}. ;
% The probability that a single component of M_{B} (say, M_{B}(j)) could be observed from T_{A}(j) is: ;
% P(M_{B}(j)|T_{A}(j)) = 1/sqrt(2\pi) \cdot 1/(\sigma/dx) \cdot \exp( -(T_{B}(j)-T_{A}(j))^{2} / (2\cdot(\sigma^{2}/dx^{2})) ). ;
%                = 1/sqrt(2\pi) \cdot dx/\sigma \cdot \exp( -(T_{B}(j)-T_{A}(j))^{2} / (2\sigma^{2}) dx^{2} ). ;
% Thus, the probability that M_{B} could be observed from T_{A}_ is: ;
% P(M_{B}|T_{A}) = \prod_{j} P(M_{B}(j)|T_{A}(j)) ;
%            = \prod_{j} 1/sqrt(2\pi) \cdot dx/\sigma \cdot \exp( -(T_{B}(j)-T_{A}(j))^{2} / (2\sigma^{2}) dx^{2} ). ;
% So the negative-log-likelihood is: ;
% -\log(P(M_{B}|T_{A})) = +0.5J\log(2\pi) + J\log(\sigma) - J\log(dx) + \frac{ \sum_{j} (T_{B}(j)-T_{A}(j))^{2} dx^{2} }{2\sigma^{2}} ;
%                       = +0.5J\log(2\pi) + J\log(\sigma) - J\log(dx) + \frac{ \int_{dx^{2}} \| T_{B}(j)-T_{A}(j) \|^{2} }{2\sigma^{2}} ;
% Now in this equation the negative-log-likelihood-ratio is: ;
% +\log(P(M_{B}|T_{A}))-\log(P(M_{B}|T_{B})) = \frac{ \int_{dx^{2}} \| T_{B}(j)-T_{A}(j) \|^{2} }{2\sigma^{2}}, ;
% where \sigma is given by M_std above, and ;
% \int_{dx^{2}} \| T_{B}(j)-T_{A}(j) \|^{2} = <T_{A},T_{A}>^2 + <T_{B},T_{B}>^2 - 2\cdot |T_{A}|\cdot |T_{B}|\cdot Correlation. ;
%%%%%%%%;
N_factor_avg = mean(N_factor_M_); 
tmp_T_l2_ = reshape(UX_T_l2_S_*N_factor_avg^2,[n_S,1]);
tmp_T_X__ = X_TT__.*(sqrt(tmp_T_l2_)*transpose(sqrt(tmp_T_l2_)));
nllr_TT__ = ( (tmp_T_l2_ + transpose(tmp_T_l2_)) - 2*tmp_T_X__ ) / (2*M_std^2) ;
flag_disp=0;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed; fig80s;
subplot(1,3,1);imagesc((tmp_T_l2_ + transpose(tmp_T_l2_))); axis image; axisnotick; colorbar; title('tmp_T_l2_','Interpreter','none');
subplot(1,3,2);imagesc(tmp_T_X__); axis image; axisnotick; colorbar; title('tmp_T_X__','Interpreter','none');
subplot(1,3,3);imagesc(nllr_TT__); axis image; axisnotick; colorbar; title('nllr_TT__','Interpreter','none');
end;%if flag_disp;

%%%%%%%%;
% Now pick out a single viewing angle nS, ;
% and visualize the maximal innerproduct with other viewing-angles. ;
%%%%%%%%;
flag_disp=0;
if flag_disp;
nS = 498; %<-- one of the equatorial viewing angles. ;
imagesc_polar_a_azimu_b_0(viewing_polar_a_all_,viewing_azimu_b_all_,nllr_TT__(:,1+nS),[0,0.25],colormap_80s,1);
end;%if flag_disp;

%%%%%%%%;
% Now look along the equator at all the templates. ;
%%%%%%%%;
flag_disp=0;
if flag_disp;
%%%%%%%%;
% First in fourier-space. ;
%%%%%%%%;
T_k_p_lim_ = max(abs(T_k_p_wkS__),[],'all')*[-1,+1]*1.25;
viewing_polar_a_mode = mode(viewing_polar_a_all_);
tmp_nS_ = efind(abs(viewing_polar_a_all_-viewing_polar_a_mode)<1e-6); n_l = numel(tmp_nS_);
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
p_row = 6; p_col = ceil(n_l/p_row); np=0;
T_k_p_lim_ = max(abs(T_k_p_wkS__),[],'all')*[-1,+1]*1.25;
for np=0:n_l-1;
tmp_nS = tmp_nS_(1+np);
T_k_p_ = T_k_p_wkS__(:,1+tmp_nS);
subplot(p_row,p_col,1+np);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_),T_k_p_lim_);
axis image; axisnotick;
title(sprintf('nS %d pi*(%0.2f,%0.2f)',tmp_nS,viewing_polar_a_all_(1+tmp_nS)/pi,viewing_azimu_b_all_(1+tmp_nS)/pi));
end;%for np=0:n_l-1;
%%%%%%%%;
viewing_polar_a_mode = mode(viewing_polar_a_all_);
tmp_nS_ = efind(abs(viewing_polar_a_all_-viewing_polar_a_mode)<1e-6); n_l = numel(tmp_nS_);
T_x_c_stack_xxw___ = zeros(n_x_u,n_x_u,n_l);
for np=0:n_l-1;
tmp_nS = tmp_nS_(1+np);
T_k_p_ = T_k_p_wkS__(:,1+tmp_nS);
T_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u^2) * n_w_sum;
T_x_c_stack_xxw___(:,:,1+np) = T_x_c_;
end;%for np=0:n_l-1;
%%%%%%%%;
% Now in real-space. ;
%%%%%%%%;
T_x_c_lim_ = max(abs(T_x_c_stack_xxw___),[],'all')*[-1,+1]*1.25;
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
p_row = 6; p_col = ceil(n_l/p_row); np=0;
for np=0:n_l-1;
subplot(p_row,p_col,1+np);
T_x_c_ = T_x_c_stack_xxw___(:,:,1+np);
%imagesc_c(n_x_u,x_c_0_,n_x_u,x_c_1_,real(T_x_c_));
imagesc(real(T_x_c_),T_x_c_lim_);
axis image; axisnotick;
title(sprintf('nS %d pi*(%0.2f,%0.2f)',tmp_nS,viewing_polar_a_all_(1+tmp_nS)/pi,viewing_azimu_b_all_(1+tmp_nS)/pi));
end;%for np=0:n_l-1;
%%%%%%%%;
end;%if flag_disp;

%%%%%%%%;
% Now use 0lsq to reconstruct the volume. ;
%%%%%%%%;
tmp_t = tic();
tmp_n_order = 5;
a_k_Y_norm_reco_yk_ = ...
cg_lsq_6( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_S ...
,S_k_p_wkS__ ...
,[] ...
,[] ...
,viewing_polar_a_all_ ...
,viewing_azimu_b_all_ ...
,zeros(n_S,1) ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_k_Y_norm_reco_yk_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ... 
  a_x_u_norm_reco_xxx_ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_norm_reco_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_x_u_norm_reco_xxx_ time %0.2fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% Now use 0lsq to reconstruct the volume, ;
% but using only templates with equatorial viewing-angles. ;
%%%%%%%%;
viewing_polar_a_mode = mode(viewing_polar_a_all_);
tmp_nS_ = efind(abs(viewing_polar_a_all_-viewing_polar_a_mode)<1e-6); n_l = numel(tmp_nS_); tmp_n_S = numel(tmp_nS_);
tmp_t = tic();
tmp_n_order = 5;
a_k_Y_norm_equa_yk_ = ...
cg_lsq_6( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,tmp_n_S ...
,S_k_p_wkS__(:,1+tmp_nS_) ...
,[] ...
,[] ...
,viewing_polar_a_all_(1+tmp_nS_) ...
,viewing_azimu_b_all_(1+tmp_nS_) ...
,zeros(tmp_n_S,1) ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_k_Y_norm_equa_yk_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ... 
  a_x_u_norm_equa_xxx_ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_norm_equa_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_x_u_norm_equa_xxx_ time %0.2fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% Now use 0lsq to reconstruct the volume, ;
% but using only templates with equatorial viewing-angles. ;
% This time we distort/squish the azimuthal viewing angles. ;
%%%%%%%%;
tmp_squi = @(azimu_b) periodize(atan2(2*sin(azimu_b),cos(azimu_b)),0,2*pi);
viewing_polar_a_mode = mode(viewing_polar_a_all_);
tmp_nS_ = efind(abs(viewing_polar_a_all_-viewing_polar_a_mode)<1e-6); n_l = numel(tmp_nS_); tmp_n_S = numel(tmp_nS_);
flag_disp=0;
if flag_disp;
tmp_azimu_b_ = sort(viewing_azimu_b_all_(1+tmp_nS_),'ascend'); 
figure(1+nf);nf=nf+1;figsml; 
plot(tmp_azimu_b_,tmp_squi(tmp_azimu_b_),'.');
xlabel('tmp_azimu_b_','Interpreter','none');
ylabel('tmp_squi(tmp_azimu_b_)','Interpreter','none');
end;%if flag_disp;
%%%%%%%%;
tmp_t = tic();
tmp_n_order = 5;
a_k_Y_norm_squi_yk_ = ...
cg_lsq_6( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,tmp_n_S ...
,S_k_p_wkS__(:,1+tmp_nS_) ...
,[] ...
,[] ...
,viewing_polar_a_all_(1+tmp_nS_) ...
,tmp_squi(viewing_azimu_b_all_(1+tmp_nS_)) ...
,zeros(tmp_n_S,1) ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_k_Y_norm_squi_yk_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ... 
  a_x_u_norm_squi_xxx_ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_norm_squi_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_x_u_norm_squi_xxx_ time %0.2fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% Now repeat with a more refined set of templates. ;
%%%%%%%%;
n_S_refine = tmp_n_S * 8;
viewing_polar_a_refine_all_ = (pi/2)*ones(n_S_refine,1);
viewing_azimu_b_refine_all_ = linspace(0,2*pi,1+n_S_refine); viewing_azimu_b_refine_all_ = transpose(viewing_azimu_b_refine_all_(1:n_S_refine));
tmp_t = tic();
[ ...
 S_refine_k_p_wkS__ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max_max ...
,n_k_p_r ...
,reshape(a_k_Y_norm_yk__,[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
,n_S_refine ...
,viewing_azimu_b_refine_all_ ...
,viewing_polar_a_refine_all_ ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% pm_template_2: %0.6fs',tmp_t)); end;
S_refine_k_p_wkS__ = reshape(S_refine_k_p_wkS__,[n_w_max*n_k_p_r,n_S_refine]);
%%%%%%%%;
% Now use 0lsq to reconstruct the volume, ;
% but using only the refined templates with equatorial viewing-angles. ;
%%%%%%%%;
tmp_t = tic();
tmp_n_order = 5;
a_k_Y_norm_refi_yk_ = ...
cg_lsq_6( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_S_refine ...
,S_refine_k_p_wkS__ ...
,[] ...
,[] ...
,viewing_polar_a_refine_all_ ...
,viewing_azimu_b_refine_all_ ...
,zeros(n_S_refine,1) ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_k_Y_norm_refi_yk_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ... 
  a_x_u_norm_refi_xxx_ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_norm_refi_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_x_u_norm_refi_xxx_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
% This time we distort/squish the refined azimuthal viewing angles. ;
%%%%%%%%;
tmp_t = tic();
tmp_n_order = 5;
a_k_Y_norm_resq_yk_ = ...
cg_lsq_6( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_S_refine ...
,S_refine_k_p_wkS__ ...
,[] ...
,[] ...
,viewing_polar_a_refine_all_ ...
,tmp_squi(viewing_azimu_b_refine_all_) ...
,zeros(n_S_refine,1) ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_k_Y_norm_resq_yk_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ... 
  a_x_u_norm_resq_xxx_ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_norm_resq_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_x_u_norm_resq_xxx_ time %0.2fs',tmp_t)); end;
%%%%%%%%;


if (verbose);
disp(sprintf(' %% a_k_Y_norm_yk_ vs a_k_Y_norm_reco_yk_: %0.16f',fnorm(a_k_Y_norm_yk_ - a_k_Y_norm_reco_yk_)/fnorm(a_k_Y_norm_yk_)));
disp(sprintf(' %% a_k_Y_norm_yk_ vs a_k_Y_norm_equa_yk_: %0.16f',fnorm(a_k_Y_norm_yk_ - a_k_Y_norm_equa_yk_)/fnorm(a_k_Y_norm_yk_)));
disp(sprintf(' %% a_k_Y_norm_yk_ vs a_k_Y_norm_refi_yk_: %0.16f',fnorm(a_k_Y_norm_yk_ - a_k_Y_norm_refi_yk_)/fnorm(a_k_Y_norm_yk_)));
disp(sprintf(' %% a_k_Y_norm_reco_yk_ vs a_k_Y_norm_equa_yk_: %0.16f',fnorm(a_k_Y_norm_reco_yk_ - a_k_Y_norm_equa_yk_)/fnorm(a_k_Y_norm_reco_yk_)));
disp(sprintf(' %% a_k_Y_norm_equa_yk_ vs a_k_Y_norm_squi_yk_: %0.16f',fnorm(a_k_Y_norm_equa_yk_ - a_k_Y_norm_squi_yk_)/fnorm(a_k_Y_norm_equa_yk_)));
disp(sprintf(' %% a_k_Y_norm_equa_yk_ vs a_k_Y_norm_resq_yk_: %0.16f',fnorm(a_k_Y_norm_equa_yk_ - a_k_Y_norm_resq_yk_)/fnorm(a_k_Y_norm_equa_yk_)));
disp(sprintf(' %% a_k_Y_norm_squi_yk_ vs a_k_Y_norm_resq_yk_: %0.16f',fnorm(a_k_Y_norm_squi_yk_ - a_k_Y_norm_resq_yk_)/fnorm(a_k_Y_norm_squi_yk_)));
end;%if (verbose);

flag_disp=0;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(2,3,1);
isosurface_f_x_u_0(a_x_u_norm_reco_xxx_,98.5);
title('a_x_u_norm_reco_xxx_','Interpreter','none');
subplot(2,3,2);
isosurface_f_x_u_0(a_x_u_norm_equa_xxx_,98.5);
title('a_x_u_norm_equa_xxx_','Interpreter','none');
subplot(2,3,3);
isosurface_f_x_u_0(a_x_u_norm_squi_xxx_,98.5);
title('a_x_u_norm_squi_xxx_','Interpreter','none');
subplot(2,3,4);
isosurface_f_x_u_0(a_x_u_norm_refi_xxx_,98.5);
title('a_x_u_norm_refi_xxx_','Interpreter','none');
subplot(2,3,5);
isosurface_f_x_u_0(a_x_u_norm_resq_xxx_,98.5);
title('a_x_u_norm_resq_xxx_','Interpreter','none');
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% Test isosurface_f_x_u_0 with a source at [+0.5,-0.1,+0.85] ;
% Note the mismatch between x_0 and x_1 in function and plot. ;
% This is because the function is defined using ndgrid, ;
% while isosurface_f_x_u_0 uses meshgrid (required for isonormals). ;
%%%%%%%%;
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_);
tmp_a_x_u_xxx_ = (abs(x_u_0___-0.50)<1e-1) & (abs(x_u_1___+0.10)<1e-1) & (abs(x_u_2___-0.85)<1e-1) ;
figure(1+nf);nf=nf+1;clf;figsml;
isosurface_f_x_u_0(tmp_a_x_u_xxx_,95);
end;%if flag_disp;

%%%%%%%%;
% test and time rotation of a_k_Y_norm_yk_. ;
% Note the transposition of x_0 and x_1 due to meshgrid in isosurface_f_x_u_0. ;
%%%%%%%%;
tmp_euler_ = pi*[+0.25 ; +0.5 ; -0.25];
tmp_t = tic();
[ ...
 b_k_Y_norm_yk_ ...
,tmp_W_beta__ ...
] = ...
rotate_spharm_to_spharm_2( ...
 verbose ...
,[] ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,a_k_Y_norm_yk_ ...
,tmp_euler_ ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% b_k_Y_norm_yk_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ... 
  b_x_u_norm_xxx_ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,b_k_Y_norm_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% b_x_u_norm_xxx_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
flag_disp=1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
isosurface_f_x_u_0(a_x_u_norm_reco_xxx_,98.5);
title('a_x_u_norm_reco_xxx_','Interpreter','none');
subplot(1,2,2);
isosurface_f_x_u_0(b_x_u_norm_xxx_,98.5);
title('b_x_u_norm_xxx_','Interpreter','none');
end;%if flag_disp;

%%%%%%%%;
% Now extract fourier components of equatorial ring of S_refine_k_p_wk_. ;
%%%%%%%%;
q_max = ceil(2*pi*k_p_r_max);
q_ = transpose(-q_max:1:+q_max);
n_q = numel(q_);
n_azimu_b = n_S_refine;
azimu_b_ = viewing_azimu_b_refine_all_;
dazimu_b = mean(diff(azimu_b_));
F_Sq__ = exp(+i*azimu_b_*transpose(q_));
F_inv_qS__ = ctranspose(F_Sq__)/n_azimu_b;
S_refine_wkq__ = S_refine_k_p_wkS__*transpose(F_inv_qS__);
%plot(q_,log10(sum(abs(S_refine_wkq__).^2,1)),'.'); xlim(q_max*[-1,+1]); set(gca,'XTick',q_); grid on;

%%%%%%%%;
% Do the same with the dominant principal-component. ;
%%%%%%%%;
tmp_A_Swk__ = transpose(bsxfun(@times,S_refine_k_p_wkS__,reshape(sqrt(weight_2d_k_all_),[n_w_sum,1])));
R2_SS__ = bsxfun(@plus,reshape(sum(abs(tmp_A_Swk__).^2,2),[n_S_refine,1]),reshape(sum(abs(tmp_A_Swk__).^2,2),[1,n_S_refine])) - 2*real(tmp_A_Swk__*ctranspose(tmp_A_Swk__));
%[tmp_B_Swk__,s_w_] = MSA_shape_speed1_1(n_w_sum,n_S_refine,[],tmp_A_Swk__,weight_2d_k_all_);
n_rank = 2;
[tmp_U_Sr__,tmp_S_rr__,tmp_V_wkr__] = svds(tmp_A_Swk__,n_rank);
tmp_US_Sr__ = tmp_A_Swk__*tmp_V_wkr__;
tmp_US_qr__ = F_inv_qS__*tmp_US_Sr__;
%plot(q_,log10(sum(abs(tmp_US_qr__).^2,2)),'.'); xlim(ceil(q_max/4)*[-1,+1]); set(gca,'XTick',q_); grid on;

%%%%%%%%;
% Now reduce the winding number from 2 to 1. ;
%%%%%%%%;
tmp_US_half_qr__ = zeros(n_q,n_rank);
S_refine_half_wkq__ = zeros(n_w_sum,n_q);
for nq1=0:n_q-1;
q1_val = q_(1+nq1);
if mod(q1_val,2)==0;
q2_val = q1_val/2;
nq2 = q2_val+q_max;
tmp_US_half_qr__(1+nq2,:) = tmp_US_qr__(1+nq1,:);
S_refine_half_wkq__(:,1+nq2) = S_refine_wkq__(:,1+nq1);
end;%if mod(q1_val,2)==0;
end;%for nq1=0:n_q-1;
tmp_US_half_Sr__ = F_Sq__*tmp_US_half_qr__;
S_refine_half_wkS__ = S_refine_half_wkq__*transpose(F_Sq__);
%%%%%%%%;
flag_disp=0;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;colormap('hsv');
subplot(2,2,1); surfline_0(real(tmp_US_Sr__(:,1)),imag(tmp_US_Sr__(:,1)),azimu_b_); title('orig');
subplot(2,2,2); surfline_0(real(tmp_US_half_Sr__(:,1)),imag(tmp_US_half_Sr__(:,1)),azimu_b_); title('half');
subplot(2,2,3); plot(q_,log10(sum(abs(tmp_US_qr__).^2,2)),'.'); xlim(ceil(q_max/4)*[-1,+1]); set(gca,'XTick',q_); grid on; title('orig');
subplot(2,2,4); plot(q_,log10(sum(abs(tmp_US_half_qr__).^2,2)),'.'); xlim(ceil(q_max/4)*[-1,+1]); set(gca,'XTick',q_); grid on; title('half');
end;%if flag_disp;

%%%%%%%%;
tmp_A_Swk__ = transpose(bsxfun(@times,S_refine_k_p_wkS__,reshape(sqrt(weight_2d_k_all_),[n_w_sum,1])));
tmp_A_half_Swk__ = transpose(bsxfun(@times,S_refine_half_wkS__,reshape(sqrt(weight_2d_k_all_),[n_w_sum,1])));
R2_SS__ = bsxfun(@plus,reshape(sum(abs(tmp_A_half_Swk__).^2,2),[n_S_refine,1]),reshape(sum(abs(tmp_A_Swk__).^2,2),[1,n_S_refine])) - 2*real(tmp_A_half_Swk__*ctranspose(tmp_A_Swk__));
%%%%%%%%;

%%%%%%%%;
% Now reconstruct halved volume. ;
%%%%%%%%;
tmp_t = tic();
tmp_n_order = 5;
a_k_Y_norm_refi_half_yk_ = ...
cg_lsq_6( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_S_refine ...
,S_refine_half_wkS__ ...
,[] ...
,[] ...
,viewing_polar_a_refine_all_ ...
,viewing_azimu_b_refine_all_ ...
,zeros(n_S_refine,1) ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_k_Y_norm_refi_yk_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ... 
  a_x_u_norm_refi_half_xxx_ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_norm_refi_half_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_x_u_norm_refi_xxx_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
flag_disp=1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
isosurface_f_x_u_0(a_x_u_norm_refi_xxx_,99.0);
title('a_x_u_norm_refi_xxx_','Interpreter','none');
subplot(1,2,2);
isosurface_f_x_u_0(a_x_u_norm_refi_half_xxx_,99.0);
title('a_x_u_norm_refi_half_xxx_','Interpreter','none');
end;%if flag_disp;
%%%%%%%%

%%%%%%%%;
% Now generate all templates from halved volume. ;
%%%%%%%%;
a_k_Y_norm_refi_half_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_norm_refi_half_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_norm_refi_half_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
template_k_eq_d = -1;
viewing_k_eq_d = 1.0; %<-- subsample just a few templates for visualization. ;
tmp_t = tic();
[ ...
 S_half_k_p_wkS__ ...
,~ ...
,n_S ...
,template_viewing_azimu_b_all_ ...
,template_viewing_polar_a_all_ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max_max ...
,n_k_p_r ...
,reshape(a_k_Y_norm_refi_half_yk__,[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% pm_template_2: %0.6fs',tmp_t)); end;
S_half_k_p_wkS__ = reshape(S_half_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
S_half_k_q_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_half_k_q_wkS__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_half_k_p_wkS__(:,1+nS));
end;%for nS=0:n_S-1;
%%%%%%%%;
% Now measure differences between template stacks. ;
%%%%%%%%;
tmp_A_Swk__ = transpose(bsxfun(@times,S_k_p_wkS__,reshape(sqrt(weight_2d_k_all_),[n_w_sum,1])));
tmp_A_half_Swk__ = transpose(bsxfun(@times,S_half_k_p_wkS__,reshape(sqrt(weight_2d_k_all_),[n_w_sum,1])));
R2_SS__ = bsxfun(@plus,reshape(sum(abs(tmp_A_half_Swk__).^2,2),[n_S,1]),reshape(sum(abs(tmp_A_Swk__).^2,2),[1,n_S])) - 2*real(tmp_A_half_Swk__*ctranspose(tmp_A_Swk__));

%%%%%%%%;
% Now do a more comprehensive comparison using ampmh_X_wSM___8. ;
% Use one of the CTF clusters for now. ;
%%%%%%%%;
ncluster = 0; %<-- pick one of the CTF clusters. ;
UX_kn__ = UX_knc___(:,:,1+ncluster);
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
T_k_p_wkS__ = reshape(bsxfun(@times,reshape(S_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_k_q_wkS__ = reshape(bsxfun(@times,reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_half_k_p_wkS__ = reshape(bsxfun(@times,reshape(S_half_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_half_k_q_wkS__ = reshape(bsxfun(@times,reshape(S_half_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
tmp_t = tic();
svd_VUXT_lwnS____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_S,T_k_q_wkS__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXT_lwnS____: %0.3fs',tmp_t)); end;
tmp_t = tic();
svd_VUXT_half_lwnS____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_S,T_half_k_q_wkS__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXT_half_lwnS____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the templates. ;
%%%%%%%%;
tmp_t = tic();
UX_T_l2_dS__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_S,n_UX_rank,svd_VUXT_lwnS____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_l2_dS__: %0.3fs',tmp_t)); end;
if (verbose); disp(sprintf(' %% average l2-norm of templates: %0.16f',mean(UX_T_l2_dS__(:))/(pi*k_p_r_max^2))); end;
tmp_TT_S_ = zeros(n_S,1);
for nS=0:n_S-1;
tmp_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_wkS__(:,1+nS),T_k_p_wkS__(:,1+nS))/(2*pi);
tmp_TT_S_(1+nS) = tmp_TT;
end;%for nS=0:n_S-1;
tmp_ij = find((FTK.delta_x_.^2 + FTK.delta_y_.^2)<1e-6);
UX_T_l2_S_ = transpose(UX_T_l2_dS__(1+tmp_ij,:));
if (verbose); disp(sprintf(' %% tmp_TT_S_ vs UX_T_l2_S_: %0.16f',fnorm(tmp_TT_S_ - UX_T_l2_S_)/fnorm(tmp_TT_S_))); end;
flag_plot=0;
if flag_plot;
tmp_ij = find((FTK.delta_x_.^2 + FTK.delta_y_.^2)<1e-6);
subplot(1,2,1); hold on; 
plot(0:n_S-1,UX_T_l2_S_/(pi*k_p_r_max^2),'rx'); xlabel('nS'); ylabel('l2');
plot(0:n_S-1,tmp_TT_S_/(pi*k_p_r_max^2),'bo'); xlabel('nS'); ylabel('l2');
hold off;
subplot(1,2,2); plot(UX_T_l2_S_/(pi*k_p_r_max^2),tmp_TT_S_/(pi*k_p_r_max^2),'g.'); xlabel('l2_A'); ylabel('l2_B');
end;%if flag_plot;
%%%%%%%%;
% repeat for half templates. ;
%%%%%%%%;
tmp_t = tic();
UX_T_half_l2_dS__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_S,n_UX_rank,svd_VUXT_half_lwnS____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_half_l2_dS__: %0.3fs',tmp_t)); end;
if (verbose); disp(sprintf(' %% average l2-norm of templates: %0.16f',mean(UX_T_half_l2_dS__(:))/(pi*k_p_r_max^2))); end;
tmp_TT_half_S_ = zeros(n_S,1);
for nS=0:n_S-1;
tmp_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_half_k_p_wkS__(:,1+nS),T_half_k_p_wkS__(:,1+nS))/(2*pi);
tmp_TT_half_S_(1+nS) = tmp_TT;
end;%for nS=0:n_S-1;
tmp_ij = find((FTK.delta_x_.^2 + FTK.delta_y_.^2)<1e-6);
UX_T_half_l2_S_ = transpose(UX_T_half_l2_dS__(1+tmp_ij,:));
if (verbose); disp(sprintf(' %% tmp_TT_half_S_ vs UX_T_half_l2_S_: %0.16f',fnorm(tmp_TT_half_S_ - UX_T_half_l2_S_)/fnorm(tmp_TT_half_S_))); end;
flag_plot=0;
if flag_plot;
tmp_ij = find((FTK.delta_x_.^2 + FTK.delta_y_.^2)<1e-6);
subplot(1,2,1); hold on; 
plot(0:n_S-1,UX_T_half_l2_S_/(pi*k_p_r_max^2),'rx'); xlabel('nS'); ylabel('l2');
plot(0:n_S-1,tmp_TT_half_S_/(pi*k_p_r_max^2),'bo'); xlabel('nS'); ylabel('l2');
hold off;
subplot(1,2,2); plot(UX_T_half_l2_S_/(pi*k_p_r_max^2),tmp_TT_half_S_/(pi*k_p_r_max^2),'g.'); xlabel('l2_A'); ylabel('l2_B');
end;%if flag_plot;
%%%%%%%%;
tmp_t = tic();
[UX_T_k_q_wnS___,UX_T_k_p_wnS___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,n_UX_rank,n_S,svd_VUXT_lwnS____,zeros(n_S,1),zeros(n_S,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_k_q_wnS___: %0.6fs',tmp_t)); end;
UX_T_k_p_wnS__ = reshape(UX_T_k_p_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
UX_T_k_q_wnS__ = reshape(UX_T_k_q_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
%%%%%%%%;
tmp_t = tic();
[UX_T_half_k_q_wnS___,UX_T_half_k_p_wnS___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,n_UX_rank,n_S,svd_VUXT_half_lwnS____,zeros(n_S,1),zeros(n_S,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_half_k_q_wnS___: %0.6fs',tmp_t)); end;
UX_T_half_k_p_wnS__ = reshape(UX_T_half_k_p_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
UX_T_half_k_q_wnS__ = reshape(UX_T_half_k_q_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
%%%%%%%%;
% Visualize: ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
figure(3);clf;figbig;fig80s;
subplot(2,2,1);imagesc(reshape(permute(log10(abs(svd_VUXT_lwnS____)),[1,3,4,2]),[FTK.n_svd_l*n_UX_rank*n_S,n_w_max]));axisnotick; colorbar;
subplot(2,2,2);plot(0:n_w_max-1,log10(sqrt(mean(reshape(permute(abs(svd_VUXT_lwnS____).^2,[1,3,4,2]),[FTK.n_svd_l*n_UX_rank*n_S,n_w_max]),1))));
subplot(2,2,3);imagesc(reshape(permute(reshape(log10(abs(UX_T_k_q_wnS__)),[n_w_max,n_UX_rank,n_S]),[2,3,1]),[n_UX_rank*n_S,n_w_max]));axisnotick;colorbar;
subplot(2,2,4);plot(0:n_w_max-1,log10(sqrt(mean(reshape(permute(reshape(abs(UX_T_k_q_wnS__).^2,[n_w_max,n_UX_rank,n_S]),[2,3,1]),[n_UX_rank*n_S,n_w_max]),1))));
end;%if flag_plot;
%%%%%%%%;
% Calculate ampmh_X_wSM___8. ;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_optimize_over_gamma_z = 1;
parameter.flag_compute_I_value = 0;
parameter.tolerance_master = tolerance_master;
parameter.pm_n_UX_rank_use = n_UX_rank;
tmp_t = tic();
[ ...
 parameter ...
,X_TT_half__ ...
,delta_x_TT_half__ ...
,delta_y_TT_half__ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK ...
,n_w_max ...
,n_UX_rank ...
,n_S ...
,UX_T_k_q_wnS__ ...
,UX_T_l2_S_ ...
,n_S ...
,svd_VUXT_half_lwnS____ ...
,UX_T_half_l2_dS__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wTT_half___: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now visualize where the good fits occur. ;
%%%%%%%%;
N_factor_avg = mean(N_factor_M_); 
tmp_T_l2_ = reshape(UX_T_l2_S_*N_factor_avg^2,[n_S,1]);
%tmp_T_X__ = X_TT__.*(sqrt(tmp_T_l2_)*transpose(sqrt(tmp_T_l2_)));
tmp_T_half_l2_ = reshape(UX_T_half_l2_S_*N_factor_avg^2,[n_S,1]);
tmp_X_TT_half__ = X_TT_half__.*(sqrt(tmp_T_l2_)*transpose(sqrt(tmp_T_half_l2_)));
nllr_TT_half__ = ( (tmp_T_l2_ + transpose(tmp_T_half_l2_)) - 2*tmp_X_TT_half__ ) / (2*M_std^2) ;
%%%%%%%%;
flag_disp=1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
imagesc_polar_a_azimu_b_0(viewing_polar_a_all_,viewing_azimu_b_all_,min(nllr_TT_half__,[],2),[0,0.0625],colormap_80s,1);
subplot(1,2,2);
imagesc_polar_a_azimu_b_0(viewing_polar_a_all_,viewing_azimu_b_all_,max(X_TT_half__,[],2),[0.5,1.0],colormap_80s,1);
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% Verdict so far: ;
% Given the signal-to-noise ratio for this molecule, ;
% as well as the N_factor_avg (from earlier). ;
% a nllr of 0.01 or less is achived for about 0.25 of the viewing-angles ;
% (i.e., those closest to the equator). ;
%%%%;
% Now to reject a null-hypothesis one might want a p-value of 0.05 or less. ;
% This means that we can only accumulate an nllr of about ;
% -log(0.05) ~ 3 across the images from the 'half' volume. ;
%%%%;
% Now we estimate how many templates can be associated with the half-volume. ;
% (This assumes one image per template). ;
%%%%%%%%;
tmp_nllr_S_ = min(nllr_TT_half__,[],2);
[tmp_nllr_sort_S_,tmp_ij_] = sort(tmp_nllr_S_,'ascend');
flag_disp=0; if flag_disp; plot(cumsum(tmp_nllr_sort_S_)); end;
tmp_n_S_use = numel(efind(cumsum(sort(tmp_nllr_S_))<-log(0.05))); %<-- 392 templates. ;
%%%%%%%%;
% Now use the remaining templates to form a new volume. ;
% We will create a new template stack that amplifies the values near the poles. ;
% This uses template_k_c_2__ from get_template_1. ;
%%%%%%%%;
tmp_ij_rem_ = tmp_ij_(tmp_n_S_use+1:end); %<-- these have nonequatorial viewing-angle. ;
template_k_r01_wkS__ = sqrt(template_k_c_0__.^2 + template_k_c_1__.^2);
template_k_r012_wkS__ = sqrt(template_k_c_0__.^2 + template_k_c_1__.^2 + template_k_c_2__.^2);
template_azimu_b_wkS__ = atan2(template_k_c_1__,template_k_c_0__);
template_polar_a_wkS__ = atan2(template_k_r01_wkS__,template_k_c_2__);
S_ampl_k_p_wkS__ = S_k_p_wkS__;
S_ampl_k_p_wkS__ = S_ampl_k_p_wkS__.*(abs(template_polar_a_wkS__-pi/2)<=0.15*pi); %<-- zero out the poles. ;
%%%%%%%%;
flag_disp=0; if flag_disp; plot(viewing_polar_a_all_(tmp_ij_rem_),'.'); end;
%%%%%%%%;
tmp_t = tic();
tmp_n_order = 5;
a_k_Y_norm_rema_yk_ = ...
cg_lsq_6( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_S ...
,S_ampl_k_p_wkS__ ...
,[] ...
,[] ...
,viewing_polar_a_all_ ...
,viewing_azimu_b_all_ ...
,zeros(n_S,1) ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_k_Y_norm_refi_yk_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ... 
  a_x_u_norm_rema_xxx_ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_norm_rema_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_x_u_norm_refi_xxx_ time %0.2fs',tmp_t)); end;
%%%%%%%%;

flag_disp=1;
if flag_disp;
prct_ = [98,98.5,99.0,99.5]; n_prct = numel(prct_);
%prct_ = [99.0,99.25,99.5,99.75]; n_prct = numel(prct_);
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 3; p_col = n_prct;
for pcol=0:p_col-1;
tmp_prct = prct_(1+pcol);
subplot(p_row,p_col,1+pcol+0*p_col);
isosurface_f_x_u_0(a_x_u_norm_refi_xxx_,tmp_prct);
title(sprintf('a_x_u_norm_refi_xxx_: %.1f%%',tmp_prct),'Interpreter','none');
subplot(p_row,p_col,1+pcol+1*p_col);
isosurface_f_x_u_0(a_x_u_norm_rema_xxx_,tmp_prct);
title(sprintf('a_x_u_norm_rema_xxx_: %.1f%%',tmp_prct),'Interpreter','none');
subplot(p_row,p_col,1+pcol+2*p_col);
isosurface_f_x_u_0(a_x_u_norm_refi_half_xxx_,tmp_prct);
title(sprintf('a_x_u_norm_refi_half_xxx_: %.1f%%',tmp_prct),'Interpreter','none');
end;%for pcol=0:p_col-1;
fname_fig_pre = sprintf('%s_jpg/test_heterogeneity_spurious_tprv1c_FIGA',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
end;%if flag_disp;
%%%%%%%%



