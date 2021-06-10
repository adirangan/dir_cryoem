% See if we can reconstruct a function using the various \psi-averaged moments of its templates. ;
% Does not consider translations or ctf. ;

clear;

%platform = 'access1';
platform = 'OptiPlex';
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;

dir_trunk = sprintf('/%s/rangan/dir_cryoem/dir_rangan_playroom/dir_principled_marching',string_root);

verbose=0;

%%%%%%%%;
% generate spherical-shell. ;
%%%%%%%%;
tmp_t=tic();
k_p_r_max = 48/(2*pi);
k_eq_d = 1.0/(2*pi)*sqrt(0.5);
[ ...
n_k_all ...
,azimu_b_all_ ...
,polar_a_all_ ...
,weight_3d_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,n_polar_a ...
,polar_a_ ...
,n_azimu_b_ ...
] = ...
sample_shell_5( ...
 k_p_r_max ...
,k_eq_d ...
,'L' ...
);
tmp_t=toc(tmp_t); if (verbose>-1); disp(sprintf(' %% sample_shell_5: %0.6fs',tmp_t)); end;
n_k_p_r = 1;
k_p_r_ = k_p_r_max;
weight_3d_k_p_r_ = (1/3)*k_p_r_max^3; %<-- sum(weight_3d_k_p_r_)*(4*pi) = (4/3)*pi*k_p_r_max^3 --> sum(weight_3d_k_p_r_) = (1/3)*k_p_r_max^3 ;
% note that sum(weight_3d_all_) = (4*pi)*k_p_r_max^2 ;

%%%%%%%%;
% set up spherical-harmonic coefficients for an iid-gaussian function. ;
%%%%%%%%;
tmp_t=tic();
sigma_init = 1.0d0;
%l_max_upb = 36;
%l_max_upb = 8;
l_max_upb = 2; %<-- given the O(l_max_upb.^8) scaling, we should be very careful here. ;
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
l_max_(1+nk_p_r) = max(0,min(l_max_upb,1+ceil(2*pi*k_p_r_(1+nk_p_r))));
end;%for nk_p_r=0:n_k_p_r-1;
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
Y_l_val_ = zeros(n_lm_sum,1);
Y_m_val_ = zeros(n_lm_sum,1);
Y_k_val_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
tmp_l_val_ = zeros(n_lm_(1+nk_p_r),1);
tmp_m_val_ = zeros(n_lm_(1+nk_p_r),1);
na=0; 
for l_val=0:l_max;
for m_val=-l_val:+l_val;
tmp_l_val_(1+na) = l_val;
tmp_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
tmp_ij_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
Y_l_val_(1+tmp_ij_) = tmp_l_val_;
Y_m_val_(1+tmp_ij_) = tmp_m_val_;
Y_k_val_(1+tmp_ij_) = k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
weight_Y_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
weight_Y_(1+tmp_ij_) = weight_3d_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;

fname_mat = sprintf('%s_mat/test_kam_0.mat',dir_trunk);
if ( exist(fname_mat,'file')); disp(sprintf(' %% %s found, not creating',fname_mat)); end;
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));

a_k_Y_true_ = sigma_init*crandn(n_lm_sum,1);

%%%%%%%%;
% generate templates. ;
%%%%%%%%;
tmp_t=tic();
%template_k_eq_d = 0.25*(2*pi*k_p_r_max) / (2*(l_max_max+1)) ; %<-- this only needs to be sufficiently large to sample each oscillation twice (i.e., nyquist). ;
template_k_eq_d = 0;
n_w_0in_ = 2*(l_max_max+1) + 6;
viewing_k_eq_d = k_eq_d*8;
[ ...
 S_k_p_true__ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,template_k_c_0__ ...
,template_k_c_1__ ...
,template_k_c_2__ ...
] = ...
get_template_1( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_true_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in_ ...
);
n_S = n_viewing_all; n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
tmp_t=toc(tmp_t); if (verbose>-1); disp(sprintf(' %% get_template_0: %0.6fs',tmp_t)); end;
%%%%%%%%;
S_k_q_true__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_q_true__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_true__(:,1+nS)); 
end;%for nS=0:n_S-1;

%%%%%%%%;
% Now reconstruct a_k_Y_0lsq_. ;
%%%%%%%%;
tmp_t=tic();
lsq_n_order = 5;
euler_polar_a_ = viewing_polar_a_all_ ;
euler_azimu_b_ = viewing_azimu_b_all_ ;
euler_gamma_z_ = zeros(n_viewing_all,1);
a_k_Y_0lsq_ = cg_lsq_1(lsq_n_order,n_k_p_r,l_max_,n_w_,n_S,S_k_p_true__,euler_polar_a_,euler_azimu_b_,euler_gamma_z_);
tmp_t=toc(tmp_t); if (verbose>1); disp(sprintf(' %% cg_lsq_1: %0.6fs',tmp_t)); end;
disp(sprintf(' %% a_k_Y_true_ vs a_k_Y_0lsq_: %0.16f',fnorm(a_k_Y_true_-a_k_Y_0lsq_)/fnorm(a_k_Y_true_)));

%%%%%%%%;
% Now create (temporary) array of templates for each basis function. ;
%%%%%%%%;
S_k_q_qlmS___ = zeros(n_w_sum,n_lm_sum,n_S);
for nlm_sum=0:n_lm_sum-1;
tmp_a_k_Y_ = zeros(n_lm_sum,1); tmp_a_k_Y_(1+nlm_sum)=1;
[ ...
 tmp_S_k_p__ ...
] = ...
get_template_1( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,tmp_a_k_Y_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in_ ...
);
tmp_S_k_q__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
tmp_S_k_q__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_S_k_p__(:,1+nS)); 
end;%for nS=0:n_S-1;
S_k_q_qlmS___(:,1+nlm_sum,:) = reshape(tmp_S_k_q__,[n_w_sum,1,n_S]);
clear tmp_S_k_q__ tmp_S_k_p__ tmp_a_k_Y_ ;
end;%for nlm_sum=0:n_lm_sum-1;

%%%%%%%%;
% Now construct the templates for each of the various basis volumes. ;
%%%%%%%%;
q_ = periodize(transpose(0:n_w_sum-1),-n_w_sum/2,+n_w_sum/2);
n_gamma = 1+2*l_max_max;
gamma_ = transpose(linspace(0,2*pi,n_gamma+1)); gamma_ = periodize(gamma_(1:n_gamma),-pi,+pi);
expiqg_qg__ = zeros(n_w_sum,n_gamma); expiqg_qg__ = exp(-i*q_*transpose(gamma_));
%for nq=0:n_w_sum-1;
%q = q_(1+nq);
%for ngamma=0:n_gamma-1;
%gamma = gamma_(1+ngamma);
%expiqg_qg__(1+nq,1+ngamma) = exp(-i*q*gamma);
%end;%for ngamma=0:n_gamma-1;
%end;%for nq=0:n_w_sum-1;
%%%%%%%%;

%%%%%%%%;
% test out moment coefficient construction, limiting the number of templates for now. ;
%%%%%%%%;
n_S = 7; S_k_q_qlmS___ = S_k_q_qlmS___(:,:,1:n_S);
%%%%%%%%;
moment_0_lm_ = zeros(n_lm_sum,1);
moment_1_glmlm___ = zeros(n_gamma,n_lm_sum,n_lm_sum);
moment_2_gglmlmlm_____ = zeros(n_gamma,n_gamma,n_lm_sum,n_lm_sum,n_lm_sum);
for nS=0:n_S-1;
for nlm0=0:n_lm_sum-1;
moment_0_lm_(1+nlm0) = moment_0_lm_(1+nlm0) + S_k_q_qlmS___(1+0,1+nlm0,1+nS);
for nlm1=0:n_lm_sum-1;
for ngamma1=0:n_gamma-1;
for nq0=0:n_w_sum-1;
nq1 = nq0;
moment_1_glmlm___(1+ngamma1,1+nlm0,1+nlm1) = moment_1_glmlm___(1+ngamma1,1+nlm0,1+nlm1) + conj(S_k_q_qlmS___(1+nq0,1+nlm0,1+nS))*S_k_q_qlmS___(1+nq1,1+nlm1,1+nS)*expiqg_qg__(1+nq1,1+ngamma1);
for nlm2=0:n_lm_sum-1;
for ngamma2=0:n_gamma-1;
for nq1=0:n_w_sum-1;
nq2 = nq0 - nq1; if (nq2<0); nq2 = nq2 + n_w_sum; end;
moment_2_gglmlmlm_____(1+ngamma1,1+ngamma2,1+nlm0,1+nlm1,1+nlm2) = moment_2_gglmlmlm_____(1+ngamma1,1+ngamma2,1+nlm0,1+nlm1,1+nlm2) + conj(S_k_q_qlmS___(1+nq0,1+nlm0,1+nS))*S_k_q_qlmS___(1+nq1,1+nlm1,1+nS)*S_k_q_qlmS___(1+nq2,1+nlm2,1+nS)*expiqg_qg__(1+nq1,1+ngamma1)*expiqg_qg__(1+nq2,1+ngamma2);
end;%for nq1=0:n_w_sum-1;
end;%for ngamma2=0:n_gamma-1;
end;%for nlm2=0:n_lm_sum-1;
end;%for nq0=0:n_w_sum-1;
end;%for ngamma1=0:n_gamma-1;
end;%for nlm1=0:n_lm_sum-1;
end;%for nlm0=0:n_lm_sum-1;
end;%for nS=0:n_S-1;
moment_0_lm_ = moment_0_lm_/n_S;
moment_1_glmlm___ =  moment_1_glmlm___/n_S;
moment_2_gglmlmlm_____ =  moment_2_gglmlmlm_____/n_S;
%%%%%%%%;
noment_0_lm_ = reshape(mean(S_k_q_qlmS___(1+0,:,:),3),[n_lm_sum,1]);
noment_1_glmlm___ = zeros(n_gamma,n_lm_sum,n_lm_sum);
for ngamma1=0:n_gamma-1;
expiqgS_qlmS___ = reshape(sparse(1:n_w_sum,1:n_w_sum,expiqg_qg__(:,1+ngamma1),n_w_sum,n_w_sum)*reshape(S_k_q_qlmS___,[n_w_sum,n_lm_sum*n_S]),[n_w_sum,n_lm_sum,n_S]);
%noment_1_glmlm___(1+ngamma1,:,:) = reshape(conj(reshape(permute(S_k_q_qlmS___,[2,1,3]),[n_lm_sum,n_w_sum*n_S]))*reshape(permute(reshape(sparse(1:n_w_sum,1:n_w_sum,expiqg_qg__(:,1+ngamma1),n_w_sum,n_w_sum)*reshape(S_k_q_qlmS___,[n_w_sum,n_lm_sum*n_S]),[n_w_sum,n_lm_sum,n_S]),[1,3,2]),[n_w_sum*n_S,n_lm_sum]),[1,n_lm_sum,n_lm_sum]);
noment_1_glmlm___(1+ngamma1,:,:) = reshape(conj(reshape(permute(S_k_q_qlmS___,[2,1,3]),[n_lm_sum,n_w_sum*n_S]))*reshape(permute(expiqgS_qlmS___,[1,3,2]),[n_w_sum*n_S,n_lm_sum]),[1,n_lm_sum,n_lm_sum]);
end;%for ngamma1=0:n_gamma-1;
noment_1_glmlm___ = noment_1_glmlm___/n_S;
noment_2_gglmlmlm_____ = zeros(n_gamma,n_gamma,n_lm_sum,n_lm_sum,n_lm_sum);
for nS=0:n_S-1;
for ngamma1=0:n_gamma-1;
expiqg1S_qlm__ = sparse(1:n_w_sum,1:n_w_sum,expiqg_qg__(:,1+ngamma1),n_w_sum,n_w_sum)*reshape(S_k_q_qlmS___(:,:,1+nS),[n_w_sum,n_lm_sum*1]);
for ngamma2=0:n_gamma-1;
expiqg2S_qlm__ = sparse(1:n_w_sum,1:n_w_sum,expiqg_qg__(:,1+ngamma2),n_w_sum,n_w_sum)*reshape(S_k_q_qlmS___(:,:,1+nS),[n_w_sum,n_lm_sum*1]);
for nq1=0:n_w_sum-1;
for nq2=0:n_w_sum-1;
nq0 = nq1 + nq2; if (nq0>n_w_sum-1); nq0 = nq0-n_w_sum; end;
expiqg1Sg2S_lmlm__ = reshape(expiqg1S_qlm__(1+nq1,:),[n_lm_sum,1])*reshape(expiqg2S_qlm__(1+nq2,:),[1,n_lm_sum]);
noment_2_gglmlmlm_____(1+ngamma1,1+ngamma2,:,:,:) = noment_2_gglmlmlm_____(1+ngamma1,1+ngamma2,:,:,:) + reshape(reshape(conj(S_k_q_qlmS___(1+nq0,:,1+nS)),[n_lm_sum,1])*reshape(expiqg1Sg2S_lmlm__,[1,n_lm_sum^2]),[1,1,n_lm_sum,n_lm_sum,n_lm_sum]);
end;%for nq2=0:n_w_sum-1;
end;%for nq1=0:n_w_sum-1;
end;%for ngamma2=0:n_gamma-1;
end;%for ngamma1=0:n_gamma-1;
end;%for nS=0:n_S-1;
noment_2_gglmlmlm_____ = noment_2_gglmlmlm_____/n_S;
%%%%%%%%;
disp(sprintf(' %% moment_0_lm_ vs noment_0_lm_: %0.16f',fnorm(moment_0_lm_-noment_0_lm_)/fnorm(moment_0_lm_)));
disp(sprintf(' %% moment_1_glmlm___ vs noment_1_glmlm___: %0.16f',fnorm(moment_1_glmlm___-noment_1_glmlm___)/fnorm(moment_1_glmlm___)));
disp(sprintf(' %% moment_2_gglmlmlm_____ vs noment_2_gglmlmlm_____: %0.16f',fnorm(moment_2_gglmlmlm_____-noment_2_gglmlmlm_____)/fnorm(moment_2_gglmlmlm_____)));
%%%%%%%%;

%%%%%%%%;
% Now test out moment calculation, using the coefficients from above. ;
%%%%%%%%;
S_true_moment_0 = 0;
for nS=0:n_S-1;
S_true_moment_0 = S_true_moment_0 + S_k_q_true__(1+0,1+nS);
end;%for nS=0:n_S-1;
S_true_moment_0 = S_true_moment_0/n_S;
S_func_moment_0 = 0;
for nlm0=0:n_lm_sum-1;
S_func_moment_0 = S_func_moment_0 + a_k_Y_true_(1+nlm0)*moment_0_lm_(1+nlm0);
end;%for nlm0=0:n_lm_sum-1;
disp(sprintf(' %% S_true_moment_0 vs S_func_moment_0: %0.16f',fnorm(S_true_moment_0-S_func_moment_0)/fnorm(S_true_moment_0)));
%%%%%%%%;
S_true_moment_1_g_ = zeros(n_gamma,1);
for nS=0:n_S-1;
for ngamma1=0:n_gamma-1;
S_true_moment_1_g_(1+ngamma1) = S_true_moment_1_g_(1+ngamma1) + sum(conj(S_k_q_true__(:,1+nS)).*S_k_q_true__(:,1+nS).*expiqg_qg__(:,1+ngamma1));
end;%for ngamma1=0:n_gamma-1;
end;%for nS=0:n_S-1;
S_true_moment_1_g_ = S_true_moment_1_g_/n_S;
S_func_moment_1_g_ = zeros(n_gamma,1);
for ngamma1=0:n_gamma-1;
S_func_moment_1_g_(1+ngamma1) = S_func_moment_1_g_(1+ngamma1) + ctranspose(a_k_Y_true_)*reshape(moment_1_glmlm___(1+ngamma1,:,:),[n_lm_sum,n_lm_sum])*a_k_Y_true_;
end;%for ngamma1=0:n_gamma-1;
disp(sprintf(' %% S_true_moment_1_g_ vs S_func_moment_1_g_: %0.16f',fnorm(S_true_moment_1_g_-S_func_moment_1_g_)/fnorm(S_true_moment_1_g_)));
%%%%%%%%;
S_true_moment_2_gg__ = zeros(n_gamma,n_gamma);
for nS=0:n_S-1;
for ngamma1=0:n_gamma-1;
for ngamma2=0:n_gamma-1;
for nq1=0:n_w_sum-1;
for nq2=0:n_w_sum-1;
nq0 = nq1 + nq2; if (nq0>n_w_sum-1); nq0 = nq0-n_w_sum; end;
S_true_moment_2_gg__(1+ngamma1,1+ngamma2) = S_true_moment_2_gg__(1+ngamma1,1+ngamma2) + conj(S_k_q_true__(1+nq0,1+nS))*S_k_q_true__(1+nq1,1+nS)*S_k_q_true__(1+nq2,1+nS)*expiqg_qg__(1+nq1,1+ngamma1)*expiqg_qg__(1+nq2,1+ngamma2);
end;%for nq2=0:n_w_sum-1;
end;%for nq1=0:n_w_sum-1;
end;%for ngamma2=0:n_gamma-1;
end;%for ngamma1=0:n_gamma-1;
end;%for nS=0:n_S-1;
S_true_moment_2_gg__ = S_true_moment_2_gg__/n_S;
S_func_moment_2_gg__ = zeros(n_gamma,n_gamma);
for ngamma1=0:n_gamma-1;
for ngamma2=0:n_gamma-1;
for nlm0=0:n_lm_sum-1;
S_func_moment_2_gg__(1+ngamma1,1+ngamma2) = S_func_moment_2_gg__(1+ngamma1,1+ngamma2) + conj(a_k_Y_true_(1+nlm0))* ( transpose(a_k_Y_true_)*reshape(moment_2_gglmlmlm_____(1+ngamma1,1+ngamma2,1+nlm0,:,:),[n_lm_sum,n_lm_sum])*a_k_Y_true_ );
end;%for nlm0=0:n_lm_sum-1;
end;%for ngamma2=0:n_gamma-1;
end;%for ngamma1=0:n_gamma-1;
disp(sprintf(' %% S_true_moment_2_gg__ vs S_func_moment_2_gg__: %0.16f',fnorm(S_true_moment_2_gg__-S_func_moment_2_gg__)/fnorm(S_true_moment_2_gg__)));



%save(fname_mat);
end;%if (~exist(fname_mat,'file'));

