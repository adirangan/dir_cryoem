function ...
[ ...
 parameter ...
,X_SM__ ...
,delta_x_SM__ ...
,delta_y_SM__ ...
,gamma_z_SM__ ...
,I_value_SM__ ...
,CTF_UX_S_l2_SM__ ...
,UX_M_l2_SM__ ...
,CTF_UX_S_l2_Sc__ ...
,index_nM_from_ncluster__ ...
,n_index_nM_from_ncluster_ ...
] = ...
ampmh_X_cluster_wrap_SM__10( ...
 parameter ...
,FTK ...
,n_w_max ...
,n_k_p_r ...
,n_S ...
,S_k_q_wkS__ ...
,n_CTF ...
,CTF_k_p_r_kC__ ...
,index_nCTF_from_nM_ ...
,n_M ...
,n_cluster ...
,index_ncluster_from_nCTF_ ...
,pm_n_UX_rank_c_ ...
,UX_knc___ ...
,X_weight_rc__ ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
,CTF_k_p_r_xavg_kc__ ...
,CTF_UX_S_k_q_wnSc___ ...
,CTF_UX_S_l2_Sc__ ...
,index_ncluster_from_nM_ ...
,index_nM_from_ncluster__ ...
,n_index_nM_from_ncluster_ ...
);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if nargin<1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% First define integral of <f,f>. ;
%%%%%%%%;
h2d_ = @(kd) 4*pi^2*(besselj(0,kd) + besselj(2,kd)); % calculates <f_j,f_k>, normalized so that <f,f> = (4*pi^2);
dh2d_ = @(kd) 4*pi^3*(besselj(-1,kd) - besselj(+3,kd));
h3d_ = @(kd) 4*pi*( sin(kd) - (kd).*cos(kd) ) ./ kd.^3 ; % calculates <f_j,f_k>, normalized so that <f,f> = 4*pi/3;
dh3d_ = @(kd) 12*pi*( (kd.^2/3 - 1) .* sin(kd) + (kd).*cos(kd) ) ./ kd.^4 ;
%%%%%%%%;
verbose=1;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/max(1e-12,k_p_r_max); TorL = 'L';
if (verbose); disp(sprintf(' %% [testing ampmh_X_cluster_wrap_SM__10.m]')); end;
%%%%%%%%;
tmp_t = tic();
[ ...
 n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,J_node_ ...
,J_weight_ ...
,J_chebfun_ ...
,J_polyval_ ...
] = ...
sample_sphere_7( ...
 verbose ...
,k_p_r_max ...
,k_eq_d ...
,TorL ...
) ;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% sample_sphere_7: %0.2fs',tmp_t)); end;
if (verbose); disp(sprintf(' %% k_p_r_max %0.2f k_eq_d %0.2f n_k_all %d n_k_p_r %d',k_p_r_max,k_eq_d,n_k_all,n_k_p_r)); end;
%%%%%%%%;
l_max_upb = 36;
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
%%%%%%%%;
delta_orig_ = [+0.12;-0.3;+0.23];
a_k_p_orig_ = exp(+2*pi*i*(k_c_0_all_*delta_orig_(1+0) + k_c_1_all_*delta_orig_(1+1) + k_c_2_all_*delta_orig_(1+2)));
tmp_t = tic;
[a_k_Y_quad_] = ...
convert_k_p_to_spharm_1( ...
 verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_p_orig_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
delta_orig_r012 = sqrt(delta_orig_(1+0)^2 + delta_orig_(1+1)^2 + delta_orig_(1+2)^2);
delta_orig_r01 = sqrt(delta_orig_(1+0)^2 + delta_orig_(1+1)^2);
delta_orig_polar_a = atan2(delta_orig_r01,delta_orig_(1+2));
delta_orig_azimu_b = atan2(delta_orig_(1+1),delta_orig_(1+0));
delta_Ylm_ = get_Ylm__(1+l_max_max,0:l_max_max,1,delta_orig_azimu_b,delta_orig_polar_a);
a_k_Y_form_ = zeros(n_lm_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
l_max = l_max_(1+nk_p_r);
for l_val=0:l_max;
tmp_x = 2*pi*k_p_r*delta_orig_r012;
tmp_jl = besselj(l_val+0.5,tmp_x)*sqrt(pi/(2*tmp_x));
for m_val=-l_val:+l_val;
a_k_Y_form_(1+na) = 4*pi * i^l_val * tmp_jl * conj(delta_Ylm_{1+l_val}(1+m_val+l_val));
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
disp(sprintf(' %% a_k_Y_form_ vs a_k_Y_quad_: %0.16f',fnorm(a_k_Y_form_ - a_k_Y_quad_)/fnorm(a_k_Y_form_)));
%%%%%%%%;
n_w_max = 2*(l_max_max+1); n_w_ = n_w_max*ones(n_k_p_r,1);
template_k_eq_d = -1;
viewing_k_eq_d = k_eq_d*128;
tmp_t = tic();
[ ...
 S_k_p__ ...
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
,a_k_Y_quad_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% get_template_1: %0.2fs',tmp_t)); end;
if (verbose); disp(sprintf(' %% n_viewing_all %d n_viewing_polar_a %d n_w_max %d',n_viewing_all,n_viewing_polar_a,max(n_w_))); end;
n_S = n_viewing_all; n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
a_k_Y_quad__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_quad__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_quad_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_t = tic();
[ ...
 tmp_S_k_p__ ...
,~ ...
,tmp_n_S ...
,~ ...
,~ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max ...
,n_k_p_r ...
,reshape(a_k_Y_quad__,[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,-1 ...
,n_w_max ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t)); end;
tmp_S_k_p__ = reshape(tmp_S_k_p__,[n_w_max*n_k_p_r,tmp_n_S]);
if (verbose); disp(sprintf(' %% S_k_p__ vs tmp_S_k_p__: %0.16f',fnorm(S_k_p__ - tmp_S_k_p__)/fnorm(S_k_p__))); end;
%%%%%%%%;

pole_k_c_0_ = zeros(n_w_sum,1);
pole_k_c_1_ = zeros(n_w_sum,1);
pole_k_c_2_ = zeros(n_w_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
for nw=0:n_w-1;
gamma_z = 2*pi*nw/n_w;
cc = cos(gamma_z); sc = sin(gamma_z);
pole_k_c_0_(1+na) = k_p_r*cc;
pole_k_c_1_(1+na) = k_p_r*sc;
pole_k_c_2_(1+na) = 0;
na=na+1;
end;%for nw=0:n_w-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
% after calling get_template_1: ;
% A template with viewing angle viewing_polar_a and viewing_azimu_b corresponds to the evaluations: ;
% let sa and ca be sin(polar_a) and cos(polar_a), respectively. ;
% let sb and cb be sin(azimu_b) and cos(azimu_b), respectively. ;
% let sc and cc be sin(gamma_z) and cos(gamma_z), respectively. ;
% [ template_k_c_0 ; template_k_c_1 ; template_k_c_2 ] = [ +cb*ca*cc - sb*sc ; +sb*ca*cc + cb*sc ; -sa*cc ] ;
% for gamma_z = 2*pi*[0:n_gamma_z-1]/n_gamma_z. ;
% Given that the original function is a plane-wave defined as: ;
% a_k_p_ = exp(+2*pi*i*( delta_orig_(1+0)*template_k_c_0 + delta_orig_(1+1)*template_k_c_1 + delta_orig_(1+2)*template_k_c_2 )) ;
% we have that the template evaluates to: ;
% S_k_p_ = exp(+2*pi*i*( delta_orig_ * Rz(azimu_b) * Ry(polar_a) * Rz(gamma_z) * [1;0;0]*k_p_r )) ;
% S_k_p_ = exp(+2*pi*i*( (Ry(-polar_a) * Rz(-azimu_b) * delta_orig_) * Rz(gamma_z) * [1;0;0]*k_p_r )) ;
%%%%%%%%;
n_test = 8;
for ntest=0:n_test-1;
nS=max(0,min(n_S-1,floor(n_S*rand())));
disp(sprintf(' %% ntest %d/%d --> nS %d',ntest,n_test,nS));
S_k_p_quad_ = S_k_p__(:,1+nS);
S_k_p_orig_ = exp(+2*pi*i*(template_k_c_0__(:,1+nS)*delta_orig_(1+0) + template_k_c_1__(:,1+nS)*delta_orig_(1+1) + template_k_c_2__(:,1+nS)*delta_orig_(1+2)));
viewing_azimu_b = viewing_azimu_b_all_(1+nS);
cb = cos(+viewing_azimu_b); sb = sin(+viewing_azimu_b);
Rz = [ +cb , -sb , 0 ; +sb , +cb , 0 ; 0 , 0 , 1 ]; %<-- rotation about the positive z-axis. ;
viewing_polar_a = viewing_polar_a_all_(1+nS);
ca = cos(+viewing_polar_a); sa = sin(+viewing_polar_a);
Ry = [ +ca , 0 , +sa ; 0 , 1 , 0 ; -sa , 0 , +ca ]; %<-- rotation about the positive y-axis. ;
delta_temp_ = transpose(Ry)*transpose(Rz)*delta_orig_;
S_k_p_form_ = exp(+2*pi*i*(pole_k_c_0_*delta_temp_(1+0) + pole_k_c_1_*delta_temp_(1+1) + pole_k_c_2_*delta_temp_(1+2)));
flag_plot=0;
if flag_plot;
figure(1);clf;
subplot(2,3,1+0); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_orig_),[-1,+1],colormap_beach()); title('real orig'); axisnotick; axis image;
subplot(2,3,1+3); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(S_k_p_orig_),[-1,+1],colormap_beach()); title('imag orig'); axisnotick; axis image;
subplot(2,3,2+0); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_quad_),[-1,+1],colormap_beach()); title('real quad'); axisnotick; axis image;
subplot(2,3,2+3); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(S_k_p_quad_),[-1,+1],colormap_beach()); title('imag quad'); axisnotick; axis image;
subplot(2,3,3+0); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_form_),[-1,+1],colormap_beach()); title('real form'); axisnotick; axis image;
subplot(2,3,3+3); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(S_k_p_form_),[-1,+1],colormap_beach()); title('imag form'); axisnotick; axis image;
figbig;
end;%if flag_plot;
error_orig_vs_quad = fnorm(S_k_p_orig_ - S_k_p_quad_)/fnorm(S_k_p_orig_);
disp(sprintf(' %% error_orig_vs_quad: %0.16f',error_orig_vs_quad));
error_orig_vs_form = fnorm(S_k_p_orig_ - S_k_p_form_)/fnorm(S_k_p_orig_);
disp(sprintf(' %% error_orig_vs_form: %0.16f',error_orig_vs_form));
end;%for ntest=0:n_test-1;
%%%%%%%%;

tmp_t = tic();
S_k_q__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_p_ = S_k_p__(:,1+nS);
S_k_q__(:,1+nS) = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,S_k_p_ ...
);
end;%for nS=0:n_S-1;
S_k_p_form__ = zeros(n_w_sum,n_S);
S_k_q_form__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
viewing_azimu_b = viewing_azimu_b_all_(1+nS);
cb = cos(+viewing_azimu_b); sb = sin(+viewing_azimu_b);
Rz = [ +cb , -sb , 0 ; +sb , +cb , 0 ; 0 , 0 , 1 ]; %<-- rotation about the positive z-axis. ;
viewing_polar_a = viewing_polar_a_all_(1+nS);
ca = cos(+viewing_polar_a); sa = sin(+viewing_polar_a);
Ry = [ +ca , 0 , +sa ; 0 , 1 , 0 ; -sa , 0 , +ca ]; %<-- rotation about the positive y-axis. ;
delta_temp_ = transpose(Ry)*transpose(Rz)*delta_orig_;
S_k_p_form_ = exp(+2*pi*i*(pole_k_c_0_*delta_temp_(1+0) + pole_k_c_1_*delta_temp_(1+1) + pole_k_c_2_*delta_temp_(1+2)));
S_k_p_form__(:,1+nS) = S_k_p_form_;
S_k_q_form_ = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,S_k_p_form_ ...
);
S_k_q_form__(:,1+nS) = S_k_q_form_;
end;%for nS=0:n_S-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% S_k_q_form__: %0.3fs',tmp_t)); end;
if (verbose); disp(sprintf(' %% S_k_q_form__ vs S_k_q__: %0.16f',fnorm(S_k_q_form__ - S_k_q__)/fnorm(S_k_q_form__))); end;
if (verbose); disp(sprintf(' %% S_k_p_form__ vs S_k_p__: %0.16f',fnorm(S_k_p_form__ - S_k_p__)/fnorm(S_k_p_form__))); end;

%%%%%%%%;
delta_r_max = 0.05; n_delta_v_requested = 64;
%delta_r_max = 0.00; n_delta_v_requested = 1;
delta_r_p = 0.05;
delta_r_s = delta_r_max/sqrt(2*log(1/delta_r_p));
delta_r_N = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2));
svd_eps = 1e-12;
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
disp(sprintf(' %% p-val %0.4f delta_r_max %0.6f sigma %0.4f N_pixel %0.4f --> FTK.n_svd_l %d, n_delta_v_requested %d',delta_r_p,delta_r_max,delta_r_s,delta_r_N,FTK.n_svd_l,n_delta_v_requested));
%%%%%%%%;
% design index_nCTF_from_nM_. ;
%%%%%%%%;

n_M = n_S - 1; %<-- just to check dimensions. ;
n_CTF = 2;
index_nCTF_from_nM_ = mod(0:n_M-1,n_CTF);
CTF_k_p_r_kC__ = zeros(n_k_p_r,n_CTF);
%CTF_k_p_r_kC__ = 1 + rand(n_k_p_r,n_CTF); %<-- invertible. ;
for nCTF=0:n_CTF-1;
CTF_k_p_r_kC__(:,1+nCTF) = 1 + sin(2*pi*(1+nCTF/n_CTF)*k_p_r_/k_p_r_max); %<-- invertible. ;
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
%  Construct CTF of same size as images. ;
%%%%%%%%;
CTF_k_p_wkC__ = zeros(n_w_sum,n_CTF);
for nCTF=0:n_CTF-1;
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
tmp_index_ = n_w_csum_(1+nk_p_r) + [0:n_w-1];
CTF_k_p_wkC__(1+tmp_index_,1+nCTF) = CTF_k_p_r_kC__(1+nk_p_r,1+nCTF);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nCTF=0:n_CTF-1;
CTF_inv_k_p_wkC__ = 1./max(1e-12,CTF_k_p_wkC__);
%%%%%%%%;
% Generate images and templates to use. ;
%%%%%%%%;
S_k_p_0use__ = S_k_p_form__;
S_k_q_0use__ = S_k_q_form__;
M_k_p_0use__ = S_k_p_form__(:,1:n_M);
M_k_q_0use__ = S_k_q_form__(:,1:n_M);
for nM=0:n_M-1;
nCTF = index_nCTF_from_nM_(1+nM);
M_k_p_0use__(:,1+nM) = M_k_p_0use__(:,1+nM).*CTF_inv_k_p_wkC__(:,1+nCTF);
M_k_q_0use__(:,1+nM) = M_k_q_0use__(:,1+nM).*CTF_inv_k_p_wkC__(:,1+nCTF);
end;%for nM=0:n_M-1;
%%%%%%%%;
n_cluster = n_CTF;
index_ncluster_from_nCTF_ = transpose(0:n_CTF-1); %<-- give each nCTF its own cluster. ;
index_ncluster_from_nM_ = index_ncluster_from_nCTF_(1+index_nCTF_from_nM_);
index_nM_from_ncluster__ = cell(n_cluster,1);
n_index_nM_from_ncluster_ = zeros(n_cluster,1);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster__{1+ncluster} = efind(index_ncluster_from_nM_==ncluster);
n_index_nM_from_ncluster_(1+ncluster) = numel(index_nM_from_ncluster__{1+ncluster});
end;%for ncluster=0:n_cluster-1;
%%%%%%%%;
% design principal-modes. ;
%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions.; 
UX_knc___ = zeros(n_k_p_r,n_UX_rank,n_cluster);
X_weight_rc__ = zeros(n_k_p_r,n_cluster);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
tmp_n_M = numel(index_nM_from_ncluster_);
assert(tmp_n_M==n_index_nM_from_ncluster_(1+ncluster));
tmp_M_k_p_0use__ = M_k_p_0use__(:,1+index_nM_from_ncluster_);
[ ...
 X_kk__ ...
,X_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,tmp_n_M ...
,tmp_M_k_p_0use__ ...
);
[tmp_UX_kn__,tmp_SX_k__,tmp_VX_kn__] = svds(X_kk__,n_UX_rank);
UX_knc___(:,1:n_UX_rank,1+ncluster) = tmp_UX_kn__(:,1:n_UX_rank);
X_weight_rc__(:,1+ncluster) = X_weight_r_;
end;%for ncluster=0:n_cluster-1;
pm_n_UX_rank_c_ = n_UX_rank*ones(n_cluster,1) - 1 - mod(transpose(0:n_cluster-1),2); %<-- just to check dimension. ;
pm_n_UX_rank_max = max(pm_n_UX_rank_c_);
%%%%%%%%;
% Now calculate the inner-products. ;
%%%%%%%%;
tmp_t = tic();
svd_VUXM_lwnM____ = zeros(FTK.n_svd_l,n_w_max,pm_n_UX_rank_max,n_M);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
tmp_n_M = numel(index_nM_from_ncluster_);
assert(tmp_n_M==n_index_nM_from_ncluster_(1+ncluster));
tmp_M_k_q_0use__ = M_k_q_0use__(:,1+index_nM_from_ncluster_);
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_weight_r_ = X_weight_rc__(:,1+ncluster);
svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+index_nM_from_ncluster_) = ...
tpmh_VUXM_lwnM____3( ...
 FTK ...
,n_k_p_r ...
,n_w_ ...
,tmp_n_M ...
,tmp_M_k_q_0use__ ...
,pm_n_UX_rank ...
,UX_kn__ ...
,X_weight_r_ ...
);
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the non-translated images. ;
%%%%%%%%;
tmp_t = tic();
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M,pm_n_UX_rank_max,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
disp(sprintf(' %% average l2-norm of images (should be 1 if CTF is unity): %0.16f',mean(UX_M_l2_dM__(:))/(pi*k_p_r_max^2)));
flag_plot=0;
if flag_plot;
plot(UX_M_l2_dM__(:)/(pi*k_p_r_max^2),'.');
end;%if flag_plot;
%%%%%%%%;
tmp_t = tic();
[UX_M_k_q_wnM___,UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank_max,n_M,svd_VUXM_lwnM____,zeros(n_M,1),zeros(n_M,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_k_q_wnM___: %0.6fs',tmp_t)); end;
UX_M_k_p_wnM__ = reshape(UX_M_k_p_wnM___(:,1:pm_n_UX_rank_max,:),[n_w_max*pm_n_UX_rank_max,n_M]);
UX_M_k_q_wnM__ = reshape(UX_M_k_q_wnM___(:,1:pm_n_UX_rank_max,:),[n_w_max*pm_n_UX_rank_max,n_M]);
%%%%%%%%;
% Visualize: ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
figure(3);clf;figbig;fig80s;
subplot(2,2,1);imagesc(reshape(permute(log10(abs(svd_VUXM_lwnM____)),[1,3,4,2]),[FTK.n_svd_l*pm_n_UX_rank_max*n_M,n_w_max]));axisnotick; colorbar;
subplot(2,2,2);plot(0:n_w_max-1,log10(sqrt(mean(reshape(permute(abs(svd_VUXM_lwnM____).^2,[1,3,4,2]),[FTK.n_svd_l*pm_n_UX_rank_max*n_M,n_w_max]),1))));
subplot(2,2,3);imagesc(reshape(permute(reshape(log10(abs(UX_M_k_q_wnM__)),[n_w_max,pm_n_UX_rank_max,n_M]),[2,3,1]),[pm_n_UX_rank_max*n_M,n_w_max]));axisnotick;colorbar;
subplot(2,2,4);plot(0:n_w_max-1,log10(sqrt(mean(reshape(permute(reshape(abs(UX_M_k_q_wnM__).^2,[n_w_max,pm_n_UX_rank_max,n_M]),[2,3,1]),[pm_n_UX_rank_max*n_M,n_w_max]),1))));
end;%if flag_plot;
%%%%%%%%;
% Calculate norm of templates. ;
%%%%%%%%;
S_k_q_wSk___ = permute(reshape(S_k_q_0use__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
tmp_n_M = numel(index_nM_from_ncluster_); assert(tmp_n_M==n_index_nM_from_ncluster_(1+ncluster));
%%%%%%%%;
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
%pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
%pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_); pm_n_lm_max = max(pm_n_lm_);
pm_weight_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
%%%%%%%%;
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_weight_r_ = X_weight_rc__(:,1+ncluster);
CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),2);
CTF_UX_S_k_q_wnS__(1:pm_n_w_sum,1:n_S) = reshape(permute(reshape(reshape(S_k_q_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_weight_r_.*CTF_k_p_r_xavg_k_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
CTF_UX_S_l2_S_ = sum(abs(CTF_UX_S_k_q_wnS__).^2,1)/n_w_max;
end;%for ncluster=0:n_cluster-1;
%%%%%%%%;
disp(sprintf(' %% average l2-norm of templates (should be 1 if CTF is unity): %0.16f',mean(CTF_UX_S_l2_S_(:))/(pi*k_p_r_max^2)));
%%%%%%%%;
% Calculate ampmh_X_cluster_wrap_SM__10. ;
%%%%%%%%;
parameter = struct('type','parameter');
tmp_t = tic();
[ ...
 parameter ...
,X_SM_calc__ ...
,delta_x_SM_calc__ ...
,delta_y_SM_calc__ ...
,gamma_z_SM_calc__ ...
,I_value_SM_calc__ ...
] = ...
ampmh_X_cluster_wrap_SM__10( ...
 parameter ...
,FTK ...
,n_w_max ...
,n_k_p_r ...
,n_S ...
,S_k_q_0use__ ...
,n_CTF ...
,CTF_k_p_r_kC__ ...
,index_nCTF_from_nM_ ...
,n_M ...
,n_cluster ...
,index_ncluster_from_nCTF_ ...
,pm_n_UX_rank_c_ ...
,UX_knc___ ...
,X_weight_rc__ ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM_calc__: %0.3fs',tmp_t)); end;

%%%%%%%%;
% Calculate true landscape of innerproducts for the same set of translations. ;
%%%%%%%%;
X_SM_form__ = zeros(n_viewing_all,n_M);
X_SM_quad__ = zeros(n_viewing_all,n_M);
for nS=0:n_viewing_all-1;
viewing_azimu_b = viewing_azimu_b_all_(1+nS);
cb = cos(+viewing_azimu_b); sb = sin(+viewing_azimu_b);
Rz = [ +cb , -sb , 0 ; +sb , +cb , 0 ; 0 , 0 , 1 ]; %<-- rotation about the positive z-axis. ;
viewing_polar_a = viewing_polar_a_all_(1+nS);
ca = cos(+viewing_polar_a); sa = sin(+viewing_polar_a);
Ry = [ +ca , 0 , +sa ; 0 , 1 , 0 ; -sa , 0 , +ca ]; %<-- rotation about the positive y-axis. ;
delta_temp_0_ = transpose(Ry)*transpose(Rz)*delta_orig_;
for nM=0:n_M-1;
nCTF = index_nCTF_from_nM_(1+nM);
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+nCTF);
CTF_inv_k_p_wk_ = CTF_inv_k_p_wkC__(:,1+nCTF);
image_viewing_azimu_b = viewing_azimu_b_all_(1+nM);
cb = cos(+image_viewing_azimu_b); sb = sin(+image_viewing_azimu_b);
Rz = [ +cb , -sb , 0 ; +sb , +cb , 0 ; 0 , 0 , 1 ]; %<-- rotation about the positive z-axis. ;
image_viewing_polar_a = viewing_polar_a_all_(1+nM);
ca = cos(+image_viewing_polar_a); sa = sin(+image_viewing_polar_a);
Ry = [ +ca , 0 , +sa ; 0 , 1 , 0 ; -sa , 0 , +ca ]; %<-- rotation about the positive y-axis. ;
delta_temp_1_ = transpose(Ry)*transpose(Rz)*delta_orig_;
%%%%;
gamma_z = gamma_z_SM_calc__(1+nS,1+nM);%gamma_z = 2*pi*nw/n_w_max;
cc = cos(+gamma_z); sc = sin(+gamma_z);
Rz = [+cc , -sc ; +sc , +cc];
delta_x = delta_x_SM_calc__(1+nS,1+nM);
delta_y = delta_y_SM_calc__(1+nS,1+nM);
delta_temp_0b_ = Rz*delta_temp_0_(1+[0:1]); %<-- rotate delta_temp_0_ by +gamma_z = rotate k by -gamma_z = rotate S_0 by +gamma_z. ;
delta_temp_1b_ = delta_temp_1_(1+[0:1]) - [delta_x;delta_y]; %<-- translate delta_temp_1_ by -delta_ = multiply S_1 by exp(-2*pi*i*dot(k_,delta_)). ;
%%%%;
S_k_p_form_ = exp(+2*pi*i*(pole_k_c_0_*delta_temp_0b_(1+0) + pole_k_c_1_*delta_temp_0b_(1+1)));
S_k_p_0use_ = S_k_p_form_.*CTF_k_p_wk_;
M_k_p_form_ = exp(+2*pi*i*(pole_k_c_0_*delta_temp_1b_(1+0) + pole_k_c_1_*delta_temp_1b_(1+1)));
M_k_p_0use_ = M_k_p_form_.*CTF_inv_k_p_wk_;
%%%%;
S_l2_quad = sqrt(sum(reshape(conj(S_k_p_0use_).*S_k_p_0use_,[n_w_max,n_k_p_r])*weight_2d_k_p_r_)/n_w_max);
M_l2_quad = sqrt(sum(reshape(conj(M_k_p_0use_).*M_k_p_0use_,[n_w_max,n_k_p_r])*weight_2d_k_p_r_)/n_w_max);
%%%%;
X_form = h2d_(2*pi*k_p_r_max*fnorm(delta_temp_0b_ - delta_temp_1b_))/(2*pi)^2 * (pi*k_p_r_max^2); %<-- note sign of translation. ;
X_SM_form__(1+nS,1+nM) = real(X_form) / (S_l2_quad*M_l2_quad) ;
%%%%;
X_quad = sum(reshape(conj(S_k_p_0use_).*M_k_p_0use_,[n_w_max,n_k_p_r])*weight_2d_k_p_r_)/n_w_max;
X_SM_quad__(1+nS,1+nM) = real(X_quad) / (S_l2_quad*M_l2_quad) ;
%%%%;
end;%for nM=0:n_M-1;
end;%for nS=0:n_viewing_all-1;
%%%%%%%%;
disp(sprintf(' %% X_SM_form__ vs X_SM_calc__: %0.16f',fnorm(X_SM_form__-X_SM_calc__)/fnorm(X_SM_form__)));
disp(sprintf(' %% X_SM_quad__ vs X_SM_calc__: %0.16f',fnorm(X_SM_quad__-X_SM_calc__)/fnorm(X_SM_quad__)));
%%%%%%%%;

disp('returning'); return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if nargin<1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

verbose = 0;
if (verbose>0); disp(sprintf(' %% [entering ampmh_X_cluster_wrap_SM__10]')); end;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); FTK=[]; end; na=na+1;
if (nargin<1+na); n_w_max=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); S_k_q_wkS__=[]; end; na=na+1;
if (nargin<1+na); n_CTF=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_r_kC__=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); n_cluster=[]; end; na=na+1;
if (nargin<1+na); index_ncluster_from_nCTF_=[]; end; na=na+1;
if (nargin<1+na); pm_n_UX_rank_c_=[]; end; na=na+1;
if (nargin<1+na); UX_knc___=[]; end; na=na+1;
if (nargin<1+na); X_weight_rc__=[]; end; na=na+1;
if (nargin<1+na); svd_VUXM_lwnM____=[]; end; na=na+1;
if (nargin<1+na); UX_M_l2_dM__=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_r_xavg_kc__=[]; end; na=na+1;
if (nargin<1+na); CTF_UX_S_k_q_wnSc___=[]; end; na=na+1;
if (nargin<1+na); CTF_UX_S_l2_Sc__=[]; end; na=na+1;
if (nargin<1+na); index_ncluster_from_nM_=[]; end; na=na+1;
if (nargin<1+na); index_nM_from_ncluster__=[]; end; na=na+1;
if (nargin<1+na); n_index_nM_from_ncluster_=[]; end; na=na+1;

%%%%%%%%;
if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;

%%%%%%%%;
% group images by cluster. ;
%%%%%%%%;
n_cluster = 1+max(index_ncluster_from_nCTF_);
%%%%;
if  isempty(index_ncluster_from_nM_);
index_ncluster_from_nM_ = index_ncluster_from_nCTF_(1+index_nCTF_from_nM_);
end;%if  isempty(index_ncluster_from_nM_);
%%%%;
if ( isempty(index_nM_from_ncluster__) | isempty(n_index_nM_from_ncluster_) );
index_nM_from_ncluster__ = cell(n_cluster,1);
n_index_nM_from_ncluster_ = zeros(n_cluster,1);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster__{1+ncluster} = efind(index_ncluster_from_nM_==ncluster);
n_index_nM_from_ncluster_(1+ncluster) = numel(index_nM_from_ncluster__{1+ncluster});
end;%for ncluster=0:n_cluster-1;
end;%if ( isempty(index_nM_from_ncluster__) | isempty(n_index_nM_from_ncluster_) );
%%%%;

X_SM__ = zeros(n_S,n_M);
delta_x_SM__ = zeros(n_S,n_M);
delta_y_SM__ = zeros(n_S,n_M);
gamma_z_SM__ = zeros(n_S,n_M);
I_value_SM__ = zeros(n_S,n_M);
CTF_UX_S_l2_SM__ = zeros(n_S,n_M);
UX_M_l2_SM__ = zeros(n_S,n_M);

if  isempty(CTF_UX_S_k_q_wnSc___);
S_k_q_wSk___ = permute(reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
end;%if  isempty(CTF_UX_S_k_q_wnSc___);

for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
tmp_n_M = numel(index_nM_from_ncluster_); assert(tmp_n_M==n_index_nM_from_ncluster_(1+ncluster));
if (tmp_n_M>0);
if (verbose); disp(sprintf(' %% ncluster %d/%d: tmp_n_M %d',ncluster,n_cluster,tmp_n_M)); end;
%%%%%%%%;
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
%pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
%pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_); pm_n_lm_max = max(pm_n_lm_);
pm_weight_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
%%%%%%%%;
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_weight_r_ = X_weight_rc__(:,1+ncluster);
%%%%%%%%;
if ~isempty(CTF_k_p_r_xavg_kc__);
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
end;%if ~isempty(CTF_k_p_r_xavg_kc__);
if  isempty(CTF_k_p_r_xavg_kc__);
CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),2);
end;%if  isempty(CTF_k_p_r_xavg_kc__);
%%%%;
if ~isempty(CTF_UX_S_k_q_wnSc___);
CTF_UX_S_k_q_wnS__ = CTF_UX_S_k_q_wnSc___(1:pm_n_w_sum,1:n_S,1+ncluster);
end;%if ~isempty(CTF_UX_S_k_q_wnSc___);
if  isempty(CTF_UX_S_k_q_wnSc___);
CTF_UX_S_k_q_wnS__(1:pm_n_w_sum,1:n_S) = reshape(permute(reshape(reshape(S_k_q_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_weight_r_.*CTF_k_p_r_xavg_k_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
end;%if  isempty(CTF_UX_S_k_q_wnSc___);
%%%%;
if ~isempty(CTF_UX_S_l2_Sc__);
CTF_UX_S_l2_S_ = CTF_UX_S_l2_Sc__(:,1+ncluster);
end;%if ~isempty(CTF_UX_S_l2_Sc__);
if  isempty(CTF_UX_S_l2_Sc__);
CTF_UX_S_l2_S_ = sum(abs(CTF_UX_S_k_q_wnS__).^2,1)/n_w_max;
end;%if  isempty(CTF_UX_S_l2_Sc__);
%%%%%%%%;
tmp_t = tic();
[ ...
 parameter ...
,X_SM__(:,1+index_nM_from_ncluster_) ...
,delta_x_SM__(:,1+index_nM_from_ncluster_) ...
,delta_y_SM__(:,1+index_nM_from_ncluster_) ...
,gamma_z_SM__(:,1+index_nM_from_ncluster_) ...
,I_value_SM__(:,1+index_nM_from_ncluster_) ...
] = ...
ampmh_X_single_cluster_SM__10( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,reshape(CTF_UX_S_k_q_wnS__(1:pm_n_w_sum,1:n_S),[pm_n_w_max,pm_n_UX_rank,n_S]) ...
,CTF_UX_S_l2_S_ ...
,tmp_n_M ...
,svd_VUXM_lwnM____(:,1:n_w_max,1:pm_n_UX_rank,1+index_nM_from_ncluster_) ...
,UX_M_l2_dM__(:,1+index_nM_from_ncluster_) ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% ampmh_X_single_cluster_SM__10: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_single_cluster_SM__10',tmp_t);
%%%%%%%%;
for tmp_nM=0:tmp_n_M-1;
nM = index_nM_from_ncluster_(1+tmp_nM);
for nS=0:n_S-1;
CTF_UX_S_l2_SM__(1+nS,1+nM) = CTF_UX_S_l2_S_(1+nS);
tmp_delta_x = delta_x_SM__(1+nS,1+nM);
tmp_delta_y = delta_y_SM__(1+nS,1+nM);
tmp_nd = efind( (FTK.delta_x_==tmp_delta_x) & (FTK.delta_y_==tmp_delta_y) );
UX_M_l2_SM__(1+nS,1+nM) = UX_M_l2_dM__(1+tmp_nd,1+nM);
end;%for nS=0:n_S-1;
end;%for tmp_nM=0:tmp_n_M-1;
%%%%%%%%;
flag_check=0;
if (flag_check);
n_test=9;
for ntest=0:n_test-1;
nl = max(0,min(tmp_n_M-1,floor(tmp_n_M*ntest/(n_test-1))));
tmp_nM = index_nM_from_ncluster_(1+nl);
disp(sprintf(' %% ntest %d/%d: nl %d tmp_nM %d',ntest,n_test,nl,tmp_nM));
tmp_parameter = struct('type','parameter');
tmp_parameter.flag_optimize_over_gamma_z = 1;
[ ...
 tmp_parameter ...
,tmp_X_Sm_ ...
,tmp_delta_x_Sm_ ...
,tmp_delta_y_Sm_ ...
,tmp_gamma_z_Sm_ ...
,tmp_I_value_Sm_ ...
] = ...
ampmh_X_wSM___8( ...
 tmp_parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_S_ ...
,1 ...
,svd_VUXM_lwnM____(:,1:n_w_max,1:pm_n_UX_rank,1+tmp_nM) ...
,UX_M_l2_dM__(:,1+tmp_nM) ...
);
disp(sprintf(' %% tmp_X_Sm_ vs X_SM__(:,1+tmp_nM): %0.16f',fnorm(tmp_X_Sm_ - X_SM__(:,1+tmp_nM))/fnorm(tmp_X_Sm_)));
disp(sprintf(' %% tmp_delta_x_Sm_ vs delta_x_SM__(:,1+tmp_nM): %0.16f',fnorm(tmp_delta_x_Sm_ - delta_x_SM__(:,1+tmp_nM))/fnorm(tmp_delta_x_Sm_)));
disp(sprintf(' %% tmp_delta_y_Sm_ vs delta_y_SM__(:,1+tmp_nM): %0.16f',fnorm(tmp_delta_y_Sm_ - delta_y_SM__(:,1+tmp_nM))/fnorm(tmp_delta_y_Sm_)));
disp(sprintf(' %% tmp_gamma_z_Sm_ vs gamma_z_SM__(:,1+tmp_nM): %0.16f',fnorm(tmp_gamma_z_Sm_ - gamma_z_SM__(:,1+tmp_nM))/fnorm(tmp_gamma_z_Sm_)));
disp(sprintf(' %% tmp_I_value_Sm_ vs I_value_SM__(:,1+tmp_nM): %0.16f',fnorm(tmp_I_value_Sm_ - I_value_SM__(:,1+tmp_nM))/fnorm(tmp_I_value_Sm_)));
end;%for ntest=0:n_test-1;
end;% if (flag_check);
%%%%%%%%;
clear CTF_UX_S_k_q_wnS__ CTF_UX_S_l2_S_;
%%%%%%%%;
end;%if (tmp_n_M>0);
end;%for ncluster=0:n_cluster-1;

if ( (nargout>5) & (isempty(I_value_SM__)) ); I_value_SM__ = ones(n_S,n_M); end;

if (verbose>0); disp(sprintf(' %% [finished ampmh_X_cluster_wrap_SM__10]')); end;
