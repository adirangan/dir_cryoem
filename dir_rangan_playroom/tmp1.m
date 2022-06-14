weight_2d_k_all_scaled_ = weight_2d_k_all_*4*pi^2;
CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_),2);
S_k_q_wkS___ = reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]);
CTF_S_k_q_wkS___ = bsxfun(@times,S_k_q_wkS___,reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1]));
CTF_S_k_q_kwS___ = permute(CTF_S_k_q_wkS___,[2,1,3]);
M_k_q_wkM___ = reshape(M_k_q_wkM__,[n_w_max,n_k_p_r,n_M]);
M_k_q_kwM___ = permute(M_k_q_wkM___,[2,1,3]);
%%%%%%%%;
set temporary array dimensions. ;
tmp_n_S = n_S;
tmp_n_M = n_M;
tmp_n_S = min(tmp_n_S,n_S);
tmp_n_M = min(tmp_n_M,n_M);
%%%%%%%%;
tot_t_fill_ = zeros(n_UX_rank,1);
tot_t_mult_ = zeros(n_UX_rank,1);
tot_t_ifft_ = zeros(n_UX_rank,1);
tot_o_fill_ = zeros(n_UX_rank,1);
tot_o_mult_ = zeros(n_UX_rank,1);
tot_o_ifft_ = zeros(n_UX_rank,1);
X_wSMe___ = zeros(n_w_max,n_S,n_M);
X_wSM0___ = zeros(n_w_max,n_S,n_M);
corr_X_X_n_ = zeros(n_UX_rank,1);
prct_X_X_n_ = zeros(n_UX_rank,1);
for nUX_rank=[n_UX_rank-1:-1:0];%for nUX_rank=n_UX_rank-1:-1:0
pm_n_UX_rank = 1+nUX_rank;
if (verbose); disp(sprintf(' %% pm_n_UX_rank %d/%d',pm_n_UX_rank,n_UX_rank)); end;
tmp_UX_weight_kn__ = diag(X_weight_r_)*UX_S_avg_kn__(:,1:pm_n_UX_rank);
tmp_SX_k_ = SX_S_avg_k_(1:pm_n_UX_rank);
%%%%%%%%;
%%%%;
n_M_batch = ceil(n_M/tmp_n_M); n_S_batch = ceil(n_S/tmp_n_S);
for nS_batch=0:n_S_batch-1;
tmp_index_nS_ = nS_batch*tmp_n_S + [0:tmp_n_S-1]; tmp_index_nS_ = intersect(tmp_index_nS_,[0:n_S-1]);
tmp_CTF_S_k_q_kwS___ = CTF_S_k_q_kwS___(:,:,1+tmp_index_nS_);
tmp_t = tic();
UX_weight_CTF_S_k_q_nwS___ = reshape(transpose(tmp_UX_weight_kn__)*reshape(tmp_CTF_S_k_q_kwS___,[n_k_p_r,n_w_max*tmp_n_S]),[pm_n_UX_rank,n_w_max,tmp_n_S]);
UX_weight_conj_CTF_S_k_q_Snw___ = conj(permute(UX_weight_CTF_S_k_q_nwS___,[3,1,2]));
tmp_o = pm_n_UX_rank*n_k_p_r*n_w_max*(tmp_n_M + tmp_n_S);
tmp_t = toc(tmp_t); if (verbose>2); disp(sprintf(' %% fill: %0.3fs --> %0.2fGHz',tmp_t,tmp_o/tmp_t/1e9)); end;
tot_o_fill_(1+nUX_rank) = tot_o_fill_(1+nUX_rank) + tmp_o;
tot_t_fill_(1+nUX_rank) = tot_t_fill_(1+nUX_rank) + tmp_t;
for nM_batch=0:n_M_batch-1;
tmp_index_nM_ = nM_batch*tmp_n_M + [0:tmp_n_M-1]; tmp_index_nM_ = intersect(tmp_index_nM_,[0:n_M-1]);
tmp_M_k_q_kwM___ = M_k_q_kwM___(:,:,1+tmp_index_nM_);
%%%%;
tmp_t = tic();
UX_weight_M_k_q_nwM___ = reshape(transpose(tmp_UX_weight_kn__)*reshape(tmp_M_k_q_kwM___,[n_k_p_r,n_w_max*tmp_n_M]),[pm_n_UX_rank,n_w_max,tmp_n_M]);
UX_weight_M_k_q_nMw___ = permute(UX_weight_M_k_q_nwM___,[1,3,2]);
tmp_o = pm_n_UX_rank*n_k_p_r*n_w_max*(tmp_n_M + tmp_n_S);
tmp_t = toc(tmp_t); if (verbose>2); disp(sprintf(' %% fill: %0.3fs --> %0.2fGHz',tmp_t,tmp_o/tmp_t/1e9)); end;
tot_o_fill_(1+nUX_rank) = tot_o_fill_(1+nUX_rank) + tmp_o;
tot_t_fill_(1+nUX_rank) = tot_t_fill_(1+nUX_rank) + tmp_t;
%%%%;
tmp_t = tic();
conj_S_CTF_weight_weight_M_k_q_SMw___ = zeros(tmp_n_S,tmp_n_M,n_w_max);
for nw=0:n_w_max-1;
conj_S_CTF_weight_weight_M_k_q_SMw___(:,:,1+nw) = UX_weight_conj_CTF_S_k_q_Snw___(:,:,1+nw) * UX_weight_M_k_q_nMw___(:,:,1+nw);
end;%for nw=0:n_w_max-1;
tmp_o = n_w_max * tmp_n_S * tmp_n_M * pm_n_UX_rank;
tmp_t = toc(tmp_t); if (verbose>2); disp(sprintf(' %% mult: %0.3fs --> %0.2fGHz',tmp_t,tmp_o/tmp_t/1e9)); end;
tot_o_mult_(1+nUX_rank) = tot_o_mult_(1+nUX_rank) + tmp_o;
tot_t_mult_(1+nUX_rank) = tot_t_mult_(1+nUX_rank) + tmp_t;
%%%%;
tmp_t = tic();
ifft_conj_S_CTF_weight_weight_M_k_q_wSM___ = ifft(permute(conj_S_CTF_weight_weight_M_k_q_SMw___,[3,1,2]),[],1);
tmp_o = tmp_n_S * tmp_n_M * n_w_max * log(n_w_max);
tmp_t = toc(tmp_t); if (verbose>2); disp(sprintf(' %% ifft: %0.3fs --> %0.2fGHz',tmp_t,tmp_o/tmp_t/1e9)); end;
tot_o_ifft_(1+nUX_rank) = tot_o_ifft_(1+nUX_rank) + tmp_o;
tot_t_ifft_(1+nUX_rank) = tot_t_ifft_(1+nUX_rank) + tmp_t;
X_wSM___(:,1+tmp_index_nS_,1+tmp_index_nM_) = real(ifft_conj_S_CTF_weight_weight_M_k_q_wSM___);
%%%%;
end;end;%for nM_batch=0:n_M_batch-1; for nS_batch=0:n_S_batch-1;
if (verbose>1); disp(sprintf(' %% tot_t_fill %0.3fs --> %0.2fGH',tot_t_fill_(1+nUX_rank),tot_o_fill_(1+nUX_rank)/tot_t_fill_(1+nUX_rank)/1e9)); end;
if (verbose>1); disp(sprintf(' %% tot_t_mult %0.3fs --> %0.2fGH',tot_t_mult_(1+nUX_rank),tot_o_mult_(1+nUX_rank)/tot_t_mult_(1+nUX_rank)/1e9)); end;
if (verbose>1); disp(sprintf(' %% tot_t_ifft %0.3fs --> %0.2fGH',tot_t_ifft_(1+nUX_rank),tot_o_ifft_(1+nUX_rank)/tot_t_ifft_(1+nUX_rank)/1e9)); end;
if (nUX_rank==n_UX_rank-1);
X_wSMe___ = X_wSM___;
end;%if (nUX_rank==n_UX_rank-1);
X_wSM0___ = X_wSM___; 
corr_X_X = corr(X_wSMe___(:),X_wSM0___(:));
corr_X_X_n_(1+nUX_rank) = corr_X_X;
[~,tmp_nw_SM_] = max(X_wSM0___,[],1); tmp_nw_SM_ = tmp_nw_SM_-1;
tmp_X_SM___ = zeros(1,n_S,n_M);
for nM=0:n_M-1; for nS=0:n_S-1;
tmp_nw = tmp_nw_SM_(1,1+nS,1+nM);
tmp_X = X_wSMe___(1+tmp_nw,1+nS,1+nM);
tmp_X_SM___(1,1+nS,1+nM) = tmp_X;
end;end;%for nM=0:n_M-1; for nS=0:n_S-1;
tmp_X_SM___ = repmat(tmp_X_SM___,[n_w_max,1,1]);
tmp_X_SM___ = X_wSMe___ > tmp_X_SM___;
tmp_p__ = sum(tmp_X_SM___,1)/n_w_max;
tmp_p = mean(tmp_p__,'all');
prct_X_X_n_(1+nUX_rank) = tmp_p;
%%%%%%%%;
tmp_t = tic();
gamma_z_ = 2*pi*[0:n_w_max-1]/n_w_max;
nS=3; nM=1;
X1_w_ = X_wSM___(:,1+nS,1+nM);
X0_w_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
gamma_z = 2*pi*nw/n_w_max;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS); M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
T_k_p_wk_ = reshape(reshape(S_k_p_wk_,[n_w_max,n_k_p_r])*diag(CTF_k_p_r_xavg_k_),[n_w_max*n_k_p_r,1]);
N_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_,-gamma_z);
X0_w_(1+nw) = real(sum(conj(T_k_p_wk_).*N_k_p_wk_.*weight_2d_k_all_scaled_));
end;%for nw=0:n_w_max-1;
tmp_o = n_w_max * n_k_p_r * n_w_max;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% test: %0.3fs --> %0.2fGHz',tmp_t,tmp_o/tmp_t/1e9)); end;
if (verbose); disp(sprintf(' %% X0_w_ vs X1_w_: %0.16f',fnorm(X0_w_-X1_w_)/fnorm(X0_w_))); end;
end;%for nUX_rank=0:n_UX_rank-1;

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
% design pricipal-modes. ;
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


%{
for i in $(squeue -u rangan -h -t PD -o %i)
do
scontrol update jobid=$i partition=ccm
done

figure(1);figmed;
x_ = 0:0.01:1;
subplot(1,2,1);
plot(x_,8*x_./(2+6*x_),'-','LineWidth',2);
xlabel('infected fraction where you live');
ylabel('probability you had covid');
title('I had a fever after vaccine!');
grid on; xlim([0,1]); ylim([0,1]); axis square;
subplot(1,2,2);
plot(x_,64*x_./(4+60*x_),'-','LineWidth',2);
xlabel('infected fraction where you live');
ylabel('probability you both had covid');
title('Me and Partner both had fever after vaccine!');
grid on; xlim([0,1]); ylim([0,1]); axis square;


% testing ampmut_wrap_wrap_3.m ;
 
fname_pre = sprintf('%s_mat/ut%.4dn%.2dr%d',dir_pm,floor(1000*delta_r_max_use),dat_n_UX_rank,dat_rseed);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);

if ( exist(fname_mat,'file'));

disp(sprintf(' %% %s found, aligning',fname_mat));
tmp_ = load(fname_mat);
if (~isfield(tmp_.parameter,'fname_pre')); tmp_.parameter.fname_pre = fname_pre; end; %<-- parameter_bookmark. ;
%%%%%%%%;
tmp_fname_pre = sprintf('%s_align_a_CTF_avg_UX_Y_',tmp_.parameter.fname_pre);
tmp_.parameter.fname_align_a_CTF_avg_UX_Y_pre = tmp_fname_pre;
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
disp(sprintf(' %% %s not found, creating',tmp_fname_mat));
ampmut_align_to_a_CTF_avg_UX_Y_0( ...
 tmp_.parameter ...
,l_max_max ...
,dat_n_UX_rank ...
,reshape(a_CTF_avg_UX_Y_quad__(:,1:dat_n_UX_rank),[n_lm_max*dat_n_UX_rank,1]) ...
,n_M ...
,tmp_euler_polar_a_true_ ...
,tmp_euler_azimu_b_true_ ...
,tmp_euler_gamma_z_true_ ...
,tmp_image_delta_x_true_ ...
,tmp_image_delta_y_true_ ...
,[] ...
,tmp_.a_CTF_avg_UX_Y__ ...
,tmp_.euler_polar_a__ ...
,tmp_.euler_azimu_b__ ...
,tmp_.euler_gamma_z__ ...
,tmp_.image_delta_x_acc__ + tmp_.image_delta_x_upd__ ...
,tmp_.image_delta_y_acc__ + tmp_.image_delta_y_upd__ ...
);
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
%%%%%%%%;
tmp_fname_pre = sprintf('%s_align_a_k_Y_',tmp_.parameter.fname_pre);
tmp_.parameter.fname_align_a_k_Y_pre = tmp_fname_pre;
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
disp(sprintf(' %% %s not found, creating',tmp_fname_mat));
ampmut_align_to_a_k_Y_0( ...
 tmp_.parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_M ...
,M_k_p__ ...
,CTF_index_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_quad_ ...
,[] ...
,tmp_.euler_polar_a__ ...
,tmp_.euler_azimu_b__ ...
,tmp_.euler_gamma_z__ ...
,tmp_.image_delta_x_acc__ + tmp_.image_delta_x_upd__ ...
,tmp_.image_delta_y_acc__ + tmp_.image_delta_y_upd__ ...
);
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
%%%%%%%%;

end;%if ( exist(fname_mat,'file'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

load('ampmut_1_debug.mat');

image_delta_x = +0.028;
image_delta_y = -0.012;
nM=0;
%%%%%%%%;
ori_M_k_p_ = ...
transf_p_to_p( ...
 n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,n_w_sum ...
,M_k_p__(:,1+nM) ...
,+image_delta_x ...
,+image_delta_y ...
);
ori_M_k_q_ = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,ori_M_k_p_ ...
);
%%%%%%%%;
ori_svd_VUXM_lwn___ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,1,ori_M_k_q_,pm_n_UX_rank,UX__,X_weight_r_);
%%%%%%%%;
[ori_UX_M_k_q_wn__,ori_UX_M_k_p_wn__] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,1,ori_svd_VUXM_lwn___,zeros(1,1),zeros(1,1));
%%%%%%%%;
pos_M_k_p_ = ...
transf_p_to_p( ...
 n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,n_w_sum ...
,M_k_p__(:,1+nM) ...
,0 ...
,0 ...
);
pos_M_k_q_ = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,pos_M_k_p_ ...
);
%%%%%%%%;
pos_svd_VUXM_lwn___ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,1,pos_M_k_q_,pm_n_UX_rank,UX__,X_weight_r_);
%%%%%%%%;
[pos_UX_M_k_q_wn__,pos_UX_M_k_p_wn__] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,1,pos_svd_VUXM_lwn___,+image_delta_x,+image_delta_y);
%%%%%%%%;
neg_M_k_p_ = ...
transf_p_to_p( ...
 n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,n_w_sum ...
,M_k_p__(:,1+nM) ...
,0 ...
,0 ...
);
neg_M_k_q_ = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,neg_M_k_p_ ...
);
%%%%%%%%;
neg_svd_VUXM_lwn___ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,1,neg_M_k_q_,pm_n_UX_rank,UX__,X_weight_r_);
%%%%%%%%;
[neg_UX_M_k_q_wn__,neg_UX_M_k_p_wn__] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,1,neg_svd_VUXM_lwn___,-image_delta_x,-image_delta_y);
%%%%%%%%;

figure(1);clf;figbig;np=0;
subplot(2,2,1+np);np=np+1;
plot(real(ori_UX_M_k_q_wn__(:)),real(pos_UX_M_k_q_wn__(:)),'kx'); xlabel('ori');ylabel('pos');title('real pos'); axis equal;
subplot(2,2,1+np);np=np+1;
plot(imag(ori_UX_M_k_q_wn__(:)),imag(pos_UX_M_k_q_wn__(:)),'kx'); xlabel('ori');ylabel('pos');title('imag pos'); axis equal;
subplot(2,2,1+np);np=np+1;
plot(real(ori_UX_M_k_q_wn__(:)),real(neg_UX_M_k_q_wn__(:)),'kx'); xlabel('ori');ylabel('neg');title('real neg'); axis equal;
subplot(2,2,1+np);np=np+1;
plot(imag(ori_UX_M_k_q_wn__(:)),imag(neg_UX_M_k_q_wn__(:)),'kx'); xlabel('ori');ylabel('neg');title('imag neg'); axis equal;
 %}



%{
sprintf('%s_jpg',dir_trunk);
delta_sigma_use = delta_sigma;
svd_eps = 1e-3;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_UX_rank = 8;
dat_n_iteration = 16;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
%delta_r_max_factor_ = [0,0.125,0.25,0.5,1,sqrt(2),2]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0]; n_delta_r_max_factor = numel(delta_r_max_factor_);
delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0.25]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('ut%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
for ndat_rseed=n_dat_rseed-1;%for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for dat_nUX_rank=0:dat_n_UX_rank-1;%for dat_nUX_rank=0:dat_n_UX_rank-1;
fname_0 = sprintf('tpmutameux_%s_n%.3d',dat_infix,dat_n_M);
fname_2 = sprintf('nUX%.3drng%.3d',dat_nUX_rank,dat_rseed);
MS_fname_mat = sprintf('%s_mat/%s_MS_%s.mat',dir_trunk,fname_0,fname_2);
SM_fname_mat = sprintf('%s_mat/%s_SM_%s.mat',dir_trunk,fname_0,fname_2);
if ( exist(MS_fname_mat,'file'));
disp(sprintf(' %% %s found',MS_fname_mat));
fname_1 = 'MS';
tpmutamf_1( ...
 sprintf('%s_mat',dir_trunk) ...
,sprintf('%s_jpg',dir_trunk) ...
,fname_0 ...
,fname_1 ...
,fname_2 ...
,n_image_sub ...
,M_k_p__ ...
,dat_nUX_rank ...
,dat_n_UX_rank ...
,dat_n_iteration ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,delta_r_max_use ...
,delta_r_max_upb ...
,svd_eps ...
,n_delta_v_requested ...
,n_CTF_rank ...
,CTF_idx_ ...
,CTF_k_p_r__ ...
,CTF_k_p__ ...
,l_max_ ...
,UX__ ...
,X_weight_r_ ...
,a_CTF_UX_Y_quad__ ...
,[] ...
,a_k_Y_quad_ ...
,euler_angle_marina_ ...
,delta_read_x_ ...
,delta_read_y_ ...
,delta_sigma ...
);
end;%if ( exist(MS_fname_mat,'file'));
if ( exist(SM_fname_mat,'file'));
disp(sprintf(' %% %s found',SM_fname_mat));
fname_1 = 'SM';
tpmutamf_1( ...
 sprintf('%s_mat',dir_trunk) ...
,sprintf('%s_jpg',dir_trunk) ...
,fname_0 ...
,fname_1 ...
,fname_2 ...
,n_image_sub ...
,M_k_p__ ...
,dat_nUX_rank ...
,dat_n_UX_rank ...
,dat_n_iteration ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,delta_r_max_use ...
,delta_r_max_upb ...
,svd_eps ...
,n_delta_v_requested ...
,n_CTF_rank ...
,CTF_idx_ ...
,CTF_k_p_r__ ...
,CTF_k_p__ ...
,l_max_ ...
,UX__ ...
,X_weight_r_ ...
,a_CTF_UX_Y_quad__ ...
,[] ...
,a_k_Y_quad_ ...
,euler_angle_marina_ ...
,delta_read_x_ ...
,delta_read_y_ ...
,delta_sigma ...
);
end;%if ( exist(SM_fname_mat,'file'));
end;%for dat_nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
 %}

%{
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_a_k_p_quad_.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_a_k_p_quad_.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_a_k_Y_quad_.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_a_k_Y_quad_.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_c_k_Y_.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_c_k_Y_.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_CTF_k_p__.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_CTF_k_p__.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_M_x_c___.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_M_x_c___.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_S_k_p__.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_S_k_p__.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_M_k_p___.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_M_k_p___.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_S_CTF_k_p___.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_S_CTF_k_p___.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_X_2d_xcor_d0047.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_X_2d_xcor_d0047.mat')
  %}
%{
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_a_k_p_quad_.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_a_k_p_quad_.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_a_k_Y_quad_A.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_a_k_Y_quad_A.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_a_k_Y_quad_.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_a_k_Y_quad_.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_CTF_k_p__.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_CTF_k_p__.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_CTF_k_p_r_xxxx__.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_CTF_k_p_r_xxxx__.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_c_x_u_.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_c_x_u_.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_delta_read_.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_delta_read_.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_euler_angle_marina_.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_euler_angle_marina_.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_M_x_c___center_.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_M_x_c___center_.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_M_x_c___sample.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_M_x_c___sample.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_S_k_p__A.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_S_k_p__A.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_S_k_p__.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_S_k_p__.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_M_k_p___.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_M_k_p___.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_M_k_q___.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_M_k_q___.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_M_k_q___spectrum.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_M_k_q___spectrum.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_S_CTF_k_p___A.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_S_CTF_k_p___A.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_S_CTF_k_p___.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_S_CTF_k_p___.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_S_CTF_k_q___.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_S_CTF_k_q___.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_S_CTF_k_q___spectrum.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_S_CTF_k_q___spectrum.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_X_2d_xcor_d0047_A.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_X_2d_xcor_d0047_A.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_X_2d_xcor_d0047_B.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_X_2d_xcor_d0047_B.jpg')
  %}

%{
  % replicates tmp_X2__ ; 
tmp_X3__ = zeros(tmp_n_delta_v,n_w_max);
%%%%%%%%;
tmp_U_d__ = zeros(tmp_n_delta_v,FTK.n_svd_l);
tmp_delta_r_ = zeros(tmp_n_delta_v,1);
for ndelta_v=0:tmp_n_delta_v-1;
tmp_U_d_ = zeros(FTK.n_svd_l,1);
tmp_delta_r_(1+ndelta_v) = fnorm(tmp_delta__(:,1+ndelta_v));
if (tmp_delta_r>FTK.svd_d_max);
disp(sprintf(' %% Warning, tmp_delta_r %0.6f > svd_d_max %0.6f',tmp_delta_r,FTK.svd_d_max));
end;%if;
tmp_svd_d = (tmp_delta_r_(1+ndelta_v) - FTK.svd_d_m)/FTK.svd_d_c;
for nl=0:FTK.n_svd_l-1;
tmp_U_d__(1+ndelta_v,1+nl) = polyval_r8_reverse_0(FTK.n_svd_d,FTK.svd_U_d_(1+0+nl*FTK.n_svd_d+(0:FTK.n_svd_d-1)),1,tmp_svd_d);
end;%for nl=0:FTK.n_svd_l-1;
end;%for ndelta_v=0:tmp_n_delta_v-1;
%%%%%%%%;
tmp_delta_w_ = atan2(tmp_delta__(1+1,:),tmp_delta__(1+0,:));
tmp_expiw__ = exp(-i*(pi/2 - reshape(tmp_delta_w_,[tmp_n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l]));
%%%%%%%%;
tmp_V_r__ = zeros(FTK.n_svd_l,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
if (2*pi*k_p_r_(1+nk_p_r)>FTK.svd_r_max);
disp(sprintf(' %% Warning, 2*pi*k_p_r_(1+nk_p_r) %0.6f > FTK.svd_r_max %0.6f',2*pi*k_p_r_(1+nk_p_r),FTK.svd_r_max));
end;%if;
tmp_svd_r = (2*pi*k_p_r_(1+nk_p_r) - FTK.svd_r_m)/FTK.svd_r_c;
for nl=0:FTK.n_svd_l-1;
tmp_V_r__(1+nl,1+nk_p_r) = polyval_r8_reverse_0(FTK.n_svd_r,FTK.svd_V_r_(1+0+nl*FTK.n_svd_r+(0:FTK.n_svd_r-1)),1,tmp_svd_r);
end;%for nl=0:FTK.n_svd_l-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_l_max = max(abs(FTK.svd_l_));
tmp_M_k_q__ = zeros(n_w_sum,1+2*tmp_l_max);
ic=0;
for nk_p_r=0:n_k_p_r-1;
for nw=0:n_w_(1+nk_p_r)-1;
tmp_n_w_t = floor(1.0d0*n_w_(1+nk_p_r)/2.0d0);
for tmp_I_l=-tmp_l_max:+tmp_l_max;
nwc = nw;
if (nwc>=tmp_n_w_t); nwc = nwc - n_w_(1+nk_p_r); end;%if;
flag_ict_overflow = 0;
nwd = nwc + tmp_I_l;
if (abs(nwd)<tmp_n_w_t); nwt = periodize(nwd,0,n_w_(1+nk_p_r));  else; nwt = 0; flag_ict_overflow = 1; end;%if;
ict = ic-nw+nwt;
if (flag_ict_overflow==0); tmp_M_k_q__(1+ic,1+tmp_l_max+tmp_I_l) = tmp_M_k_q_(1+ict); end;%if;
end;%for tmp_I_l=-tmp_l_max:+tmp_l_max;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_M_k_q_rwl___ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,weight_2d_k_p_r_/(2*pi),n_w_,1+2*tmp_l_max,tmp_M_k_q__);
tmp_S_k_q_rw__ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,weight_2d_k_p_r_/(2*pi),n_w_,1,tmp_S_k_q_);
%%%%%%%%;
tmp_VSM__ = zeros(FTK.n_svd_l,n_w_max);
for nl=0:FTK.n_svd_l-1;
tmp_VSM__(1+nl,:) = ifft(tmp_V_r__(1+nl,:)*(conj(tmp_S_k_q_rw__(:,:)).*tmp_M_k_q_rwl___(:,:,1+tmp_l_max+FTK.svd_l_(1+nl))))*n_w_max;
end;%for nl=0:FTK.n_svd_l-1;
tmp_USEVSM__ = (tmp_U_d__.*tmp_expiw__)*diag(FTK.svd_s_)*tmp_VSM__;
%%%%%%%%;
tmp_X3__ = tmp_USEVSM__;
disp(sprintf(' %% tmp_X0__ vs tmp_X3__: %0.16f',fnorm(tmp_X0__-tmp_X3__)/fnorm(tmp_X0__)));
 %}

%{
  % replicates tmp_X2__ ; 
tmp_X3__ = zeros(tmp_n_delta_v,n_w_max);
%%%%%%%%;
tmp_U_d__ = zeros(tmp_n_delta_v,FTK.n_svd_l);
tmp_delta_r_ = zeros(tmp_n_delta_v,1);
for ndelta_v=0:tmp_n_delta_v-1;
tmp_U_d_ = zeros(FTK.n_svd_l,1);
tmp_delta_r_(1+ndelta_v) = fnorm(tmp_delta__(:,1+ndelta_v));
if (tmp_delta_r>FTK.svd_d_max);
disp(sprintf(' %% Warning, tmp_delta_r %0.6f > svd_d_max %0.6f',tmp_delta_r,FTK.svd_d_max));
end;%if;
tmp_svd_d = (tmp_delta_r_(1+ndelta_v) - FTK.svd_d_m)/FTK.svd_d_c;
for nl=0:FTK.n_svd_l-1;
tmp_U_d__(1+ndelta_v,1+nl) = polyval_r8_reverse_0(FTK.n_svd_d,FTK.svd_U_d_(1+0+nl*FTK.n_svd_d+(0:FTK.n_svd_d-1)),1,tmp_svd_d);
end;%for nl=0:FTK.n_svd_l-1;
end;%for ndelta_v=0:tmp_n_delta_v-1;
%%%%%%%%;
tmp_delta_w_ = atan2(tmp_delta__(1+1,:),tmp_delta__(1+0,:));
tmp_expiw__ = exp(-i*(pi/2 - reshape(tmp_delta_w_,[tmp_n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l]));
%%%%%%%%;
tmp_V_r__ = zeros(FTK.n_svd_l,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
if (2*pi*k_p_r_(1+nk_p_r)>FTK.svd_r_max);
disp(sprintf(' %% Warning, 2*pi*k_p_r_(1+nk_p_r) %0.6f > FTK.svd_r_max %0.6f',2*pi*k_p_r_(1+nk_p_r),FTK.svd_r_max));
end;%if;
tmp_svd_r = (2*pi*k_p_r_(1+nk_p_r) - FTK.svd_r_m)/FTK.svd_r_c;
for nl=0:FTK.n_svd_l-1;
tmp_V_r__(1+nl,1+nk_p_r) = polyval_r8_reverse_0(FTK.n_svd_r,FTK.svd_V_r_(1+0+nl*FTK.n_svd_r+(0:FTK.n_svd_r-1)),1,tmp_svd_r);
end;%for nl=0:FTK.n_svd_l-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_l_max = max(abs(FTK.svd_l_));
tmp_M_k_q__ = zeros(n_w_sum,1+2*tmp_l_max);
ic=0;
for nk_p_r=0:n_k_p_r-1;
for nw=0:n_w_(1+nk_p_r)-1;
tmp_n_w_t = floor(1.0d0*n_w_(1+nk_p_r)/2.0d0);
for tmp_I_l=-tmp_l_max:+tmp_l_max;
nwc = nw;
if (nwc>=tmp_n_w_t); nwc = nwc - n_w_(1+nk_p_r); end;%if;
flag_ict_overflow = 0;
nwd = nwc + tmp_I_l;
if (abs(nwd)<tmp_n_w_t); nwt = periodize(nwd,0,n_w_(1+nk_p_r));  else; nwt = 0; flag_ict_overflow = 1; end;%if;
ict = ic-nw+nwt;
if (flag_ict_overflow==0); tmp_M_k_q__(1+ic,1+tmp_l_max+tmp_I_l) = tmp_M_k_q_(1+ict); end;%if;
end;%for tmp_I_l=-tmp_l_max:+tmp_l_max;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_M_k_q_rwl___ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,weight_2d_k_p_r_/(2*pi),n_w_,1+2*tmp_l_max,tmp_M_k_q__);
tmp_S_k_q_rw__ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,weight_2d_k_p_r_/(2*pi),n_w_,1,tmp_S_k_q_);
%%%%%%%%;
for ndelta_v=0:tmp_n_delta_v-1;
tmp_Z_q_ = zeros(n_w_sum,1);
tmp_C_q_ = zeros(n_w_(1+n_k_p_r-1),1);
for nw=0:n_w_max-1;
for nk_p_r=0:n_k_p_r-1;
for nl=0:FTK.n_svd_l-1;
tmp_C_q_(1+nw) = tmp_C_q_(1+nw) + ...
  tmp_U_d__(1+ndelta_v,1+nl) * ...
  FTK.svd_s_(1+nl) * ...
  tmp_V_r__(1+nl,1+nk_p_r) * ...
  tmp_expiw__(1+ndelta_v,1+nl) * ...
  tmp_M_k_q_rwl___(1+nk_p_r,1+nw,1+tmp_l_max+FTK.svd_l_(1+nl)) * ...
  conj(tmp_S_k_q_rw__(1+nk_p_r,1+nw));  
end;%for nl=0:FTK.n_svd_l-1;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nw=0:n_w_max-1;
%%%%%%%%;
tmp_X3__(1+ndelta_v,:) = ifft(tmp_C_q_)*n_w_max;
end;%for ndelta_v=0:tmp_n_delta_v-1;
disp(sprintf(' %% tmp_X0__ vs tmp_X3__: %0.16f',fnorm(tmp_X0__-tmp_X3__)/fnorm(tmp_X0__)));
 %}

%{
  % replicates tmp_X2__ ; 
tmp_X3__ = zeros(tmp_n_delta_v,n_w_max);
%%%%%%%%;
tmp_U_d__ = zeros(tmp_n_delta_v,FTK.n_svd_l);
tmp_delta_r_ = zeros(tmp_n_delta_v,1);
for ndelta_v=0:tmp_n_delta_v-1;
tmp_U_d_ = zeros(FTK.n_svd_l,1);
tmp_delta_r_(1+ndelta_v) = fnorm(tmp_delta__(:,1+ndelta_v));
if (tmp_delta_r>FTK.svd_d_max);
disp(sprintf(' %% Warning, tmp_delta_r %0.6f > svd_d_max %0.6f',tmp_delta_r,FTK.svd_d_max));
end;%if;
tmp_svd_d = (tmp_delta_r_(1+ndelta_v) - FTK.svd_d_m)/FTK.svd_d_c;
for nl=0:FTK.n_svd_l-1;
tmp_U_d__(1+ndelta_v,1+nl) = polyval_r8_reverse_0(FTK.n_svd_d,FTK.svd_U_d_(1+0+nl*FTK.n_svd_d+(0:FTK.n_svd_d-1)),1,tmp_svd_d);
end;%for nl=0:FTK.n_svd_l-1;
end;%for ndelta_v=0:tmp_n_delta_v-1;
%%%%%%%%;
tmp_delta_w_ = atan2(tmp_delta__(1+1,:),tmp_delta__(1+0,:));
tmp_expiw__ = exp(-i*(pi/2 - reshape(tmp_delta_w_,[tmp_n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l]));
%%%%%%%%;
tmp_V_r__ = zeros(FTK.n_svd_l,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
if (2*pi*k_p_r_(1+nk_p_r)>FTK.svd_r_max);
disp(sprintf(' %% Warning, 2*pi*k_p_r_(1+nk_p_r) %0.6f > FTK.svd_r_max %0.6f',2*pi*k_p_r_(1+nk_p_r),FTK.svd_r_max));
end;%if;
tmp_svd_r = (2*pi*k_p_r_(1+nk_p_r) - FTK.svd_r_m)/FTK.svd_r_c;
for nl=0:FTK.n_svd_l-1;
tmp_V_r__(1+nl,1+nk_p_r) = polyval_r8_reverse_0(FTK.n_svd_r,FTK.svd_V_r_(1+0+nl*FTK.n_svd_r+(0:FTK.n_svd_r-1)),1,tmp_svd_r);
end;%for nl=0:FTK.n_svd_l-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_l_max = max(abs(FTK.svd_l_));
tmp_M_k_q__ = zeros(n_w_sum,1+2*tmp_l_max);
ic=0;
for nk_p_r=0:n_k_p_r-1;
for nw=0:n_w_(1+nk_p_r)-1;
tmp_n_w_t = floor(1.0d0*n_w_(1+nk_p_r)/2.0d0);
for tmp_I_l=-tmp_l_max:+tmp_l_max;
nwc = nw;
if (nwc>=tmp_n_w_t); nwc = nwc - n_w_(1+nk_p_r); end;%if;
flag_ict_overflow = 0;
nwd = nwc + tmp_I_l;
if (abs(nwd)<tmp_n_w_t); nwt = periodize(nwd,0,n_w_(1+nk_p_r));  else; nwt = 0; flag_ict_overflow = 1; end;%if;
ict = ic-nw+nwt;
if (flag_ict_overflow==0); tmp_M_k_q__(1+ic,1+tmp_l_max+tmp_I_l) = tmp_M_k_q_(1+ict); end;%if;
end;%for tmp_I_l=-tmp_l_max:+tmp_l_max;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
for ndelta_v=0:tmp_n_delta_v-1;
tmp_Z_q_ = zeros(n_w_sum,1);
tmp_C_q_ = zeros(n_w_(1+n_k_p_r-1),1);
ic=0;
for nk_p_r=0:n_k_p_r-1;
tmp_dw = 2*pi/(1.0d0*max(1,n_w_(1+nk_p_r))); tmp_dA = weight_2d_k_p_r_(1+nk_p_r)/(2*pi); tmp_dAn = tmp_dA*tmp_dw;
for nw=0:n_w_(1+nk_p_r)-1;
n_w_t = floor(n_w_(1+nk_p_r)/2);
tmp_Z_q_(1+ic) = 0;
for nl=0:FTK.n_svd_l-1;
tmp_Z_q_(1+ic) = ...
  tmp_Z_q_(1+ic) + ...
  tmp_U_d__(1+ndelta_v,1+nl) * ...
  FTK.svd_s_(1+nl) * ...
  tmp_V_r__(1+nl,1+nk_p_r) * ...
  tmp_expiw__(1+ndelta_v,1+nl) * ...
  tmp_M_k_q__(1+ic,1+tmp_l_max+FTK.svd_l_(1+nl)) ;
end;%for nl=0:FTK.n_svd_l-1;
if (nw>n_w_t); nw_fix = nw - n_w_(1+nk_p_r) + n_w_max; end; 
if (nw<n_w_t); nw_fix = nw; end;
if (nw~=n_w_t); tmp_C_q_(1+nw_fix) = tmp_C_q_(1+nw_fix) + conj(tmp_S_k_q_(1+ic)) * tmp_Z_q_(1+ic) * tmp_dAn ; end;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_X3__(1+ndelta_v,:) = ifft(tmp_C_q_)*n_w_max;
end;%for ndelta_v=0:tmp_n_delta_v-1;
disp(sprintf(' %% tmp_X0__ vs tmp_X3__: %0.16f',fnorm(tmp_X0__-tmp_X3__)/fnorm(tmp_X0__)));
 %}

%{
  % replicates tmp_X2__ ; 
tmp_X3__ = zeros(tmp_n_delta_v,n_w_max);
%%%%%%%%;
tmp_U_d__ = zeros(tmp_n_delta_v,FTK.n_svd_l);
tmp_delta_r_ = zeros(tmp_n_delta_v,1);
for ndelta_v=0:tmp_n_delta_v-1;
tmp_U_d_ = zeros(FTK.n_svd_l,1);
tmp_delta_r_(1+ndelta_v) = fnorm(tmp_delta__(:,1+ndelta_v));
if (tmp_delta_r>FTK.svd_d_max);
disp(sprintf(' %% Warning, tmp_delta_r %0.6f > svd_d_max %0.6f',tmp_delta_r,FTK.svd_d_max));
end;%if;
tmp_svd_d = (tmp_delta_r_(1+ndelta_v) - FTK.svd_d_m)/FTK.svd_d_c;
for nl=0:FTK.n_svd_l-1;
tmp_U_d__(1+ndelta_v,1+nl) = polyval_r8_reverse_0(FTK.n_svd_d,FTK.svd_U_d_(1+0+nl*FTK.n_svd_d+(0:FTK.n_svd_d-1)),1,tmp_svd_d);
end;%for nl=0:FTK.n_svd_l-1;
end;%for ndelta_v=0:tmp_n_delta_v-1;
%%%%%%%%;
tmp_delta_w_ = atan2(tmp_delta__(1+1,:),tmp_delta__(1+0,:));
tmp_expiw__ = exp(-i*(pi/2 - reshape(tmp_delta_w_,[tmp_n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l]));
%%%%%%%%;
tmp_V_r__ = zeros(FTK.n_svd_l,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
if (2*pi*k_p_r_(1+nk_p_r)>FTK.svd_r_max);
disp(sprintf(' %% Warning, 2*pi*k_p_r_(1+nk_p_r) %0.6f > FTK.svd_r_max %0.6f',2*pi*k_p_r_(1+nk_p_r),FTK.svd_r_max));
end;%if;
tmp_svd_r = (2*pi*k_p_r_(1+nk_p_r) - FTK.svd_r_m)/FTK.svd_r_c;
for nl=0:FTK.n_svd_l-1;
tmp_V_r__(1+nl,1+nk_p_r) = polyval_r8_reverse_0(FTK.n_svd_r,FTK.svd_V_r_(1+0+nl*FTK.n_svd_r+(0:FTK.n_svd_r-1)),1,tmp_svd_r);
end;%for nl=0:FTK.n_svd_l-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_l_max = max(abs(FTK.svd_l_));
tmp_M_k_q__ = zeros(n_w_sum,1+2*tmp_l_max);
ic=0;
for nk_p_r=0:n_k_p_r-1;
for nw=0:n_w_(1+nk_p_r)-1;
tmp_n_w_t = floor(1.0d0*n_w_(1+nk_p_r)/2.0d0);
for tmp_I_l=-tmp_l_max:+tmp_l_max;
nwc = nw;
if (nwc>=tmp_n_w_t); nwc = nwc - n_w_(1+nk_p_r); end;%if;
flag_ict_overflow = 0;
nwd = nwc + tmp_I_l;
if (abs(nwd)<tmp_n_w_t); nwt = periodize(nwd,0,n_w_(1+nk_p_r));  else; nwt = 0; flag_ict_overflow = 1; end;%if;
ict = ic-nw+nwt;
if (flag_ict_overflow==0); tmp_M_k_q__(1+ic,1+tmp_l_max+tmp_I_l) = tmp_M_k_q_(1+ict); end;%if;
end;%for tmp_I_l=-tmp_l_max:+tmp_l_max;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
for ndelta_v=0:tmp_n_delta_v-1;
tmp_Z_q_ = zeros(n_w_sum,1);
tmp_C_w_ = zeros(FTK.n_svd_l,1);
ic=0;
for nk_p_r=0:n_k_p_r-1;
for nw=0:n_w_(1+nk_p_r)-1;
tmp_Z_q_(1+ic) = 0;
for nl=0:FTK.n_svd_l-1;
tmp_Z_q_(1+ic) = tmp_Z_q_(1+ic) + ...
                 tmp_U_d__(1+ndelta_v,1+nl) * ...
                 FTK.svd_s_(1+nl) * ...
                 tmp_V_r__(1+nl,1+nk_p_r) * ...
                 tmp_expiw__(1+ndelta_v,1+nl) * ...
                 tmp_M_k_q__(1+ic,1+tmp_l_max+FTK.svd_l_(1+nl)) ;
end;%for nl=0:FTK.n_svd_l-1;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_T_k_q_ = tmp_Z_q_;
%%%%%%%%;
tmp_C_q_ = zeros(n_w_(1+n_k_p_r-1),1);
ic = 0;
for nk_p_r=0:n_k_p_r-1;
tmp_dw = 2*pi/(1.0d0*max(1,n_w_(1+nk_p_r)));
tmp_dA = weight_2d_k_p_r_(1+nk_p_r)/(2*pi);
% We assume that the fourier basis is orthonormal (not merely orthogonal);
tmp_dAn = tmp_dA*tmp_dw;
for nw=0:n_w_(1+nk_p_r)-1;
if (nw>n_w_(1+nk_p_r)/2);
nw_fix = nw - n_w_(1+nk_p_r) + n_w_max;
tmp_C_q = conj(tmp_S_k_q_(1+ic))*tmp_T_k_q_(1+ic);
nw_C = nw_fix;
tmp_C_q_(1+nw_C) = tmp_C_q_(1+nw_C) + tmp_C_q*tmp_dAn;
end;%if (nw>n_w_(1+nk_p_r)/2);
if (nw<n_w_(1+nk_p_r)/2);
nw_fix = nw;
tmp_C_q = conj(tmp_S_k_q_(1+ic))*tmp_T_k_q_(1+ic);
nw_C = nw_fix;
tmp_C_q_(1+nw_C) = tmp_C_q_(1+nw_C) + tmp_C_q*tmp_dAn;
end%if;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_X3__(1+ndelta_v,:) = ifft(tmp_C_q_)*n_w_max;
end;%for ndelta_v=0:tmp_n_delta_v-1;
disp(sprintf(' %% tmp_X0__ vs tmp_X3__: %0.16f',fnorm(tmp_X0__-tmp_X3__)/fnorm(tmp_X0__)));
 %}

%{
  % replicates tmp_X2__ ; 
tmp_X3__ = zeros(tmp_n_delta_v,n_w_max);
for ndelta_v=0:tmp_n_delta_v-1;
%%%%%%%%;
tmp_Z_q_ = zeros(n_w_sum,1);
tmp_U_d_ = zeros(FTK.n_svd_l,1);
tmp_delta_r = fnorm(tmp_delta__(:,1+ndelta_v));
if (tmp_delta_r>FTK.svd_d_max);
disp(sprintf(' %% Warning, tmp_delta_r %0.6f > svd_d_max %0.6f',tmp_delta_r,FTK.svd_d_max));
end;%if;
tmp_delta_w = atan2(tmp_delta__(1+1,1+ndelta_v),tmp_delta__(1+0,1+ndelta_v));
tmp_svd_d = (tmp_delta_r - FTK.svd_d_m)/FTK.svd_d_c;
for nl=0:FTK.n_svd_l-1;
tmp_U_d_(1+nl) = polyval_r8_reverse_0(FTK.n_svd_d,FTK.svd_U_d_(1+0+nl*FTK.n_svd_d+(0:FTK.n_svd_d-1)),1,tmp_svd_d);
end;%for nl=0:FTK.n_svd_l-1;
tmp_V_r_ = zeros(FTK.n_svd_l*n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
if (2*pi*k_p_r_(1+nk_p_r)>FTK.svd_r_max);
disp(sprintf(' %% Warning, 2*pi*k_p_r_(1+nk_p_r) %0.6f > FTK.svd_r_max %0.6f',2*pi*k_p_r_(1+nk_p_r),FTK.svd_r_max));
end;%if;
tmp_svd_r = (2*pi*k_p_r_(1+nk_p_r) - FTK.svd_r_m)/FTK.svd_r_c;
for nl=0:FTK.n_svd_l-1;
tmp_V_r_(1+nl+nk_p_r*FTK.n_svd_l) = polyval_r8_reverse_0(FTK.n_svd_r,FTK.svd_V_r_(1+0+nl*FTK.n_svd_r+(0:FTK.n_svd_r-1)),1,tmp_svd_r);
end;%for nl=0:FTK.n_svd_l-1;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_C_w_ = zeros(FTK.n_svd_l,1);
ic=0;
for nk_p_r=0:n_k_p_r-1;
for nl=0:FTK.n_svd_l-1;
tmp_theta = FTK.svd_l_(1+nl)*(pi/2 - tmp_delta_w);
tmp_C_w = +cos(tmp_theta) - i*sin(tmp_theta);
tmp_D_V_r = tmp_V_r_(1+nl+nk_p_r*FTK.n_svd_l);
tmp_D_U_d = tmp_U_d_(1+nl);
tmp_D_s = FTK.svd_s_(1+nl);
tmp_C_w_(1+nl) = (tmp_D_U_d * tmp_D_s * tmp_D_V_r) * tmp_C_w;
end;%for nl=0:FTK.n_svd_l-1;
for nw=0:n_w_(1+nk_p_r)-1;
tmp_Z_q_(1+ic) = 0;
tmp_n_w_t = floor(1.0d0*n_w_(1+nk_p_r)/2.0d0);
for nl=0:FTK.n_svd_l-1;
tmp_I_l = FTK.svd_l_(1+nl);
tmp_C_q = tmp_C_w_(1+nl);
nwc = nw;
if (nwc>=tmp_n_w_t);
nwc = nwc - n_w_(1+nk_p_r);
end;%if;
flag_ict_overflow = 0;
nwd = nwc + tmp_I_l;
if (abs(nwd)<tmp_n_w_t);
nwt = periodize(nwd,0,n_w_(1+nk_p_r));
 else;
nwt = 0;
flag_ict_overflow = 1;
end;%if;
ict = ic-nw+nwt;
if (flag_ict_overflow==0);
tmp_Z_q_(1+ic) = tmp_Z_q_(1+ic) + tmp_C_q*tmp_M_k_q_(1+ict);
end;%if;
end;%for nl=0:FTK.n_svd_l-1;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_T_k_q_ = tmp_Z_q_;
%%%%%%%%;
tmp_C_q_ = zeros(n_w_(1+n_k_p_r-1),1);
ic = 0;
for nk_p_r=0:n_k_p_r-1;
tmp_dw = 2*pi/(1.0d0*max(1,n_w_(1+nk_p_r)));
tmp_dA = weight_2d_k_p_r_(1+nk_p_r)/(2*pi);
% We assume that the fourier basis is orthonormal (not merely orthogonal);
tmp_dAn = tmp_dA*tmp_dw;
for nw=0:n_w_(1+nk_p_r)-1;
if (nw>n_w_(1+nk_p_r)/2);
nw_fix = nw - n_w_(1+nk_p_r) + n_w_max;
tmp_C_q = conj(tmp_S_k_q_(1+ic))*tmp_T_k_q_(1+ic);
nw_C = nw_fix;
tmp_C_q_(1+nw_C) = tmp_C_q_(1+nw_C) + tmp_C_q*tmp_dAn;
end;%if (nw>n_w_(1+nk_p_r)/2);
if (nw<n_w_(1+nk_p_r)/2);
nw_fix = nw;
tmp_C_q = conj(tmp_S_k_q_(1+ic))*tmp_T_k_q_(1+ic);
nw_C = nw_fix;
tmp_C_q_(1+nw_C) = tmp_C_q_(1+nw_C) + tmp_C_q*tmp_dAn;
end%if;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_X3__(1+ndelta_v,:) = ifft(tmp_C_q_)*n_w_max;
end;%for ndelta_v=0:tmp_n_delta_v-1;
 %}

%{
tmp_U__ = 2*pi * reshape(FTK.svd_polyval_U_d_,[FTK.n_delta_v,FTK.n_svd_l]) .* exp(-i*atan2(FTK.delta_y_(:),FTK.delta_x_(:))*transpose(FTK.svd_l_(:))) * diag((-i).^FTK.svd_l_);
tmp_S__ = diag(FTK.svd_s_);
tmp_V__ = reshape(FTK.svd_polyval_V_r_,[n_k_p_r,FTK.n_svd_l]);
tmp_s__ = zeros(n_w_max,n_k_p_r); %<-- new zero-mode is at n_w_max/2 instead of 0;
tmp_m__ = zeros(n_w_max,n_k_p_r);
tmp_n_w_max_2 = round(n_w_max/2);
for nk_p_r=0:n_k_p_r-1;
tmp_n_w_2 = round(n_w_(1+nk_p_r)/2)-1; 
tmp_ij_sub_ = 0:tmp_n_w_2-1;
tmp_ij_set_ = tmp_ij_sub_ + tmp_n_w_max_2;
tmp_ij_ = n_w_csum_(1+nk_p_r) + tmp_ij_sub_;
tmp_s__(1+tmp_ij_set_,1+nk_p_r) = tmp_S_k_q_(1+tmp_ij_);
tmp_m__(1+tmp_ij_set_,1+nk_p_r) = tmp_M_k_q_(1+tmp_ij_);
tmp_ij_sub_ = n_w_(1+nk_p_r)-tmp_n_w_2:n_w_(1+nk_p_r)-1;
tmp_ij_set_ = tmp_ij_sub_ + tmp_n_w_max_2 - n_w_(1+nk_p_r);
tmp_ij_ = n_w_csum_(1+nk_p_r) + tmp_ij_sub_;
tmp_s__(1+tmp_ij_set_,1+nk_p_r) = tmp_S_k_q_(1+tmp_ij_);
tmp_m__(1+tmp_ij_set_,1+nk_p_r) = tmp_M_k_q_(1+tmp_ij_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_svwm__ = zeros(FTK.n_svd_l,n_w_max);
for nsvd_l=0:FTK.n_svd_l-1;
tmp_l = FTK.svd_l_(1+nsvd_l);
if (tmp_l==0);
tmp_ij_ = 0:n_w_max-1; 
tmp_ = sum(conj(tmp_s__).*repmat(transpose(weight_2d_k_p_r_(:).*tmp_V__(:,1+nsvd_l)),[n_w_max,1]).*tmp_m__(1+tmp_ij_,:),2);
tmp_ij_ = [[tmp_n_w_max_2:n_w_max-1],[0:tmp_n_w_max_2-1]];
tmp_svwm__(1+nsvd_l,:) = ifft(tmp_(1+tmp_ij_))*n_w_max;
end;%if (tmp_l==0);
if (tmp_l> 0); 
tmp_ij_ = 0:n_w_max-1-tmp_l; 
tmp_ = sum(conj(tmp_s__).*repmat(transpose(weight_2d_k_p_r_(:).*tmp_V__(:,1+nsvd_l)),[n_w_max,1]).*[zeros(tmp_l,n_k_p_r);tmp_m__(1+tmp_ij_,:)],2);
tmp_ij_ = [[tmp_n_w_max_2:n_w_max-1],[0:tmp_n_w_max_2-1]];
tmp_svwm__(1+nsvd_l,:) = ifft(tmp_(1+tmp_ij_))*n_w_max;
end;%if (tmp_l> 0); 
if (tmp_l< 0); 
tmp_ij_ = -tmp_l:n_w_max-1;
tmp_ = sum(conj(tmp_s__).*repmat(transpose(weight_2d_k_p_r_(:).*tmp_V__(:,1+nsvd_l)),[n_w_max,1]).*[tmp_m__(1+tmp_ij_,:);zeros(-tmp_l,n_k_p_r)],2);
tmp_ij_ = [[tmp_n_w_max_2:n_w_max-1],[0:tmp_n_w_max_2-1]];
tmp_svwm__(1+nsvd_l,:) = ifft(tmp_(1+tmp_ij_))*n_w_max;
end;%if (tmp_l< 0); 
end;%for nsvd_l=0:FTK.n_svd_l-1;
%%%%%%%%;
tmp_X3__ = tmp_U__*tmp_S__*tmp_svwm__;
%%%%%%%%;
flag_plot=0;
if flag_plot;
clf;
colormap(colormap_beach());
subplot(2,2,1+0); imagesc(real(tmp_X0__)); axisnotick; title('X0');
subplot(2,2,1+1); imagesc(real(tmp_X1__)); axisnotick; title('X1');
subplot(2,2,1+2); imagesc(real(tmp_X2__)); axisnotick; title('X2');
subplot(2,2,1+3); imagesc(real(tmp_X3__)); axisnotick; title('X3');
figbig;
end;%if flag_plot;
 %}


%{
function M_q_ = transf_svd_q_to_q_FTK_5(FTK,n_r,grid_p_,n_w_,n_A,S_q_,delta_x,delta_y);

warning_flag = 1;

svd_r_max = FTK.svd_r_max;
n_svd_r = FTK.n_svd_r;
svd_r_ = FTK.svd_r_;
svd_d_max = FTK.svd_d_max;
n_svd_d = FTK.n_svd_d;
svd_d_ = FTK.svd_d_;
n_svd_l = FTK.n_svd_l;
svd_l_ = FTK.svd_l_;
svd_U_d_ = FTK.svd_U_d_;
svd_s_ = FTK.svd_s_;
svd_V_r_ = FTK.svd_V_r_;

Z_q_ = zeros(n_A,1);
U_d_ = zeros(n_svd_l,1);
svd_d_m = svd_d_max / 2.0;
svd_d_c = svd_d_m;
delta = sqrt(delta_x^2 + delta_y^2);
omega = atan2(delta_y,delta_x);
if (delta>svd_d_max & warning_flag);
disp(sprintf(' %% Warning, delta %0.6f > svd_d_max %0.6f',delta,svd_d_max));
end;%if;
svd_d = (delta - svd_d_m)/svd_d_c;
for nl=0:n_svd_l-1;
U_d_(1+nl) = polyval_r8_reverse_0(n_svd_d,svd_U_d_(1+0+nl*n_svd_d+(0:n_svd_d-1)),1,svd_d);
end;%for nl=0:n_svd_l-1;
V_r_ = zeros(n_svd_l*n_r,1);
svd_r_m = svd_r_max / 2.0;
svd_r_c = svd_r_m;
for nr=0:n_r-1;
if (grid_p_(1+nr)>svd_r_max & warning_flag);
disp(sprintf(' %% Warning, grid_p_(1+nr) %0.6f > svd_r_max %0.6f',grid_p_(1+nr),svd_r_max));
end;%if;
svd_r = (grid_p_(1+nr) - svd_r_m)/svd_r_c;
for nl=0:n_svd_l-1;
V_r_(1+nl+nr*n_svd_l) = polyval_r8_reverse_0(n_svd_r,svd_V_r_(1+0+nl*n_svd_r+(0:n_svd_r-1)),1,svd_r);
end;%for nl=0:n_svd_l-1;
end;%for nr=0:n_r-1;
C_w_ = zeros(n_svd_l,1);
ic=0;
for nr=0:n_r-1;
for nl=0:n_svd_l-1;
theta = svd_l_(1+nl)*(pi/2 - omega);
C_w = +cos(theta) - i*sin(theta);
D_V_r = V_r_(1+nl+nr*n_svd_l);
D_U_d = U_d_(1+nl);
D_s = svd_s_(1+nl);
C_w_(1+nl) = (D_U_d * D_s * D_V_r) * C_w;
end;%for nl=0:n_svd_l-1;
for nw=0:n_w_(1+nr)-1;
Z_q_(1+ic) = 0;
n_w_t = floor(1.0d0*n_w_(1+nr)/2.0d0);
for nl=0:n_svd_l-1;
I_l = svd_l_(1+nl);
C_q = C_w_(1+nl);
nwc = nw;
if (nwc>=n_w_t);
nwc = nwc - n_w_(1+nr);
end;%if;
flag_ict_overflow = 0;
nwd = nwc + I_l;
if (abs(nwd)<n_w_t);
nwt = periodize(nwd,0,n_w_(1+nr));
 else;
nwt = 0;
flag_ict_overflow = 1;
end;%if;
ict = ic-nw+nwt;
if (flag_ict_overflow==0);
Z_q_(1+ic) = Z_q_(1+ic) + C_q*S_q_(1+ict);
end;%if;
end;%for nl=0:n_svd_l-1;
ic = ic + 1;
end;%for nw=0:n_w_(1+nr)-1;
end;%for nr=0:n_r-1;
M_q_ = Z_q_;
 %}

%{
function C_q_ = innerproduct_q_k_stretch_quad_0(n_r,grid_p_,weight_p_,n_w_,n_A,T_q_,M_q_) ;
% Assumes that M_q_ is the same size and dimensions as T_q_. ;
% Assumes quasi-uniform polar-grid. ;
% Assumes that C_q_ is large enough to hold all n_w_(1+n_r-1) modes ;
% (assuming of course that n_w_(1+n_r-1) is the largest value within n_w_). ;
% Stores C_q_ in fft ordering (mode 0 , mode 1 , ... , mode -1) ;
% Upgraded to ignore frequencies of magnitude n_w_(1+nr)/2 or larger. ;
% Upgraded to include radial weight. ;
verbose=0;
C_q_ = zeros(n_w_(1+n_r-1),1);
if (verbose>0); disp(sprintf(' %% [entering innerproduct_q_k_stretch_quad_0] n_r %d',n_r)); end%if;
n_w_max = n_w_(1+n_r-1);
if (verbose>0); disp(sprintf(' %% n_w_max %d',n_w_max)); end%if;
for nw=0:n_w_max-1; C_q_(1+nw) = 0.0; end;%for nw=0:n_w_max-1; C_q = 0.0;
ic = 0;
for nr=0:n_r-1;
dw = 2*pi/(1.0d0*max(1,n_w_(1+nr)));
dA = weight_p_(1+nr);
% We assume that the fourier basis is orthonormal (not merely orthogonal);
dAn = dA*dw;
for nw=0:n_w_(1+nr)-1;
if (nw>n_w_(1+nr)/2);
nw_fix = nw - n_w_(1+nr) + n_w_max;
if (verbose>2 & nr<5); disp(sprintf(' %% nw %d nw_fix %d (full loop)',nw,nw_fix)); end%if;
C_q = conj(T_q_(1+ic))*M_q_(1+ic);
nw_C = nw_fix;
C_q_(1+nw_C) = C_q_(1+nw_C) + C_q*dAn;
 elseif (nw==n_w_(1+nr)/2);
nw_fix = nw;
if (verbose>2 & nr<5); disp(sprintf(' %% nw %d nw_fix %d (first orig)',nw,nw_fix)); end%if;
C_q = 0.0d0*0.5d0*conj(T_q_(1+ic))*M_q_(1+ic);
nw_C = nw_fix;
C_q_(1+nw_C) = C_q_(1+nw_C) + C_q*dAn;
nw_fix = nw - n_w_(1+nr) + n_w_max;
if (verbose>2 & nr<5); disp(sprintf(' %% nw %d nw_fix %d (then loop)',nw,nw_fix)); end%if;
C_q = 0.0d0*0.5d0*conj(T_q_(1+ic))*M_q_(1+ic);
nw_C = nw_fix;
C_q_(1+nw_C) = C_q_(1+nw_C) + C_q*dAn;
 else;
nw_fix = nw;
if (verbose>2 & nr<5); disp(sprintf(' %% nw %d nw_fix %d (full orig)',nw,nw_fix)); end;%if;
C_q = conj(T_q_(1+ic))*M_q_(1+ic);
nw_C = nw_fix;
C_q_(1+nw_C) = C_q_(1+nw_C) + C_q*dAn;
end%if;
ic = ic + 1;
end;%for nw=0:n_w_(1+nr)-1;
end;%for nr=0:n_r-1;
if (verbose>0); disp(sprintf(' %% [finished innerproduct_q_k_stretch_quad_0]')); end;%if;
 %}

%{


 %tmp_ms_ = load('dir_principled_marching_mat/tpmhtameux_2d_xcor_d0047_le3v128_n1024_MS_nUX002rng000.mat');
 %tmp_sm_ = load('dir_principled_marching_mat/tpmhtameux_2d_xcor_d0047_le3v128_n1024_SM_nUX002rng000.mat');
 tmp_ms_ = load('dir_principled_marching_mat/tpmhtamsux_2d_xcor_d0047_le7v128_n1024s0125_MS_nUX003rng000.mat');
 tmp_sm_ = load('dir_principled_marching_mat/tpmhtamsux_2d_xcor_d0047_le7v128_n1024s0125_SM_nUX003rng000.mat');
  % plots histograms of image parameters. ;
   figure(1);clf;
  subplot(2,3,1); hist(tmp_ms_.euler_polar_a_MS__(:,end),linspace(0,2*pi,64)); title('ms polar');
  subplot(2,3,2); hist(tmp_ms_.euler_azimu_b_MS__(:,end),linspace(0,2*pi,64)); title('ms azimu');
  subplot(2,3,3); hist(tmp_ms_.euler_gamma_z_MS__(:,end),linspace(0,2*pi,64)); title('ms gamma');
  subplot(2,3,4); hist(tmp_sm_.euler_polar_a_SM__(:,end),linspace(0,2*pi,64)); title('sm polar');
  subplot(2,3,5); hist(tmp_sm_.euler_azimu_b_SM__(:,end),linspace(0,2*pi,64)); title('sm azimu');
  subplot(2,3,6); hist(tmp_sm_.euler_gamma_z_SM__(:,end),linspace(0,2*pi,64)); title('sm gamma');
  figbig;
  % plots histograms of image parameters. ;
   figure(2);clf;
  subplot(2,2,1); hist(tmp_ms_.image_delta_x_MS__(:,end),linspace(-0.11,+0.11,16)); title('ms polar');
  subplot(2,2,2); hist(tmp_ms_.image_delta_y_MS__(:,end),linspace(-0.11,+0.11,16)); title('ms azimu');
  subplot(2,2,3); hist(tmp_sm_.image_delta_x_SM__(:,end),linspace(-0.11,+0.11,16)); title('sm polar');
  subplot(2,2,4); hist(tmp_sm_.image_delta_y_SM__(:,end),linspace(-0.11,+0.11,16)); title('sm azimu');
  figbig;
  % plots scatterplots of image parameters. ;
   figure(3);clf;
   ms_ij = 32;
   sm_ij = 32;
   subplot(2,3,1); plot(tmp_ms_.euler_polar_a_MS__(:,ms_ij),tmp_sm_.euler_polar_a_SM__(:,sm_ij),'x'); axis equal;
   subplot(2,3,2); plot(tmp_ms_.euler_azimu_b_MS__(:,ms_ij),tmp_sm_.euler_azimu_b_SM__(:,sm_ij),'x'); axis equal;
   subplot(2,3,3); plot(tmp_ms_.euler_gamma_z_MS__(:,ms_ij),tmp_sm_.euler_gamma_z_SM__(:,sm_ij),'x'); axis equal;
   subplot(2,3,4); plot(tmp_ms_.image_delta_x_MS__(:,ms_ij),tmp_sm_.image_delta_x_SM__(:,sm_ij),'x'); axis equal;
   subplot(2,3,5); plot(tmp_ms_.image_delta_y_MS__(:,ms_ij),tmp_sm_.image_delta_y_SM__(:,sm_ij),'x'); axis equal;
  figbig;


  %}
