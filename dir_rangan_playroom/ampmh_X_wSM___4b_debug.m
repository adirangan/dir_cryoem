function ...
[ ...
 X_wSM___ ...
,delta_j_wSM___ ...
,I_value_wSM___ ...
] = ...
ampmh_X_wSM___4b_debug(...
 FTK...
,n_w_max...
,pm_n_UX_rank...
,n_S...
,n_S_per_Sbatch...
,CTF_UX_S_k_q_wS__ ...
,UX_S_l2_...
,n_M...
,n_M_per_Mbatch...
,svd_VUXM_lwnM____...
,UX_M_l2_dM__...
);
%%%%%%%%;
% uses raw templates. ;
%%%%%%%%;

if nargin<1;

%%%%%%%%;
% First define integral of <f,f>. ;
%%%%%%%%;
h2d_ = @(kd) 4*pi^2*(besselj(0,kd) + besselj(2,kd)); % calculates <f_j,f_k>, normalized so that <f,f> = (4*pi^2);
dh2d_ = @(kd) 4*pi^3*(besselj(-1,kd) - besselj(+3,kd));
h3d_ = @(kd) 4*pi*( sin(kd) - (kd).*cos(kd) ) ./ kd.^3 ; % calculates <f_j,f_k>, normalized so that <f,f> = 4*pi/3;
dh3d_ = @(kd) 12*pi*( (kd.^2/3 - 1) .* sin(kd) + (kd).*cos(kd) ) ./ kd.^4 ;

verbose=1;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi); TorL = 'L';
if (verbose); disp(sprintf(' %% [testing ampmh_X_wSM___4b_debug.m]')); end;
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
template_k_eq_d = k_eq_d*2;
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
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% get_template_1: %0.2fs',tmp_t)); end;
if (verbose); disp(sprintf(' %% n_viewing_all %d n_viewing_polar_a %d n_w_max %d',n_viewing_all,n_viewing_polar_a,max(n_w_))); end;
n_S = n_viewing_all; n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
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
nS=max(0,min(n_S-1,floor(n_S*rand())));
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
%%%%%%%%;
n_w_uni_ = n_w_max*ones(n_k_p_r,1);
n_w_uni_csum_ = [0;cumsum(n_w_uni_)];
n_w_uni_sum = sum(n_w_uni_);
weight_uni_2d_k_all_ = zeros(n_w_uni_sum,1);
weight_uni_2d_k_all_ = reshape(ones(n_w_max,1)*transpose(weight_2d_k_p_r_)/n_w_max,[n_w_uni_sum,1]);
%%%%%%%%;
pole_uni_k_c_0_ = zeros(n_w_uni_sum,1);
pole_uni_k_c_1_ = zeros(n_w_uni_sum,1);
pole_uni_k_c_2_ = zeros(n_w_uni_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_w_uni = n_w_uni_(1+nk_p_r);
for nw=0:n_w_uni-1;
gamma_z = 2*pi*nw/n_w_uni;
cc = cos(gamma_z); sc = sin(gamma_z);
pole_uni_k_c_0_(1+na) = k_p_r*cc;
pole_uni_k_c_1_(1+na) = k_p_r*sc;
pole_uni_k_c_2_(1+na) = 0;
na=na+1;
end;%for nw=0:n_w_uni-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_t = tic();
S_uni_k_q__ = zeros(n_w_uni_sum,n_S);
for nS=0:n_S-1;
viewing_azimu_b = viewing_azimu_b_all_(1+nS);
cb = cos(+viewing_azimu_b); sb = sin(+viewing_azimu_b);
Rz = [ +cb , -sb , 0 ; +sb , +cb , 0 ; 0 , 0 , 1 ]; %<-- rotation about the positive z-axis. ;
viewing_polar_a = viewing_polar_a_all_(1+nS);
ca = cos(+viewing_polar_a); sa = sin(+viewing_polar_a);
Ry = [ +ca , 0 , +sa ; 0 , 1 , 0 ; -sa , 0 , +ca ]; %<-- rotation about the positive y-axis. ;
delta_temp_ = transpose(Ry)*transpose(Rz)*delta_orig_;
S_uni_k_p_form_ = exp(+2*pi*i*(pole_uni_k_c_0_*delta_temp_(1+0) + pole_uni_k_c_1_*delta_temp_(1+1) + pole_uni_k_c_2_*delta_temp_(1+2)));
S_uni_k_q__(:,1+nS) = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_uni_ ...
,n_w_uni_sum ...
,S_uni_k_p_form_ ...
);
end;%for nS=0:n_S-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% S_uni_k_q__: %0.3fs',tmp_t)); end;
%%%%%%%%;
delta_r_max = 0.05; n_delta_v_requested = 16;
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
% Now calculate the inner-products. ;
%%%%%%%%;
n_UX_rank = n_k_p_r;
tmp_t = tic();
svd_VUXS_lwnS____ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_uni_,n_S,S_uni_k_q__,n_UX_rank,eye(n_UX_rank,n_UX_rank),sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXS_lwnS____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the non-translated images. ;
%%%%%%%%;
tmp_t = tic();
UX_S_l2_dS__ = ampmh_UX_M_l2_dM__1(FTK,n_w_uni_,n_S,n_UX_rank,svd_VUXS_lwnS____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_S_l2_dS__: %0.3fs',tmp_t)); end;
disp(sprintf(' %% average l2-norm of templates (should be 1): %0.16f',mean(UX_S_l2_dS__(:))/(pi*k_p_r_max^2)));
flag_plot=0;
if flag_plot;
plot(UX_S_l2_dS__(:)/(pi*k_p_r_max^2),'.');
end;%if flag_plot;
%%%%%%%%;
tmp_t = tic();
[UX_S_k_q_wnS___,UX_S_k_p_wnS___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_uni_,n_UX_rank,n_S,svd_VUXS_lwnS____,zeros(n_S,1),zeros(n_S,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_S_k_q_wnS___: %0.6fs',tmp_t)); end;
UX_S_k_p_wnS__ = reshape(UX_S_k_p_wnS___(:,1:n_k_p_r,:),[n_w_max*n_k_p_r,n_S]);
UX_S_k_q_wnS__ = reshape(UX_S_k_q_wnS___(:,1:n_k_p_r,:),[n_w_max*n_k_p_r,n_S]);
%%%%%%%%;
% Calculate ampmh_X_wSS___4b_debug. ;
%%%%%%%%;
tmp_t = tic();
SS_k_q_ = svd(UX_S_k_q_wnS__);
n_S_rank = min(find(SS_k_q_/SS_k_q_(1)<1e-3)); if isempty(n_S_rank); n_S_rank = min(size(UX_S_k_q_wnS__)); end;
[US_k_q__,SS_k_q__,VS_k_q__] = svds(UX_S_k_q_wnS__,n_S_rank);
if (verbose>1); disp(sprintf(' %% n_S %d --> n_S_rank %d',n_S,n_S_rank)); end;
%%%%%%%%;
n_S_per_Sbatch = 24;
tmp_t = tic();
[ ...
 X_wSS___ ...
,delta_j_wSS___ ...
,I_value_wSS___ ...
] = ...
ampmh_X_wSM___4b_debug( ...
 FTK ...
,n_w_max ...
,n_UX_rank ...
,n_S ...
,n_S_per_Sbatch ...
,UX_S_k_q_wnS__ ...
,UX_S_l2_dS__ ...
,n_S ...
,n_S_per_Sbatch ...
,svd_VUXS_lwnS____ ...
,UX_S_l2_dS__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wSS___: %0.3fs',tmp_t)); end;
X_wSS___ = real(X_wSS___);
%%%%%%%%;
% Finally check ampmh_X_wSM___4b_debug to see how the rotations and translations are implemented. ;
%%%%%%%%;
error_form_vs_quad_max = 0;
error_form_vs_calc_max = 0;
for nS0=0:n_S-1;
viewing_azimu_b = viewing_azimu_b_all_(1+nS0);
cb = cos(+viewing_azimu_b); sb = sin(+viewing_azimu_b);
Rz = [ +cb , -sb , 0 ; +sb , +cb , 0 ; 0 , 0 , 1 ]; %<-- rotation about the positive z-axis. ;
viewing_polar_a = viewing_polar_a_all_(1+nS0);
ca = cos(+viewing_polar_a); sa = sin(+viewing_polar_a);
Ry = [ +ca , 0 , +sa ; 0 , 1 , 0 ; -sa , 0 , +ca ]; %<-- rotation about the positive y-axis. ;
delta_temp_0_ = transpose(Ry)*transpose(Rz)*delta_orig_;
S_0_uni_k_p_form_ = exp(+2*pi*i*(pole_uni_k_c_0_*delta_temp_0_(1+0) + pole_uni_k_c_1_*delta_temp_0_(1+1) + pole_uni_k_c_2_*delta_temp_0_(1+2)));
for nS1=0:n_S-1;
viewing_azimu_b = viewing_azimu_b_all_(1+nS1);
cb = cos(+viewing_azimu_b); sb = sin(+viewing_azimu_b);
Rz = [ +cb , -sb , 0 ; +sb , +cb , 0 ; 0 , 0 , 1 ]; %<-- rotation about the positive z-axis. ;
viewing_polar_a = viewing_polar_a_all_(1+nS1);
ca = cos(+viewing_polar_a); sa = sin(+viewing_polar_a);
Ry = [ +ca , 0 , +sa ; 0 , 1 , 0 ; -sa , 0 , +ca ]; %<-- rotation about the positive y-axis. ;
delta_temp_1_ = transpose(Ry)*transpose(Rz)*delta_orig_;
S_1_uni_k_p_form_ = exp(+2*pi*i*(pole_uni_k_c_0_*delta_temp_1_(1+0) + pole_uni_k_c_1_*delta_temp_1_(1+1) + pole_uni_k_c_2_*delta_temp_1_(1+2)));
%%%%;
X_form_ = zeros(n_w_max,1);
X_quad_ = zeros(n_w_max,1);
X_calc_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
gamma_z = 2*pi*nw/n_w_max;
cc = cos(+gamma_z); sc = sin(+gamma_z);
Rz = [+cc , -sc ; +sc , +cc];
delta_j = delta_j_wSS___(1+nw,1+nS0,1+nS1);
delta_x = FTK.delta_x_(1+delta_j);
delta_y = FTK.delta_y_(1+delta_j);
delta_temp_0b_ = Rz*delta_temp_0_(1+[0:1]); %<-- rotate delta_temp_0_ by +gamma_z = rotate k by -gamma_z = rotate S_0 by +gamma_z. ;
delta_temp_1b_ = delta_temp_1_(1+[0:1]) - [delta_x;delta_y]; %<-- translate delta_temp_1_ by -delta_ = multiply S_1 by exp(-2*pi*i*dot(k_,delta_)). ;
X_form_(1+nw) = h2d_(2*pi*k_p_r_max*fnorm(delta_temp_0b_ - delta_temp_1b_))/(2*pi)^2 * (pi*k_p_r_max^2); %<-- note sign of translation. ;
T_0_uni_k_p_form_ = rotate_p_to_p_fftw(n_k_p_r,n_w_uni_,n_w_uni_sum,S_0_uni_k_p_form_,+gamma_z); %<-- rotate S_0 by +gamma_z. ;
T_1_uni_k_p_form_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_uni_,n_w_uni_sum,S_1_uni_k_p_form_,+delta_x,+delta_y); %<-- multiply S_1 by exp(-2*pi*i*dot(k_,delta_)). ;
X_quad_(1+nw) = real(sum(conj(T_0_uni_k_p_form_).*T_1_uni_k_p_form_.*weight_uni_2d_k_all_));
X_calc_(1+nw) = X_wSS___(1+nw,1+nS0,1+nS1)*(pi*k_p_r_max^2);
end;%for nw=0:n_w_max-1;
flag_plot=0;
if flag_plot;
gamma_z_ = 2*pi*transpose([0:n_w_max-1])/n_w_max;
plot(gamma_z_,X_form_,'k.-',gamma_z_,X_quad_,'ro-',gamma_z_,X_calc_,'bx'); xlim([0,2*pi]); legend({'form','quad','calc'});
end;%if flag_plot;
error_form_vs_quad = fnorm(X_form_ - X_quad_)/fnorm(X_form_);
if (verbose>0); disp(sprintf(' %% nS0 %.2d nS1 %.2d / %.2d error_form_vs_quad: %0.16f',nS0,nS1,n_S,error_form_vs_quad)); end;
error_form_vs_calc = fnorm(X_form_ - X_calc_)/fnorm(X_form_);
if (verbose>0); disp(sprintf(' %% nS0 %.2d nS1 %.2d / %.2d error_form_vs_calc: %0.16f',nS0,nS1,n_S,error_form_vs_calc)); end;
error_form_vs_quad_max = max(error_form_vs_quad_max,error_form_vs_quad);
error_form_vs_calc_max = max(error_form_vs_calc_max,error_form_vs_calc);
%%%%;
end;%for nS1=0:n_S-1;
end;%for nS0=0:n_S-1;
if (verbose); disp(sprintf(' %% error_form_vs_calc_max: %0.16f',error_form_vs_calc_max)); end;
if (verbose); disp(sprintf(' %% error_form_vs_quad_max: %0.16f',error_form_vs_quad_max)); end;

disp('returning'); return;
end;% if nargin<1;

verbose = 1;
if (verbose>0); disp(sprintf(' %% [entering ampmh_X_wSM___4b_debug]')); end;
%n_w_max = max(n_w_);
X_wSM___ = zeros(n_w_max,n_S,n_M);
delta_j_wSM___ = zeros(n_w_max,n_S,n_M);
I_value_wSM___ = zeros(n_w_max,n_S,n_M);
CTF_UX_S_k_q_nSw___ = permute(reshape(CTF_UX_S_k_q_wS__,[n_w_max,pm_n_UX_rank,n_S]),[2,3,1]);
if (verbose>0); 
tmp_str = 'X_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'delta_j_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'I_value_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'CTF_UX_S_k_q_nSw___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if (verbose>0); 
n_Mbatch = ceil(n_M/n_M_per_Mbatch);
if (verbose>1); disp(sprintf(' %% n_Mbatch %d',n_Mbatch)); end;
n_Sbatch = ceil(n_S/n_S_per_Sbatch);
if (verbose>1); disp(sprintf(' %% n_Sbatch %d',n_Sbatch)); end;
for nMbatch=0:n_Mbatch-1;
M_ij_ = nMbatch*n_M_per_Mbatch + (0:n_M_per_Mbatch-1);
M_ij_ = M_ij_(find(M_ij_<n_M)); n_M_sub = numel(M_ij_);
if (verbose>1); disp(sprintf(' %% nMbatch %d/%d M_ij_ %d-->%d',nMbatch,n_Mbatch,M_ij_(1+0),M_ij_(1+n_M_sub-1))); end;
if (verbose>0 & mod(nMbatch,1)==0); disp(sprintf(' %% nMbatch %d/%d M_ij_ %d-->%d',nMbatch,n_Mbatch,M_ij_(1+0),M_ij_(1+n_M_sub-1))); end;
if (n_M_sub>0);
tmp_t = tic();
svd_VUXM_nMwl____ = zeros(pm_n_UX_rank,n_M_sub,n_w_max,FTK.n_svd_l);
svd_VUXM_nMwl____ = permute(svd_VUXM_lwnM____(:,:,:,1+M_ij_),[3,4,2,1]);
svd_SVUXM_SMwl____ = zeros(n_S,n_M_sub,n_w_max,FTK.n_svd_l);
for nl=0:FTK.n_svd_l-1;
for nw=0:n_w_max-1;
svd_SVUXM_SMwl____(:,:,1+nw,1+nl) = ctranspose(CTF_UX_S_k_q_nSw___(:,:,1+nw))*svd_VUXM_nMwl____(:,:,1+nw,1+nl);
end;%for nw=0:n_w_max-1;
end;%for nl=0:FTK.n_svd_l-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_SVUXM_SMwl____: %0.6f',tmp_t)); end;
tmp_t = tic();
svd_SVUXM_SMwl____ = permute(ifft(permute(svd_SVUXM_SMwl____,[3,4,1,2]),[],1)*n_w_max,[3,4,1,2]);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_SVUXM_SMwl____: %0.6f',tmp_t)); end;
for nSbatch=0:n_Sbatch-1;
S_ij_ = nSbatch*n_S_per_Sbatch + (0:n_S_per_Sbatch-1);
S_ij_ = S_ij_(find(S_ij_<n_S)); n_S_sub = numel(S_ij_);
if (verbose>1); disp(sprintf(' %% nSbatch %d/%d S_ij_ %d-->%d',nSbatch,n_Sbatch,S_ij_(1+0),S_ij_(1+n_S_sub-1))); end;
if (verbose>0 & mod(nSbatch,32)==0); disp(sprintf(' %% nSbatch %d/%d S_ij_ %d-->%d',nSbatch,n_Sbatch,S_ij_(1+0),S_ij_(1+n_S_sub-1))); end;
if (n_S_sub>0);
tmp_t = tic();
%%%%%%%%;
svd_SVUXM_lwSM____ = permute(svd_SVUXM_SMwl____(1+S_ij_,:,:,:),[4,3,1,2]);
%%%%%%%%;
svd_USESVUXM_dwSM____ = reshape(FTK.svd_U_d_expiw_s__*reshape(svd_SVUXM_lwSM____,[FTK.n_svd_l,n_w_max*n_S_sub*n_M_sub]),[FTK.n_delta_v,n_w_max,n_S_sub,n_M_sub]);
%%%%%%%%;
l2_dSM___ = permute(reshape(reshape(sqrt(UX_S_l2_(1+S_ij_)),[n_S_sub,1])*reshape(sqrt(UX_M_l2_dM__(:,1+M_ij_)),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
n2_dSM___ = 1./max(1e-14,l2_dSM___);
f2_dSM___ = permute(reshape(reshape(sqrt(UX_S_l2_(1+S_ij_)),[n_S_sub,1])*reshape(1./max(1e-14,sqrt(UX_M_l2_dM__(:,1+M_ij_))),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
ss_S_ = reshape(UX_S_l2_(1+S_ij_),[n_S_sub,1]);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_USESVUXM_dwSM____: %0.6f',tmp_t)); end;
if (nMbatch==0 && nSbatch==0 && verbose>0); 
tmp_str = 'svd_VUXM_nMwl____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_SVUXM_SMwl____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_SVUXM_lwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_USESVUXM_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'l2_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'n2_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'f2_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if (verbose>1); 
tmp_t = tic();
X_sub_dwSM____ = repmat(reshape(n2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* svd_USESVUXM_dwSM____; %<-- correlation. ;
%%%%%%%%;
%tmp_s_ = sparse(1:FTK.n_delta_v*1*n_S_sub*n_M_sub,1:FTK.n_delta_v*1*n_S_sub*n_M_sub,f2_dSM___(:),FTK.n_delta_v*1*n_S_sub*n_M_sub,FTK.n_delta_v*1*n_S_sub*n_M_sub);
%tmp_0_ = reshape(permute(X_sub_dwSM____,[1,3,4,2]),[FTK.n_delta_v*n_S_sub*n_M_sub,n_w_max]);
%tmp_1_ = tmp_s_*tmp_0_;
%tmp_2_ = permute(reshape(tmp_1_,[FTK.n_delta_v,n_S_sub,n_M_sub,n_w_max]),[1,4,2,3]);
%tmp_3_ = max(0,tmp_2_);
%%%%;
I_value_sub_dwSM____ = repmat(reshape(f2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* X_sub_dwSM____; %<-- I_value. ;
I_value_use_dwSM____ = max(0,real(I_value_sub_dwSM____));
%%%%%%%%;
% L = -I*I*<M,M> + 2*I*<M,S> - <S,S> ;
% L_sub_dwSM____ = -abs(I_value_use_dwSM____).^2.*repmat(reshape(real(UX_M_l2_dM__(:,1+M_ij_)),[FTK.n_delta_v,1,1,n_M_sub]),[1,n_w_max,n_S_sub,1]) + 2*I_value_use_dwSM____.*real(svd_USESVUXM_dwSM____) - repmat(reshape(real(ss_S_),[1,1,n_S_sub,1]),[FTK.n_delta_v,n_w_max,1,n_M_sub]); %<-- log-likelihood. ;
if (nMbatch==0 && nSbatch==0 && verbose>0); 
%tmp_str = 'L_sub_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if (verbose>1); 
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% L_sub_dwSM____: %0.6f',tmp_t)); end;
tmp_t = tic();
%[tmp_L_,tmp_delta_j_] = max(reshape(L_sub_dwSM____,[FTK.n_delta_v,n_w_max*n_S_sub*n_M_sub]),[],1); %<-- maximize log-likelihood. ;
[tmp_X_,tmp_delta_j_] = max(reshape(X_sub_dwSM____,[FTK.n_delta_v,n_w_max*n_S_sub*n_M_sub]),[],1); %<-- maximize correlation. ;
tmp_delta_j_ = tmp_delta_j_ - 1; assert(min(tmp_delta_j_)>=0); assert(max(tmp_delta_j_)<FTK.n_delta_v);
tmp_I_value_use_dwSM__ = reshape(I_value_use_dwSM____,[FTK.n_delta_v,n_w_max*n_S_sub*n_M_sub]);
tmp_I_value_use_wSM_ = zeros(n_w_max*n_S_sub*n_M_sub,1);
tmp_t2=tic();
for nl=0:n_w_max*n_S_sub*n_M_sub-1;
tmp_I_value_use_wSM_(1+nl) = tmp_I_value_use_dwSM__(1+tmp_delta_j_(1+nl),1+nl);
end;%for nl=0:n_w_max*n_S_sub*n_M_sub-1;
tmp_t2 = toc(tmp_t2); if (verbose>1); disp(sprintf(' %% tmp_I_value_use_wSM_ %0.6fs',tmp_t2)); end;
X_wSM___(:,1+S_ij_,1+M_ij_) = reshape(tmp_X_,[n_w_max,n_S_sub,n_M_sub]);
delta_j_wSM___(:,1+S_ij_,1+M_ij_) = reshape(tmp_delta_j_,[n_w_max,n_S_sub,n_M_sub]);
I_value_wSM___(:,1+S_ij_,1+M_ij_) = reshape(tmp_I_value_use_wSM_,[n_w_max,n_S_sub,n_M_sub]);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wSM___: %0.6f',tmp_t)); end;
end;%if (n_S_sub>0);
end;%for nSbatch=0:n_Sbatch-1;
end;%if (n_M_sub>0);
end;%for nMbatch=0:n_Mbatch-1;
if (verbose>1); disp(sprintf(' %% [finished ampmh_X_wSM___4b_debug]')); end;


