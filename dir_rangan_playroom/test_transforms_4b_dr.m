%%%%%%%%;
% An addendum to test_transforms_4_dr.m ;
% sets up a simple template-image comparison ;
% tests ampmh_X_dwSM____8.m ;
%%%%%%%%;

flag_verbose = 1;
flag_disp = 1; nf=0;

%%%%%%%%;
% spatial grid. ;
%%%%%%%%;
n_x = 1+128;
diameter_x_c = 2.0;
half_diameter_x_c = diameter_x_c/2.0;
x_c_0_ = linspace(-half_diameter_x_c,+half_diameter_x_c,n_x);
x_c_1_ = linspace(-half_diameter_x_c,+half_diameter_x_c,n_x);
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_);

%%%%%%%%;
% Set up polar-quadrature-weights. ;
%%%%%%%%;
k_p_r_max = 48/(2*pi); k_eq_d = 0.5/(2*pi); TorL = 'L';
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,TorL ...
);
%%%%;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
n_w_max = 2*2*(l_max_max+1); n_w_0in_ = n_w_max*ones(n_k_p_r,1);
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_0in_ ...
);
%%%%%%%%;
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);

%%%%%%%%;
% Use ampmh_X_wSM___8 to estimate integral of plane-wave times bessel times plane-wave. ;
% Note that ampmh_X_wSM___8 assumes an isotropic CTF. ;
% I.e., the inputs CTF_UX_S_k_q_wnS__ and CTF_UX_S_l2_ ;
% combine an angularly-independent CTF with the templates S_k_q_wnS__. ;
%%%%%%%%;
rng(0);
flag_CTF_use = 1;
flag_delta_use = 1;
n_S = 7; n_M = 5;
if (flag_verbose); disp(sprintf(' %% testing ampmh_X_wSM___8.m using n_S %d templates and n_M %d images. ',n_S,n_M)); end;
delta_max = 0.1;
delta_S_S_ = rand(n_S,1)*delta_max;
omega_S_S_ = 2*pi*rand(n_S,1);
delta_S_2S__ = zeros(2,n_S);
S_k_p_wkS__ = zeros(n_w_sum,n_S);
S_k_q_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
tmp_delta_S = delta_S_S_(1+nS);
tmp_omega_S = omega_S_S_(1+nS);
tmp_delta_S_ = tmp_delta_S * [cos(tmp_omega_S) ; sin(tmp_omega_S)];
delta_S_2S__(:,1+nS) = tmp_delta_S_;
S_k_p_wkS__(:,1+nS) = exp(2*pi*i*k_p_r_wk_.*tmp_delta_S.*cos(k_p_w_wk_-tmp_omega_S));
S_k_q_wkS__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS));
end;%for nS=0:n_S-1;
delta_M_M_ = rand(n_M,1)*delta_max;
omega_M_M_ = 2*pi*rand(n_M,1);
delta_M_M_(1+0) = delta_S_S_(1+0); omega_M_M_(1+0) = omega_S_S_(1+0); %<-- ensure that the first image matches the first template. ;
delta_M_2M__ = zeros(2,n_M);
M_k_p_wkM__ = zeros(n_w_sum,n_M);
M_k_q_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
tmp_delta_M = delta_M_M_(1+nM);
tmp_omega_M = omega_M_M_(1+nM);
tmp_delta_M_ = tmp_delta_M * [cos(tmp_omega_M) ; sin(tmp_omega_M)];
delta_M_2M__(:,1+nM) = tmp_delta_M_;
M_k_p_wkM__(:,1+nM) = exp(2*pi*i*k_p_r_wk_.*tmp_delta_M.*cos(k_p_w_wk_-tmp_omega_M));
M_k_q_wkM__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wkM__(:,1+nM));
end;%for nM=0:n_M-1;
CTF_alpha = flag_CTF_use*1.5;
CTF_k_p_r_k_ = besselj(0,CTF_alpha*k_p_r_);
CTF_k_p_r_wk_ = besselj(0,CTF_alpha*k_p_r_wk_); CTF_wk_ = CTF_k_p_r_wk_;
%%%%%%%%;
% Prepare FTK. ;
%%%%%%%%;
delta_r_max = flag_delta_use*0.05; n_delta_v_requested = flag_delta_use*16;
delta_r_p = 0.05;
delta_r_s = delta_r_max/sqrt(2*log(1/delta_r_p));
delta_r_N = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2));
svd_eps = 1e-6;
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested,[],[],32,64,65);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
disp(sprintf(' %% p-val %0.4f delta_r_max %0.6f sigma %0.4f N_pixel %0.4f --> FTK.n_svd_l %d, n_delta_v_requested %d',delta_r_p,delta_r_max,delta_r_s,delta_r_N,FTK.n_svd_l,n_delta_v_requested));
%%%%%%%%;
% Prepare principal-modes. ;
%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions.; 
UX_kn__ = zeros(n_k_p_r,n_UX_rank);
X_weight_r_ = zeros(n_k_p_r,1);
[ ...
 X_kk__ ...
,X_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
);
[tmp_UX_kn__,tmp_SX_k__,tmp_VX_kn__] = svds(X_kk__,n_UX_rank);
UX_kn__(:,1:n_UX_rank) = tmp_UX_kn__(:,1:n_UX_rank);
pm_n_UX_rank = n_UX_rank; %<-- just to check dimension. ;
%%%%%%%%;
% Prepare quasi-images. ;
%%%%%%%%;
n_UX_rank = n_k_p_r;
tmp_t = tic();
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M,M_k_q_wkM__,pm_n_UX_rank,UX_kn__,X_weight_r_);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the non-translated images. ;
%%%%%%%%;
tmp_t = tic();
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M,pm_n_UX_rank,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
disp(sprintf(' %% average l2-norm of images: %0.16f %%<-- should be 1 ',mean(UX_M_l2_dM__(:))/(pi*k_p_r_max^2)));
%%%%%%%%;
% Prepare CTF_UX_S_k_q_wnS__. ;
%%%%%%%%;
S_k_q_wSk___ = permute(reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
CTF_UX_S_k_q_wnS__ = reshape(permute(reshape(reshape(S_k_q_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_weight_r_.*CTF_k_p_r_k_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
CTF_UX_S_l2_S_ = sum(abs(CTF_UX_S_k_q_wnS__).^2,1)/n_w_max;
%%%%%%%%;
% Calculate ampmh_X_wSM___8. ;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_compute_I_value = 0;
tmp_t = tic();
[ ...
 parameter ...
,X_wSM_ampm___ ...
,delta_x_wSM___ ...
,delta_y_wSM___ ...
,gamma_z_wSM___ ...
,I_value_wSM___ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_S_ ...
,n_M ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% X_wSM___: %0.3fs',tmp_t)); end;
X_wSM_ampm___ = real(X_wSM_ampm___);
%%%%%%%%;
% compare against a la carte calculation. ;
%%%%%%%%;
nM=2;
nS=3;
nw=71;
tmp_delta_S = delta_S_S_(1+nS); tmp_omega_S = omega_S_S_(1+nS);
tmp_delta_S_ = tmp_delta_S * [cos(tmp_omega_S) ; sin(tmp_omega_S)];
tmp_delta_M = delta_M_M_(1+nM); tmp_omega_M = omega_M_M_(1+nM);
tmp_delta_M_ = tmp_delta_M * [cos(tmp_omega_M) ; sin(tmp_omega_M)];
gamma_z = (2*pi*nw)/n_w_max; assert(gamma_z==gamma_z_wSM___(1+nw,1+nS,1+nM));
cc = cos(+gamma_z); sc = sin(+gamma_z);
Rz = [+cc , -sc ; +sc , +cc];
delta_x = delta_x_wSM___(1+nw,1+nS,1+nM);
delta_y = delta_y_wSM___(1+nw,1+nS,1+nM);
if (flag_verbose>0); disp(sprintf(' %% nM %.2d nS %.2d nw %.3d delta [%0.4f,%0.4f]',nM,nS,nw,delta_x,delta_y)); end;
tmp_delta_d_ = [delta_x;delta_y];
tmp_delta_d = fnorm(tmp_delta_d_);
tmp_omega_d = atan2(delta_y,delta_x);
tmp_delta_TM_ = tmp_delta_M_ - tmp_delta_d_; %<-- multiply image by exp(-2*pi*i*dot(k_,delta_)). ;
tmp_delta_TM = fnorm(tmp_delta_TM_);
tmp_omega_TM = atan2(tmp_delta_TM_(1+1),tmp_delta_TM_(1+0));
tmp_delta_RS_ = Rz*tmp_delta_S_; %<-- rotate template by +gamma_z. ;
tmp_delta_RS = fnorm(tmp_delta_RS_);
tmp_omega_RS = atan2(tmp_delta_RS_(1+1),tmp_delta_RS_(1+0));
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
RS_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,+gamma_z);
RS_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,RS_k_p_wk_);
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
M_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_);
TM_k_q_wk_ = transf_svd_q_to_q_FTK_6(FTK,n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_q_wk_,delta_x,delta_y);
tmp_I = plane_bessel_plane_integral_0(k_p_r_max,tmp_delta_RS,tmp_omega_RS,tmp_delta_TM,tmp_omega_TM,CTF_alpha);
TM_k_p_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,TM_k_q_wk_);
tmp_Q = sum(conj(CTF_wk_.*RS_k_p_wk_).*TM_k_p_wk_.*weight_2d_wk_,'all')*(2*pi)^2;
if (flag_verbose>0); disp(sprintf(' %% tmp_I vs tmp_Q: %0.16f %%<-- should be <1e-2 ',fnorm(tmp_I-tmp_Q)/max(1e-12,fnorm(tmp_I)))); end;
%%%%;
if flag_disp>1;
M_x_c_xx__ = interp_k_p_to_x_c_xxnufft(n_x,diameter_x_c,n_x,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_k_p_wk_);
TM_x_c_xx__ = interp_k_p_to_x_c_xxnufft(n_x,diameter_x_c,n_x,diameter_x_c,n_k_p_r,k_p_r_,n_w_,TM_k_p_wk_);
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
imagesc_c(n_x,x_c_0_,n_x,x_c_1_,real(M_x_c_xx__)); axis image; axisnotick; title('original M');
subplot(1,2,2);
imagesc_c(n_x,x_c_0_,n_x,x_c_1_,real(TM_x_c_xx__)); axis image; axisnotick; title('translated M');
end;%if flag_disp>1;
%%%%;

%%%%%%%%;
% Calculate true landscape of innerproducts for the same set of translations. ;
%%%%%%%%;
X_wSM_form___ = zeros(n_w_max,n_S,n_M);
X_wSM_quad___ = zeros(n_w_max,n_S,n_M);
for nS=0:n_S-1;
if (flag_verbose>0); disp(sprintf(' %% nS %d/%d',nS,n_S)); end;
tmp_delta_S = delta_S_S_(1+nS);
tmp_omega_S = omega_S_S_(1+nS);
tmp_delta_S_ = tmp_delta_S * [cos(tmp_omega_S) ; sin(tmp_omega_S)];
for nM=0:n_M-1;
tmp_delta_M = delta_M_M_(1+nM);
tmp_omega_M = omega_M_M_(1+nM);
tmp_delta_M_ = tmp_delta_M * [cos(tmp_omega_M) ; sin(tmp_omega_M)];
for nw=0:n_w_max-1;
gamma_z = (2*pi*nw)/n_w_max; assert(gamma_z==gamma_z_wSM___(1+nw,1+nS,1+nM));
cc = cos(+gamma_z); sc = sin(+gamma_z);
Rz = [+cc , -sc ; +sc , +cc];
delta_x = delta_x_wSM___(1+nw,1+nS,1+nM);
delta_y = delta_y_wSM___(1+nw,1+nS,1+nM);
tmp_delta_d_ = [delta_x;delta_y];
tmp_delta_d = fnorm(tmp_delta_d_);
tmp_omega_d = atan2(delta_y,delta_x);
tmp_delta_TM_ = tmp_delta_M_ - tmp_delta_d_; %<-- multiply image by exp(-2*pi*i*dot(k_,delta_)). ;
tmp_delta_TM = fnorm(tmp_delta_TM_);
tmp_omega_TM = atan2(tmp_delta_TM_(1+1),tmp_delta_TM_(1+0));
tmp_delta_RS_ = Rz*tmp_delta_S_; %<-- rotate template by +gamma_z. ;
tmp_delta_RS = fnorm(tmp_delta_RS_);
tmp_omega_RS = atan2(tmp_delta_RS_(1+1),tmp_delta_RS_(1+0));
%%%%;
%{
tmp_I = integral2( ...
@(k_p_r,psi) k_p_r ...
 .* conj(besselj(0,CTF_alpha*k_p_r).*exp(2*pi*i*k_p_r.*tmp_delta_RS.*cos(psi-tmp_omega_RS))) ...
 .* exp(2*pi*i*k_p_r.*tmp_delta_TM.*cos(psi-tmp_omega_TM)) ...
,0,k_p_r_max,0,2*pi);
%}
%%%%;
RS_wk_ = exp(2*pi*i*k_p_r_wk_.*tmp_delta_RS.*cos(k_p_w_wk_-tmp_omega_RS));
TM_wk_ = exp(2*pi*i*k_p_r_wk_.*tmp_delta_TM.*cos(k_p_w_wk_-tmp_omega_TM));
tmp_Q = sum(conj(CTF_wk_.*RS_wk_).*TM_wk_.*weight_2d_wk_,'all')*(2*pi)^2;
tmp_RS_l2 = sum(conj(CTF_wk_.*RS_wk_).*(CTF_wk_.*RS_wk_).*weight_2d_wk_,'all')*(2*pi)^2;
tmp_TM_l2 = sum(conj(TM_wk_).*(TM_wk_).*weight_2d_wk_,'all')*(2*pi)^2;
%%%%;
tmp_I = plane_bessel_plane_integral_0(k_p_r_max,tmp_delta_RS,tmp_omega_RS,tmp_delta_TM,tmp_omega_TM,CTF_alpha);
%%%%;
X_wSM_form___(1+nw,1+nS,1+nM) = tmp_I/max(1e-12,sqrt(tmp_RS_l2*tmp_TM_l2));
X_wSM_quad___(1+nw,1+nS,1+nM) = tmp_Q/max(1e-12,sqrt(tmp_RS_l2*tmp_TM_l2));
end;%for nw=0:n_w_max-1;
end;%for nM=0:n_M-1;
end;%for nS=0:n_S-1;
if (flag_verbose>0); disp(sprintf(' %% X_wSM_form___ vs X_wSM_quad___: %0.16f %%<-- should be <1e-2 ',fnorm(X_wSM_form___-X_wSM_quad___)/max(1e-12,fnorm(X_wSM_form___)))); end;
if (flag_verbose>0); disp(sprintf(' %% X_wSM_form___ vs X_wSM_ampm___: %0.16f %%<-- should be <1e-2 ',fnorm(X_wSM_form___-X_wSM_ampm___)/max(1e-12,fnorm(X_wSM_form___)))); end;
%%%%%%%%;

%%%%%%%%;
% Now prepare a single image and translate by a particular [delta_x,delta_y];
%%%%%%%%;
delta_r_max = 0.15; svd_eps = 1e-6; n_delta_v_requested = 2;
omega = -pi/6;
delta_x_ = [0;0.95*delta_r_max*cos(omega)];
delta_y_ = [0;0.95*delta_r_max*sin(omega)];
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested,delta_x_,delta_y_,32,64,65);
%%%%;
n_source = 7;
source_0_ = 2*rand(n_source,1)-1; source_0_(1+0) = 0.0;
source_1_ = 2*rand(n_source,1)-1; source_1_(1+0) = 0.0;
source_s_ = 0.25*rand(n_source,1); source_s_(1+0) = 1.0;
source_sigma = 1/128;
N_x_c_xx__ = zeros(n_x,n_x);
for nsource=0:n_source-1;
source_0 = source_0_(1+nsource);
source_1 = source_1_(1+nsource);
source_s = source_s_(1+nsource);
N_x_c_xx__ = N_x_c_xx__ + source_s*exp(-( (x_c_0__-source_0).^2 + (x_c_1__-source_1).^2 )/(2*source_sigma^2));
end;%for nsource=0:n_source-1;
M_k_p_wk_ = interp_x_c_to_k_p_xxnufft(n_x,diameter_x_c,n_x,diameter_x_c,N_x_c_xx__,n_k_p_r,k_p_r_,n_w_);
M_x_c_xx__ = interp_k_p_to_x_c_xxnufft(n_x,diameter_x_c,n_x,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_k_p_wk_);
Mlim_ = prctile(real(M_x_c_xx__(:)),[ 5,95]);
Mlim_ = mean(Mlim_) + 1.25*0.5*diff(Mlim_)*[-1,+1];
%%%%;
delta_x = delta_x_(1+1);
delta_y = delta_y_(1+1);
M_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_);
TM_k_q_wk_ = transf_svd_q_to_q_FTK_6(FTK,n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_q_wk_,delta_x,delta_y);
TM_k_p_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,TM_k_q_wk_);
TM_x_c_xx__ = interp_k_p_to_x_c_xxnufft(n_x,diameter_x_c,n_x,diameter_x_c,n_k_p_r,k_p_r_,n_w_,TM_k_p_wk_);
%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figmed;
markersize_use = 6;
subplot(1,2,1);
hold on;
imagesc_c(n_x,x_c_0_,n_x,x_c_1_,real(M_x_c_xx__),Mlim_); axis image; axisnotick; title('original M (source at triangle)');
for nsource=0:n_source-1;
source_0 = source_0_(1+nsource);
source_1 = source_1_(1+nsource);
plot(source_0,source_1,'k^','MarkerSize',markersize_use,'MarkerFaceColor','b');
plot(source_0+delta_x,source_1+delta_y,'kd','MarkerSize',markersize_use,'MarkerFaceColor','c');
end;%for nsource=0:n_source-1;
hold off;
subplot(1,2,2);
hold on;
imagesc_c(n_x,x_c_0_,n_x,x_c_1_,real(TM_x_c_xx__),Mlim_); axis image; axisnotick; title('translated M (source at diamond)');
for nsource=0:n_source-1;
source_0 = source_0_(1+nsource);
source_1 = source_1_(1+nsource);
plot(source_0,source_1,'k^','MarkerSize',markersize_use,'MarkerFaceColor','b');
plot(source_0+delta_x,source_1+delta_y,'kd','MarkerSize',markersize_use,'MarkerFaceColor','c');
end;%for nsource=0:n_source-1;
hold off;
sgtitle(sprintf('delta = [%+0.6f,%+0.6f]',delta_x,delta_y));
end;%if flag_disp>0;
%%%%;


