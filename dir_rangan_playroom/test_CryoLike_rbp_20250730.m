function nf = test_CryoLike_rbp_20250730(k_eq_d_double,k_int,viewing_k_eq_d_double,nf);
%%%%%%%%;
% Tests raw back-propagation onto spherical-grid. ;
%%%%%%%%;
% Input: ;
% k_eq_d_double: double resolution-parameter for spherical-grid. ;
%     k_eq_d_double is defined as (2*pi)*k_eq_d. ;
%     k_eq_d is the distance (along the equator) between sample-point in k-space on the largest k-shell. ;
% k_int: int largest wavenumber. ;
%     k_p_r_max is defined as k_int/(2*pi). ;
% viewing_k_eq_d_double: double viewing-angle resolution-parameter. ;
%     viewing_k_eq_d = viewing_k_eq_d_double*1.0/max(1e-12,k_p_r_max). ;
% nf: int figure number for plotting. ;
%%%%%%%%;
% Additional internal parameters. ;
% n_source: number of plane-wave sources. ;
% n_w_int: int prefactor for determining n_w_max (default 1.0). ;
%%%%%%%%;
% defaults: k_eq_d_double = 1.0; k_int = 48; nf=0;

str_thisfunction = 'test_CryoLike_rbp_20250730';

if (nargin<1);
viewing_k_eq_d_double_ = [1.0,0.5]; n_viewing_k_eq_d_double = numel(viewing_k_eq_d_double_);
nf=0;
for nviewing_k_eq_d_double=0:n_viewing_k_eq_d_double-1;
viewing_k_eq_d_double = viewing_k_eq_d_double_(1+nviewing_k_eq_d_double);
nf = test_CryoLike_rbp_20250730([],[],viewing_k_eq_d_double,nf);
end;%for nviewing_k_eq_d_double=0:n_viewing_k_eq_d_double-1;
disp(sprintf(' %% returning')); return;
end;%if (nargin<1);

na=0;
if nargin<1+na; k_eq_d_double = []; end; na=na+1;
if nargin<1+na; k_int = []; end; na=na+1;
if nargin<1+na; viewing_k_eq_d_double = []; end; na=na+1;
if nargin<1+na; nf = []; end; na=na+1;

% defaults for testing: ;
% clear; str_thisfunction = 'test_CryoLike_rbp_20250730'; k_eq_d_double = []; k_int = []; viewing_k_eq_d_double = []; nf = [];

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

rng(0);
flag_recalc=0; flag_replot=0;
flag_verbose=1; flag_disp=1;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(k_eq_d_double); k_eq_d_double = 1.0; end;
if isempty(k_int); k_int = 48; end;
if isempty(viewing_k_eq_d_double); viewing_k_eq_d_double = 1.0; end;
if isempty(nf); nf = 0; end;
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% setting k_eq_d_double = %0.6f %%<-- should be roughly 1.0 or less for accurate integration. ',k_eq_d_double)); end;
if (flag_verbose>0); disp(sprintf(' %% setting k_int = %d ',k_int)); end;
if (flag_verbose>0); disp(sprintf(' %% setting viewing_k_eq_d_double = %0.6f %%<-- should be roughly 1.0 or less for accurate integration. ',viewing_k_eq_d_double)); end;

%%%%%%%%;
% Now set up k-quadrature on sphere. ;
%%%%%%%%;
k_p_r_max = k_int/(2*pi); k_eq_d = k_eq_d_double/(2*pi); str_T_vs_L = 'L';
flag_uniform_over_n_k_p_r = 1; flag_uniform_over_polar_a = 0; %<-- This is set to match test_ssnll_from_a_k_Y_12 ;
[ ...
 n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_qk_ ...
,k_c_1_qk_ ...
,k_c_2_qk_ ...
,~ ...
,~ ...
,~ ...
,~ ...
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
] = ...
sample_sphere_7( ...
 0*flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_T_vs_L ...
,flag_uniform_over_n_k_p_r ...
,flag_uniform_over_polar_a ...
) ; %<-- sum(weight_3d_k_p_r_)*(4*pi) = (4/3)*pi*k_p_r_max^3 --> sum(weight_3d_k_p_r_) = (1/3)*k_p_r_max^3 ;
%%%%%%%%;
% extract outer shell for visualization. ;
%%%%%%%%;
index_outer_shell_ = n_qk_csum_(1+n_k_p_r-1):n_qk_csum_(1+n_k_p_r-1+1)-1;
n_outer_shell = numel(index_outer_shell_);
k_outer_shell_r_q_ = k_p_r_qk_(1+index_outer_shell_);
assert(std(k_outer_shell_r_q_,1)< 1e-6);
k_outer_shell_azimu_b_q_ = k_p_azimu_b_qk_(1+index_outer_shell_);
k_outer_shell_polar_a_q_ = k_p_polar_a_qk_(1+index_outer_shell_);
%%%%%%%%;
% Set up k-quadrature on disc. ;
%%%%%%%%;
n_w_int = 1;
l_max_upb = round(2*pi*k_p_r_max);
n_w_max = n_w_int*2*(l_max_upb+1);
template_k_eq_d = -1;
n_w_0in_ = n_w_max*ones(n_k_p_r,1);
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_p_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 0*flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
,weight_3d_k_p_r_ ...
);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;

%%%%%%%%;
% Now set up volume with plane-wave sources. ;
%%%%%%%%;
delta_a_c_3s__ = [  ...
 +1.5 , -0.5 ...
;-0.5 , -1.5 ...
;+0.3 , +2.0 ...
] / 2 / max(1e-12,k_p_r_max) ;
n_source = size(delta_a_c_3s__,2);
a_k_p_form_ = zeros(n_qk,1);
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c_3s__(:,1+nsource);
a_k_p_form_ = a_k_p_form_ + exp(+i*2*pi*(k_c_0_qk_*delta_a_c_(1+0) + k_c_1_qk_*delta_a_c_(1+1) + k_c_2_qk_*delta_a_c_(1+2)));
end;%for nsource=0:n_source-1;
%%%%%%%%;

%%%%%%%%;
% select some viewing-angles. ;
%%%%%%%%;
viewing_k_eq_d = viewing_k_eq_d_double*1.0/max(1e-12,k_p_r_max);
flag_uniform_over_polar_a = 0;
[ ...
 n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
] = ...
sample_shell_6( ...
 1.0 ...
,viewing_k_eq_d ...
,'L' ...
,flag_uniform_over_polar_a ...
) ; %<-- obtain viewing angles on outer shell. ;
n_S = n_viewing_S;
if (flag_verbose>0); disp(sprintf(' %% viewing_k_eq_d_double: %0.2f, n_S: %d',viewing_k_eq_d_double,n_S)); end;

%%%%%%%%;
% define rotations in 2d and 3d. ;
%%%%%%%%;
R2 = @(gamma_z) ...
[ +cos(gamma_z) -sin(gamma_z) ; ...
  +sin(gamma_z) +cos(gamma_z) ; ...
] ;
%%%%%%%%;
Rz = @(azimu_b) ...
[ +cos(azimu_b) -sin(azimu_b) 0 ; ...
  +sin(azimu_b) +cos(azimu_b) 0 ; ...
   0             0            1 ; ...
] ;
%%%%%%%%;
Ry = @(polar_a) ...
[ +cos(polar_a) 0 +sin(polar_a) ; ...
   0            1  0            ; ...
  -sin(polar_a) 0 +cos(polar_a) ; ...
];
%%%%%%%%;

%%%%%%%%;
% Now generate the templates analytically. ;
%%%%%%%%;
S_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
tmp_azimu_b = viewing_azimu_b_S_(1+nS);
tmp_polar_a = viewing_polar_a_S_(1+nS);
tmp_gamma_z = 0.0;
tmp_R__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
S_k_p_wk_ = zeros(n_w_sum,1);
for nsource=0:n_source-1;
tmp_delta_ = tmp_R__*delta_a_c_3s__(:,1+nsource);
S_k_p_wk_ = S_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_(1+0) + k_c_1_wk_*tmp_delta_(1+1)));
end;%for nsource=0:n_source-1;
S_k_p_wkS__(:,1+nS) = S_k_p_wk_;
end;%for nS=0:n_S-1;
%%%%%%%%;

%%%%;
% Now we generate a selection of synthetic-images, copied from the templates, ;
% as well as a collection of anisotropic planar-CTF-functions. ;
%%%%;
n_M = n_S*2;
index_nS_from_nM_ = max(0,min(n_S-1,floor(mod(transpose(0:n_M-1),n_S))));
euler_polar_a_M_ = viewing_polar_a_S_(1+index_nS_from_nM_);
euler_azimu_b_M_ = viewing_azimu_b_S_(1+index_nS_from_nM_);
euler_gamma_z_M_ = 2*pi*transpose(0:n_M-1)/max(1,n_M) + pi/2;
n_CTF = 8; %<-- some number of CTF-functions. ;
CTF_phi_C_ = zeros(n_CTF,1);
CTF_k_p_wkC__ = zeros(n_w_sum,n_CTF);
for nCTF=0:n_CTF-1;
CTF_phi = 2*pi*nCTF/max(1,n_CTF) + pi/2;
CTF_phi_C_(1+nCTF) = CTF_phi;
CTF_k_p_wk_ = 2*k_p_r_wk_.*cos(k_p_w_wk_ - CTF_phi);
CTF_k_p_wkC__(:,1+nCTF) = CTF_k_p_wk_;
end;%for nCTF=0:n_CTF-1;
index_nCTF_from_nM_ = mod(transpose([0:n_M-1]),n_CTF);
M_k_p_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
nCTF = index_nCTF_from_nM_(1+nM);
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+nCTF);
nS = index_nS_from_nM_(1+nM);
gamma_z = euler_gamma_z_M_(1+nM);
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
M_k_p_wk_ = CTF_k_p_wk_.*rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,+gamma_z);
M_k_p_wkM__(:,1+nM) = M_k_p_wk_;
end;%for nM=0:n_M-1;
%%%%;
% alternatively, since n_w_ is uniform, we can perform the rotations all at once. ;
%%%%;
tmp_q_ = periodize(transpose(0:n_w_max-1),-n_w_max/2,+n_w_max/2);
N_k_p_wkM__ = zeros(n_w_sum,n_M);
N_k_p_wkM__ = ...
CTF_k_p_wkC__(:,1+index_nCTF_from_nM_) ...
.* ...
reshape( ...
	 ifft( ...
	       bsxfun(@times ...
		      ,fft(reshape(S_k_p_wkS__(:,1+index_nS_from_nM_),[n_w_max,n_k_p_r,n_M]),[],1+0) ...
		      ,reshape(exp(-i*bsxfun(@times,reshape(tmp_q_,[n_w_max,1]),reshape(euler_gamma_z_M_,[1,n_M]))),[n_w_max,1,n_M]) ...
		     ) ...
	       ,[],1+0 ...
	     ) ...
	 ,[n_w_sum,n_M] ...
       );
%%%%%%%%;
fnorm_disp(flag_verbose,'N_k_p_wkM__',N_k_p_wkM__,'M_k_p_wkM__',M_k_p_wkM__,' %%<-- should be zero');

%%%%%%%%;
% Now we perform a quick and dirty 'raw' back-propagation. ;
% Note that this assumes that n_w_ is uniform (i.e., that n_w_ == n_w_max*ones(n_k_p_r,1)). ;
% This also assumes that the spherical-grid is uniform (i.e., that flag_uniform_over_n_k_p_r==1). ;
%%%%%%%%;
CTF_k_p_wkM__ = CTF_k_p_wkC__(:,1+index_nCTF_from_nM_);
M0C2_k_p_wkM__ = abs(CTF_k_p_wkM__).^2;
M1C1_k_p_wkM__ = M_k_p_wkM__.*CTF_k_p_wkM__;
M0C2_k_p_wMk__ = reshape(permute(reshape(M0C2_k_p_wkM__,[n_w_max,n_k_p_r,n_M]),[1,3,2]),[n_w_max*n_M,n_k_p_r]);
M1C1_k_p_wMk__ = reshape(permute(reshape(M1C1_k_p_wkM__,[n_w_max,n_k_p_r,n_M]),[1,3,2]),[n_w_max*n_M,n_k_p_r]);
%%%%;
n_3 = 3;
[ ...
 k_p_polar_a_wM__ ...
,k_p_azimu_b_wM__ ...
,k_c_0_wM__ ...
,k_c_1_wM__ ...
,k_c_2_wM__ ...
] = ...
cg_rhs_2( ...
 n_M ...
,n_w_max ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
);
data_wM3___ = cat(3,k_c_0_wM__,k_c_1_wM__,k_c_2_wM__);
%%%%;
k_outer_shell_c_0_q_ = k_outer_shell_r_q_.*cos(k_outer_shell_azimu_b_q_).*sin(k_outer_shell_polar_a_q_);
k_outer_shell_c_1_q_ = k_outer_shell_r_q_.*sin(k_outer_shell_azimu_b_q_).*sin(k_outer_shell_polar_a_q_);
k_outer_shell_c_2_q_ = k_outer_shell_r_q_.*cos(k_outer_shell_polar_a_q_);
outer_shell_q3__ = cat(2,k_outer_shell_c_0_q_,k_outer_shell_c_1_q_,k_outer_shell_c_2_q_);
%%;
% we can use knnsearch: ;
%%;
tmp_t = tic();
ij_outer_shell_from_data_wM_ = knnsearch(outer_shell_q3__,reshape(data_wM3___,n_w_max*n_M,n_3));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% knnsearch: time %0.6fs',tmp_t)); end;;
index_outer_shell_from_data_wM_ = ij_outer_shell_from_data_wM_ - 1;
outer_shell_from_data_qwM__ = sparse(1+index_outer_shell_from_data_wM_,1:n_w_max*n_M,1,n_outer_shell,n_w_max*n_M);
%%;
% or a simple angular-binning (possible because str_T_vs_L == 'T', rather than 'L');
%%;
tmp_t = tic();
polar_a_a_ = polar_a_ka__{1+0}; n_polar_a = numel(polar_a_a_);
assert(std(diff(polar_a_a_),1)<1e-12);
dpolar_a = mean(diff(polar_a_a_)); %assert(abs(pi/(n_polar_a-1)+dpolar_a)<=1e-12); %<-- assertion only true for str_T_vs_L = 'C' or 'C2'. ;
n_azimu_b_a_ = n_azimu_b_ka__{1+0}; dazimu_b_a_ = (2*pi)./max(1,n_azimu_b_a_);
n_azimu_b_csum_ = cumsum([0;n_azimu_b_a_]); assert(n_outer_shell==n_azimu_b_csum_(1+n_polar_a));
npolar_a_from_data_wM_ = reshape(max(0,min(n_polar_a-1,round((pi-k_p_polar_a_wM__)/max(1e-12,-dpolar_a)))),[n_w_max*n_M,1]); %<-- note that here k_p_polar_a_ ranges from pi down to 0. ;
nazimu_b_from_data_wMa__ = zeros(n_w_max*n_M,n_polar_a);
for npolar_a=0:n_polar_a-1;
n_azimu_b = n_azimu_b_a_(1+npolar_a); dazimu_b = dazimu_b_a_(1+npolar_a);
nazimu_b_from_data_wM_ = reshape(periodize(round(k_p_azimu_b_wM__/max(1e-12,dazimu_b)),0,n_azimu_b),[n_w_max*n_M,1]); %<-- periodize conflicts with trivial vectorization. ;
nazimu_b_from_data_wMa__(:,1+npolar_a) = nazimu_b_from_data_wM_;
end;%for npolar_a=0:n_polar_a-1;
index_babin_shell_from_data_wM_ = zeros(n_w_max*n_M,1);
for npolar_a=0:n_polar_a-1;
tmp_index_ = efind(npolar_a_from_data_wM_==npolar_a);
index_babin_shell_from_data_wM_(1+tmp_index_) = n_azimu_b_csum_(1+npolar_a) + nazimu_b_from_data_wMa__(1+tmp_index_,1+npolar_a);
end;%for npolar_a=0:n_polar_a-1;
babin_shell_from_data_qwM__ = sparse(1+index_babin_shell_from_data_wM_,1:n_w_max*n_M,1,n_outer_shell,n_w_max*n_M);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% angular-binning: time %0.6fs',tmp_t)); end;;
%%%%;
a_k_p_outer_ = zeros(n_qk,1);
for nk_p_r=0:n_k_p_r-1;
index_single_shell_ = n_qk_csum_(1+nk_p_r):n_qk_csum_(1+nk_p_r+1)-1;
n_single_shell = numel(index_single_shell_); assert(n_single_shell==n_outer_shell);
M0C2_single_shell_q_ = outer_shell_from_data_qwM__*M0C2_k_p_wMk__(:,1+nk_p_r);
M1C1_single_shell_q_ = outer_shell_from_data_qwM__*M1C1_k_p_wMk__(:,1+nk_p_r);
a_single_shell_q_ = M1C1_single_shell_q_./max(1e-12,M0C2_single_shell_q_);
a_k_p_outer_(1+index_single_shell_) = a_single_shell_q_;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_outer_',a_k_p_outer_,' %<-- total error can be 1e-1 or so');
%%%%%%%%;
% Note that, because flag_uniform_over_n_k_p_r==1 and n_w_ is uniform, ;
% a_k_p_outer_ can be calculated using matrix-multiplication: ;
%%%%%%%%;
a_k_p_outer_ = reshape ((outer_shell_from_data_qwM__*M1C1_k_p_wMk__)./max(1e-12,outer_shell_from_data_qwM__*M0C2_k_p_wMk__),[n_qk,1]);
fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_outer_',a_k_p_outer_,' %<-- total error can be 1e-1 or so');
%%%%%%%%;
a_k_p_babin_ = reshape ((babin_shell_from_data_qwM__*M1C1_k_p_wMk__)./max(1e-12,babin_shell_from_data_qwM__*M0C2_k_p_wMk__),[n_qk,1]);
fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_babin_',a_k_p_babin_,' %<-- total error can be 1e-1 or so');
%%%%%%%%;

%%%%%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 1; p_col = 3; np=0;
alim_ = max(abs(a_k_p_form_))*[-1,+1]; flag_2d_vs_3d = 0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc_polar_a_azimu_b_0(k_outer_shell_polar_a_q_,k_outer_shell_azimu_b_q_,real(a_k_p_form_(1+index_outer_shell_)),alim_,colormap_beach(),flag_2d_vs_3d,k_p_r_max);
axis image; axisnotick3d; title('a_k_p_form_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc_polar_a_azimu_b_0(k_outer_shell_polar_a_q_,k_outer_shell_azimu_b_q_,real(a_k_p_outer_(1+index_outer_shell_)),alim_,colormap_beach(),flag_2d_vs_3d,k_p_r_max);
axis image; axisnotick3d; title('a_k_p_outer_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc_polar_a_azimu_b_0(k_outer_shell_polar_a_q_,k_outer_shell_azimu_b_q_,real(a_k_p_babin_(1+index_outer_shell_)),alim_,colormap_beach(),flag_2d_vs_3d,k_p_r_max);
axis image; axisnotick3d; title('a_k_p_babin_','Interpreter','none');
%%%%;
end;%if flag_disp>0;
%%%%%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;




