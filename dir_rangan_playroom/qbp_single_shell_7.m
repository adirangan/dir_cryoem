function ...
[ ...
 parameter ...
,a_k_Y_lm_ ...
,n_quad_from_data_q_ ...
] = ...
qbp_single_shell_7( ...
 parameter ...
,l_max ...
,n_w ...
,n_M ...
,M_k_p_wM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wC__ ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
);
%%%%%%%%;
% Applies quadrature-back-propagation to solve for a single shell of a_k_Y_. ;
% Associates CTF_k_p_wC__(:,1+index_nCTF_from_nM_(1+nM)) with image M_k_p_wM__(:,1+nM);
% ;
% Input: ;
% l_max: integer order used for a_k_Y_ on shell. ;
% n_w: integer number of inplane_gamma_z values recorded at that ring. ;
% n_M: integer number of images. ;
% M_k_p_wM__: complex array of size (n_w,n_M). stack of images on ring in k_p_ format. ;
% index_nCTF_from_nM_: integer array of size n_M. index_nCTF_from_nM_(1+nM) is the (base 0) CTF_index used for image M_k_p_wM__(:,1+nM). ;
%             This can be empty or set to 1, in which case the same CTF_k_p_ will be used for each image. ;
% CTF_k_p_wC__: complex array of size(n_w,n_CTF). stack of ctf-functions in k_p_ format. ;
%            If index_nCTF_from_nM_ is empty or set to 1, then we assume this contains only a single CTF_k_p_, ;
%            which will then be used for all images. ;
% euler_polar_a_: real array of size n_M. polar_a used for each image ;
% euler_azimu_b_: real array of size n_M. azimu_b used for each image ;
% euler_gamma_z_: real array of size n_M. gamma_z used for each image ;
% ;
% Output: ;
% a_k_Y_lm_: complex array of size n_lm. output functions in k_Y_ format. ;
%%%%%%%%;

str_thisfunction = 'qbp_single_shell_7';

%%%%%%%%;
if (nargin<1);
%%%%%%%%;
flag_verbose=1; nf=0;
if (flag_verbose); disp(sprintf(' %% testing %s',str_thisfunction)); end;
k_p_r_max = 48.0/(2*pi); k_eq_d = 1.0/(2*pi);
n_k_p_r = 1; k_p_r_1 = 1.0; k_p_r_ = k_p_r_1;
[ ...
 n_shell ...
,azimu_b_shell_ ...
,polar_a_shell_ ...
,weight_shell_ ...
,k_c_0_shell_ ...
,k_c_1_shell_ ...
,k_c_2_shell_ ...
] = ...
sample_shell_6( ...
 k_p_r_1 ...
,k_eq_d/k_p_r_max ...
) ;
k_p_r_shell_ = k_p_r_(1+0)*ones(n_shell,1);
%%%%;
l_max_upb = 96;
l_max = max(0,min(l_max_upb,1+ceil(2*pi*k_p_r_max)));
n_lm = (l_max+1).^2;
m_max_ = -l_max : +l_max;
n_m_max = length(m_max_);
Y_l_val_ = zeros(n_lm,1);
Y_m_val_ = zeros(n_lm,1);
tmp_l_val_ = zeros(n_lm,1);
tmp_m_val_ = zeros(n_lm,1);
na=0; 
for l_val=0:l_max;
for m_val=-l_val:+l_val;
tmp_l_val_(1+na) = l_val;
tmp_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
tmp_index_ = 0:n_lm-1;
Y_l_val_(1+tmp_index_) = tmp_l_val_;
Y_m_val_(1+tmp_index_) = tmp_m_val_;
weight_Y_ = ones(n_lm,1);
weight_3d_k_p_r_ = 4*pi;
%%%%;
a_k_Y_form_ = zeros(n_lm,1);
Y_l_use = +28; Y_m_use = -18;
a_k_Y_form_ = +1.0.*(Y_l_val_==Y_l_use).*(Y_m_val_==Y_m_use);
if ~exist('Ylm_uklma___','var'); Ylm_uklma___=[]; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
 a_k_p_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_spharm_to_k_p_4( ...
 flag_verbose ...
,n_shell ...
,[0,n_shell] ...
,k_p_r_shell_ ...
,azimu_b_shell_ ...
,polar_a_shell_ ...
,weight_3d_k_p_r_ ...
,weight_shell_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max ...
,a_k_Y_form_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
%%%%;
[ ...
 a_k_Y_quad_ ...
] = ...
convert_k_p_to_spharm_4( ...
 flag_verbose ...
,n_shell ...
,[0,n_shell] ...
,k_p_r_shell_ ...
,azimu_b_shell_ ...
,polar_a_shell_ ...
,weight_3d_k_p_r_ ...
,weight_shell_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max ...
,a_k_p_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
if (flag_verbose); disp(sprintf(' %% a_k_Y_form_ vs a_k_Y_quad_: %0.16f',fnorm(a_k_Y_form_-a_k_Y_quad_)/fnorm(a_k_Y_form_))); end;
%%%%;
n_w = 2*2*l_max;
viewing_k_eq_d = 0.5*1.0/k_p_r_max;
[ ...
 S_k_p_wS__ ...
,~ ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
] = ...
pm_template_2( ...
 flag_verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_form_ ...
,viewing_k_eq_d ...
,-1 ...
,n_w ...
);
n_S = n_viewing_S;
if (flag_verbose); disp(sprintf(' %% n_w %d, viewing_k_eq_d %0.6f, n_S %d',n_w,viewing_k_eq_d,n_S)); end;
S_k_p_wS__ = reshape(S_k_p_wS__,[n_w,n_S]);
%%%%;
parameter_qbp = struct('type','parameter');
parameter_qbp.flag_verbose = flag_verbose;
[ ...
 parameter_qbp ...
,a_k_Y_0qbp_ ...
,n_quad_from_data_q_ ...
] = ...
qbp_single_shell_7( ...
 parameter_qbp ...
,l_max ...
,n_w ...
,n_S ...
,S_k_p_wS__ ...
,[] ...
,[] ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
);
%%%%;
[ ...
 a_k_p_0qbp_ ...
] = ...
convert_spharm_to_k_p_4( ...
 flag_verbose ...
,n_shell ...
,[0,n_shell] ...
,k_p_r_shell_ ...
,azimu_b_shell_ ...
,polar_a_shell_ ...
,weight_3d_k_p_r_ ...
,weight_shell_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max ...
,a_k_Y_0qbp_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
%%%%;
flag_2d_vs_3d = 0; %<-- 3d;
alim_ = prctile(real(a_k_p_quad_),[5,95]); alim_ = mean(alim_) + 0.5*1.25*diff(alim_)*[-1,+1];
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
imagesc_polar_a_azimu_b_0( ...
 polar_a_shell_ ... 
,azimu_b_shell_ ... 
,real(a_k_p_quad_) ... 
,alim_ ... 
,colormap_81s ... 
,flag_2d_vs_3d ...
,k_p_r_1 ...
);
title('real(a_k_p_quad_)','Interpreter','none');
subplot(1,2,2);
imagesc_polar_a_azimu_b_0( ...
 polar_a_shell_ ... 
,azimu_b_shell_ ... 
,real(a_k_p_0qbp_) ... 
,alim_ ... 
,colormap_81s ... 
,flag_2d_vs_3d ...
,k_p_r_1 ...
);
title('real(a_k_p_0qbp_)','Interpreter','none');
%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,3,1); plot(Y_l_val_,abs(a_k_Y_form_-a_k_Y_0qbp_),'.');
xlabel('Y_l_val_','Interpreter','none'); ylabel('abs(a_k_Y_form_-a_k_Y_0qbp_)','Interpreter','none'); grid on;
subplot(1,3,2); plot(Y_m_val_,abs(a_k_Y_form_-a_k_Y_0qbp_),'.');
xlabel('Y_m_val_','Interpreter','none'); ylabel('abs(a_k_Y_form_-a_k_Y_0qbp_)','Interpreter','none'); grid on;
disp(sprintf(' %% a_k_Y_form_ vs a_k_Y_0qbp_: %0.16f',fnorm(a_k_Y_form_-a_k_Y_0qbp_)/fnorm(a_k_Y_form_)));
subplot(1,3,3); plot(abs(a_k_Y_form_),abs(a_k_Y_0qbp_),'o');
%%%%;
disp(sprintf(' %% returning')); return;
%%%%%%%%;
end;%if (nargin<1);
%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); l_max=[]; end; na=na+1;
if (nargin<1+na); n_w=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wM__=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_wC__=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_=[]; end; na=na+1;
if (nargin<1+na); euler_gamma_z_=[]; end; na=na+1;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
flag_verbose = parameter.flag_verbose;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
tolerance_master = parameter.tolerance_master;
if (~isfield(parameter,'factor_max')); parameter.factor_max = 8; end; %<-- parameter_bookmark. ;
factor_max = parameter.factor_max;
%%%%%%%%;

if isempty(euler_polar_a_); euler_polar_a_ = zeros(n_M,1); end;
if isempty(euler_azimu_b_); euler_azimu_b_ = zeros(n_M,1); end;
if isempty(euler_gamma_z_); euler_gamma_z_ = zeros(n_M,1); end;
if isempty(index_nCTF_from_nM_); index_nCTF_from_nM_ = zeros(n_M,1); end;
if isempty(CTF_k_p_wC__); CTF_k_p_wC__ = ones(n_w,1); end;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
n_lm = (1+l_max).^2; %<-- total number of spherical-harmonic coefficients. ;
CTF_k_p_wM__ = CTF_k_p_wC__(:,1+index_nCTF_from_nM_); %<-- CTF for each image. ;
%%%%%%%%;
% set resolution of quadrature grid. ;
% Mollification involves surface-diffusion: 1/(2*pi)/sigma_use^2 .* exp(-d2/(2*sigma_use^2)). ;
% This corresponds to a diffusion-time of t_diffuse:=(sigma_use^2)/2. ;
% eventual deconvolution involves: exp(+(Y_l_val_).*(1+Y_l_val_).*t_diffuse.^2). ;
% factor_max = exp(+(l_max)*(1+l_max)*t_diffuse);
% t_diffuse_max = log(factor_max)/max(1,l_max*(1+l_max));
% sigma_max = sqrt(log(factor_max)*2/max(1,(l_max*(1+l_max))));
% Note that the first neglected nearest-neighbor will be roughly: ;
% d2_max = (sqrt(3)/2 * n_ring*quad_k_eq_d)^2 away. ;
% This should corresond to a (relatively) negligible mollified quantity: ;
% epsilon_check = exp(-d2_max/(2*sigma_use^2)) = exp(-(n_ring*quad_k_eq_d)^2/(2*sigma_use^2)) < tolerance_master. ;
% n_ring = ceil(sqrt(max(1,-log(tolerance_master)*(2*sigma_use^2)/quad_k_eq_d^2)));
%%%%%%%%;
factor_max = 8;
t_diffuse_max = log(factor_max)/max(1,l_max*(1+l_max));
sigma_max = sqrt(2*t_diffuse_max);
%quad_k_eq_d = sqrt(4*pi./max(1,n_lm)); %<-- older setting from qbp_6. consider changing. ;
quad_k_eq_d = 0.5*sqrt(4*pi./max(1,n_lm)); %<-- older setting from qbp_6. consider changing. ;
sigma_use = min(sigma_max,quad_k_eq_d); %<-- ensure maximum amplification. ;
t_diffuse = (sigma_use^2)/2;
n_ring = ceil(sqrt(max(1,-log(tolerance_master)*(2*sigma_use^2)/quad_k_eq_d^2))); %<-- number of nearest neighbors requested for each point. ;
d2_max = (sqrt(3)/2 * n_ring*quad_k_eq_d)^2;
epsilon_check = exp(-d2_max/(2*sigma_use^2));
assert(epsilon_check<=tolerance_master);
n_nearest = 1+6*n_ring*(n_ring+1)/2; %<-- rough number of neighbors on hexagonal grid (i.e., 1+6+12+18+...). ;
if (flag_verbose); disp(sprintf(' %% quad_k_eq_d %0.6f n_nearest_k %d',quad_k_eq_d,n_nearest)); end;
%%%%;
tmp_t = tic();
[ ...
 quad_n_all ...
,quad_azimu_b_all_ ...
,quad_polar_a_all_ ...
,quad_weight_all_ ...
,quad_k_c_0_all_ ...
,quad_k_c_1_all_ ...
,quad_k_c_2_all_ ...
,~ ...
,~ ...
,~ ...
] = ...
sample_shell_5( ...
 1.0 ...
,quad_k_eq_d ...
,'L' ...
) ;
quad_k_c_qd__ = [ quad_k_c_0_all_ , quad_k_c_1_all_ , quad_k_c_2_all_ ];
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% sample_shell_5 (should be a precomputation): %0.2fs',tmp_t)); end;
if (flag_verbose); disp(sprintf(' %% quad_n_all: %.d',quad_n_all)); end;
%%%%;
tmp_t = tic();
Ylm__ = get_Ylm__(1+l_max,0:l_max,quad_n_all,quad_azimu_b_all_,quad_polar_a_all_);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% get_Ylm__ (should be a precomputation): %0.2fs',tmp_t)); end;
tmp_t = tic();
Ylm_yq__ = zeros(n_lm,quad_n_all);
Y_l_val_ = zeros(n_lm,1);
Y_m_val_ = zeros(n_lm,1);
nml=0;
for l_val=0:l_max;
for m_val=-l_val:+l_val;
Y_l_val_(1+nml) = l_val;
Y_m_val_(1+nml) = m_val;
Ylm_yq__(1+nml,:) = Ylm__{1+l_val}(1+l_val+m_val,:);
nml=nml+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
Ylm_weight_yq__ = Ylm_yq__ * sparse(1:quad_n_all,1:quad_n_all,quad_weight_all_,quad_n_all,quad_n_all);
deconvolve_lm_ = exp(+(Y_l_val_).*(1+Y_l_val_).*t_diffuse);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% Ylm_weight_yq__ (should be a precomputation): %0.2fs',tmp_t)); end;
%%%%;
tmp_t = tic();
[ ...
 data_k_p_polar_a__ ...
,data_k_p_azimu_b__ ...
,data_k_c_0__ ...
,data_k_c_1__ ...
,data_k_c_2__ ...
] = ...
cg_rhs_1( ...
 n_M ...
,n_w ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,+euler_gamma_z_ ...
);
data_k_c_wMd__ = [ data_k_c_0__(:) , data_k_c_1__(:) , data_k_c_2__(:) ];
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% cg_rhs_1: %0.2fs',tmp_t)); end;
%%%%;
tmp_t = tic();
[ij_quad_from_data_wMn__,d1_quad_from_data_wMn__] = knnsearch(quad_k_c_qd__,data_k_c_wMd__,'K',n_nearest);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% knnsearch (technically unnecessary): %0.2fs',tmp_t)); end;
tmp_t = tic();
index_quad_from_data_wMn__ = ij_quad_from_data_wMn__ - 1;
index_data_wMn__ = repmat(transpose(0:n_w*n_M-1),[1,n_nearest]);
d2_quad_from_data_wMn__ = d1_quad_from_data_wMn__.^2;
mollify_quad_from_data_wMn__ = 1/(2*pi)/sigma_use.^2 .* exp(-d2_quad_from_data_wMn__/(2*sigma_use^2));
quad_from_data_qwM__ = sparse(1+index_quad_from_data_wMn__(:),1+index_data_wMn__(:),mollify_quad_from_data_wMn__(:),quad_n_all,n_w*n_M);
n_quad_from_data_q_ = quad_from_data_qwM__*ones(n_w*n_M,1); %<-- effective number of data-points per quadrature-point. ;
CTF2_k_p_q_ = quad_from_data_qwM__*abs(CTF_k_p_wM__(:)).^2;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% quad_from_data_qwM__: %0.2fs',tmp_t)); end;
%%%%;
tmp_t = tic();
a_k_Y_lm_ = zeros(n_lm,1);
M_CTF_k_p_wM__ = M_k_p_wM__.*CTF_k_p_wM__;
quad_from_data_M_CTF_k_p_normalized_q_ = (quad_from_data_qwM__ * M_CTF_k_p_wM__(:))./max(1e-12,CTF2_k_p_q_);
a_k_Y_lm_ = conj(Ylm_weight_yq__)*quad_from_data_M_CTF_k_p_normalized_q_;
a_k_Y_lm_ = a_k_Y_lm_.*deconvolve_lm_; %<-- deconvolve (i.e., unmollify). ;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% a_k_Y_lm_: %0.2fs',tmp_t)); end;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
