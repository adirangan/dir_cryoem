function ...
[ ...
 parameter ...
,a_qbp_k_Y_lm_ ...
,a_evi_k_Y_lm_ ...
,a_var_k_Y_lm_ ...
,a_num_k_Y_lm_ ...
,a_0_k_Y_lm_ ...
,a_1_k_Y_lm_ ...
,a_2_k_Y_lm_ ...
,perturbation_a_qbp_k_Y_lm_ ...
,perturbation_a_evi_k_Y_lm_ ...
,perturbation_a_var_k_Y_lm_ ...
,perturbation_a_num_k_Y_lm_ ...
,perturbation_a_0_k_Y_lm_ ...
,perturbation_a_1_k_Y_lm_ ...
,perturbation_a_2_k_Y_lm_ ...
] = ...
qbp_single_shell_9( ...
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
,perturbation_euler_polar_a_ ...
,perturbation_euler_azimu_b_ ...
,perturbation_euler_gamma_z_ ...
);
%%%%%%%%;
% fix get_ylm__ ; %<-- done. ;
% fix knnsearch__ ; %<-- done (box3d does not work well). ;
% allow precomputation ;
% test using nonuniform mu(tau) ;
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
% perturbation_euler_polar_a_: real array of size n_M. perturbation to polar_a used for each image ;
% perturbation_euler_azimu_b_: real array of size n_M. perturbation to azimu_b used for each image ;
% perturbation_euler_gamma_z_: real array of size n_M. perturbation to gamma_z used for each image ;
% ;
% Output: ;
% a_qbp_k_Y_lm_: complex array of size n_lm. output function for back-propagated-volume in k_Y_ format. ;
% a_evi_k_Y_lm_: complex array of size n_lm. output function for back-propagated-evidence in k_Y_ format. ;
% a_var_k_Y_lm_: complex array of size n_lm. output function for back-propagated-variance in k_Y_ format. ;
% a_num_k_Y_lm_: complex array of size n_lm. output function for back-propagated-rawcount in k_Y_ format. ;
% a_0_k_Y_lm_: complex array of size n_lm. output function for back-propagated 0th-moment in k_Y_ format. ;
% a_1_k_Y_lm_: complex array of size n_lm. output function for back-propagated 1st-moment in k_Y_ format. ;
% a_2_k_Y_lm_: complex array of size n_lm. output function for back-propagated 2nd-moment in k_Y_ format. ;
% perturbation_a_qbp_k_Y_lm_: complex array of size n_lm. output function for perturbation to back-propagated-volume in k_Y_ format. ;
% perturbation_a_evi_k_Y_lm_: complex array of size n_lm. output function for perturbation to back-propagated-evidence in k_Y_ format. ;
% perturbation_a_var_k_Y_lm_: complex array of size n_lm. output function for perturbation to back-propagated-variance in k_Y_ format. ;
% perturbation_a_num_k_Y_lm_: complex array of size n_lm. output function for perturbation to back-propagated-rawcount in k_Y_ format. ;
% perturbation_a_0_k_Y_lm_: complex array of size n_lm. output function for perturbation to back-propagated 0th-moment in k_Y_ format. ;
% perturbation_a_1_k_Y_lm_: complex array of size n_lm. output function for perturbation to back-propagated 1st-moment in k_Y_ format. ;
% perturbation_a_2_k_Y_lm_: complex array of size n_lm. output function for perturbation to back-propagated 2nd-moment in k_Y_ format. ;
%%%%%%%%;

str_thisfunction = 'qbp_single_shell_9';

%%%%%%%%;
if (nargin<1);
%%%%%%%%;
tmp12;
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
if (nargin<1+na); perturbation_euler_polar_a_=[]; end; na=na+1;
if (nargin<1+na); perturbation_euler_azimu_b_=[]; end; na=na+1;
if (nargin<1+na); perturbation_euler_gamma_z_=[]; end; na=na+1;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
flag_verbose = parameter.flag_verbose;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
tolerance_master = parameter.tolerance_master;
if (~isfield(parameter,'qbp_deconvolution_factor_max')); parameter.qbp_deconvolution_factor_max = 8; end; %<-- parameter_bookmark. ;
qbp_deconvolution_factor_max = parameter.qbp_deconvolution_factor_max;
if (~isfield(parameter,'flag_qbp_Ylm_1_vs_0')); parameter.flag_qbp_Ylm_1_vs_0 = 1; end; %<-- parameter_bookmark. ;
flag_qbp_Ylm_1_vs_0 = parameter.flag_qbp_Ylm_1_vs_0;
if (~isfield(parameter,'flag_qbp_knn_vs_box3d')); parameter.flag_qbp_knn_vs_box3d = 1; end; %<-- parameter_bookmark. ;
flag_qbp_knn_vs_box3d = parameter.flag_qbp_knn_vs_box3d;
%%%%%%%%;

if isempty(euler_polar_a_); euler_polar_a_ = zeros(n_M,1); end;
if isempty(euler_azimu_b_); euler_azimu_b_ = zeros(n_M,1); end;
if isempty(euler_gamma_z_); euler_gamma_z_ = zeros(n_M,1); end;
if isempty(index_nCTF_from_nM_); index_nCTF_from_nM_ = zeros(n_M,1); end;
if isempty(CTF_k_p_wC__); CTF_k_p_wC__ = ones(n_w,1); end;
if isempty(perturbation_euler_polar_a_); perturbation_euler_polar_a_ = zeros(n_M,1); end;
if isempty(perturbation_euler_azimu_b_); perturbation_euler_azimu_b_ = zeros(n_M,1); end;
if isempty(perturbation_euler_gamma_z_); perturbation_euler_gamma_z_ = zeros(n_M,1); end;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

flag_d0 = 1;
flag_d1 = (nargin>=1+1*8); %<-- first-derivative requested. ;

%%%%%%%%;
n_lm = (1+l_max).^2; %<-- total number of spherical-harmonic coefficients. ;
%%%%%%%%;
% Rotate each of the CTF-functions. ;
%%%%%%%%;
CTF_k_p_wM__ = zeros(n_w,n_M); %<-- CTF for each image. ;
for nM=0:n_M-1;
euler_gamma_z = euler_gamma_z_(1+nM);
nCTF = index_nCTF_from_nM_(1+nM);
CTF_k_p_w_ = CTF_k_p_wC__(:,1+nCTF);
CTF_k_p_wM__(:,1+nM) = rotate_p_to_p_fftw(1,n_w,n_w,CTF_k_p_w_,+euler_gamma_z);
end;%for nM=0:n_M-1;
%%%%%%%%;
% set resolution of quadrature grid. ;
% Mollification involves surface-diffusion: 1/(2*pi)/sigma_use^2 .* exp(-d2/(2*sigma_use^2)). ;
% This corresponds to a diffusion-time of t_diffuse:=(sigma_use^2)/2. ;
% eventual deconvolution involves: exp(+(Y_l_val_).*(1+Y_l_val_).*t_diffuse). ;
% qbp_deconvolution_factor_max = exp(+(l_max)*(1+l_max)*t_diffuse);
% t_diffuse_max = log(qbp_deconvolution_factor_max)/max(1,l_max*(1+l_max));
% sigma_max = sqrt(log(qbp_deconvolution_factor_max)*2/max(1,(l_max*(1+l_max))));
% Note that the first neglected nearest-neighbor will be roughly: ;
% d2_max = (sqrt(3)/2 * n_ring*quad_k_eq_d)^2 away. ;
% This should corresond to a (relatively) negligible mollified quantity: ;
% epsilon_check = exp(-d2_max/(2*sigma_use^2)) = exp(-(n_ring*quad_k_eq_d)^2/(2*sigma_use^2)) < tolerance_master. ;
% n_ring = ceil(sqrt(max(1,-log(tolerance_master)*(2*sigma_use^2)/quad_k_eq_d^2)));
%%%%%%%%;
qbp_deconvolution_factor_max = 8;
t_diffuse_max = log(qbp_deconvolution_factor_max)/max(1,l_max*(1+l_max));
sigma_max = sqrt(2*t_diffuse_max);
%quad_k_eq_d = sqrt(4*pi./max(1,n_lm)); %<-- older setting from qbp_6. consider changing. ;
quad_k_eq_d = 0.5*sqrt(4*pi./max(1,n_lm)); %<-- one half the older setting from qbp_6. increased density of quadrature points. ;
sigma_use = min(sigma_max,quad_k_eq_d); %<-- ensure maximum amplification. ;
t_diffuse = (sigma_use^2)/2;
n_ring = ceil(sqrt(max(1,-log(tolerance_master)*(2*sigma_use^2)/quad_k_eq_d^2))); %<-- number of nearest neighbors requested for each point. ;
d1_max = sqrt(3)/2 * n_ring*quad_k_eq_d ;
d2_max = d1_max^2;
epsilon_check = exp(-d2_max/(2*sigma_use^2));
assert(epsilon_check<=tolerance_master);
n_nearest = 1+6*n_ring*(n_ring+1)/2; %<-- rough number of neighbors on hexagonal grid (i.e., 1+6+12+18+...). ;
if (flag_verbose>0); disp(sprintf(' %% tolerance_master %0.6f: d1_max %0.6f d2_max %0.6f quad_k_eq_d %0.6f n_nearest_k %d',tolerance_master,d1_max,d2_max,quad_k_eq_d,n_nearest)); end;
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
,quad_n_polar_a ...
,quad_polar_a_ ...
,quad_n_azimu_b_ ...
] = ...
sample_shell_5( ...
 1.0 ...
,quad_k_eq_d ...
,'L' ...
) ;
quad_k_c_qd__ = [ quad_k_c_0_all_ , quad_k_c_1_all_ , quad_k_c_2_all_ ];
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% sample_shell_5 (should be a precomputation): %0.2fs',tmp_t)); end;
if (flag_verbose>0); disp(sprintf(' %% quad_n_all: %.d',quad_n_all)); end;
%%%%;
tmp_t = tic();
if flag_qbp_Ylm_1_vs_0==0;
Ylm__ = get_Ylm__(1+l_max,0:l_max,quad_n_all,quad_azimu_b_all_,quad_polar_a_all_); %<-- potentially unstable for l_max>=88. ;
end;%if flag_qbp_Ylm_1_vs_0==0;
if flag_qbp_Ylm_1_vs_0==1;
Ylm__ = get_Ylm__1(1+l_max,0:l_max,quad_n_all,quad_azimu_b_all_,quad_polar_a_all_); %<-- more stable. ;
end;%if flag_qbp_Ylm_1_vs_0==1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% get_Ylm__1 (should be a precomputation): %0.2fs',tmp_t)); end;
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
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% Ylm_weight_yq__ (should be a precomputation): %0.2fs',tmp_t)); end;
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
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% cg_rhs_1: %0.2fs',tmp_t)); end;
%%%%;
if flag_qbp_knn_vs_box3d==0;
if (flag_verbose>0); disp(sprintf(' %% quad_n_all %d n_w*n_M %d',quad_n_all,n_w*n_M)); end;
if (flag_verbose>-1); disp(sprintf(' %% Warning, box3d_S2_search_0 option not implemented, returning early')); return; end;
tmp_t = tic();
d_req = d1_max;
[ ...
 index_quad_from_data_wMn__ ...
,n_index_quad_from_data_wM_ ...
,d2_quad_from_data_wMn__ ...
] = ...
box3d_S2_search_0( ...
 quad_n_all ...
,quad_k_c_qd__ ...
,n_w ...
,n_M ...
,data_k_c_wMd__ ...
,d_req ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% box3d_S2_search_0 (should be a precomputation): %0.2fs',tmp_t)); end;
end;%if flag_qbp_knn_vs_box3d==0;
%%%%;
if flag_qbp_knn_vs_box3d==1;
tmp_t = tic();
[ij_quad_from_data_wMn__,d1_quad_from_data_wMn__] = knnsearch(quad_k_c_qd__,data_k_c_wMd__,'K',n_nearest);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% knnsearch (technically unnecessary): %0.2fs',tmp_t)); end;
tmp_t = tic();
index_quad_from_data_wMn__ = ij_quad_from_data_wMn__ - 1;
d2_quad_from_data_wMn__ = d1_quad_from_data_wMn__.^2;
end;%if flag_qbp_knn_vs_box3d==1;
%%%%;
index_data_wMn__ = repmat(transpose(0:n_w*n_M-1),[1,n_nearest]);
mollify_quad_from_data_wMn__ = 1/(2*pi)/sigma_use.^2 .* exp(-d2_quad_from_data_wMn__/(2*sigma_use^2));
if flag_qbp_knn_vs_box3d==1; tmp_index_ = 0:n_w*n_M*n_nearest-1; end;
if flag_qbp_knn_vs_box3d==0; tmp_index_ = efind(index_quad_from_data_wMn__>=0); end;
quad_from_data_qwM__ = sparse(1+index_quad_from_data_wMn__(1+tmp_index_),1+index_data_wMn__(1+tmp_index_),mollify_quad_from_data_wMn__(1+tmp_index_),quad_n_all,n_w*n_M);
tmp_sml = nnz(quad_from_data_qwM__);
tmp_big = prod(size(quad_from_data_qwM__));
tmp_big_G = tmp_big*8/1e9;
tmp_sml_G = tmp_sml*8/1e9;
tmp_emp_G = memory_GB('quad_from_data_qwM__');
if (flag_verbose>0); disp(sprintf(' %% quad_from_data_qwM__ (%d,%d) <-- sparsity %d/%d = %0.6f <-- (%0.2fG,%0.2fG,%0.2fG) ',size(quad_from_data_qwM__),tmp_sml,tmp_big,tmp_sml/max(1,tmp_big),tmp_sml_G,tmp_emp_G,tmp_big_G)); end;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% quad_from_data_qwM__: %0.2fs',tmp_t)); end;
%%%%;
tmp_t = tic();
num_k_p_q_ = quad_from_data_qwM__*ones(n_w*n_M,1); %<-- effective number of data-points per quadrature-point. ;
C_C_k_p_q_ = quad_from_data_qwM__*reshape(abs(CTF_k_p_wM__).^2,[n_w*n_M,1]); %<-- accumulated CTF-squared per quadrature-point. ;
M_C_k_p_q_ = quad_from_data_qwM__*reshape(M_k_p_wM__.*CTF_k_p_wM__,[n_w*n_M,1]); %<-- accumulated M-times-CTF per quadrature-point. ;
M_M_k_p_q_ = quad_from_data_qwM__*reshape(abs(M_k_p_wM__).^2,[n_w*n_M,1]); %<-- accumulated M-times-M per quadrature-point. ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% apply quad_from_data_qwM__: %0.2fs',tmp_t)); end;
%%%%;
a_qbp_k_Y_lm_ = zeros(n_lm,1);
a_evi_k_Y_lm_ = zeros(n_lm,1);
a_var_k_Y_lm_ = zeros(n_lm,1);
a_num_k_Y_lm_ = zeros(n_lm,1);
a_0_k_Y_lm_ = zeros(n_lm,1);
a_1_k_Y_lm_ = zeros(n_lm,1);
a_2_k_Y_lm_ = zeros(n_lm,1);
%M_C_k_p_wM__ = M_k_p_wM__.*CTF_k_p_wM__;
%quad_from_data_M_C_k_p_normalized_q_ = (quad_from_data_qwM__ * M_C_k_p_wM__(:))./max(1e-12,C_C_k_p_q_);
%a_qbp_k_Y_lm_ = conj(Ylm_weight_yq__)*quad_from_data_M_C_k_p_normalized_q_;
%a_qbp_k_Y_lm_ = a_qbp_k_Y_lm_.*deconvolve_lm_; %<-- deconvolve (i.e., unmollify). ;
tmp_t = tic();
a_qbp_k_Y_lm_ = (conj(Ylm_weight_yq__)*( M_C_k_p_q_./max(1e-12,C_C_k_p_q_) )).*deconvolve_lm_;
a_evi_k_Y_lm_ = (conj(Ylm_weight_yq__)*( C_C_k_p_q_ )).*deconvolve_lm_;
a_var_k_Y_lm_ = (conj(Ylm_weight_yq__)*( M_M_k_p_q_./max(1e-12,C_C_k_p_q_) - abs(M_C_k_p_q_./max(1e-12,C_C_k_p_q_)).^2 )).*deconvolve_lm_;
a_num_k_Y_lm_ = (conj(Ylm_weight_yq__)*( num_k_p_q_ )).*deconvolve_lm_;
a_0_k_Y_lm_ = (conj(Ylm_weight_yq__)*( C_C_k_p_q_ )).*deconvolve_lm_;
a_1_k_Y_lm_ = (conj(Ylm_weight_yq__)*( M_C_k_p_q_ )).*deconvolve_lm_;
a_2_k_Y_lm_ = (conj(Ylm_weight_yq__)*( M_M_k_p_q_ )).*deconvolve_lm_;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% apply conj(Ylm_weight_yq__) and deconvolve_lm_: %0.2fs',tmp_t)); end;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
