function ...
[ ...
 parameter ...
,ssnll ...
] = ...
ssnll_from_S_k_p_single_shell_9( ...
 parameter ...
,l_max ...
,n_w ...
,n_S ...
,S_k_p_wS__ ...
,euler_polar_a_S_ ...
,euler_azimu_b_S_ ...
,n_M ...
,M_k_p_wM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wC__ ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
);
%%%%%%%%;
% Calculates ssnll (sigma*sigma*log_unlikelihood). ;
% Applies to a single shell of a_k_Y_. ;
% Associates CTF_k_p_wC__(:,1+index_nCTF_from_nM_(1+nM)) with image M_k_p_wM__(:,1+nM);
% ;
% Input: ;
% l_max: integer order used for a_k_Y_ on shell. ;
% n_w: integer number of inplane_gamma_z values recorded at that ring. ;
% n_S: integer number of templates. ;
% S_k_p_wS__: complex array of size (n_w,n_S). stack of templates on ring in k_p_ format. ;
% euler_polar_a_S_: real array of size n_S. polar_a used for each template ;
% euler_azimu_b_S_: real array of size n_S. azimu_b used for each template ;
% n_M: integer number of images. ;
% M_k_p_wM__: complex array of size (n_w,n_M). stack of images on ring in k_p_ format. ;
% index_nCTF_from_nM_: integer array of size n_M. index_nCTF_from_nM_(1+nM) is the (base 0) CTF_index used for image M_k_p_wM__(:,1+nM). ;
%             This can be empty or set to 1, in which case the same CTF_k_p_ will be used for each image. ;
% CTF_k_p_wC__: complex array of size(n_w,n_CTF). stack of ctf-functions in k_p_ format. ;
%            If index_nCTF_from_nM_ is empty or set to 1, then we assume this contains only a single CTF_k_p_, ;
%            which will then be used for all images. ;
% euler_polar_a_M_: real array of size n_M. polar_a used for each image ;
% euler_azimu_b_M_: real array of size n_M. azimu_b used for each image ;
% euler_gamma_z_M_: real array of size n_M. gamma_z used for each image ;
% ;
% Output: ;
% ssnll: double (sigma*sigma*log_unlikelihood). ;
%%%%%%%%;

str_thisfunction = 'ssnll_from_S_k_p_single_shell_9';

%%%%%%%%;
if (nargin<1);
%%%%%%%%;
disp(sprintf(' %% returning')); return;
%%%%%%%%;
end;%if (nargin<1);
%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); l_max=[]; end; na=na+1;
if (nargin<1+na); n_w=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); S_k_p_wS__=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_S_=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_S_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wM__=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_wC__=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_M_=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_M_=[]; end; na=na+1;
if (nargin<1+na); euler_gamma_z_M_=[]; end; na=na+1;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
flag_verbose = parameter.flag_verbose;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
tolerance_master = parameter.tolerance_master;
%%%%%%%%;

if isempty(euler_polar_a_M_); euler_polar_a_M_ = zeros(n_M,1); end;
if isempty(euler_azimu_b_M_); euler_azimu_b_M_ = zeros(n_M,1); end;
if isempty(euler_gamma_z_M_); euler_gamma_z_M_ = zeros(n_M,1); end;
if isempty(index_nCTF_from_nM_); index_nCTF_from_nM_ = zeros(n_M,1); end;
if isempty(CTF_k_p_wC__); CTF_k_p_wC__ = ones(n_w,1); end;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
n_lm = (1+l_max).^2; %<-- total number of spherical-harmonic coefficients. ;
%%%%%%%%;
% Not rotate each of the CTF-functions. ;
%%%%%%%%;
tmp_t = tic();
CTF_k_p_wM__ = zeros(n_w,n_M); %<-- CTF for each image. ;
for nM=0:n_M-1;
euler_gamma_z = euler_gamma_z_M_(1+nM);
nCTF = index_nCTF_from_nM_(1+nM);
CTF_k_p_w_ = CTF_k_p_wC__(:,1+nCTF);
CTF_k_p_wM__(:,1+nM) = rotate_p_to_p_fftw(1,n_w,n_w,CTF_k_p_w_,+0*euler_gamma_z); %<-- Not rotate the CTF-functions. ;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% CTF_k_p_wM__: %0.2fs',tmp_t)); end;
%%%%%%%%;
% Yes rotate each of the images. ;
%%%%%%%%;
tmp_t = tic();
N_k_p_wM__ = zeros(n_w,n_M); %<-- aligned images. ;
for nM=0:n_M-1;
euler_gamma_z = euler_gamma_z_M_(1+nM);
M_k_p_w_ = M_k_p_wM__(:,1+nM);
N_k_p_w_ = rotate_p_to_p_fftw(1,n_w,n_w,M_k_p_w_,-1*euler_gamma_z); %<-- Yes rotate the images. ;
N_k_p_wM__(:,1+nM) = N_k_p_w_;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% N_k_p_wM__: %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
index_nS_from_nM_ = zeros(n_M,1);
for nM=0:n_M-1;
euler_polar_a = euler_polar_a_M_(1+nM);
euler_azimu_b = euler_azimu_b_M_(1+nM);
[tmp_val,tmp_ij] = min(abs(euler_polar_a_S_ - euler_polar_a) + abs(euler_azimu_b_S_ - euler_azimu_b)); %<-- assuming a perfect match exists. ;
if (tmp_val>1e-12); disp(sprintf(' %% Warning, imperfect match in %s',str_thisfunction)); end;
nS = tmp_ij-1;
index_nS_from_nM_(1+nM) = nS;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% index_nS_from_nM_: %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
ssnll = 0.5d0 * sum(abs(S_k_p_wS__(:,1+index_nS_from_nM_).*CTF_k_p_wM__ - N_k_p_wM__).^2,'all') / max(1,n_w);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% ssnll: %0.2fs',tmp_t)); end;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
