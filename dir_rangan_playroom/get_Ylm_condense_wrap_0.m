function ...
[ ...
 Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
get_Ylm_condense_wrap_0( ...
 flag_verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,n_k_p_r ...
,l_max_ ...
);
% Gets array of Ylm_klma___ ... ;
% Condenses adjacent k if the evaluations are the same. ;
% Assumes that the evaluations k_p_azimu_b_ and k_p_polar_a_ on each shell are a function of n_k_per_shell. ;
% ;
% inputs: ;
% ;
% flag_verbose = integer verbosity_level. ;
% n_k_all = integer total number of points. ;
% n_k_all_csum_ = integer array of starting indices associated with each k-value. ;
% k_p_azimu_b_all_ = real array of azimu_b-values for each point. ;
% k_p_polar_a_all_ = real array of polar_a-values for each point. ;
% n_k_p_r = integer maximum k. ;
% l_max_ = integer array of length n_k_p_r; l_max_(nk_p_r) = spherical harmonic order on shell nk_p_r; l_max_(nk_p_r) corresponds to n_lm_(nk_p_r) = (l_max_(nk_p_r)+1)^2 coefficients. ;
% ;
% outputs: ;
% ;
% Ylm_klma___ = complex cell array of length n_k_p_r. ; Ylm_klma___{1+nk_p_r} holds the array Ylm_lma__. ;
% Ylm_lma__ = complex array of size (n_lm,n_sub), where n_lm record (l,m) pairs (m varying quickly), and n_sub records the ordered list of points on the shell at nk_p_r. ;
% Ylm_uklma___ = complex cell array of length n_u_n_k_per_shell. ; unique sets of Ylm_lma__, one for each nu_n_k_per_k_shell. Stores l_val up to l_max_uk_(1+nu_k_per_k_shell). ;
% index_k_per_shell_uka__ = integer cell array of length n_u_n_k_per_shell. ; index sets (of nk_p_r) associated with each nu_n_k_per_k_shell. ;
% index_nu_n_k_per_shell_from_nk_p_r_ = integer array of length n_k_p_r. ; cross reference. ;
% l_max_uk_ = integer array of length n_u_n_k_per_shell. ; largest l_max across the shells with nk_p_r in index_k_per_shell_uka__{1+nu_k_per_k_shell}. ;
% k_p_azimu_b_sub_uka__ = double array of length n_u_n_k_per_shell. ; k_p_azimu_b_ for the shells with nk_p_r in index_k_per_shell_uka__{1+nu_k_per_k_shell}. ;
% k_p_polar_a_sub_uka__ = double array of length n_u_n_k_per_shell. ; k_p_polar_a_ for the shells with nk_p_r in index_k_per_shell_uka__{1+nu_k_per_k_shell}. ;

str_thisfunction = 'get_Ylm_condense_wrap_0';

if nargin<1;
flag_verbose=1;
if (flag_verbose>0); disp(sprintf(' %% testing %s',str_thisfunction)); end;
for k_int = [ 8];%for k_int = [ 8,16,32,48,64];
for k_eq_d_double = [1.0];%for k_eq_d_double = [1.0,0.5];
%%%%;
k_p_r_max = k_int/(2*pi); k_eq_d = k_eq_d_double/(2*pi); str_L = 'L';
flag_uniform_over_n_k_p_r = 1;
flag_uniform_over_polar_a = 0; %<-- This is set to match test_ssnll_from_a_k_p_14 ;
[ ...
 n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
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
,str_L ...
,flag_uniform_over_n_k_p_r ...
,flag_uniform_over_polar_a ...
) ;
%%%%;
l_max_upb = round(2*pi*k_p_r_max); %<-- typically sufficient for 2-3 digits of precision. ;
l_max_ = zeros(n_k_p_r,1); for nk_p_r=0:n_k_p_r-1; l_max_(1+nk_p_r) = max(0,min(l_max_upb,1+ceil(2*pi*k_p_r_(1+nk_p_r)))); end;%for nk_p_r=0:n_k_p_r-1;
[ ...
 Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
get_Ylm_condense_wrap_0( ...
 0*flag_verbose ...
,n_qk ...
,n_qk_csum_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,n_k_p_r ...
,l_max_ ...
);
n_u = numel(Ylm_uklma___);
if (flag_verbose>0);
n_uklma=0; for nu=0:n_u-1; n_uklma = n_uklma + numel(Ylm_uklma___{1+nu}); end;%for nu=0:n_u-1;
disp(sprintf(' %% k_int %.2d --> n_k_p_r %.3d l_max_upb %.3d n_u %.2d n_uklma %0.6d --> Ylm_uklma___: %0.6fGB',k_int,n_k_p_r,l_max_upb,n_u,n_uklma,(n_uklma*16/1e9)));
end;%if (flag_verbose>0);
%%%%;
% compare to ../dir_rangan_python/Ylm_0lma__.ascii. ;
%%%%;
dir_ascii = '../dir_rangan_python/dir_ascii';
n_d = 5;
for nd=0:n_d-1;
na=0;
if nd==na; tmp_str_infix = 'Ylm_0lma__'; tmp_array = Ylm_uklma___{1}; end; na=na+1;
if nd==na; tmp_str_infix = 'k_p_azimu_b_sub_0a_'; tmp_array = k_p_azimu_b_sub_uka__{1}; end; na=na+1;
if nd==na; tmp_str_infix = 'k_p_polar_a_sub_0a_'; tmp_array = k_p_polar_a_sub_uka__{1}; end; na=na+1;
if nd==na; tmp_str_infix = 'index_nu_n_k_per_shell_from_nk_p_r_'; tmp_array = index_nu_n_k_per_shell_from_nk_p_r_; end; na=na+1;
if nd==na; tmp_str_infix = 'index_k_per_shell_0a_'; tmp_array = index_k_per_shell_uka__{1}; end; na=na+1;
tmp_error=0;
fname_ascii = sprintf('%s/%s.ascii',dir_ascii,tmp_str_infix);
if ~exist(fname_ascii,'file'); disp(sprintf(' %% Warning, %s not found',fname_ascii)); end;
if  exist(fname_ascii,'file');
if (flag_verbose>1); disp(sprintf(' %% %s found, loading',fname_ascii)); end;
fid = fopen(fname_ascii,'r'); tmp_load_ascii = textscan(fid,'(%f)'); fclose(fid);
if isempty(tmp_load_ascii{1});
fid = fopen(fname_ascii,'r'); tmp_load_ascii = textscan(fid,'%f'); fclose(fid);
end;%if isempty(tmp_load_ascii{1});
tmp_array_ascii = reshape(cell2mat(tmp_load_ascii),size(tmp_array));
tmp_error = tmp_error + fnorm_disp(max(0,flag_verbose-1),tmp_str_infix,tmp_array,'ascii',tmp_array_ascii,' %% <-- should be <1e-6');
end;%if  exist(fname_ascii,'file');
if (flag_verbose>0); disp(sprintf(' %% nd %d/%d: %s: tmp_error: %0.16f; %%<-- should be <1e-6 per entry',nd,n_d,tmp_str_infix,tmp_error)); end;
end;%for nd=0:n_d-1;
%%%%;
end;%for k_int = [ 8,16,32,48,64];
end;%for k_eq_d_double = [1.0,0.5];
%%%%;
disp('returning');return;
end;%if nargin<1;

na=0;
if (nargin<1+na); flag_verbose=[]; end; na=na+1;
if (nargin<1+na); n_k_all=[]; end; na=na+1;
if (nargin<1+na); n_k_all_csum_=[]; end; na=na+1;
if (nargin<1+na); k_p_azimu_b_all_=[]; end; na=na+1;
if (nargin<1+na); k_p_polar_a_all_=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;

n_lm_k_ = (l_max_+1).^2;
n_k_per_shell_k_ = zeros(n_k_p_r,1);
index_sub_ka__ = cell(n_k_p_r,1);
k_p_azimu_b_sub_ka__ = cell(n_k_p_r,1);
k_p_polar_a_sub_ka__ = cell(n_k_p_r,1);

for nk_p_r=0:n_k_p_r-1;
n_k_all_csum = n_k_all_csum_(1+nk_p_r);
if (flag_verbose>1); disp(sprintf(' %% nk_p_r %d/%d, n_k_all_csum %d --> %0.2f%%',nk_p_r,n_k_p_r,n_k_all_csum,n_k_all_csum/n_k_all)); end;
if (nk_p_r<n_k_p_r-1); n_k_per_shell = n_k_all_csum_(1+nk_p_r+1) - n_k_all_csum_(1+nk_p_r); else n_k_per_shell = n_k_all - n_k_all_csum ; end;
index_sub_ = n_k_all_csum + (0:n_k_per_shell-1);
k_p_azimu_b_sub_ = k_p_azimu_b_all_(1+index_sub_);
k_p_polar_a_sub_ = k_p_polar_a_all_(1+index_sub_);
n_k_per_shell_k_(1+nk_p_r) = n_k_per_shell;
index_sub_ka__{1+nk_p_r} = index_sub_;
k_p_azimu_b_sub_ka__{1+nk_p_r} = k_p_azimu_b_sub_;
k_p_polar_a_sub_ka__{1+nk_p_r} = k_p_polar_a_sub_;
end;%for nk_p_r=0:n_k_p_r-1;

[u_n_k_per_shell_,~,ij_nu_n_k_per_shell_from_nk_p_r_] = unique(n_k_per_shell_k_);
n_u_n_k_per_shell = numel(u_n_k_per_shell_);
index_nu_n_k_per_shell_from_nk_p_r_ = ij_nu_n_k_per_shell_from_nk_p_r_-1;
index_k_per_shell_uka__ = cell(n_u_n_k_per_shell,1);
l_max_uk_ = zeros(n_u_n_k_per_shell,1);
for nu_n_k_per_shell=0:n_u_n_k_per_shell-1;
u_n_k_per_shell = u_n_k_per_shell_(1+nu_n_k_per_shell);
index_k_per_shell_uka__{1+nu_n_k_per_shell} = efind(n_k_per_shell_k_==u_n_k_per_shell);
l_max_uk_(1+nu_n_k_per_shell) = max(l_max_(1+index_k_per_shell_uka__{1+nu_n_k_per_shell}));
end;%for nu_n_k_per_shell=0:n_u_n_k_per_shell-1;

n_lm_uk_ = (l_max_uk_+1).^2;
n_k_per_shell_uk_ = zeros(n_u_n_k_per_shell,1);
k_p_azimu_b_sub_uka__ = cell(n_u_n_k_per_shell,1);
k_p_polar_a_sub_uka__ = cell(n_u_n_k_per_shell,1);
for nu_n_k_per_shell=0:n_u_n_k_per_shell-1;
u_n_k_per_shell = u_n_k_per_shell_(1+nu_n_k_per_shell);
index_k_per_shell_uka__{1+nu_n_k_per_shell} = efind(n_k_per_shell_k_==u_n_k_per_shell);
nk_p_r_0 = index_k_per_shell_uka__{1+nu_n_k_per_shell}(1+0);
k_p_azimu_b_sub_uka__{1+nu_n_k_per_shell} = k_p_azimu_b_sub_ka__{1+nk_p_r_0};
k_p_polar_a_sub_uka__{1+nu_n_k_per_shell} = k_p_polar_a_sub_ka__{1+nk_p_r_0};
end;%for nu_n_k_per_shell=0:n_u_n_k_per_shell-1;

if (flag_verbose); disp(sprintf(' %% [entering %s] n_k_all %d, n_lm_sum %d, n_u_n_k_per_shell %d',str_thisfunction,n_k_all,sum(n_lm_k_),n_u_n_k_per_shell)); end;

Ylm_uklma___ = cell(n_u_n_k_per_shell,1);
for nu_n_k_per_shell=0:n_u_n_k_per_shell-1;
k_p_azimu_b_sub_ = k_p_azimu_b_sub_uka__{1+nu_n_k_per_shell};
k_p_polar_a_sub_ = k_p_polar_a_sub_uka__{1+nu_n_k_per_shell};
n_k_per_shell = numel(k_p_polar_a_sub_);
l_max = l_max_uk_(1+nu_n_k_per_shell);
n_l_max = l_max+1; flag_flip=0;
n_lm = (l_max+1).^2;
%%%%;
if (flag_verbose); disp(sprintf(' %% nu_n_k_per_shell %d/%d: Ylm_lma__ %.2dGB',nu_n_k_per_shell,n_u_n_k_per_shell,n_lm*n_k_per_shell*8/1e9)); end;
tmp_t = tic();
[Ylm_sub__] = get_Ylm__(n_l_max,0:l_max,n_k_per_shell,k_p_azimu_b_sub_,k_p_polar_a_sub_,flag_flip);
Ylm_lma__ = zeros(n_lm,n_k_per_shell);
ix1=0;
for l_val=0:l_max;
n_m = 2*l_val+1;
for nm=0:n_m-1;
m_val = nm - l_val;
Ylm_lma__(1+ix1,:) = Ylm_sub__{1+l_val}(1+nm,:);
ix1 = ix1+1;
end;%for nm=0:n_m-1;
end;%for l_val=0:l_max;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% nu_n_k_per_shell %d/%d Ylm_lma__: %0.6fs',nu_n_k_per_shell,n_u_n_k_per_shell,tmp_t)); end;
Ylm_uklma___{1+nu_n_k_per_shell} = Ylm_lma__;
clear Ylm_sub__;
end;%for nu_n_k_per_shell=0:n_u_n_k_per_shell-1;

if (flag_verbose); disp(sprintf(' %% [finished %s] n_k_all %d, n_lm_sum %d, n_u_n_k_per_shell %d',str_thisfunction,n_k_all,sum(n_lm_k_),n_u_n_k_per_shell)); end;
