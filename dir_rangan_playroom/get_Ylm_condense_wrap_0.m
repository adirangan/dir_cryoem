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
 verbose ...
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
% verbose = integer verbosity_level. ;
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

na=0;
if (nargin<1+na); verbose=[]; end; na=na+1;
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
if (verbose>1); disp(sprintf(' %% nk_p_r %d/%d k_p_r %0.2f, n_k_all_csum %d --> %0.2f%%',nk_p_r,n_k_p_r,k_p_r,n_k_all_csum,n_k_all_csum/n_k_all)); end;
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

if (verbose); disp(sprintf(' %% [entering %s] n_k_all %d, n_lm_sum %d, n_u_n_k_per_shell %d',str_thisfunction,n_k_all,sum(n_lm_k_),n_u_n_k_per_shell)); end;

Ylm_uklma___ = cell(n_u_n_k_per_shell,1);
for nu_n_k_per_shell=0:n_u_n_k_per_shell-1;
k_p_azimu_b_sub_ = k_p_azimu_b_sub_uka__{1+nu_n_k_per_shell};
k_p_polar_a_sub_ = k_p_polar_a_sub_uka__{1+nu_n_k_per_shell};
n_k_per_shell = numel(k_p_polar_a_sub_);
l_max = l_max_uk_(1+nu_n_k_per_shell);
n_l_max = l_max+1; flag_flip=0;
n_lm = (l_max+1).^2;
%%%%;
if (verbose); disp(sprintf(' %% nu_n_k_per_shell %d/%d: Ylm_lma__ %.2dGB',nu_n_k_per_shell,n_u_n_k_per_shell,n_lm*n_k_per_shell*8/1e9)); end;
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
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% nu_n_k_per_shell %d/%d Ylm_lma__: %0.6fs',nu_n_k_per_shell,n_u_n_k_per_shell,tmp_t)); end;
Ylm_uklma___{1+nu_n_k_per_shell} = Ylm_lma__;
clear Ylm_sub__;
end;%for nu_n_k_per_shell=0:n_u_n_k_per_shell-1;

if (verbose); disp(sprintf(' %% [finished %s] n_k_all %d, n_lm_sum %d, n_u_n_k_per_shell %d',str_thisfunction,n_k_all,sum(n_lm_k_),n_u_n_k_per_shell)); end;
