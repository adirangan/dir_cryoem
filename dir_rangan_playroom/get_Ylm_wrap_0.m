function ...
[ ...
 Ylm_klma___ ...
] = ...
get_Ylm_wrap_0( ...
 verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,n_k_p_r ...
,l_max_ ...
);
% Gets array of Ylm_klma___ ...
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

str_thisfunction = 'get_Ylm_wrap_0';

na=0;
if (nargin<1+na); verbose=[]; end; na=na+1;
if (nargin<1+na); n_k_all=[]; end; na=na+1;
if (nargin<1+na); n_k_all_csum_=[]; end; na=na+1;
if (nargin<1+na); k_p_azimu_b_all_=[]; end; na=na+1;
if (nargin<1+na); k_p_polar_a_all_=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;

n_lm_ = (l_max_+1).^2;
if (verbose); disp(sprintf(' %% [entering %s] n_k_all %d, n_lm_sum %d',str_thisfunction,n_k_all,sum(n_lm_))); end;

Ylm_klma___ = cell(n_k_p_r,1);
ix0=0;
for nk_p_r=0:n_k_p_r-1;
n_k_all_csum = n_k_all_csum_(1+nk_p_r);
if (verbose>1); disp(sprintf(' %% nk_p_r %d/%d k_p_r %0.2f, n_k_all_csum %d --> %0.2f%%',nk_p_r,n_k_p_r,k_p_r,n_k_all_csum,n_k_all_csum/n_k_all)); end;
if (nk_p_r<n_k_p_r-1); n_sub = n_k_all_csum_(1+nk_p_r+1) - n_k_all_csum_(1+nk_p_r); else n_sub = n_k_all - n_k_all_csum ; end;
index_sub_ = n_k_all_csum + (0:n_sub-1);
k_p_azimu_b_sub_ = k_p_azimu_b_all_(1+index_sub_);
k_p_polar_a_sub_ = k_p_polar_a_all_(1+index_sub_);
l_max = l_max_(1+nk_p_r);
n_l_max = l_max+1; flag_flip=0;
n_lm = n_lm_(1+nk_p_r);
%%%%;
tmp_t = tic();
[Ylm_sub__] = get_Ylm__(n_l_max,0:n_l_max,n_sub,k_p_azimu_b_sub_,k_p_polar_a_sub_,flag_flip);
Ylm_lma__ = zeros(n_lm,n_sub);
ix1=0;
for l_val=0:l_max;
n_m = 2*l_val+1;
for nm=0:n_m-1;
m_val = nm - l_val;
Ylm_lma__(1+ix1,:) = Ylm_sub__{1+l_val}(1+nm,:);
ix1 = ix1+1;
end;%for nm=0:n_m-1;
end;%for l_val=0:l_max;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% nk_p_r %d/%d Ylm_lma__: %0.6fs',nk_p_r,n_k_p_r,tmp_t)); end;
Ylm_klma___{1+nk_p_r} = Ylm_lma__;
clear Ylm_sub__;
end;%for nk_p_r=0:n_k_p_r-1;

if (verbose); disp(sprintf(' %% [finished %s] n_k_all %d, n_lm_sum %d',str_thisfunction,n_k_all,sum(n_lm_))); end;
