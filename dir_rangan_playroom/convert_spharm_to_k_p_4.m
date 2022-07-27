function ...
[ ...
 a_k_all_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_spharm_to_k_p_4( ...
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
,a_k_Y_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
% uses spherical-harmonic-expansion a_k_Y_ to evaluate a_k_p_ on a collection of points on spherical shells determined by k_p_r_. ;
% We assume that the polar-representation and quadrature weights associated with these points have been previously calculated. ; 
% ;
% inputs: ;
% ;
% verbose = integer verbosity_level. ;
% n_k_all = integer total number of points. ;
% n_k_all_csum_ = integer array of starting indices associated with each k-value. ;
% k_p_r_all_ = real array of k-values for each point. ;
% k_p_azimu_b_all_ = real array of azimu_b-values for each point. ;
% k_p_polar_a_all_ = real array of polar_a-values for each point. ;
% weight_3d_k_all_ = real array of quadrature weights for volume integral (for each point) (unused). ;
% weight_shell_k_ = real array of quadrature weights for shell integral (for each point). ;
% n_k_p_r = integer maximum k. ;
% k_p_r_ = real array of length n_k_p_r; k_p_r_(nk_p_r) = k_p_r_value for shell nk_p_r. ;
% weight_3d_k_p_r_ = real array of length n_k_p_r; radial quadrature weights (already assumed to be a factor of weight_3d_k_all_) (unused). ;
% l_max_ = integer array of length n_k_p_r; l_max_(nk_p_r) = spherical harmonic order on shell nk_p_r; l_max_(nk_p_r) corresponds to n_lm_(nk_p_r) = (l_max_(nk_p_r)+1)^2 coefficients. ;
% a_k_Y_ = complex array of length \sum_{nk_p_r} (n_lm_(nk_p_r)+1)^2 ; coefficients are ordered in a row, with m varying quickly, l varying slowly and k varying most slowly. ;
% Ylm_klma___ = complex cell array of length n_k_p_r. ; Ylm_klma___{1+nk_p_r} holds the array Ylm_lma__. ;
% Ylm_lma__ = complex array of size (n_lm,n_sub), where n_lm record (l,m) pairs (m varying quickly), and n_sub records the ordered list of points on the shell at nk_p_r. ;
% See get_Ylm_condense_wrap_0 for additional outputs. ;
% ;
% outputs: ;
% ;
% a_k_all_ = complex array of a-values for each point. ;
% Ylm_klma___ = complex cell array of length n_k_p_r. ; Ylm_klma___{1+nk_p_r} holds the array Ylm_lma__. ;
% See get_Ylm_condense_wrap_0 for additional inputs. ;

str_thisfunction = 'convert_spharm_to_k_p_4';

na=0;
if (nargin<1+na); verbose=[]; end; na=na+1;
if (nargin<1+na); n_k_all=[]; end; na=na+1;
if (nargin<1+na); n_k_all_csum_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_all_=[]; end; na=na+1;
if (nargin<1+na); k_p_azimu_b_all_=[]; end; na=na+1;
if (nargin<1+na); k_p_polar_a_all_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_all_=[]; end; na=na+1;
if (nargin<1+na); weight_shell_k_=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_=[]; end; na=na+1;
if (nargin<1+na); Ylm_uklma___=[]; end; na=na+1;

flag_Ylm_create = 0; 
if isempty(Ylm_uklma___);
flag_Ylm_create=1;
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
end;%if isempty(Ylm_uklma___);

n_lm_ = (l_max_+1).^2;
if (verbose); disp(sprintf(' %% [entering %s] n_k_all %d, n_lm_sum %d',str_thisfunction,n_k_all,sum(n_lm_))); end;
a_k_all_ = zeros(n_k_all,1);
ix0=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_k_all_csum = n_k_all_csum_(1+nk_p_r);
if (verbose>1); disp(sprintf(' %% nk_p_r %d/%d k_p_r %0.2f, n_k_all_csum %d --> %0.2f%%',nk_p_r,n_k_p_r,k_p_r,n_k_all_csum,n_k_all_csum/n_k_all)); end;
if (nk_p_r<n_k_p_r-1); n_sub = n_k_all_csum_(1+nk_p_r+1) - n_k_all_csum_(1+nk_p_r); else n_sub = n_k_all - n_k_all_csum ; end;
index_sub_ = n_k_all_csum + (0:n_sub-1);
k_p_r_sub_ = k_p_r_all_(1+index_sub_); 
assert(sum(k_p_r_sub_==k_p_r)==n_sub);
k_p_azimu_b_sub_ = k_p_azimu_b_all_(1+index_sub_);
k_p_polar_a_sub_ = k_p_polar_a_all_(1+index_sub_);
weight_3d_k_sub_ = weight_shell_k_(1+index_sub_)/k_p_r^2;
l_max = l_max_(1+nk_p_r);
n_l_max = l_max+1; flag_flip=0;
n_lm = n_lm_(1+nk_p_r);
%%%%;
nu_n_k_per_shell = index_nu_n_k_per_shell_from_nk_p_r_(1+nk_p_r);
tmp_k_p_azimu_b_sub_ = k_p_azimu_b_sub_uka__{1+nu_n_k_per_shell};
tmp_k_p_polar_a_sub_ = k_p_polar_a_sub_uka__{1+nu_n_k_per_shell};
if (fnorm(k_p_azimu_b_sub_-tmp_k_p_azimu_b_sub_)>1e-6);
disp(sprintf(' %% Warning, fnorm(k_p_azimu_b_sub_-tmp_k_p_azimu_b_sub_) %0.16f in %s',fnorm(k_p_azimu_b_sub_-tmp_k_p_azimu_b_sub_),str_thisfunction));
end;%if (fnorm(k_p_azimu_b_sub_-tmp_k_p_azimu_b_sub_)>1e-6);
if (fnorm(k_p_polar_a_sub_-tmp_k_p_polar_a_sub_)>1e-6);
disp(sprintf(' %% Warning, fnorm(k_p_polar_a_sub_-tmp_k_p_polar_a_sub_) %0.16f in %s',fnorm(k_p_polar_a_sub_-tmp_k_p_polar_a_sub_),str_thisfunction));
end;%if (fnorm(k_p_polar_a_sub_-tmp_k_p_polar_a_sub_)>1e-6);
Ylm_lma__ = Ylm_uklma___{1+nu_n_k_per_shell}(1:n_lm,:);
%%%%;
tmp_t = tic();
a_k_all_sub_ = zeros(1,n_sub);
a_k_all_sub_ = reshape(a_k_Y_(1+ix0 + [0:n_lm-1]),[1,n_lm])*Ylm_lma__;
ix0 = ix0 + n_lm;
a_k_all_(1+index_sub_) = a_k_all_sub_;
clear Ylm_lma__;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% nk_p_r %d/%d a_k_all_sub_: %0.6fs',nk_p_r,n_k_p_r,tmp_t)); end;
end;%for nk_p_r=0:n_k_p_r-1;
if (verbose); disp(sprintf(' %% [finished %s] n_k_all %d, n_lm_sum %d',str_thisfunction,n_k_all,sum(n_lm_))); end;

%{
  
%%%%%%%%;
% First run test_pm_trpv1c_9c.m ;
%%%%%%%%;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi);
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
,'L' ...
,1 ...
) ; 
%%%%%%%%;
tmp_t = tic();
[ ...
 a_k_p_reco1_  ...
] =  ...
convert_spharm_to_k_p_1( ...
 1 ...
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
,a_k_Y_quad_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% convert_spharm_to_k_p_1: %0.6fs',tmp_t));
%%%%%%%%;
[ ...
 Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
get_Ylm_condense_wrap_0( ...
 1 ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,n_k_p_r ...
,l_max_ ...
);
%%%%%%%%;
tmp_t = tic();
[ ...
 a_k_p_reco2_  ...
] =  ...
convert_spharm_to_k_p_4( ...
 1 ...
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
,a_k_Y_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% convert_spharm_to_k_p_4: %0.6fs',tmp_t));
%%%%%%%%;
disp(sprintf(' %% a_k_p_reco1_ vs a_k_p_reco2_: %0.16f',fnorm(a_k_p_reco1_ - a_k_p_reco2_)/fnorm(a_k_p_reco1_)));
%%%%%%%%;

%%%%%%%%;
tmp_t = tic();
[ ...
 a_k_Y_reco1_  ...
] =  ...
convert_k_p_to_spharm_1( ...
 1 ...
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
,a_k_p_reco1_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% convert_k_p_to_spharm_1: %0.6fs',tmp_t));
%%%%%%%%;
tmp_t = tic();
[ ...
 a_k_Y_reco2_  ...
] =  ...
convert_k_p_to_spharm_4( ...
 1 ...
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
,a_k_p_reco2_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% convert_k_p_to_spharm_4: %0.6fs',tmp_t));
%%%%%%%%;
disp(sprintf(' %% a_k_Y_reco1_ vs a_k_Y_reco2_: %0.16f',fnorm(a_k_Y_reco1_ - a_k_Y_reco2_)/fnorm(a_k_Y_reco1_)));
%%%%%%%%;


 %}
