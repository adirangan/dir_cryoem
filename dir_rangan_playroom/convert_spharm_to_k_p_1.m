function [a_k_all_] = ...
convert_spharm_to_k_p_1( ...
 verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_k_p_r_ ...
,l_max_ ...
,a_k_Y_ ...
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
% weight_k_all_ = real array of quadrature weights for volume integral (for each point) (unused). ;
% weight_shell_k_ = real array of quadrature weights for shell integral (for each point). ;
% n_k_p_r = integer maximum k. ;
% k_p_r_ = real array of length n_k_p_r; k_p_r_(nk_p_r) = k_p_r_value for shell nk_p_r. ;
% weight_k_p_r_ = real array of length n_k_p_r; radial quadrature weights (already assumed to be a factor of weight_k_all_) (unused). ;
% l_max_ = integer array of length n_k_p_r; l_max_(nk_p_r) = spherical harmonic order on shell nk_p_r; l_max_(nk_p_r) corresponds to n_lm_(nk_p_r) = (l_max_(nk_p_r)+1)^2 coefficients. ;
% a_k_Y_ = complex array of length \sum_{nk_p_r} (n_lm_(nk_p_r)+1)^2 ; coefficients are ordered in a row, with m varying quickly, l varying slowly and k varying most slowly. ;
% ;
% outputs: ;
% ;
% a_k_all_ = complex array of a-values for each point. ;

n_lm_ = (l_max_+1).^2;
if (verbose); disp(sprintf(' %% [entering convert_spharm_to_k_p_1] n_k_all %d, n_lm_sum %d',n_k_all,sum(n_lm_))); end;
a_k_all_ = zeros(n_k_all,1);
ix=0;
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
weight_k_sub_ = weight_shell_k_(1+index_sub_)/k_p_r^2;
l_max = l_max_(1+nk_p_r);
n_l_max = l_max+1; flag_flip=0;
[Ylm_sub__] = get_Ylm__(n_l_max,0:n_l_max,n_sub,k_p_azimu_b_sub_,k_p_polar_a_sub_,flag_flip);
a_k_sub_ = zeros(1,n_sub);
for l_val=0:l_max;
n_m = 2*l_val+1;
for nm=0:n_m-1;
m_val = nm - l_val;
a_k_sub_ = a_k_sub_ + a_k_Y_(1+ix)*Ylm_sub__{1+l_val}(1+nm,:);
ix = ix+1;
end;%for nm=0:n_m-1;
end;%for l_val=0:l_max;
a_k_all_(1+index_sub_) = a_k_sub_;
end;%for nk_p_r=0:n_k_p_r-1;
if (verbose); disp(sprintf(' %% [finished convert_spharm_to_k_p_1] n_k_all %d, n_lm_sum %d',n_k_all,sum(n_lm_))); end;



