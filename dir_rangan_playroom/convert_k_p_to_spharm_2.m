function ...
[a_k_Y_] = ...
convert_k_p_to_spharm_2( ...
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
,a_k_all_ ...
);
% reconstitutes spherical-harmonic-expansion a_k_Y_ from a collection of points on spherical shells determined by k_p_r_. ;
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
% a_k_all_ = complex array of a-values for each point. ;
% ;
% outputs: ;
% ;
% a_k_Y_ = complex array of length \sum_{nk_p_r} (n_lm_(nk_p_r)+1)^2 ; coefficients are ordered in a row, with m varying quickly, l varying slowly and k varying most slowly. ;

if (nargin<1);
%%%%%%%%;
verbose=0; k_p_r_max = 15.0d0; eq_d = 0.125*k_p_r_max;
[n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,k_c_0_all_,k_c_1_all_,k_c_2_all_] = sample_sphere_6(verbose,k_p_r_max,eq_d,'L') ;
disp(sprintf(' %% Setting shell number nk_p_r to have spherical harmonic (nk_p_r%%3,nk_p_r%%5-2), provided that k_p_r is sufficiently large.'));
%%%%%%%%;
a_k_all_ = zeros(n_k_all,1);
disp(sprintf(' %% Setting l_max_ = 1+ceil(2*pi*k_p_r_/k_p_r_max). Note that this should be limited by the sampling density at each k_p_r!'));
l_max_ = 1+ceil(2*pi*k_p_r_/k_p_r_max); n_lm_ = (l_max_+1).^2; n_lm_sum = sum(n_lm_); a_k_Y_ = zeros(n_lm_sum,1);
ix=0;
for nk_p_r=0:n_k_p_r-1;
disp(sprintf(' %% nk_p_r %d/%d; l_max_max %d',nk_p_r,n_k_p_r,l_max_(1+nk_p_r)));
k_p_r = k_p_r_(1+nk_p_r);
n_k_all_csum = n_k_all_csum_(1+nk_p_r);
if (nk_p_r<n_k_p_r-1); n_sub = n_k_all_csum_(1+nk_p_r+1) - n_k_all_csum_(1+nk_p_r); else n_sub = n_k_all - n_k_all_csum ; end;
index_sub_ = n_k_all_csum + (0:n_sub-1);
k_p_r_sub_ = k_p_r_all_(1+index_sub_); 
assert(sum(k_p_r_sub_==k_p_r)==n_sub);
k_p_azimu_b_sub_ = k_p_azimu_b_all_(1+index_sub_);
k_p_polar_a_sub_ = k_p_polar_a_all_(1+index_sub_);
weight_k_sub_ = weight_shell_k_(1+index_sub_)/k_p_r^2;
l_max_max = l_max_(1+nk_p_r);
n_l_max = l_max_max+1; flag_flip=0;
[Ylm_sub__] = get_Ylm__(n_l_max,0:n_l_max,n_sub,k_p_azimu_b_sub_,k_p_polar_a_sub_,flag_flip);
l_max_set = mod(nk_p_r,3)-0;
m_val_set = mod(nk_p_r,5)-2;
if (l_max_set>l_max_max); l_max_set = l_max_max; end;
if (abs(m_val_set)>l_max_max); m_val_set = sign(m_val_set)*l_max_max; end;
a_k_sub_ = zeros(n_sub,1);
if (abs(m_val_set)<=l_max_set); a_k_sub_ = Ylm_sub__{1+l_max_set}(1+l_max_set + m_val_set,:); end;
for l_val=0:l_max_max;
n_m = 2*l_val+1;
for nm=0:n_m-1;
m_val = nm - l_val;
a_k_Y_(1+ix) = (l_val==l_max_set) & (m_val==m_val_set);
ix = ix+1;
end;%for nm=0:n_m-1;
end;%for l_val=0:l_max_max;
a_k_all_(1+index_sub_) = a_k_sub_;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
[b_k_all_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,a_k_Y_);
disp(sprintf(' %% a_k_all_ vs b_k_all_ error %0.16f',fnorm(a_k_all_-b_k_all_)/fnorm(a_k_all_)));
%%%%%%%%;
[b_k_Y_] = convert_k_p_to_spharm_2(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,a_k_all_);
disp(sprintf(' %% a_k_Y_ vs b_k_Y_ error %0.16f',fnorm(a_k_Y_-b_k_Y_)/fnorm(a_k_Y_)));
%%%%%%%%;
ix=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
l_max_max = l_max_(1+nk_p_r);
n_l_max = l_max_max+1; flag_flip=0;
l_max_set = mod(nk_p_r,3)-0;
m_val_set = mod(nk_p_r,5)-2;
if (l_max_set>l_max_max); l_max_set = l_max_max; end;
if (abs(m_val_set)>l_max_max); m_val_set = sign(m_val_set)*l_max_max; end;
for l_val=0:l_max_max;
n_m = 2*l_val+1;
for nm=0:n_m-1;
m_val = nm - l_val;
if ((a_k_Y_(1+ix)~=0) & (a_k_Y_(1+ix)~=b_k_Y_(1+ix))); 
disp(sprintf(' %% nk_p_r %d; l_max_max %d; k_p_r %0.2f; (l,m) = (%d,%d); (l_set,m_set) = (%d,%d); a_k_Y_ %0.16f b_k_Y_ %0.16f',nk_p_r,l_max_max,k_p_r,l_val,m_val,l_max_set,m_val_set,a_k_Y_(1+ix),b_k_Y_(1+ix))); 
end;%if ((a_k_Y_(1+ix)~=0) & (a_k_Y_(1+ix)~=b_k_Y_(1+ix))); 
ix = ix+1;
end;%for nm=0:n_m-1;
end;%for l_val=0:l_max_max;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
verbose=1;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi); TorL = 'L';
%%%%%%%%;
tmp_t = tic();
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
,TorL ...
) ;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% sample_sphere_7: %0.2fs',tmp_t)); end;
if (verbose); disp(sprintf(' %% k_p_r_max %0.2f k_eq_d %0.2f n_k_all %d n_k_p_r %d',k_p_r_max,k_eq_d,n_k_all,n_k_p_r)); end;
%%%%%%%%;
l_max_upb = 36;
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
l_max_(1+nk_p_r) = max(0,min(l_max_upb,1+ceil(2*pi*k_p_r_(1+nk_p_r))));
end;%for nk_p_r=0:n_k_p_r-1;
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
%%%%%%%%;
delta_orig_ = [+0.12;-0.3;+0.23];
a_k_p_orig_ = exp(+2*pi*i*(k_c_0_all_*delta_orig_(1+0) + k_c_1_all_*delta_orig_(1+1) + k_c_2_all_*delta_orig_(1+2)));
tmp_t = tic;
[a_k_Y_quad_] = ...
convert_k_p_to_spharm_2( ...
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
,a_k_p_orig_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
delta_orig_r012 = sqrt(delta_orig_(1+0)^2 + delta_orig_(1+1)^2 + delta_orig_(1+2)^2);
delta_orig_r01 = sqrt(delta_orig_(1+0)^2 + delta_orig_(1+1)^2);
delta_orig_polar_a = atan2(delta_orig_r01,delta_orig_(1+2));
delta_orig_azimu_b = atan2(delta_orig_(1+1),delta_orig_(1+0));
delta_Ylm_ = get_Ylm__(1+l_max_max,0:l_max_max,1,delta_orig_azimu_b,delta_orig_polar_a);
a_k_Y_form_ = zeros(n_lm_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
l_max = l_max_(1+nk_p_r);
for l_val=0:l_max;
tmp_x = 2*pi*k_p_r*delta_orig_r012;
tmp_jl = besselj(l_val+0.5,tmp_x)*sqrt(pi/(2*tmp_x));
for m_val=-l_val:+l_val;
a_k_Y_form_(1+na) = 4*pi * i^l_val * tmp_jl * conj(delta_Ylm_{1+l_val}(1+m_val+l_val));
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
disp(sprintf(' %% a_k_Y_form_ vs a_k_Y_quad_: %0.16f',fnorm(a_k_Y_form_ - a_k_Y_quad_)/fnorm(a_k_Y_form_)));
%%%%%%%%;
disp('returning'); return;
end;%if (nargin<1);

str_thisfunction = 'convert_k_p_to_spharm_2';

n_lm_ = (l_max_+1).^2;
if (verbose); disp(sprintf(' %% [entering %s] n_k_all %d, n_lm_sum %d',str_thisfunction,n_k_all,sum(n_lm_))); end;
a_k_Y_ = zeros(sum(n_lm_),1);
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
tmp_t = tic();
[Ylm_sub__] = get_Ylm__(n_l_max,0:n_l_max,n_sub,k_p_azimu_b_sub_,k_p_polar_a_sub_,flag_flip);
disp(size(Ylm_sub__)),;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% nk_p_r %d/%d Ylm_sub__: %0.6fs',nk_p_r,n_k_p_r,tmp_t)); end;
tmp_t = tic();
a_k_sub_ = a_k_all_(1+index_sub_);
for l_val=0:l_max;
n_m = 2*l_val+1;
for nm=0:n_m-1;
m_val = nm - l_val;
a_k_Y_(1+ix) = dot(Ylm_sub__{1+l_val}(1+nm,:),a_k_sub_(:).*weight_k_sub_(:));
ix = ix+1;
end;%for nm=0:n_m-1;
end;%for l_val=0:l_max;
clear k_p_r_sub_ k_p_azimu_b_sub_ k_p_polar_a_sub_ Ylm_sub__ a_k_sub_ ;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% nk_p_r %d/%d a_k_sub_: %0.6fs',nk_p_r,n_k_p_r,tmp_t)); end;
end;%for nk_p_r=0:n_k_p_r-1;
if (verbose); disp(sprintf(' %% [finished %s] n_k_all %d, n_lm_sum %d',str_thisfunction,n_k_all,sum(n_lm_))); end;




