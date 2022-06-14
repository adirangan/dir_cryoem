function ...
a_k_Y_ ...
= ...
cg_lsq_6( ...
 n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
,image_delta_x_ ...
,image_delta_y_ ...
,image_I_value_...
);
%%%%%%%%;
% simply applies conjugate-gradient to the least-squares problem to solve for a_k_Y_. ;
% Associates CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+nM)) with image M_k_p_wkM__(:,1+nM);
% ;
% Input: ;
% n_order: integer polynomial interpolation order used for conjugate-gradient interpolation operator. ;
% n_k_p_r: integer number of shells. ;
% k_p_r_: real array of size n_k_p_r. k-values for each shell. ;
% l_max_: integer array of size n_k_p_r. l_max_(1+nk_p_r) is the order used for a_k_Y_ on shell nk_p_r. ;
% n_w_: integer array of size n_k_p_r. n_w_(1+nk_p_r) is the number of inplane_gamma_z values recorded at that ring. ;
% n_M: integer number of images. ;
% M_k_p_wkM__: complex array of size (n_w_sum,n_M). stack of images in k_p_ format. ;
% index_nCTF_from_nM_: integer array of size n_M. index_nCTF_from_nM_(1+nM) is the (base 0) CTF_index used for image M_k_p_wkM__(:,1+nM). ;
%             This can be empty or set to 1, in which case the same CTF_k_p_ will be used for each image. ;
% CTF_k_p_wkC__: complex array of size(n_w_sum,n_CTF). stack of ctf-functions in k_p_ format. ;
%            If index_nCTF_from_nM_ is empty or set to 1, then we assume this contains only a single CTF_k_p_, ;
%            which will then be used for all images. ;
% euler_polar_a_: real array of size n_M. polar_a used for each image ;
% euler_azimu_b_: real array of size n_M. azimu_b used for each image ;
% euler_gamma_z_: real array of size n_M. gamma_z used for each image ;
% image_delta_x_: real array of size n_M. delta_x used for each image ;
% image_delta_y_: real array of size n_M. delta_y used for each image ;
% image_I_value_: real array of size n_M. I_value used for each image ;
% ;
% Output: ;
% a_k_Y_: complex array of size n_lm_sum. output functions in k_Y_ format. ;
%%%%%%%%;

na=0;
if (nargin<1+na); n_order=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wkM__=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_wkC__=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_=[]; end; na=na+1;
if (nargin<1+na); euler_gamma_z_=[]; end; na=na+1;
if (nargin<1+na); image_delta_x_=[]; end; na=na+1;
if (nargin<1+na); image_delta_y_=[]; end; na=na+1;
if (nargin<1+na); image_I_value_=[]; end; na=na+1;

if (isempty(image_delta_x_)); image_delta_x_ = zeros(n_M,1); end;
if (isempty(image_delta_y_)); image_delta_y_ = zeros(n_M,1); end;
if (isempty(image_I_value_)); image_I_value_ = ones(n_M,1); end;

n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
% construct CTF of same size as images. ;
%%%%%%%%
if (size(CTF_k_p_wkC__,1)==n_k_p_r);
n_CTF = size(CTF_k_p_wkC__,2);
CTF_k_p_r_kC__ = CTF_k_p_wkC__;
CTF_k_p_wkC__ = zeros(n_w_sum,n_CTF);
for nCTF=0:n_CTF-1;
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
tmp_index_ = n_w_csum_(1+nk_p_r) + [0:n_w-1];
CTF_k_p_wkC__(1+tmp_index_,1+nCTF) = CTF_k_p_r_kC__(1+nk_p_r,1+nCTF);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nCTF=0:n_CTF-1;
end;%if (size(CTF_k_p_wkC__,1)==n_k_p_r);

T_k_p_wkM__ = M_k_p_wkM__;
for nM=0:n_M-1;
T_k_p_wkM__(:,1+nM) = image_I_value_(1+nM) * transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wkM__(:,1+nM),+image_delta_x_(1+nM),+image_delta_y_(1+nM));
end;%for nM=0:n_M-1;

a_k_Y_ = cg_lsq_3(n_order,n_k_p_r,l_max_,n_w_,n_M,T_k_p_wkM__,index_nCTF_from_nM_,CTF_k_p_wkC__,euler_polar_a_,euler_azimu_b_,euler_gamma_z_);
