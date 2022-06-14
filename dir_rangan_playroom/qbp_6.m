function ...
[ ...
 a_k_Y_ ...
,n_quad_from_data_q_ ...
] = ...
qbp_6( ...
 qbp_eps ...
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
% Applies quadrature-back-propagation to solve for a_k_Y_. ;
% Associates CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+nM)) with image M_k_p_wkM__(:,1+nM);
% ;
% Input: ;
% qbp_eps: pseudoinverse threshold used for solving local least-squares for ctf. ;
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
if (nargin<1+na); qbp_eps=[]; end; na=na+1;
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

if isempty(qbp_eps); qbp_eps = 1e-3; end;
if isempty(image_delta_x_); image_delta_x_ = zeros(n_M,1); end;
if isempty(image_delta_y_); image_delta_y_ = zeros(n_M,1); end;
if isempty(image_I_value_); image_I_value_ = ones(n_M,1); end;
if isempty(index_nCTF_from_nM_); index_nCTF_from_nM_ = zeros(n_M,1); end;
if isempty(CTF_k_p_wkC__); CTF_k_p_wkC__ = ones(n_w_sum,1); end;
if (qbp_eps>1); qbp_eps = max(1e-12,0.1^(qbp_eps)); end;

verbose=0;
if (verbose>0); disp(sprintf(' %% [entering qbp_6]')); end;

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
if ( (fnorm(image_delta_x_)>0) | (fnorm(image_delta_y_)>0) | (fnorm(image_I_value_-ones(n_M,1))>0) );
for nM=0:n_M-1;
T_k_p_wkM__(:,1+nM) = image_I_value_(1+nM) * transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wkM__(:,1+nM),+image_delta_x_(1+nM),+image_delta_y_(1+nM));
end;%for nM=0:n_M-1;
end;%if ( (fnorm(image_delta_x_)>0) | (fnorm(image_delta_y_)>0) | (fnorm(image_I_value_-ones(n_M,1))>0) );

l_max_ = l_max_(1:n_k_p_r);
n_lm_ = (1+l_max_).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
n_w_max = max(n_w_(1:n_k_p_r));
n_w_sum = sum(n_w_(1:n_k_p_r));
n_w_csum_ = cumsum([0;n_w_(1:n_k_p_r)]);
if (size(CTF_k_p_wkC__,1)==n_k_p_r);
%%%%%%%%;
% construct CTF of same size as images. ;
%%%%%%%%
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

quad_k_eq_d_ = sqrt(4*pi./n_lm_);
if (verbose); for nk_p_r=0:n_k_p_r-1; disp(sprintf(' %% nk_p_r %d/%d: quad_k_eq_d %0.6f',nk_p_r,n_k_p_r,quad_k_eq_d_(1+nk_p_r))); end; end;

flag_unique_n = 0;
if (numel(unique(l_max_))==1 & numel(unique(n_lm_))==1 & numel(unique(n_w_))==1);
flag_unique_n = 1;
quad_k_eq_d = quad_k_eq_d_(1+0);
l_max = l_max_(1+0);
n_lm = n_lm_(1+0);
n_w = n_w_(1+0);
%%%%;
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
%%%%;
Ylm__ = get_Ylm__(1+l_max,0:l_max,quad_n_all,quad_azimu_b_all_,quad_polar_a_all_);
Ylm_yq__ = zeros(n_lm,quad_n_all);
nml=0;
for l_val=0:l_max;
for m_val=-l_val:+l_val;
Ylm_yq__(1+nml,:) = Ylm__{1+l_val}(1+l_val+m_val,:);
nml=nml+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
Ylm_w_yq__ = Ylm_yq__ * sparse(1:quad_n_all,1:quad_n_all,quad_weight_all_,quad_n_all,quad_n_all);
%%%%;
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
%%%%;
index_quad_from_data_ = knnsearch(quad_k_c_qd__,data_k_c_wMd__,'K',1); index_quad_from_data_ = index_quad_from_data_ - 1;
quad_from_data_qwM__ = sparse(1+index_quad_from_data_,1:n_w*n_M,1,quad_n_all,n_w*n_M);
n_quad_from_data_q_ = quad_from_data_qwM__*ones(n_w*n_M,1); %<-- number of data-points per quadrature-point. ;
data_from_quad_wMq__ = bsxfun(@rdivide,transpose(quad_from_data_qwM__),max(1,transpose(n_quad_from_data_q_)));
CTF_wMn__ = reshape(permute(reshape(CTF_k_p_wkC__(:,1+index_nCTF_from_nM_),[n_w,n_k_p_r,n_M]),[1,3,2]),[n_w*n_M,n_k_p_r]);
CTF2_qk__ = quad_from_data_qwM__*abs(CTF_wMn__).^2;
%%%%;
end;%if (numel(unique(l_max_))==1 & numel(unique(n_lm_))==1 & numel(unique(n_w_))==1);

a_k_Y_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
T_wM__ = T_k_p_wkM__(1+index_nw_,:);
%%%%%%%%;
if flag_unique_n==1;
CTF_wM_ = CTF_wMn__(:,1+nk_p_r);
CTF2_q_ = CTF2_qk__(:,1+nk_p_r);
end;%if flag_unique_n==1;
if flag_unique_n==0;
% redefine workspace. ;
%%%%;
quad_k_eq_d = quad_k_eq_d_(1+nk_p_r);
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
%%%%;
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
%%%%;
Ylm__ = get_Ylm__(1+l_max,0:l_max,quad_n_all,quad_azimu_b_all_,quad_polar_a_all_);
Ylm_yq__ = zeros(n_lm,quad_n_all);
nml=0;
for l_val=0:l_max;
for m_val=-l_val:+l_val;
Ylm_yq__(1+nml,:) = Ylm__{1+l_val}(1+l_val+m_val,:);
nml=nml+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
Ylm_w_yq__ = Ylm_yq__ * sparse(1:quad_n_all,1:quad_n_all,quad_weight_all_,quad_n_all,quad_n_all);
%%%%;
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
%%%%;
index_quad_from_data_ = knnsearch(quad_k_c_qd__,data_k_c_wMd__,'K',1); index_quad_from_data_ = index_quad_from_data_ - 1;
quad_from_data_qwM__ = sparse(1+index_quad_from_data_,1:n_w*n_M,1,quad_n_all,n_w*n_M);
n_quad_from_data_q_ = quad_from_data_qwM__*ones(n_w*n_M,1); %<-- number of data-points per quadrature-point. ;
data_from_quad_wMq__ = bsxfun(@rdivide,transpose(quad_from_data_qwM__),max(1,transpose(n_quad_from_data_q_)));
CTF_wM_ = reshape(CTF_k_p_wkC__(1+index_nw_,1+index_nCTF_from_nM_),[n_w*n_M,1]);
CTF2_q_ = quad_from_data_qwM__*abs(CTF_wM_).^2;
%%%%;
end;%if flag_unique_n==0;
%%%%%%%%;
T_CTF_wM__ = T_wM__.*reshape(CTF_wM_,[n_w,n_M]);
quad_from_data_T_CTF_normalized_q_ = (quad_from_data_qwM__ * T_CTF_wM__(:))./max(1e-12,CTF2_q_);
a_k_Y_(1+index_Y_) = conj(Ylm_w_yq__)*quad_from_data_T_CTF_normalized_q_;
%%%%%%%%;
end;%for nk_p_r=0:n_k_p_r-1;

if (verbose>0); disp(sprintf(' %% [finished qbp_6]')); end;
