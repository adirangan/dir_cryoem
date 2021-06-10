function ...
a_k_Y_ = ...
cg_lsq_3( ...
 n_order ...
,n_k_p_r ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p__ ...
,CTF_index_ ...
,CTF_k_p__ ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
);
%%%%%%%%;
% simply applies conjugate-gradient to the least-squares problem to solve for a_k_Y_. ;
% Associates CTF_k_p__(:,1+CTF_index_(1+nM)) with image M_k_p__(:,1+nM);
% ;
% Input: ;
% n_order: integer polynomial interpolation order used for conjugate-gradient interpolation operator. ;
% n_k_p_r: integer number of shells. ;
% l_max_: integer array of size n_k_p_r. l_max_(1+nk_p_r) is the order used for a_k_Y_ on shell nk_p_r. ;
% n_w_: integer array of size n_k_p_r. n_w_(1+nk_p_r) is the number of inplane_gamma_z values recorded at that ring. ;
% n_M: integer number of images. ;
% M_k_p__: complex array of size (n_w_sum,n_M). stack of images in k_p_ format. ;
% CTF_index_: integer array of size n_M. CTF_index_(1+nM) is the (base 0) CTF_index used for image M_k_p__(:,1+nM). ;
%             This can be empty or set to 1, in which case the same CTF_k_p_ will be used for each image. ;
% CTF_k_p__: complex array of size(n_w_sum,n_CTF). stack of ctf-functions in k_p_ format. ;
%            If CTF_index_ is empty or set to 1, then we assume this contains only a single CTF_k_p_, ;
%            which will then be used for all images. ;
% euler_polar_a_: real array of size n_M. polar_a used for each image ;
% euler_azimu_b_: real array of size n_M. azimu_b used for each image ;
% euler_gamma_z_: real array of size n_M. gamma_z used for each image ;
% ;
% Output: ;
% a_k_Y_: complex array of size n_lm_sum. output functions in k_Y_ format. ;
%%%%%%%%;

n_lm_ = (1+l_max_).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);

n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);

if isempty(CTF_k_p__); CTF_index_ = 0; CTF_k_p__ = ones(n_w_sum,1); end;

a_k_Y_ = zeros(n_lm_sum,1);
l_max_pre = 0;
n_w_pre = 0;
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
tmp_Y_ij_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
n_w = n_w_(1+nk_p_r);
tmp_M_ij_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
if (isempty(CTF_index_) | numel(CTF_index_)==1);
tmp_CTF_k_p_ = reshape(repmat(CTF_k_p__(1+tmp_M_ij_),[1,n_M]),[n_w*n_M,1]);
 else;
tmp_CTF_k_p_ = reshape(CTF_k_p__(1+tmp_M_ij_,1+CTF_index_(1:n_M)),[n_w*n_M,1]);
end;%if (isempty(CTF_index_) | numel(CTF_index_)==1);
%%%%%%%%;
% Reuse previous operators if appropriate. ;
%%%%%%%%;
if ( (n_w~=n_w_pre) | (l_max~=l_max_pre) );
[k_p_polar_a__,k_p_azimu_b__] = cg_rhs_1(n_M,n_w,euler_polar_a_,euler_azimu_b_,+euler_gamma_z_);
n_polar_a = ceil(n_w/2);
n_azimu_b = max(1+2*l_max,2*n_polar_a);
[legendre_evaluate_ljm___,legendre_evaluate_mlj___,expil__,expi__] = legendre_evaluate_ljm___0(l_max,cos(linspace(0,pi,n_polar_a)),n_azimu_b);
tensor_to_scatter__ = cg_interpolate_n_1(n_order,n_polar_a,n_azimu_b,n_w*n_M,k_p_polar_a__(:),k_p_azimu_b__(:));
scatter_to_tensor__ = ctranspose(tensor_to_scatter__); %<-- this conjugation is not necessary, since the matrix should be real. ;
end;%if ( (n_w~=n_w_pre) | (l_max~=l_max_pre) );
%%%%%%%%;
An__ = @(a_k_Y_) tmp_CTF_k_p_.*(tensor_to_scatter__*reshape(cg_evaluate_n_1(l_max,convert_spharm_to_spharm__0(l_max,a_k_Y_),n_polar_a,n_azimu_b,legendre_evaluate_ljm___),[n_polar_a*n_azimu_b,1]));
At__ = @(a_k_X_) convert_spharm__to_spharm_0(l_max,cg_evaluate_t_1(n_polar_a,n_azimu_b,reshape(scatter_to_tensor__*(conj(tmp_CTF_k_p_).*a_k_X_),[n_polar_a,n_azimu_b]),l_max,legendre_evaluate_mlj___,expil__,expi__));
AtAn__ = @(a_k_Y_) At__(An__(a_k_Y_));
[a_k_Y_(1+tmp_Y_ij_),~] = pcg(AtAn__,At__(reshape(M_k_p__(1+tmp_M_ij_,:),[n_w*n_M,1])));
l_max_pre = l_max;
n_w_pre = n_w;
end;%for nk_p_r=0:n_k_p_r-1;

