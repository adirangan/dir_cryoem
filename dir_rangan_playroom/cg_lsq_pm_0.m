function ...
a_UCTF_UX_Y_ync__ = ...
cg_lsq_pm_0( ...
 n_order ...
,pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_l_max_ ...
,pm_n_w_ ...
,n_M ...
,UX_M_k_p_wnM__ ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
,image_delta_x_ ...
,image_delta_y_ ...
,image_I_value_...
);
%%%%%%%%;
% Applies conjugate-gradient to the least-squares problem to solve for a_UX_Y_. ;
% Associates CTF_k_p_r__(:,1+CTF_index_(1+nM)) with image M_k_p__(:,1+nM);
% ;
% Input: ;
% n_order: integer polynomial interpolation order used for conjugate-gradient interpolation operator. ;
% pm_n_UX_rank = pm_n_k_p_r: integer number of principal-volume shells retained. ;
% (unused) pm_k_p_r_: real array of size pm_n_k_p_r. k-values for each principal-volume shell (assumed to be all ones). ;
% pm_l_max_: integer array of size pm_n_k_p_r. pm_l_max_(1+pm_nk_p_r) is the order used for a_UCTF_UX_Y_ync__ on principal-volume-shell pm_nk_p_r. ;
% pm_n_w_: integer array of size pm_n_k_p_r. pm_n_w_(1+pm_nk_p_r) is the number of inplane_gamma_z values recorded at that principal-image-ring. ;
% n_M: integer number of images. ;
% UX_M_k_p_wnM__: complex array of size (pm_n_w_sum,n_M). stack of principal-images in k_p_ format. ;
% We assume that each column of UX_M_k_p_wnM__ corresponds to a single principal-image, which is itself a stack of principal-image-rings. ;
% n_CTF_rank: integer number of CTF ranks to consider. ;
% (implicit) UCTF_kc__: real array of size (n_k_p_r,n_CTF_rank). ;
% VSCTF_Mc__: real array of size (n_M,n_CTF_rank). ;
% We assume that the full (isotropic) CTF-function CTF_k_p_r__, given by: ;
% (implicit) CTF_k_p_r__: real array of size (n_k_p_r,n_M). ;
% can be approximated via: ;
% CTF_k_p_r__ = UCTF_kc__*transpose(VSCTF_Mc__);
% which is a low-rank approximation with n_CTF_rank terms. ;
% If n_CTF_rank<=0 or VSCTF_Mc__ is empty, we assume that: ;
% n_CTF_rank=1; UCTF_kc__ = ones(pm_n_k_p_r,1); VSCTF_Mc__ = ones(n_M,1);
% euler_polar_a_: real array of size n_M. polar_a used for each image ;
% euler_azimu_b_: real array of size n_M. azimu_b used for each image ;
% euler_gamma_z_: real array of size n_M. gamma_z used for each image ;
% (unused) image_delta_x_: real array of size n_M. delta_x used for each image ;
% (unused) image_delta_y_: real array of size n_M. delta_y used for each image ;
% image_I_value_: real array of size n_M. I_value used for each image ;
% ;
% Output: ;
% a_UCTF_UX_Y_ync__: complex array of size (pm_n_lm_sum,n_UCTF_rank). output functions in k_Y_ format. ;
% This output function satisfies the least-square problem: ;
% \sum_{nCTF_rank=0}^{n_CTF_rank-1} S * [ \tau_{1+nM} * VSCTF_Mc__(1+nM,1+nCTF_rank) ] * a_UCTF_UX_Y_ync___(:,1+pm_nUX_rank,1+nCTF_rank) = [ UX_M_k_p_wnM___(:,1+pm_nUX_rank,1+nM) ] \forall nM \in [0,\ldots,n_M-1] and \forall pm_nUX_rank \in [0,pm_n_UX_rank-1]. ;
% where : ;
% \tau_{1+nM} corresponds to rotation by the viewing-angle associated with image nM, and ;
% S is the template-operator (i.e., equatorial-evaluation), and ;
% UX_M_k_p_wnM___(:,1+pm_nUX_rank,1+nM) = UX_M_k_p_wnM__(1+pm_n_w_csum_(1+pm_nUX_rank) + (0:pm_n_w_(1+pm_nUX_rank)-1),1+nM), and ;
% a_UCTF_UX_Y_ync___(:,1+pm_nUX_rank,1+nCTF_rank) = \sum_{nk_p_r=0}^{n_k_p_r} UCTF_(1+nk_p_r,1+nCTF_rank) * UX_(1+nk_p_r,1+pm_nUX_rank) * a_UX_Y__(:,1+nk_p_r). ;
%%%%%%%%;

if (isempty(image_I_value_)); image_I_value_ = ones(n_M,1); end;

pm_n_UX_rank = pm_n_k_p_r;
pm_l_max_ = pm_l_max_(1:pm_n_UX_rank);
pm_n_lm_ = (1+pm_l_max_).^2;
pm_n_lm_max = max(pm_n_lm_);
pm_n_lm_sum = sum(pm_n_lm_);
pm_n_lm_csum_ = cumsum([0;pm_n_lm_]);

pm_n_w_max = max(pm_n_w_(1:pm_n_k_p_r));
pm_n_w_sum = sum(pm_n_w_(1:pm_n_k_p_r));
pm_n_w_csum_ = cumsum([0;pm_n_w_(1:pm_n_k_p_r)]);

if (n_CTF_rank<=0); n_CTF_rank = 1; VSCTF_Mc__ = ones(n_M,1); end;
if (isempty(VSCTF_Mc__)); n_CTF_rank = 1; VSCTF_Mc__ = ones(n_M,1); end;

a_UCTF_UX_Y_ync__ = zeros(pm_n_lm_sum,n_CTF_rank);
for pm_nk_p_r=0:pm_n_k_p_r-1;
pm_l_max = pm_l_max_(1+pm_nk_p_r);
pm_n_lm = pm_n_lm_(1+pm_nk_p_r);
index_Y_ = pm_n_lm_csum_(1+pm_nk_p_r) + (0:pm_n_lm-1);
pm_n_w = pm_n_w_(1+pm_nk_p_r);
index_M_ = pm_n_w_csum_(1+pm_nk_p_r) + (0:pm_n_w-1);
VSCTF_wMc__ = reshape(repmat(reshape(VSCTF_Mc__,[1,n_M,n_CTF_rank]),[pm_n_w,1,1]),[pm_n_w*n_M,n_CTF_rank]);
[k_p_polar_a__,k_p_azimu_b__] = cg_rhs_1(n_M,pm_n_w,euler_polar_a_,euler_azimu_b_,+euler_gamma_z_);
n_polar_a = ceil(pm_n_w/2);
n_azimu_b = max(1+2*pm_l_max,2*n_polar_a);
[legendre_evaluate_ljm___,legendre_evaluate_mlj___,expil__,expi__] = legendre_evaluate_ljm___0(pm_l_max,cos(linspace(0,pi,n_polar_a)),n_azimu_b);
tensor_to_scatter__ = cg_interpolate_n_1(n_order,n_polar_a,n_azimu_b,pm_n_w*n_M,k_p_polar_a__(:),k_p_azimu_b__(:));
scatter_to_tensor__ = ctranspose(tensor_to_scatter__); %<-- this conjugation is not necessary, since the matrix should be real. ;
An__ = @(a_UX_Y_yc_) An(n_CTF_rank,n_M,VSCTF_wMc__,pm_n_w,pm_l_max,a_UX_Y_yc_,tensor_to_scatter__,n_polar_a,n_azimu_b,legendre_evaluate_ljm___);
At__ = @(UX_M_) At(n_CTF_rank,n_M,VSCTF_wMc__,pm_n_w,pm_l_max,UX_M_,scatter_to_tensor__,n_polar_a,n_azimu_b,legendre_evaluate_mlj___,expil__,expi__);
AtAn__ = @(a_UX_Y_yc_) At__(An__(a_UX_Y_yc_));
[a_UX_Y_yc_,~] = pcg(AtAn__,At__(reshape(UX_M_k_p_wnM__(1+index_M_,:)*sparse(1:n_M,1:n_M,image_I_value_,n_M,n_M),[pm_n_w*n_M,1])));
for nCTF_rank=0:n_CTF_rank-1;
tmp_index_ = nCTF_rank*pm_n_lm + (0:pm_n_lm-1);
a_UCTF_UX_Y_ync__(1+index_Y_,1+nCTF_rank) = a_UX_Y_yc_(1+tmp_index_);
end;%for nCTF_rank=0:n_CTF_rank-1;
end;%for pm_nk_p_r=0:pm_n_k_p_r-1;

function An_a_UX_Y_yc_ = An(n_CTF_rank,n_M,VSCTF_wMc__,pm_n_w,pm_l_max,a_UX_Y_yc_,tensor_to_scatter__,n_polar_a,n_azimu_b,legendre_evaluate_ljm___);
pm_n_lm = (1+pm_l_max).^2;
An_a_UX_Y_yc_ = zeros(pm_n_w*n_M,1);
for nCTF_rank=0:n_CTF_rank-1;
tmp_index_ = nCTF_rank*pm_n_lm + (0:pm_n_lm-1);
An_a_UX_Y_yc_ = An_a_UX_Y_yc_ + VSCTF_wMc__(:,1+nCTF_rank) .* (tensor_to_scatter__*reshape(cg_evaluate_n_1(pm_l_max,convert_spharm_to_spharm__0(pm_l_max,a_UX_Y_yc_(1+tmp_index_)),n_polar_a,n_azimu_b,legendre_evaluate_ljm___),[n_polar_a*n_azimu_b,1]));
end;%for nCTF_rank=0:n_CTF_rank-1;

function At_UX_M_ = At(n_CTF_rank,n_M,VSCTF_wMc__,pm_n_w,pm_l_max,UX_M_,scatter_to_tensor__,n_polar_a,n_azimu_b,legendre_evaluate_mlj___,expil__,expi__);
pm_n_lm = (1+pm_l_max).^2;
At_UX_M_ = zeros(pm_n_lm*n_CTF_rank,1);
for nCTF_rank=0:n_CTF_rank-1;
tmp_index_ = nCTF_rank*pm_n_lm + (0:pm_n_lm-1);
At_UX_M_(1+tmp_index_) = convert_spharm__to_spharm_0(pm_l_max,cg_evaluate_t_1(n_polar_a,n_azimu_b,reshape(scatter_to_tensor__*(conj(VSCTF_wMc__(:,1+nCTF_rank)).*UX_M_),[n_polar_a,n_azimu_b]),pm_l_max,legendre_evaluate_mlj___,expil__,expi__));
end;%for nCTF_rank=0:n_CTF_rank-1;
  

  
