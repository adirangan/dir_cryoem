function ...
[ ...
 parameter ...
,a_UCTF_UX_Y_0qbp_yncx___ ...
] = ...
a_UCTF_UX_Y_0qbp_yncx___0( ...
 parameter ...
,pm_n_k_p_r ...
,pm_l_max_ ...
,pm_n_w_ ...
,n_M ...
,UX_M_k_p_wnM__ ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
,image_I_value_...
,n_neighborhood ...
,index_neighborhood_MP__ ...
);
%%%%%%%%;
% Applies quadrature-based back-propagation to solve for a_UX_Y_. ;
% Associates CTF_k_p_r__(:,1+CTF_index_(1+nM)) with image M_k_p__(:,1+nM);
% Constructs a different principal-volume for each neighborhood ;
% listed in index_neighborhood_MP__. ;
% ;
% Input: ;
% quad_k_eq_d: real equatorial distance used for determining quadrature nodes on sphere (radius assumed to be 1). ;
% pm_n_UX_rank = pm_n_k_p_r: integer number of principal-volume shells retained. ;
% pm_l_max_: integer array of size pm_n_k_p_r. pm_l_max_(1+pm_nk_p_r) is the order used for a_UCTF_UX_Y_0qbp_ync__ on principal-volume-shell pm_nk_p_r. ;
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
% image_I_value_: real array of size n_M. I_value used for each image ;
% n_neighborhood: integer number of neighborhoods to consider. ;
% index_neighborhood_MP__: cell array of size n_neighborhood. ;
%                          each cell contains a list of image-indices in a corresponding neighborhood. ;
% ;
% Output: ;
% a_UCTF_UX_Y_0qbp_yncx___: complex array of size (pm_n_lm_sum,n_UCTF_rank,n_neighborhood). output functions in k_Y_ format. ;
% If we define, for a particular nneighborhood, the index-set: ;
% index_nM_ = index_neighborhood_exp_MP__{1+nneighborhood}, ;
% Then for that neighborhood the function: ;
% a_UCTF_UX_Y_ync__ = a_UCTF_UX_Y_0qbp_yncx___(:,:,1+nneighborhood) ;
% should approximately satisfy the least-square problem: ;
% \sum_{nCTF_rank=0}^{n_CTF_rank-1} S * [ \tau_{1+nM} * VSCTF_Mc__(1+nM,1+nCTF_rank) ] * a_UCTF_UX_Y_ync___(:,1+pm_nUX_rank,1+nCTF_rank) = [ UX_M_k_p_wnM___(:,1+pm_nUX_rank,1+nM) ] \forall nM \in index_nM_ and \forall pm_nUX_rank \in [0,pm_n_UX_rank-1]. ;
% where : ;
% \tau_{1+nM} corresponds to rotation by the viewing-angle associated with image nM, and ;
% S is the template-operator (i.e., equatorial-evaluation), and ;
% UX_M_k_p_wnM___(:,1+pm_nUX_rank,1+nM) = UX_M_k_p_wnM__(1+pm_n_w_csum_(1+pm_nUX_rank) + (0:pm_n_w_(1+pm_nUX_rank)-1),1+nM), and ;
% a_UCTF_UX_Y_ync___(:,1+pm_nUX_rank,1+nCTF_rank) = \sum_{nk_p_r=0}^{n_k_p_r} UCTF_(1+nk_p_r,1+nCTF_rank) * UX_(1+nk_p_r,1+pm_nUX_rank) * a_UX_Y__(:,1+nk_p_r). ;
% ;
% Rather than using least-squares to solve this problem, we instead generate a quadrature-grid on the sphere, ;
% map each data-point to its closest quadrature-gridpoint, and then numerically integrate to recover a_UCTF_UX_Y_ync__. ;
% We solve for the dominant CTF-modes of a first, then use an approximate residual to solve for the second (and subsequent) CTF-modes. ;
% This approximate residual is obtained simply by evaluating a_ at the quadrature-grid, and comparing the result to the aggregated data-points. ;
%%%%%%%%;

verbose = 0;
if (verbose>0); disp(sprintf(' %% [entering a_UCTF_UX_Y_0qbp_yncx___0]')); end;
pm_n_UX_rank = pm_n_k_p_r;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'quad_k_eq_d')); %<-- parameter_bookmark. ;
if (~isfield(parameter,'template_viewing_k_eq_d_min'));
parameter.quad_k_eq_d = 1/(2*pi)/2; %<-- parameter_bookmark. ;
end;%if (~isfield(parameter,'template_viewing_k_eq_d_min'));
if ( isfield(parameter,'template_viewing_k_eq_d_min'));
parameter.quad_k_eq_d = parameter.template_viewing_k_eq_d_min/2; %<-- parameter_bookmark. ;
end;%if ( isfield(parameter,'template_viewing_k_eq_d_min'));
end;%if (~isfield(parameter,'quad_k_eq_d'));
%%%%%%%%;
quad_k_eq_d = parameter.quad_k_eq_d;

[ ...
 parameter ...
,pm_n_UX_rank ...
,pm_n_w_ ...
,n_M ...
,UX_M_k_p_wnM__ ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
] = ...
cg_lsq_pm_reduce_1( ...
 parameter ...
,pm_n_UX_rank ...
,pm_n_w_ ...
,n_M ...
,UX_M_k_p_wnM__ ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
);

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

flag_unique_pm_n = 0;
if (numel(unique(pm_l_max_))==1 & numel(unique(pm_n_lm_))==1 & numel(unique(pm_n_w_))==1);
flag_unique_pm_n = 1;
pm_l_max = pm_l_max_(1+0);
Ylm__ = get_Ylm__(1+pm_l_max,0:pm_l_max,quad_n_all,quad_azimu_b_all_,quad_polar_a_all_);
pm_n_lm = pm_n_lm_(1+0);
Ylm_yq__ = zeros(pm_n_lm,quad_n_all);
nml=0;
for l_val=0:pm_l_max;
for m_val=-l_val:+l_val;
Ylm_yq__(1+nml,:) = Ylm__{1+l_val}(1+l_val+m_val,:);
nml=nml+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:pm_l_max;
Ylm_w_yq__ = Ylm_yq__ * sparse(1:quad_n_all,1:quad_n_all,quad_weight_all_,quad_n_all,quad_n_all);
pm_n_w = pm_n_w_(1+0);
[ ...
 data_k_p_polar_a__ ...
,data_k_p_azimu_b__ ...
,data_k_c_0__ ...
,data_k_c_1__ ...
,data_k_c_2__ ...
] = ...
cg_rhs_1( ...
 n_M ...
,pm_n_w ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,+euler_gamma_z_ ...
);
data_k_c_wMd__ = [ data_k_c_0__(:) , data_k_c_1__(:) , data_k_c_2__(:) ];
index_quad_from_data_ = knnsearch(quad_k_c_qd__,data_k_c_wMd__,'K',1); index_quad_from_data_ = index_quad_from_data_ - 1;
index_quad_from_data_wM__ = reshape(index_quad_from_data_,[pm_n_w_max,n_M]);
quad_from_data_qwMx___ = cell(n_neighborhood,1);
n_quad_from_data_qx__ = cell(n_neighborhood,1);
data_from_quad_wMqx___ = cell(n_neighborhood,1);
VSCTF_Mcx___ = cell(n_neighborhood,1);
VSCTF_wMcx___ = cell(n_neighborhood,1);
VSCTF2_qcx___ = cell(n_neighborhood,1);
for nneighborhood=0:n_neighborhood-1;
tmp_index_M_ = index_neighborhood_MP__{1+nneighborhood};
tmp_n_M = numel(tmp_index_M_);
tmp_index_quad_from_data_ = reshape(index_quad_from_data_wM__(:,1+tmp_index_M_),[pm_n_w_max*tmp_n_M,1]);
quad_from_data_qwM__ = sparse(1+tmp_index_quad_from_data_,1:pm_n_w*tmp_n_M,1,quad_n_all,pm_n_w*tmp_n_M);
n_quad_from_data_q_ = quad_from_data_qwM__*ones(pm_n_w*tmp_n_M,1); %<-- number of data-points per quadrature-point. ;
data_from_quad_wMq__ = bsxfun(@rdivide,transpose(quad_from_data_qwM__),max(1,transpose(n_quad_from_data_q_)));
VSCTF_wMc__ = reshape(repmat(reshape(VSCTF_Mc__(1+tmp_index_M_,:),[1,tmp_n_M,n_CTF_rank]),[pm_n_w,1,1]),[pm_n_w*tmp_n_M,n_CTF_rank]);
VSCTF2_qc__ = quad_from_data_qwM__*abs(VSCTF_wMc__).^2;
quad_from_data_qwMx___{1+nneighborhood} = quad_from_data_qwM__;
n_quad_from_data_qx__{1+nneighborhood} = n_quad_from_data_q_;
data_from_quad_wMqx___{1+nneighborhood} = data_from_quad_wMq__;
VSCTF_Mcx___{1+nneighborhood} = VSCTF_Mc__(1+tmp_index_M_,:);
VSCTF_wMcx___{1+nneighborhood} = VSCTF_wMc__;
VSCTF2_qcx___{1+nneighborhood} = VSCTF2_qc__;
end;%for nneighborhood=0:n_neighborhood-1;
end;%if (numel(unique(pm_l_max_))==1 & numel(unique(pm_n_lm_))==1 & numel(unique(pm_n_w_))==1);

if (~flag_unique_pm_n);
disp(sprintf(' %% Warning, consider setting all values of pm_l_max, pm_n_lm and pm_n_w to be the same in cg_lsq_pm_1.'));
end;%if (~flag_unique_pm_n);

a_UCTF_UX_Y_0qbp_yncx___ = zeros(pm_n_lm_sum,n_CTF_rank,n_neighborhood);
for nneighborhood=0:n_neighborhood-1;
tmp_index_M_ = index_neighborhood_MP__{1+nneighborhood};
tmp_n_M = numel(tmp_index_M_);
quad_from_data_qwM__ = quad_from_data_qwMx___{1+nneighborhood};
n_quad_from_data_q_ = n_quad_from_data_qx__{1+nneighborhood};
data_from_quad_wMq__ = data_from_quad_wMqx___{1+nneighborhood};
tmp_VSCTF_Mc__ = VSCTF_Mcx___{1+nneighborhood};
VSCTF_wMc__ = VSCTF_wMcx___{1+nneighborhood};
VSCTF2_qc__ = VSCTF2_qcx___{1+nneighborhood};
%%%%%%%%;
a_UCTF_UX_Y_0qbp_ync__ = zeros(pm_n_lm_sum,n_CTF_rank);
for pm_nk_p_r=0:pm_n_k_p_r-1;
if (~flag_unique_pm_n);
error('Error: non-unique values of pm_l_max, pm_n_lm or pm_n_w in a_UCTF_UX_Y_0qbp_yncx___0');
end;%if (~flag_unique_pm_n);
index_Y_ = pm_n_lm_csum_(1+pm_nk_p_r) + (0:pm_n_lm-1);
index_pm_nw = pm_n_w_csum_(1+pm_nk_p_r) + (0:pm_n_w-1);
UX_M_I_wM__ = bsxfun(@times,UX_M_k_p_wnM__(1+index_pm_nw,1+tmp_index_M_),transpose(image_I_value_(1+tmp_index_M_)));
for nCTF_rank=0:n_CTF_rank-1;
UX_M_I_VSCTF_wM__ = bsxfun(@times,UX_M_I_wM__,transpose(tmp_VSCTF_Mc__(:,1+nCTF_rank)));
quad_from_data_UX_M_I_VSCTF_normalized_q_ = (quad_from_data_qwM__ * UX_M_I_VSCTF_wM__(:))./max(1e-12,VSCTF2_qc__(:,1+nCTF_rank));
a_UCTF_UX_Y_0qbp_ync__(1+index_Y_,1+nCTF_rank) = conj(Ylm_w_yq__)*quad_from_data_UX_M_I_VSCTF_normalized_q_;
UX_M_I_wM__ = UX_M_I_wM__ - bsxfun(@times,reshape(data_from_quad_wMq__*(transpose(Ylm_yq__)*a_UCTF_UX_Y_0qbp_ync__(1+index_Y_,1+nCTF_rank)),[pm_n_w,tmp_n_M]),transpose(tmp_VSCTF_Mc__(:,1+nCTF_rank))); %<-- replace with approximate residual. ;
end;%for nCTF_rank=0:n_CTF_rank-1;
end;%for pm_nk_p_r=0:pm_n_k_p_r-1;
%%%%%%%%;
a_UCTF_UX_Y_0qbp_yncx___(:,:,1+nneighborhood) = a_UCTF_UX_Y_0qbp_ync__;
end;%for nneighborhood=0:n_neighborhood-1;

if (verbose>0); disp(sprintf(' %% [finished a_UCTF_UX_Y_0qbp_yncx___0]')); end;  
