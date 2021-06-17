function ...
[ ...
 parameter ...
,SV_loading_Ml__ ...
] = ...
get_loading_qbp_0( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p_wnM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wnc__ ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
,image_delta_x_ ...
,image_delta_y_ ...
,image_I_value_ ...
);

verbose = 1;
if (verbose>1); disp(sprintf(' %% [entering get_loading_qbp_0]')); end;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'rseed')); parameter.rseed = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_loading')); parameter.n_loading = 3; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_loading_iteration')); parameter.n_loading_iteration = 32; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'cg_lsq_n_order')); parameter.cg_lsq_n_order = 5; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'qbp_eps')); parameter.qbp_eps = 1e-3; end; %<-- parameter_bookmark. ;
%%%%%%%%;
tolerance_master = parameter.tolerance_master;
rseed = parameter.rseed;
n_loading = parameter.n_loading;
n_loading_iteration = parameter.n_loading_iteration;
cg_lsq_n_order = parameter.cg_lsq_n_order;
qbp_eps = parameter.qbp_eps;

if (nargin<14); image_delta_x_ = []; end;
if (nargin<15); image_delta_y_ = []; end;
if (nargin<16); image_I_value_ = []; end;

if isempty(image_delta_x_); image_delta_x_ = zeros(n_M,1); end;
if isempty(image_delta_y_); image_delta_y_ = zeros(n_M,1); end;
if isempty(image_I_value_); image_I_value_ = ones(n_M,1); end;
if isempty(index_nCTF_from_nM_); index_nCTF_from_nM_ = zeros(n_M,1); end;
if isempty(CTF_k_p_wnc__); CTF_k_p_wnc__ = ones(n_w_sum,1); end;
if (qbp_eps>1); qbp_eps = max(1e-12,0.1^(qbp_eps)); end;

verbose=1;
%%%%%%%%;
l_max_ = l_max_(1:n_k_p_r);
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); 
n_w_max = max(n_w_(1:n_k_p_r));
n_w_sum = sum(n_w_(1:n_k_p_r));
n_w_csum_ = cumsum([0;n_w_(1:n_k_p_r)]);
%%%%%%%%;
% Solve qbp. ;
%%%%%%%%;
tmp_t = tic();
a_k_Y_0qbp_yn_ ...
= ...
qbp_5(...
 qbp_eps ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p_wnM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wnc__ ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
,image_delta_x_ ...
,image_delta_y_ ...
,image_I_value_...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_k_Y_0qbp_yn_: %0.3fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% Center and scale images. ;
%%%%%%%%;
T_k_p_wnM__ = M_k_p_wnM__;
if ( (fnorm(image_delta_x_)>0) | (fnorm(image_delta_y_)>0) | (fnorm(image_I_value_-ones(n_M,1))>0) );
for nM=0:n_M-1;
T_k_p_wnM__(:,1+nM) = image_I_value_(1+nM) * transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wnM__(:,1+nM),+image_delta_x_(1+nM),+image_delta_y_(1+nM));
end;%for nM=0:n_M-1;
end;%if ( (fnorm(image_delta_x_)>0) | (fnorm(image_delta_y_)>0) | (fnorm(image_I_value_-ones(n_M,1))>0) );

quad_k_eq_d_ = sqrt(4*pi./n_lm_);
if (verbose>1); for nk_p_r=0:n_k_p_r-1; disp(sprintf(' %% nk_p_r %d/%d: quad_k_eq_d %0.6f',nk_p_r,n_k_p_r,quad_k_eq_d_(1+nk_p_r))); end; end;

CTF_k_p_wM__ = cell(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
CTF_k_p_wM__{1+nk_p_r} = reshape(CTF_k_p_wnc__(1+index_nw_,1+index_nCTF_from_nM_),[n_w*n_M,1]);
end;%for nk_p_r=0:n_k_p_r-1;

%%%%%%%%;
% precomputation. ;
%%%%%%%%;
flag_unique_n = 0;
if (numel(unique(l_max_))==1 & numel(unique(n_lm_))==1 & numel(unique(n_w_))==1);
flag_unique_n = 1;
quad_k_eq_d = quad_k_eq_d_(1+0);
l_max = l_max_(1+0);
n_lm = n_lm_(1+0);
n_w = n_w_(1+0);
weight_3d_yn_ = reshape(repmat(reshape(weight_3d_k_p_r_,[1,n_k_p_r]),[n_lm,1]),[n_lm*n_k_p_r,1]);
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
quad_from_data_Mqw___ = cell(n_M,1);
tmp_t = tic();
for nM=0:n_M-1;
index_nw_ = n_w*nM + (0:n_w-1);
quad_from_data_Mqw___{1+nM} = quad_from_data_qwM__(:,1+index_nw_);
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% quad_from_data_Mqw___: %0.3fs',tmp_t)); end;
n_quad_from_data_q_ = quad_from_data_qwM__*ones(n_w*n_M,1); %<-- number of data-points per quadrature-point. ;
data_from_quad_wMq__ = bsxfun(@rdivide,transpose(quad_from_data_qwM__),max(1,transpose(n_quad_from_data_q_)));
CTF_wMn__ = reshape(permute(reshape(CTF_k_p_wnc__(:,1+index_nCTF_from_nM_),[n_w,n_k_p_r,n_M]),[1,3,2]),[n_w*n_M,n_k_p_r]);
CTF2_qn__ = quad_from_data_qwM__*abs(CTF_wMn__).^2;
%%%%;
[k_p_polar_a__,k_p_azimu_b__] = cg_rhs_1(n_M,n_w,euler_polar_a_,euler_azimu_b_,+euler_gamma_z_);
n_polar_a = ceil(n_w/2);
n_azimu_b = max(1+2*l_max,2*n_polar_a);
[legendre_evaluate_ljm___,legendre_evaluate_mlj___,expil__,expi__] = legendre_evaluate_ljm___0(l_max,cos(linspace(0,pi,n_polar_a)),n_azimu_b);
tensor_to_scatter__ = cg_interpolate_n_1(cg_lsq_n_order,n_polar_a,n_azimu_b,n_w*n_M,k_p_polar_a__(:),k_p_azimu_b__(:));
scatter_to_tensor__ = ctranspose(tensor_to_scatter__); %<-- this conjugation is not necessary, since the matrix should be real. ;
%%%%;
An___ = cell(n_k_p_r,1);
At___ = cell(n_k_p_r,1);
AtAn___ = cell(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
An___{1+nk_p_r} = @(a_k_Y_) CTF_k_p_wM__{1+nk_p_r}.*(tensor_to_scatter__*reshape(cg_evaluate_n_1(l_max,convert_spharm_to_spharm__0(l_max,a_k_Y_),n_polar_a,n_azimu_b,legendre_evaluate_ljm___),[n_polar_a*n_azimu_b,1]));
At___{1+nk_p_r} = @(M_k_p_) convert_spharm__to_spharm_0(l_max,cg_evaluate_t_1(n_polar_a,n_azimu_b,reshape(scatter_to_tensor__*(conj(CTF_k_p_wM__{1+nk_p_r}).*M_k_p_),[n_polar_a,n_azimu_b]),l_max,legendre_evaluate_mlj___,expil__,expi__));
AtAn___{1+nk_p_r} = @(a_k_Y_) At___{1+nk_p_r}(An___{1+nk_p_r}(a_k_Y_));
end;%for nk_p_r=0:n_k_p_r-1;
end;%if (numel(unique(l_max_))==1 & numel(unique(n_lm_))==1 & numel(unique(n_w_))==1);

if (flag_unique_n==0);
error(sprintf(' %% Error, consider setting l_max_ and n_w_ to be uniform in get_loading_qbp_0'));
end;%if (flag_unique_n==0);

%%%%%%%%;
% calculate residual from qbp. ;
%%%%%%%%;
tmp_t = tic();
S_k_p_wnM__ = zeros(n_w_sum,n_M);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
An__ = An___{1+nk_p_r}; At__ = At___{1+nk_p_r};
S_k_p_wnM__(1+index_nw_,:) = reshape(An__(a_k_Y_0qbp_yn_(1+index_Y_)),[n_w,n_M]);
end;%for nk_p_r=0:n_k_p_r-1;
R_k_p_wnM__ = T_k_p_wnM__ - S_k_p_wnM__;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% R_k_p_wnM__: %0.3fs',tmp_t)); end;
if (verbose); disp(sprintf(' %% fnorm(R_k_p_wnM__)/fnorm(T_k_p_wnM__): %0.6f/%0.6f = %0.6f',fnorm(R_k_p_wnM__),fnorm(T_k_p_wnM__),fnorm(R_k_p_wnM__)/fnorm(T_k_p_wnM__))); end;

%%%%%%%%;
% Initialize residual clustering. ;
%%%%%%%%;
tmp_t = tic();
H_ynM__ = zeros(n_lm_sum,n_M);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
R_wM__ = R_k_p_wnM__(1+index_nw_,:);
CTF_wM_ = CTF_wMn__(:,1+nk_p_r);
CTF2_q_ = CTF2_qn__(:,1+nk_p_r);
R_CTF_wM__ = R_wM__.*reshape(CTF_wM_,[n_w,n_M]);
quad_from_data_R_CTF_normalized_qM__ = zeros(quad_n_all,n_M);
for nM=0:n_M-1;
quad_from_data_R_CTF_normalized_qM__(:,1+nM) = (quad_from_data_Mqw___{1+nM} * R_CTF_wM__(:,1+nM))./max(1e-12,CTF2_q_);
end;%for nM=0:n_M-1;
H_ynM__(1+index_Y_,:) = conj(Ylm_w_yq__)*quad_from_data_R_CTF_normalized_qM__;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% H_ynM__: %0.3fs',tmp_t)); end;
%%%%%%%%;
HH_M_ = sum(bsxfun(@times,abs(H_ynM__).^2,reshape(weight_3d_yn_,[n_lm_sum,1])),1);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Iterate to cluster residuals. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_t = tic();
rng(rseed);
SV_loading_Ml__ = randn(n_M,n_loading);
U_ynl__ = randn(n_lm_sum,n_loading); 
U_l2_l_ = zeros(n_loading,1);
for nloading=0:n_loading-1;
U_l2_l_(1+nloading) = sqrt(ctranspose(U_ynl__(:,1+nloading))*(weight_3d_yn_.*U_ynl__(:,1+nloading))); 
U_ynl__(:,1+nloading) = U_ynl__(:,1+nloading)/U_l2_l_(1+nloading);
end;%for nloading=0:n_loading-1; 
for niteration=0:n_loading_iteration-1;
%%%%;
if (niteration==0);
% do nothing;
end;%if (niteration==0);
if (niteration>0);
U_ynl__ = H_ynM__*SV_loading_Ml__;
for nloading_0=0:n_loading-1;
U_l2_l_(1+nloading_0) = sqrt(ctranspose(U_ynl__(:,1+nloading_0))*(weight_3d_yn_.*U_ynl__(:,1+nloading_0))); 
U_ynl__(:,1+nloading_0) = U_ynl__(:,1+nloading_0)/U_l2_l_(1+nloading_0);
for nloading_1=nloading_0+1:n_loading-1;
U_ynl__(:,1+nloading_1) = U_ynl__(:,1+nloading_1) - (ctranspose(U_ynl__(:,1+nloading_0))*(weight_3d_yn_.*U_ynl__(:,1+nloading_1)))*U_ynl__(:,1+nloading_0);
end;%for nloading_1=nloading_0+1:n_loading-1;
end;%for nloading_0=0:n_loading-1;
end;%if (niteration>0);
%%%%;
for nM=0:n_M-1;
HH = HH_M_(1+nM);
for nloading=0:n_loading-1;
HU = ctranspose(H_ynM__(:,1+nM))*(weight_3d_yn_.*U_ynl__(:,1+nloading));
SV_loading_Ml__(1+nM,1+nloading) = real(HU)/max(1e-12,real(HH));
end;%for nloading=0:n_loading-1;
end;%for nM=0:n_M-1;
%%%%;
end;%for niteration=0:n_loading_iteration-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% SV_loading_Ml__: %0.3fs',tmp_t)); end;
%%%%;

if (verbose>1); disp(sprintf(' %% [finished get_loading_qbp_0]')); end;
