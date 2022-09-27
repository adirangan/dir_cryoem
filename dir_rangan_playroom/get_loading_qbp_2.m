function ...
[ ...
 parameter ...
,SV_loading_Ml__ ...
,a_k_Y_0qbp_yn_ ...
,T_k_p_wnM__ ...
,CTF_k_p_rwM__ ...
,S_k_p_wnM__ ...
,R_k_p_wnM__ ...
,AUSV_wnM__ ...
,H_ynM__ ...
,U_ynl__ ...
,S_l_ ...
,V_Ml__ ...
] = ...
get_loading_qbp_2( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
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

flag_verbose = 1; flag_disp = flag_verbose; nf=0;
if (flag_verbose>1); disp(sprintf(' %% [entering get_loading_qbp_2]')); end;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wnM__=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_wnc__=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_=[]; end; na=na+1;
if (nargin<1+na); euler_gamma_z_=[]; end; na=na+1;
if (nargin<1+na); image_delta_x_=[]; end; na=na+1;
if (nargin<1+na); image_delta_y_=[]; end; na=na+1;
if (nargin<1+na); image_I_value_=[]; end; na=na+1;

SV_loading_Ml__ = [];
a_k_Y_0qbp_yn_ = [];
T_k_p_wnM__ = [];
CTF_k_p_rwM__ = [];
S_k_p_wnM__ = [];
R_k_p_wnM__ = [];
AUSV_wnM__ = [];
H_ynM__ = [];
U_ynl__ = [];
S_l_ = [];
V_Ml__ = [];

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'rseed')); parameter.rseed = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_loading')); parameter.n_loading = 3; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_loading_iteration')); parameter.n_loading_iteration = 64; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'cg_lsq_n_order')); parameter.cg_lsq_n_order = 5; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'qbp_eps')); parameter.qbp_eps = 1e-3; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_loading_svd_vs_iterate')); parameter.flag_loading_svd_vs_iterate = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_loading_skip_loading')); parameter.flag_loading_skip_loading = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_SV_lsq_vs_dot')); parameter.flag_SV_lsq_vs_dot = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_U_lsq_vs_dot')); parameter.flag_U_lsq_vs_dot = 0; end; %<-- parameter_bookmark. ;
%%%%%%%%;
tolerance_master = parameter.tolerance_master;
rseed = parameter.rseed;
n_loading = parameter.n_loading;
n_loading_iteration = parameter.n_loading_iteration;
cg_lsq_n_order = parameter.cg_lsq_n_order;
qbp_eps = parameter.qbp_eps;
flag_loading_svd_vs_iterate = parameter.flag_loading_svd_vs_iterate;
flag_loading_skip_loading = parameter.flag_loading_skip_loading;
flag_SV_lsq_vs_dot =  parameter.flag_SV_lsq_vs_dot;
flag_U_lsq_vs_dot =  parameter.flag_U_lsq_vs_dot;

if isempty(image_delta_x_); image_delta_x_ = zeros(n_M,1); end;
if isempty(image_delta_y_); image_delta_y_ = zeros(n_M,1); end;
if isempty(image_I_value_); image_I_value_ = ones(n_M,1); end;
if isempty(index_nCTF_from_nM_); index_nCTF_from_nM_ = zeros(n_M,1); end;
if isempty(CTF_k_p_wnc__); CTF_k_p_wnc__ = ones(n_w_sum,1); end;
if (qbp_eps>1); qbp_eps = max(1e-12,0.1^(qbp_eps)); end;

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
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% a_k_Y_0qbp_yn_: %0.3fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% Center and scale images. ;
% T_k_p_wnM__ hold the centered and scaled images. ;
%%%%%%%%;
T_k_p_wnM__ = M_k_p_wnM__;
if ( (fnorm(image_delta_x_)>0) | (fnorm(image_delta_y_)>0) | (fnorm(image_I_value_-ones(n_M,1))>0) );
for nM=0:n_M-1;
T_k_p_wnM__(:,1+nM) = image_I_value_(1+nM) * transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wnM__(:,1+nM),+image_delta_x_(1+nM),+image_delta_y_(1+nM));
end;%for nM=0:n_M-1;
end;%if ( (fnorm(image_delta_x_)>0) | (fnorm(image_delta_y_)>0) | (fnorm(image_I_value_-ones(n_M,1))>0) );

quad_k_eq_d_ = sqrt(4*pi./n_lm_);
if (flag_verbose>1); for nk_p_r=0:n_k_p_r-1; disp(sprintf(' %% nk_p_r %d/%d: quad_k_eq_d %0.6f',nk_p_r,n_k_p_r,quad_k_eq_d_(1+nk_p_r))); end; end;

%%%%%%%%;
% CTF_k_p_rwM__ hold the CTF_k_p_wr__ for each image. ;
% Specifically: ;
% CTF_k_p_rwM__{1+nk_p_r}(1+nw+nM*n_w_(1+nk_p_r)) holds CTF_k_p_wr__ for image nM. ;
%%%%%%%%;
CTF_k_p_rwM__ = cell(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
CTF_k_p_rwM__{1+nk_p_r} = reshape(CTF_k_p_wnc__(1+index_nw_,1+index_nCTF_from_nM_),[n_w*n_M,1]);
end;%for nk_p_r=0:n_k_p_r-1;

%%%%%%%%;
% precomputation. ;
% We use sample_shell_5 to construct a quadrature grid on the shell, ;
% referenced by quad_n_all through quad_k_c_2_all_. ;
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
% The spherical-harmonics on the shell are referenced by Ylm__ through Ylm_yq__. ;
% The quadrature-weighted spherical-harmonics (for integration) are referenced by Ylm_w_yq__. ;
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
% The point-locations of the images (oriented according to their euler-angles) on the shell ;
% are referenced by data_k_p_polar_a__ through data_k_c_wMd__. ;
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
% We construct a sparse matrix linking the point-locations on the shell to the quadrature-nodes. ;
% This is referenced by quad_from_data_qwM__ and quad_from_data_Mqw___. ;
%%%%;
index_quad_from_data_ = knnsearch(quad_k_c_qd__,data_k_c_wMd__,'K',1); index_quad_from_data_ = index_quad_from_data_ - 1;
quad_from_data_qwM__ = sparse(1+index_quad_from_data_,1:n_w*n_M,1,quad_n_all,n_w*n_M);
quad_from_data_Mqw___ = cell(n_M,1);
tmp_t = tic();
for nM=0:n_M-1;
index_nw_ = n_w*nM + (0:n_w-1);
quad_from_data_Mqw___{1+nM} = quad_from_data_qwM__(:,1+index_nw_);
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% quad_from_data_Mqw___: %0.3fs',tmp_t)); end;
%%%%;
% We then calculate the number of data-points sent to each quadrature-node, ;
% (i.e., the row-sums of quad_from_data_qwM__). ;
% and use these values to normalize the transpose of the matrix above. ;
% This creates a matrix which links the quadrature-nodes to the point-locations on the shell, ;
% and has column-sums equal to one. ;
% This matrix is denoted by data_from_quad_wMq__, but is unused below. ;
%%%%;
n_quad_from_data_q_ = quad_from_data_qwM__*ones(n_w*n_M,1); %<-- number of data-points per quadrature-point. ;
data_from_quad_wMq__ = bsxfun(@rdivide,transpose(quad_from_data_qwM__),max(1,transpose(n_quad_from_data_q_)));
%%%%;
% The sum of the squares of the CTF-coefficients mapped to each quadrature-node ;
% is referenced by CTF2_qn__. ;
% We will use this below to solve a local-least-squares problem for each quadrature-node. ;
%%%%;
CTF_wMn__ = reshape(permute(reshape(CTF_k_p_wnc__(:,1+index_nCTF_from_nM_),[n_w,n_k_p_r,n_M]),[1,3,2]),[n_w*n_M,n_k_p_r]);
CTF2_qn__ = quad_from_data_qwM__*abs(CTF_wMn__).^2;
%%%%;
% We set up the operators An___ and At___ for each of the nk_p_r, ;
% representing the Y-to-data evaluation-operator (composed with CTF-multiplication) ;
% along with its transpose. ;
% These will not necessarily be used in calculating the loadings. ;
%%%%;
%[k_p_polar_a__,k_p_azimu_b__] = cg_rhs_1(n_M,n_w,euler_polar_a_,euler_azimu_b_,+euler_gamma_z_);
k_p_polar_a__ = data_k_p_polar_a__; k_p_azimu_b__ = data_k_p_azimu_b__;
n_polar_a = ceil(n_w/2);
n_azimu_b = max(1+2*l_max,2*n_polar_a);
[legendre_evaluate_ljm___,legendre_evaluate_mlj___,expil__,expi__] = legendre_evaluate_ljm___0(l_max,cos(linspace(0,pi,n_polar_a)),n_azimu_b);
tensor_to_scatter__ = cg_interpolate_n_1(cg_lsq_n_order,n_polar_a,n_azimu_b,n_w*n_M,k_p_polar_a__(:),k_p_azimu_b__(:));
scatter_to_tensor__ = ctranspose(tensor_to_scatter__); %<-- this conjugation is not necessary, since the matrix should be real. ;
%%%%;
% Note that, even though nk_p_r, tensor_to_scatter, l_max, n_polar_a, n_azimu_b, legendre_evaluate_ljm___, etc. will all be overwritten, ;
% the functional definitions of An___ and At___ will remain static. ;
%%%%;
An___ = cell(n_k_p_r,1);
At___ = cell(n_k_p_r,1);
AtAn___ = cell(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
An___{1+nk_p_r} = @(a_k_Y_) CTF_k_p_rwM__{1+nk_p_r}.*(tensor_to_scatter__*reshape(cg_evaluate_n_1(l_max,convert_spharm_to_spharm__0(l_max,a_k_Y_),n_polar_a,n_azimu_b,legendre_evaluate_ljm___),[n_polar_a*n_azimu_b,1]));
At___{1+nk_p_r} = @(M_k_p_) convert_spharm__to_spharm_0(l_max,cg_evaluate_t_1(n_polar_a,n_azimu_b,reshape(scatter_to_tensor__*(conj(CTF_k_p_rwM__{1+nk_p_r}).*M_k_p_),[n_polar_a,n_azimu_b]),l_max,legendre_evaluate_mlj___,expil__,expi__));
AtAn___{1+nk_p_r} = @(a_k_Y_) At___{1+nk_p_r}(An___{1+nk_p_r}(a_k_Y_));
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
% We set up qbp_An and qbp_Ai, the qbp evaluation and inverse operators. ;
% Again, these will not necessarily be used in calculating the loadings. ;
%%%%;
qbp_An = @(a_k_Y_yn_) ...
get_loading_qbp_An( ...
 n_k_p_r ...
,l_max_ ...
,n_w_ ...
,a_k_Y_yn_ ...
,legendre_evaluate_ljm___ ...
,tensor_to_scatter__ ...
,n_M ...
,CTF_k_p_rwM__ ...
);
qbp_Ai = @(P_k_p_wnM__) ...
get_loading_qbp_Ai( ...
 n_k_p_r ...
,n_w_ ...
,n_M ...
,P_k_p_wnM__ ...
,CTF_wMn__ ...
,CTF2_qn__ ...
,quad_from_data_qwM__ ...
,Ylm_w_yq__ ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if (numel(unique(l_max_))==1 & numel(unique(n_lm_))==1 & numel(unique(n_w_))==1);

if (flag_unique_n==0);
error(sprintf(' %% Error, consider setting l_max_ and n_w_ to be uniform in get_loading_qbp_2'));
end;%if (flag_unique_n==0);

%%%%%%%%;
% calculate residual from qbp. ;
% The S_k_p_wnM__ represent the templates, and R_k_p_wnM__ the residuals. ;
%%%%%%%%;
tmp_t = tic();
S_k_p_wnM__ = qbp_An(a_k_Y_0qbp_yn_);
R_k_p_wnM__ = T_k_p_wnM__ - S_k_p_wnM__;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% R_k_p_wnM__: %0.3fs',tmp_t)); end;
if (flag_verbose); disp(sprintf(' %% fnorm(R_k_p_wnM__)/fnorm(T_k_p_wnM__): %0.6f/%0.6f = %0.6f',fnorm(R_k_p_wnM__),fnorm(T_k_p_wnM__),fnorm(R_k_p_wnM__)/fnorm(T_k_p_wnM__))); end;
%%%%%%%%;
flag_check=0;
if flag_check;
tmp_t = tic();
tmp_S_k_p_wnM__ = zeros(n_w_sum,n_M);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
An__ = An___{1+nk_p_r}; At__ = At___{1+nk_p_r};
tmp_S_k_p_wnM__(1+index_nw_,:) = reshape(An__(a_k_Y_0qbp_yn_(1+index_Y_)),[n_w,n_M]);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_R_k_p_wnM__ = T_k_p_wnM__ - tmp_S_k_p_wnM__;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% tmp_R_k_p_wnM__: %0.3fs',tmp_t)); end;
if (flag_verbose); disp(sprintf(' %% S_k_p_wnM__ vs tmp_S_k_p_wnM__: %0.16f',fnorm(S_k_p_wnM__ - tmp_S_k_p_wnM__)/fnorm(S_k_p_wnM__))); end;
if (flag_verbose); disp(sprintf(' %% R_k_p_wnM__ vs tmp_R_k_p_wnM__: %0.16f',fnorm(R_k_p_wnM__ - tmp_R_k_p_wnM__)/fnorm(R_k_p_wnM__))); end;
end;%if flag_check;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
flag_check=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_check;
%%%%%%%%;
%{
if (flag_verbose); disp(sprintf(' %% Testing flag_U_lsq_vs_dot with basic reconstruction. ;')); end;
tmp_t = tic();
a_k_Y_0lsq_ = ...
cg_lsq_3( ...
 cg_lsq_n_order ...
,n_k_p_r ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p_wnM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wnc__ ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,+euler_gamma_z_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% a_k_Y_0lsq_: %0.3fs',tmp_t)); end;
tmp_t = tic();
a_k_Y_1lsq_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
An__ = An___{1+nk_p_r}; At__ = At___{1+nk_p_r}; AtAn__ = AtAn___{1+nk_p_r};
tmp_LHS_y_ = zeros(n_lm,1);
tmp_M_k_p_wM_ = reshape(M_k_p_wnM__(1+index_nw_,:),[n_w*n_M,1]);
tmp_LHS_y_ = VSAt__(n_lm,1,n_w,n_M,At__,ones(n_M,1),tmp_M_k_p_wM_);
tmp_RHS_yy__ = @(a_k_Y_yl_) VSAtAnSV__(n_lm,1,n_w,n_M,At__,An__,ones(n_M,1),a_k_Y_yl_);
[tmp_a_y_,~] = pcg(tmp_RHS_yy__,tmp_LHS_y_);
a_k_Y_1lsq_(1+index_Y_,1) = reshape(tmp_a_y_,[n_lm,1]);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% a_k_Y_1lsq_: %0.3fs',tmp_t)); end;
if (flag_verbose); disp(sprintf(' %% a_k_Y_0lsq_ vs a_k_Y_1lsq_: %0.16f',fnorm(a_k_Y_0lsq_ - a_k_Y_1lsq_)/fnorm(a_k_Y_0lsq_))); end;
if (flag_verbose); disp(sprintf(' %% real(corr(a_k_Y_0lsq_,a_k_Y_1lsq_)): %0.16f',real(corr(a_k_Y_0lsq_,a_k_Y_1lsq_)))); end;
Ana_k_p_wnM__ = qbp_An(a_k_Y_0lsq_);
n_lin = 1+8; tmp_eps_ = linspace(-1,+1,n_lin); tmp_eps_scale = 0.1*fnorm(a_k_Y_0lsq_);
n_test = 16; l2_et_ = zeros(n_lin,n_test);
for ntest=0:n_test-1;
if (flag_verbose); disp(sprintf(' %% ntest %d/%d',ntest,n_test)); end;
tmp_da_yn_ = (randn(n_lm_sum,1) + i*randn(n_lm_sum,1));
for nlin=0:n_lin-1;
tmp_eps = tmp_eps_scale*tmp_eps_(1+nlin);
b_k_Y_0lsq_ = a_k_Y_0lsq_ + tmp_eps*tmp_da_yn_;
Anb_k_p_wnM__ = qbp_An(b_k_Y_0lsq_); Anc_k_p_wnM__ = Anb_k_p_wnM__ - M_k_p_wnM__;
l2_et__(1+nlin,1+ntest) = real(sum(innerproduct_p_quad_(n_k_p_r,weight_2d_k_p_r_,n_w_,n_M,Anc_k_p_wnM__,Anc_k_p_wnM__)));
end;%for nlin=0:n_lin-1;
end;%for ntest=0:n_test-1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
hold on;
plot(0,l2_et__(1+(n_lin-1)/2,:),'ko','MarkerFaceColor','r');
plot(tmp_eps_,l2_et__,'kx'); xlabel('tmp_eps_','Interpreter','none'); ylabel('l2 error');
hold off;
end;%if flag_disp;
return;
 %}
%%%%%%%%;
%{
if (flag_verbose); disp(sprintf(' %% Testing flag_U_lsq_vs_dot with residual reconstruction. ;')); end;
tmp_t = tic();
a_k_Y_0lsq_ = ...
cg_lsq_3( ...
 cg_lsq_n_order ...
,n_k_p_r ...
,l_max_ ...
,n_w_ ...
,n_M ...
,R_k_p_wnM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wnc__ ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,+euler_gamma_z_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% a_k_Y_0lsq_: %0.3fs',tmp_t)); end;
tmp_t = tic();
a_k_Y_1lsq_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
An__ = An___{1+nk_p_r}; At__ = At___{1+nk_p_r}; AtAn__ = AtAn___{1+nk_p_r};
tmp_LHS_y_ = zeros(n_lm,1);
tmp_R_k_p_wM_ = reshape(R_k_p_wnM__(1+index_nw_,:),[n_w*n_M,1]);
tmp_LHS_y_ = VSAt__(n_lm,1,n_w,n_M,At__,ones(n_M,1),tmp_R_k_p_wM_);
tmp_RHS_yy__ = @(a_k_Y_yl_) VSAtAnSV__(n_lm,1,n_w,n_M,At__,An__,ones(n_M,1),a_k_Y_yl_);
[tmp_a_y_,~] = pcg(tmp_RHS_yy__,tmp_LHS_y_);
a_k_Y_1lsq_(1+index_Y_,1) = reshape(tmp_a_y_,[n_lm,1]);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% a_k_Y_1lsq_: %0.3fs',tmp_t)); end;
if (flag_verbose); disp(sprintf(' %% a_k_Y_0lsq_ vs a_k_Y_1lsq_: %0.16f',fnorm(a_k_Y_0lsq_ - a_k_Y_1lsq_)/fnorm(a_k_Y_0lsq_))); end;
if (flag_verbose); disp(sprintf(' %% real(corr(a_k_Y_0lsq_,a_k_Y_1lsq_)): %0.16f',real(corr(a_k_Y_0lsq_,a_k_Y_1lsq_)))); end;
Ana_k_p_wnM__ = qbp_An(a_k_Y_0lsq_);
n_lin = 1+8; tmp_eps_ = linspace(-1,+1,n_lin); tmp_eps_scale = 0.1*fnorm(a_k_Y_0lsq_);
n_test = 8; l2_et_ = zeros(n_lin,n_test);
for ntest=0:n_test-1;
if (flag_verbose); disp(sprintf(' %% ntest %d/%d',ntest,n_test)); end;
tmp_da_yn_ = (randn(n_lm_sum,1) + i*randn(n_lm_sum,1));
for nlin=0:n_lin-1;
tmp_eps = tmp_eps_scale*tmp_eps_(1+nlin);
b_k_Y_0lsq_ = a_k_Y_0lsq_ + tmp_eps*tmp_da_yn_;
Anb_k_p_wnM__ = qbp_An(b_k_Y_0lsq_); Anc_k_p_wnM__ = Anb_k_p_wnM__ - R_k_p_wnM__;
l2_et__(1+nlin,1+ntest) = real(sum(innerproduct_p_quad_(n_k_p_r,weight_2d_k_p_r_,n_w_,n_M,Anc_k_p_wnM__,Anc_k_p_wnM__)));
end;%for nlin=0:n_lin-1;
end;%for ntest=0:n_test-1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
hold on;
plot(0,l2_et__(1+(n_lin-1)/2,:),'ko','MarkerFaceColor','r');
plot(tmp_eps_,l2_et__,'kx'); xlabel('tmp_eps_','Interpreter','none'); ylabel('l2 error');
hold off;
end;%if flag_disp;
return;
 %}
%%%%%%%%;
%{
if (flag_verbose); disp(sprintf(' %% Testing flag_U_lsq_vs_dot with residual reconstruction and single loading. ;')); end;
tmp_t = tic();
tmp_n_loading = 1; rng(0);
tmp_SV_loading_Ml__ = randn(n_M,tmp_n_loading);
a_k_Y_1lsq_ynl__ = zeros(n_lm_sum,tmp_n_loading);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
An__ = An___{1+nk_p_r}; At__ = At___{1+nk_p_r}; AtAn__ = AtAn___{1+nk_p_r};
tmp_LHS_yl_ = zeros(n_lm*tmp_n_loading,1);
tmp_R_k_p_wM_ = reshape(R_k_p_wnM__(1+index_nw_,:),[n_w*n_M,1]);
tmp_LHS_yl_ = VSAt__(n_lm,tmp_n_loading,n_w,n_M,At__,tmp_SV_loading_Ml__,tmp_R_k_p_wM_);
tmp_RHS_ylyl__ = @(a_k_Y_yl_) VSAtAnSV__(n_lm,tmp_n_loading,n_w,n_M,At__,An__,tmp_SV_loading_Ml__,a_k_Y_yl_);
[tmp_a_yl_,~] = pcg(tmp_RHS_ylyl__,tmp_LHS_yl_);
a_k_Y_1lsq_ynl__(1+index_Y_,:) = reshape(tmp_a_yl_,[n_lm,tmp_n_loading]);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% a_k_Y_1lsq_ynl__: %0.3fs',tmp_t)); end;
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;figmed;
subplot(1,1,1); plot(log(abs(a_k_Y_1lsq_ynl__(:,1+0))),'.'); title(sprintf('1lsq: nloading %d',0));
end;%if flag_disp;
%%%%;
tmp_SV_loading_0_Ml__ = SV_lsq(n_k_p_r,weight_2d_k_p_r_,n_w_,n_M,R_k_p_wnM__,l_max_,tmp_n_loading,qbp_An,a_k_Y_1lsq_ynl__);
n_lin = 1+14; tmp_eps_ = linspace(-1,+1,n_lin); tmp_eps_scale = 0.01*fnorm(tmp_SV_loading_0_Ml__);
n_test = 14; l2_et_ = zeros(n_lin,n_test);
for ntest=0:n_test-1;
if (flag_verbose); disp(sprintf(' %% ntest %d/%d',ntest,n_test)); end;
tmp_dSV_Ml__ = randn(n_M,tmp_n_loading);
for nlin=0:n_lin-1;
tmp_eps = tmp_eps_scale*tmp_eps_(1+nlin);
tmp_SV_loading_1_Ml__ = tmp_SV_loading_0_Ml__ + tmp_eps*tmp_dSV_Ml__;
dR_k_p_wnM__ = AnSV_all__(n_k_p_r,n_w_,n_M,l_max_,tmp_n_loading,An___,tmp_SV_loading_1_Ml__,a_k_Y_1lsq_ynl__) - R_k_p_wnM__;
l2_et__(1+nlin,1+ntest) = real(sum(innerproduct_p_quad_(n_k_p_r,weight_2d_k_p_r_,n_w_,n_M,dR_k_p_wnM__,dR_k_p_wnM__)));
end;%for nlin=0:n_lin-1;
end;%for ntest=0:n_test-1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
hold on;
plot(0,l2_et__(1+(n_lin-1)/2,:),'ko','MarkerFaceColor','r');
plot(tmp_eps_,l2_et__,'kx'); xlabel('tmp_eps_','Interpreter','none'); ylabel('l2 error');
hold off;
end;%if flag_disp;
%%%%;
return;
 %}
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% Testing flag_U_lsq_vs_dot with residual reconstruction and two loadings. ;')); end;
tmp_t = tic();
tmp_n_loading = 2; rng(0);
tmp_SV_loading_Ml__ = randn(n_M,tmp_n_loading);
a_k_Y_2lsq_ynl__ = zeros(n_lm_sum,tmp_n_loading);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
An__ = An___{1+nk_p_r}; At__ = At___{1+nk_p_r}; AtAn__ = AtAn___{1+nk_p_r};
tmp_LHS_yl_ = zeros(n_lm*tmp_n_loading,1);
tmp_R_k_p_wM_ = reshape(R_k_p_wnM__(1+index_nw_,:),[n_w*n_M,1]);
tmp_LHS_yl_ = VSAt__(n_lm,tmp_n_loading,n_w,n_M,At__,tmp_SV_loading_Ml__,tmp_R_k_p_wM_);
tmp_RHS_ylyl__ = @(a_k_Y_yl_) VSAtAnSV__(n_lm,tmp_n_loading,n_w,n_M,At__,An__,tmp_SV_loading_Ml__,a_k_Y_yl_);
[tmp_a_yl_,~] = pcg(tmp_RHS_ylyl__,tmp_LHS_yl_);
a_k_Y_2lsq_ynl__(1+index_Y_,:) = reshape(tmp_a_yl_,[n_lm,tmp_n_loading]);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% a_k_Y_2lsq_ynl__: %0.3fs',tmp_t)); end;
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;figmed;
subplot(1,2,1); plot(log(abs(a_k_Y_2lsq_ynl__(:,1+0))),'.'); title(sprintf('2lsq: nloading %d',0));
subplot(1,2,2); plot(log(abs(a_k_Y_2lsq_ynl__(:,1+1))),'.'); title(sprintf('2lsq: nloading %d',1));
drawnow();
end;%if flag_disp;
%%%%;
tmp_SV_loading_0_Ml__ = SV_lsq(n_k_p_r,weight_2d_k_p_r_,n_w_,n_M,R_k_p_wnM__,l_max_,tmp_n_loading,qbp_An,a_k_Y_2lsq_ynl__);
n_lin = 1+12; tmp_eps_ = linspace(-1,+1,n_lin); tmp_eps_scale = 0.01*fnorm(tmp_SV_loading_0_Ml__);
n_test = 6; l2_et_ = zeros(n_lin,n_test);
for ntest=0:n_test-1;
if (flag_verbose); disp(sprintf(' %% ntest %d/%d',ntest,n_test)); end;
tmp_dSV_Ml__ = randn(n_M,tmp_n_loading);
for nlin=0:n_lin-1;
tmp_eps = tmp_eps_scale*tmp_eps_(1+nlin);
tmp_SV_loading_1_Ml__ = tmp_SV_loading_0_Ml__ + tmp_eps*tmp_dSV_Ml__;
dR_k_p_wnM__ = AnSV_all__(n_k_p_r,n_w_,n_M,l_max_,tmp_n_loading,An___,tmp_SV_loading_1_Ml__,a_k_Y_2lsq_ynl__) - R_k_p_wnM__;
l2_et__(1+nlin,1+ntest) = real(sum(innerproduct_p_quad_(n_k_p_r,weight_2d_k_p_r_,n_w_,n_M,dR_k_p_wnM__,dR_k_p_wnM__)));
end;%for nlin=0:n_lin-1;
end;%for ntest=0:n_test-1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
hold on;
plot(0,l2_et__(1+(n_lin-1)/2,:),'ko','MarkerFaceColor','r');
plot(tmp_eps_,l2_et__,'kx'); xlabel('tmp_eps_','Interpreter','none'); ylabel('l2 error');
hold off;
title(sprintf('SV_loading_Ml__'),'Interpreter','none');
drawnow();
end;%if flag_disp;
%%%%;
AnSVa_k_p_wnM__ = AnSV_all__(n_k_p_r,n_w_,n_M,l_max_,tmp_n_loading,An___,tmp_SV_loading_Ml__,a_k_Y_2lsq_ynl__);
n_lin = 1+12; tmp_eps_ = linspace(-1,+1,n_lin); tmp_eps_scale = 0.01*fnorm(a_k_Y_2lsq_ynl__);
n_test = 8; l2_et_ = zeros(n_lin,n_test);
for ntest=0:n_test-1;
if (flag_verbose); disp(sprintf(' %% ntest %d/%d',ntest,n_test)); end;
tmp_da_ynl__ = (randn(n_lm_sum,tmp_n_loading) + i*randn(n_lm_sum,tmp_n_loading));
for nlin=0:n_lin-1;
tmp_eps = tmp_eps_scale*tmp_eps_(1+nlin);
b_k_Y_2lsq_ynl__ = a_k_Y_2lsq_ynl__ + tmp_eps*tmp_da_ynl__;
AnSVb_k_p_wnM__ = AnSV_all__(n_k_p_r,n_w_,n_M,l_max_,tmp_n_loading,An___,tmp_SV_loading_Ml__,b_k_Y_2lsq_ynl__);
AnSVc_k_p_wnM__ = AnSVb_k_p_wnM__ - R_k_p_wnM__;
l2_et__(1+nlin,1+ntest) = real(sum(innerproduct_p_quad_(n_k_p_r,weight_2d_k_p_r_,n_w_,n_M,AnSVc_k_p_wnM__,AnSVc_k_p_wnM__)));
end;%for nlin=0:n_lin-1;
end;%for ntest=0:n_test-1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
hold on;
plot(0,l2_et__(1+(n_lin-1)/2,:),'ko','MarkerFaceColor','r');
plot(tmp_eps_,l2_et__,'kx'); xlabel('tmp_eps_','Interpreter','none'); ylabel('l2 error');
hold off;
title(sprintf('a_k_Y_2lsq_ynl__'),'Interpreter','none');
drawnow();
end;%if flag_disp;
%%%%%%%%;
if (flag_verbose>-1); disp(sprintf(' %% flag_check==1; returning')); end;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_check;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ~flag_loading_skip_loading;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_loading_svd_vs_iterate==0) | (flag_loading_svd_vs_iterate==1) ;
%%%%%%%%;
% Initialize residual clustering. ;
% For each nk_p_r, we load the residual R_wM__ for that shell, ;
% and weight it by CTF_wM_ to form R_CTF_wM__. ;
% We then map the residuals for each image (individually) onto the quadrature-nodes, ;
% assuming the other images are all zero (hence normalizing by CTF2_q_). ;
% Finally, we integrate to obtain an (image- and residual-specific) volume H_ynM__(:,1+nM). ;
% (Note that we account for quadrature-weights in Ylm_w_yq__ below). ;
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
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% H_ynM__: %0.3fs',tmp_t)); end;
%%%%;
% The norm-squared of each of the image-specific H-volumes is denoted by HH_M_(1+nM). ;
%%%%;
HH_M_ = sum(bsxfun(@times,abs(H_ynM__).^2,reshape(weight_3d_yn_,[n_lm_sum,1])),1);
%%%%%%%%;
end;%if (flag_loading_svd_vs_iterate==0) | (flag_loading_svd_vs_iterate==1) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_loading_svd_vs_iterate==-1; %<-- sequential iteration, designed for a single loading. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Iterate to determine a single loading. ;
% Our model is: ;
% R_k_p_wnM__(:,1+nM) * SV_loading_Ml__(1+nM,1+nloading) = ... ;
% An_ ( U_ynl__(:,1+nloading) ), ;
% meaning that, for a single loading (constrained to be norm 1), ;
% we try and choose SV_loading_Ml__ so that the output of the evaluation operator ;
% (applied to U_ynl__(:,1+nloading)) aligns with the respective residuals R_k_p_wnM__. ;
%%%%%%%%;
SV_loading_Ml__ = zeros(n_M,n_loading);
U_ynl__ = zeros(n_lm_sum,n_loading);
Q_k_p_wnM__ = R_k_p_wnM__;
nloading=0;
tmp_t = tic();
rng(rseed + nloading);
SV_loading_M_ = randn(n_M,1); SV_loading_M_ = SV_loading_M_/max(1e-12,fnorm(SV_loading_M_));
U_yn_ = zeros(n_lm_sum,1);
flag_continue=1; niteration=0;
while flag_continue;
U_yn_ = qbp_Ai(bsxfun(@times,Q_k_p_wnM__,reshape(SV_loading_M_,[1,n_M])));
P_k_p_wnM__ = qbp_An(U_yn_);
for nM=0:n_M-1;
tmp_PP = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,P_k_p_wnM__(:,1+nM),P_k_p_wnM__(:,1+nM));
tmp_QQ = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,Q_k_p_wnM__(:,1+nM),Q_k_p_wnM__(:,1+nM));
tmp_PQ = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,P_k_p_wnM__(:,1+nM),Q_k_p_wnM__(:,1+nM));
SV_loading_new_M_(1+nM) = real(tmp_PQ/max(1e-12,tmp_QQ));
end;%for nM=0:n_M-1;
SV_loading_new_M_ = SV_loading_new_M_/max(1e-12,fnorm(SV_loading_new_M_));
SV_loading_upd_M_ = SV_loading_new_M_ - SV_loading_M_;
SV_loading_upd_l2 = fnorm(SV_loading_upd_M_).^2;
SV_loading_upd_f = sqrt(SV_loading_upd_l2)/max(1e-12,fnorm(SV_loading_M_));
if (flag_verbose); disp(sprintf(' %% niteration %d/%d, SV_loading_upd_f %0.6f',niteration,n_loading_iteration,SV_loading_upd_f)); end;
niteration = niteration+1;
flag_continue = (niteration< n_loading_iteration) & (SV_loading_upd_f> tolerance_master);
SV_loading_M_ = SV_loading_new_M_;
end;%while flag_continue;
SV_loading_Ml__(:,1+nloading) = SV_loading_M_;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% nloading %d/%d, SV_loading_M_: %0.3fs',nloading,n_loading,tmp_t)); end;
%%%%;
S_l_ = zeros(n_loading,1);
V_Ml__ = SV_loading_Ml__;
for nloading=0:n_loading-1;
S_l_(1+nloading) = fnorm(V_Ml__(:,1+nloading));
V_Ml__(:,1+nloading) = V_Ml__(:,1+nloading)/max(1e-12,S_l_(1+nloading));
end;%for nloading=0:n_loading-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_loading_svd_vs_iterate==-1; %<-- sequential iteration. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_loading_svd_vs_iterate==0; %<-- simultaneous iteration. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Iterate to cluster residuals. ;
% Our model is: ;
% R_k_p_wnM__(:,1+nM) = ... ;
% An_ ( \sum_{nloading} U_ynl__(:,1+nloading) * SV_loading_Ml__(1+nM,1+nloading) ), ;
% Or equivalently that: ;
% H_ynM__(:,1+nM) \approx ... ;
% \sum_{nloading} U_ynl__(:,1+nloading) * SV_loading_Ml__(1+nM,1+nloading). ;
%%%%%%%%;
% Our iteration proceeds in three steps: ;
% Step 1: ;
% We define: ;
% U_ynl__(:,1+nloading) = \sum_{nM} H_ynM__(:,1+nM) * SV_loading_Ml__(1+nM,1+nloading). ;
% Step 2: ;
% We ensure that each of the U_ynl__(:,1+nloading) are of unit-norm ;
% and orthonormal to one another. ;
% Step 3 (crude): ;
% We approximate SV_loading_Ml__(1+nM,1+nloading) by taking the dot-product: ;
% SV_loading_Ml__(1+nM,1+nloading) \approx ... ;
% dot( H_ynM__(:,1+nM) , U_ynl__(:,1+nloading) ). ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_t = tic();
rng(rseed);
SV_loading_Ml__ = randn(n_M,n_loading);
U_ynl__ = randn(n_lm_sum,n_loading); 
U_l2_l_ = zeros(n_loading,1);
for nloading=0:n_loading-1;
U_l2_l_(1+nloading) = sqrt(ctranspose(U_ynl__(:,1+nloading))*(weight_3d_yn_.*U_ynl__(:,1+nloading))); 
U_ynl__(:,1+nloading) = U_ynl__(:,1+nloading)/max(1e-12,U_l2_l_(1+nloading));
end;%for nloading=0:n_loading-1; 
niteration=0; flag_continue=1;
while flag_continue;
%%%%;
if (niteration==0);
% do nothing;
end;%if (niteration==0);
if (niteration>0);
%%%%;
if flag_U_lsq_vs_dot==0;
U_ynl__ = H_ynM__*SV_loading_Ml__;
end;%end;%if flag_U_lsq_vs_dot==0;
%%%%;
if flag_U_lsq_vs_dot==1;
%%;
% Solving LHS = RHS*U, ;
% The LHS is: ;
% \sum_{j} V_{h,j} \ctranspose(A_{j}) \odot R_{j} ; %<-- one row for each h. ;
% Thee matrix on the RHS is: ;
% \sum_{j} V_{h,j} \ctranspose(A_{j}) \odot A_{j} V_{h^{\prime},j} ; %<-- one row for each h and one column for each h^{\prime}. ;
%%;
U_ynl__ = zeros(n_lm_sum,n_loading);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
An__ = An___{1+nk_p_r}; At__ = At___{1+nk_p_r}; AtAn__ = AtAn___{1+nk_p_r};
tmp_LHS_yl__ = zeros(n_lm,n_loading);
tmp_R_k_p_wM_ = reshape(R_k_p_wnM__(1+index_nw_,:),[n_w*n_M,1]);
tmp_LHS_yl_ = VSAt__(n_lm,n_loading,n_w,n_M,At__,SV_loading_Ml__,tmp_R_k_p_wM_);
tmp_RHS_ylyl__ = @(a_k_Y_yl_) VSAtAnSV__(n_lm,n_loading,n_w,n_M,At__,An__,SV_loading_Ml__,a_k_Y_yl_);
tmp_pgc_U_n_iteration = max(1,min(6,floor(niteration/2)));
[tmp_U_yl_,~] = pcg(tmp_RHS_ylyl__,tmp_LHS_yl_,[],tmp_pgc_U_n_iteration);
U_ynl__(1+index_Y_,:) = reshape(tmp_U_yl_,[n_lm,n_loading]);
end;%for nk_p_r=0:n_k_p_r-1;
end;%if flag_U_lsq_vs_dot==1;
%%%%;
% now orthonormalize. ;
%%%%;
for nloading_0=0:n_loading-1;
U_l2_l_(1+nloading_0) = sqrt(ctranspose(U_ynl__(:,1+nloading_0))*(weight_3d_yn_.*U_ynl__(:,1+nloading_0))); 
U_ynl__(:,1+nloading_0) = U_ynl__(:,1+nloading_0)/max(1e-12,U_l2_l_(1+nloading_0));
for nloading_1=nloading_0+1:n_loading-1;
U_ynl__(:,1+nloading_1) = U_ynl__(:,1+nloading_1) - (ctranspose(U_ynl__(:,1+nloading_0))*(weight_3d_yn_.*U_ynl__(:,1+nloading_1)))*U_ynl__(:,1+nloading_0);
end;%for nloading_1=nloading_0+1:n_loading-1;
end;%for nloading_0=0:n_loading-1;
end;%if (niteration>0);
%%%%;
if flag_SV_lsq_vs_dot==0;
SV_loading_new_Ml__ = zeros(n_M,n_loading);
for nM=0:n_M-1;
HH = HH_M_(1+nM);
for nloading=0:n_loading-1;
HU = ctranspose(H_ynM__(:,1+nM))*(weight_3d_yn_.*U_ynl__(:,1+nloading));
SV_loading_new_Ml__(1+nM,1+nloading) = real(HU)/max(1e-12,real(HH));
end;%for nloading=0:n_loading-1;
end;%for nM=0:n_M-1;
end;%if flag_SV_lsq_vs_dot==0;
%%;
if flag_SV_lsq_vs_dot==1;
SV_loading_new_Ml__ = SV_lsq(n_k_p_r,weight_2d_k_p_r_,n_w_,n_M,R_k_p_wnM__,l_max_,n_loading,qbp_An,U_ynl__);
end;%if flag_SV_lsq_vs_dot==1;
%%%%;
SV_loading_upd_Ml__ = SV_loading_new_Ml__ - SV_loading_Ml__;
SV_loading_upd_rn = fnorm(SV_loading_upd_Ml__)/max(1e-12,fnorm(SV_loading_new_Ml__));
if (flag_verbose); disp(sprintf(' %% niteration %d/%d, SV_loading_upd_rn: %0.2f',niteration,n_loading_iteration,SV_loading_upd_rn)); end;
SV_loading_Ml__ = SV_loading_new_Ml__;
niteration = niteration+1;
flag_continue = (niteration< n_loading_iteration) & (SV_loading_upd_rn> tolerance_master);
%%%%;
end;%while flag_continue;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% SV_loading_Ml__: %0.3fs',tmp_t)); end;
%%%%;
S_l_ = zeros(n_loading,1);
V_Ml__ = SV_loading_Ml__;
for nloading=0:n_loading-1;
S_l_(1+nloading) = fnorm(V_Ml__(:,1+nloading));
V_Ml__(:,1+nloading) = V_Ml__(:,1+nloading)/max(1e-12,S_l_(1+nloading));
end;%for nloading=0:n_loading-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%flag_loading_svd_vs_iterate==0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_loading_svd_vs_iterate==1; %<-- svd. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Call svd to cluster residuals. ;
% Our model is: ;
% R_k_p_wnM__(:,1+nM) = ... ;
% An_ ( \sum_{nloading} U_ynl__(:,1+nloading) * SV_loading_Ml__(1+nM,1+nloading) ), ;
% Or equivalently that: ;
% H_ynM__(:,1+nM) \approx ... ;
% \sum_{nloading} U_ynl__(:,1+nloading) * S_l_(1+nloading) * V_Ml__(1+nM,1+nloading). ;
% Note that we must reweight the entries so that vector-orthonormality ;
% corresponds to l2-orthonormality. ;
%%%%;
% e.g., this is a direct svd: ;
%[U_ynl__,S_l__,V_Ml__] = svds(H_ynM__),n_loading); S_l_ = diag(S_l__);
%SV_loading_Ml__ = real(V_Ml__*diag(S_l_));
%%%%;
% instead, we will perform: ;
tmp_t = tic();
[U_ynl__,S_l__,V_Ml__] = svds(bsxfun(@times,H_ynM__,reshape(sqrt(weight_3d_yn_),[n_lm_sum,1])),n_loading); S_l_ = diag(S_l__);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% SV_loading_Ml__: %0.3fs',tmp_t)); end;
SV_loading_Ml__ = real(ctranspose(H_ynM__)*U_ynl__);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_loading_svd_vs_iterate==1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_loading_svd_vs_iterate==2; %<-- rank-by-rank svd. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Our model is: ;
% R_k_p_wnM__(:,1+nM) = ... ;
% An_ ( \sum_{nloading} U_ynl__(:,1+nloading) \times SV_loading_Ml__(1+nM,1+nloading) ), ;
% which we will write as: ;
% R_{j} = \sum_{h} A_{j} U_{h} V_{h,j} ;
% with : ;
% j = image-index. ;
% h = rank. ;
% R_{j} = residual for image-j. ;
% A_{j} = evaluation-operator for image-j. ;
% U_{h} = volume for rank-h. ;
% V_{h,j} = loading for image-j and rank-h. ;
%%%%;
% Using this notation, we can immediately see that: ;
% if the U_{h} are given, then the V_{h,j} for each fixed j solve a lsq: ;
% ctranspose( An_{j} U_{h} ) \times R_{j} = \sum_{h^{\prime}} ctranspose( An_{j} U_{h} ) \times An_{j} U_{h^{\prime}} \times V_{h^{\prime} , j} ; %<-- this holds for each h. ;
% Here the relevant quantities are the LHS: ;
% [ dot(An_{j}U_{h},R_{j}) ] ; %<-- one row for each h. ;
% And the matrix on the RHS. ;
% [ dot(An_{j}U_{h},An_{j}U_{h^{\prime}}) ] ; %<-- one row for each h and one column for each h^{\prime}. ;
%%%%;
% Similarly, if the V_{h,j} are given then the U_{h} solve a lsq: ;
% \sum_{j} V_{h,j} \ctranspose(A_{j}) \odot R_{j} = \sum_{j,h^{\prime}} V_{h,j} \ctranspose(A_{j}) \odot A_{j} V_{h^{\prime},j} U_{h^[\prime}} ; %<-- this holds for each h. ;
% So the LHS is: ;
% \sum_{j} V_{h,j} \ctranspose(A_{j}) \odot R_{j} ; %<-- one row for each h. ;
% And the matrix on the RHS is: ;
% \sum_{j} V_{h,j} \ctranspose(A_{j}) \odot A_{j} V_{h^{\prime},j} ; %<-- one row for each h and one column for each h^{\prime}. ;
% Note that, just as for a standard volume reconstruction, this lsq for U_{h} is blocked by radius k. ;
%%%%;
% While one could in principle solve the lsq for U_{h} for multiple ranks simultaneously (e.g., using conjugate-gradient), ;
% a faster (and more robust) alternative is to solve for the ranks sequentially. ;
% In this case we can use the approximate-inverse of qbp. ;
%%%%;
% For a single rank (h==0), we can write: ;
% R_{j} = A_{j} U V_{j}, ;
% So given V_{j}, U solves a reconstruction problem: ;
% R_{j}/V_{j} = A_{j} U. ;
% And given U, V_{j} is simply a dot-product: ;
% dot(R_{j},A_{j}U) / |A_{j}U|. ;
%%%%;
% Once we have converged to a particular U and V_{j}, ;
% we can form a new (corrected) residual as: ;
% R_{j} := R_{j} - A_{j} U V_{j} ;
% and then find the next U. ;
% Along the way we can ensure that the U are a well conditioned (e.g., orthonormal) basis. ;
%%%%;
SV_loading_Ml__ = zeros(n_M,n_loading);
U_ynl__ = zeros(n_lm_sum,n_loading);
Q_k_p_wnM__ = R_k_p_wnM__;
for nloading_0=0:n_loading-1;
tmp_t = tic();
rng(rseed + nloading_0);
SV_loading_M_ = randn(n_M,1);
SV_loading_l2 = fnorm(SV_loading_M_);
U_yn_ = zeros(n_lm_sum,1);
flag_continue=1; niteration=0;
while flag_continue;
%%%%;
% define U_yn_ via solving approximate inverse. ;
%%%%;
U_yn_ = qbp_Ai(bsxfun(@times,Q_k_p_wnM__,reshape(SV_loading_M_./(1e-12 + SV_loading_M_.^2),[1,n_M]))); %<-- pseudo-division. ;
%%%%;
% orthogonalize. ;
%%%%;
for nloading_1=0:nloading_0-1;
pre_U_yn_ = U_ynl__(:,1+nloading_1);
U_yn_ = U_yn_ - (ctranspose(pre_U_yn_)*(weight_3d_yn_.*U_yn_))*pre_U_yn_;
end;%for nloading_1=0:nloading_0-1;
%%%%;
% normalize. ;
%%%%;
U_l2 = sqrt(ctranspose(U_yn_)*(weight_3d_yn_.*U_yn_)); 
U_yn_ = U_yn_/max(1e-12,U_l2);
%%%%;
% define new loadings. ;
%%%%;
P_k_p_wnM__ = qbp_An(U_yn_);
Q_l2_M_ = zeros(n_M,1);
for nM=0:n_M-1;
tmp_PP = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,P_k_p_wnM__(:,1+nM),P_k_p_wnM__(:,1+nM));
tmp_QQ = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,Q_k_p_wnM__(:,1+nM),Q_k_p_wnM__(:,1+nM));
tmp_PQ = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,P_k_p_wnM__(:,1+nM),Q_k_p_wnM__(:,1+nM));
SV_loading_new_M_(1+nM) = real(tmp_PQ/max(1e-12,tmp_PP));
Q_l2_M_(1+nM) = tmp_QQ;
end;%for nM=0:n_M-1;
SV_loading_new_l2 = fnorm(SV_loading_new_M_);
%%%%;
% determine flag_continue. ;
%%%%;
SV_loading_upd_M_ = SV_loading_new_M_ - SV_loading_M_;
SV_loading_upd_l2 = fnorm(SV_loading_upd_M_);
SV_loading_upd_f = SV_loading_upd_l2/max(1e-12,SV_loading_l2);
if (flag_verbose); disp(sprintf(' %% nloading_0 %d/%d niteration %d/%d, SV_loading_upd_f %0.6f',nloading_0,n_loading,niteration,n_loading_iteration,SV_loading_upd_f)); end;
niteration = niteration+1;
flag_continue = (niteration< n_loading_iteration) & (SV_loading_upd_f> tolerance_master);
SV_loading_M_ = SV_loading_new_M_; SV_loading_l2 = SV_loading_new_l2;
end;%while flag_continue;
SV_loading_Ml__(:,1+nloading_0) = SV_loading_M_;
U_ynl__(:,1+nloading_0) = U_yn_;
Q_k_p_wnM__ = Q_k_p_wnM__ - bsxfun(@times,P_k_p_wnM__,reshape(SV_loading_M_,[1,n_M]));
%%%%;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% nloading_0 %d/%d, Q_l2 %0.6f, SV_loading_M_: %0.3fs',nloading_0,n_loading,sum(Q_l2_M_),tmp_t)); end;
end;%for nloading_0=0:n_loading-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_loading_svd_vs_iterate==2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now use the loadings and volumes to approximate the residuals. ;
%%%%%%%%;
tmp_t = tic();
AUSV_wnM__ = zeros(n_w_sum,n_M);
for nloading=0:n_loading-1;
tmp_AU_wnM__ = zeros(n_w_sum,n_M);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
An__ = An___{1+nk_p_r}; At__ = At___{1+nk_p_r};
tmp_AU_wnM__(1+index_nw_,:) = reshape(An__(U_ynl__(1+index_Y_,1+nloading)),[n_w,n_M]);
end;%for nk_p_r=0:n_k_p_r-1;
AUSV_wnM__ = AUSV_wnM__ + tmp_AU_wnM__*diag(SV_loading_Ml__(:,1+nloading));
end;%for nloading=0:n_loading-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% AUSV_wnM__: %0.3fs',tmp_t)); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ~flag_loading_skip_loading;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose>1); disp(sprintf(' %% [finished get_loading_qbp_2]')); end;

%function a_k_Y_0qbp_yn_ = qbp_Ai(n_k_p_r,n_w_,n_M,P_k_p_wnM__,CTF_wMn__,CTF2_qn__,quad_from_data_qwM__,Ylm_w_yq__);
%n_w = max(n_w_); assert(min(n_w_)==n_w);
%n_lm_sum = n_k_p_r*size(Ylm_w_yq__,1);
%a_k_Y_0qbp_yn_ = zeros(n_lm_sum,1);
%P_wMn___ = permute(reshape(P_k_p_wnM__,[n_w,n_k_p_r,n_M]),[1,3,2]);
%P_CTF_wMn___ = bsxfun(@times,P_wMn___,reshape(CTF_wMn__,[n_w,n_M,n_k_p_r]));
%quad_from_data_P_CTF_normalized_qn__ = (quad_from_data_qwM__*reshape(P_CTF_wMn___,[n_w*n_M,n_k_p_r]))./max(1e-12,CTF2_qn__);
%a_k_Y_0qbp_yn_ = reshape(conj(Ylm_w_yq__)*quad_from_data_P_CTF_normalized_qn__,[n_lm_sum,1]);

function a_k_Y_yl_ = VSAt__(n_lm,n_loading,n_w,n_M,At__,SV_loading_Ml__,M_k_p_wM_);
a_k_Y_yl__ = zeros(n_lm,n_loading);
M_k_p_wM__ = reshape(M_k_p_wM_,[n_w,n_M]);
for nloading=0:n_loading-1;
tmp_SV_loading_M_ = SV_loading_Ml__(:,1+nloading);
a_k_Y_yl__(:,1+nloading) = At__(reshape(bsxfun(@times,M_k_p_wM__,reshape(tmp_SV_loading_M_,[1,n_M])),[n_w*n_M,1]));
end;%for nloading=0:n_loading-1;
a_k_Y_yl_ = reshape(a_k_Y_yl__,[n_lm*n_loading,1]);

function M_k_p_wM_ = AnSV__(n_w,n_M,n_lm,n_loading,An__,SV_loading_Ml__,a_k_Y_yl_);
a_k_Y_yl__ = reshape(a_k_Y_yl_,[n_lm,n_loading]);
M_k_p_wM__ = zeros(n_w,n_M);
for nloading=0:n_loading-1;
tmp_SV_loading_M_ = SV_loading_Ml__(:,1+nloading);
M_k_p_wM__ = M_k_p_wM__ + bsxfun(@times,reshape(An__(a_k_Y_yl__(:,1+nloading)),[n_w,n_M]),reshape(tmp_SV_loading_M_,[1,n_M]));
end;%for nloading=0:n_loading-1;
M_k_p_wM_ = reshape(M_k_p_wM__,[n_w*n_M,1]);

function a_k_Y_out_yl_ = VSAtAnSV__(n_lm,n_loading,n_w,n_M,At__,An__,SV_loading_Ml__,a_k_Y_0in_yl_);
M_k_p_wM_ = AnSV__(n_w,n_M,n_lm,n_loading,An__,SV_loading_Ml__,a_k_Y_0in_yl_);
a_k_Y_out_yl_ = VSAt__(n_lm,n_loading,n_w,n_M,At__,SV_loading_Ml__,M_k_p_wM_);

function M_k_p_wnM__ = AnSV_all__(n_k_p_r,n_w_,n_M,l_max_,n_loading,An___,SV_loading_Ml__,a_k_Y_ynl__);
%%%%;
n_w_ = reshape(n_w_(1:n_k_p_r),[n_k_p_r,1]);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
l_max_ = reshape(l_max_(1:n_k_p_r),[n_k_p_r,1]);
l_max_max = max(l_max_);
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
%%%%;
M_k_p_wnM__ = zeros(n_w_sum,n_M);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
An__ = An___{1+nk_p_r};
a_k_Y_yl_ = reshape(a_k_Y_ynl__(1+index_Y_,:),[n_lm*n_loading,1]);
M_k_p_wM_ = AnSV__(n_w,n_M,n_lm,n_loading,An__,SV_loading_Ml__,a_k_Y_yl_);
M_k_p_wnM__(1+index_nw_,:) = reshape(M_k_p_wM_,[n_w,n_M]);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;

function SV_loading_Ml__ = SV_lsq(n_k_p_r,weight_2d_k_p_r_,n_w_,n_M,R_k_p_wnM__,l_max_,n_loading,qbp_An,U_k_Y_ynl__);
%%%%;
n_w_ = reshape(n_w_(1:n_k_p_r),[n_k_p_r,1]);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
l_max_ = reshape(l_max_(1:n_k_p_r),[n_k_p_r,1]);
l_max_max = max(l_max_);
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
%%%%;
P_k_p_wnMl___ = zeros(n_w_sum,n_M,n_loading);
for nloading=0:n_loading-1;
P_k_p_wnMl___(:,:,1+nloading) = qbp_An(U_k_Y_ynl__(:,1+nloading));
end;%for nloading=0:n_loading-1;
PP_llM___ = zeros(n_loading,n_loading,n_M);
PR_lM__ = zeros(n_loading,n_M);
RR_M_ = innerproduct_p_quad_(n_k_p_r,weight_2d_k_p_r_,n_w_,n_M,R_k_p_wnM__,R_k_p_wnM__);
for nloading_0=0:n_loading-1;
P_k_p_wnM0__ = P_k_p_wnMl___(:,:,1+nloading_0);
PR_lM__(1+nloading_0,:) = innerproduct_p_quad_(n_k_p_r,weight_2d_k_p_r_,n_w_,n_M,P_k_p_wnM0__,R_k_p_wnM__);
for nloading_1=0:n_loading-1;
P_k_p_wnM1__ = P_k_p_wnMl___(:,:,1+nloading_1);
PP_llM___(1+nloading_0,1+nloading_1,:) = innerproduct_p_quad_(n_k_p_r,weight_2d_k_p_r_,n_w_,n_M,P_k_p_wnM0__,P_k_p_wnM1__);
end;%for nloading_1=0:n_loading-1;
end;%for nloading_0=0:n_loading-1;
for nM=0:n_M-1;
Atb_l_ = reshape(PR_lM__(:,1+nM),[n_loading,1]);
AtAn_ll__ = reshape(PP_llM___(:,:,1+nM),[n_loading,n_loading]);
x_l_ = pinv(AtAn_ll__,1e-6)*Atb_l_;
SV_loading_Ml__(1+nM,:) = real(x_l_);
end;%for nM=0:n_M-1;
%%%%;
