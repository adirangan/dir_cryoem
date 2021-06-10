function M_loading__ = get_loading_lsq(n_residual_loading,n_residual_iteration,n_order,n_k_p_r,weight_k_p_r_,weight_2d_k_p_r_,l_max_,n_w_,n_M,M_k_p__,euler_polar_a_,euler_azimu_b_,euler_gamma_z_,label_);
verbose=1;
%%%%%%%%;
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); 
%%%%%%%%;
% Solve lsq with estimated viewing angles. ;
%%%%%%%%;
n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
n_lm_ = (1+l_max_).^2; n_lm_sum = sum(n_lm_);
tmp_t = tic();
a_k_Y_0lsq_ = cg_lsq_3(n_order,n_k_p_r,l_max_,n_w_,n_M,M_k_p__,1,ones(n_w_sum,1),euler_polar_a_,euler_azimu_b_,euler_gamma_z_);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_k_Y_0lsq_: %0.3fs',tmp_t)); end;
%%%%%%%%;
% calculate residual from lsq. ;
%%%%%%%%;
tmp_t = tic();
CTF_index_ = 1; CTF_k_p__ = ones(n_w_sum,1); %<-- no CTF necessary. ;
S_k_p__ = zeros(n_w_sum,n_M);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
Y_ij_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
n_w = n_w_(1+nk_p_r);
M_ij_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
if (isempty(CTF_index_) | numel(CTF_index_)==1);
CTF_k_p_ = reshape(repmat(CTF_k_p__(1+M_ij_),[1,n_M]),[n_w*n_M,1]);
 else;
CTF_k_p_ = reshape(CTF_k_p__(1+M_ij_,CTF_index_(1:n_M)),[n_w*n_M,1]);
end;%if (isempty(CTF_index_) | numel(CTF_index_)==1);
[k_p_polar_a__,k_p_azimu_b__] = cg_rhs_1(n_M,n_w,euler_polar_a_,euler_azimu_b_,+euler_gamma_z_);
n_polar_a = ceil(n_w/2);
n_azimu_b = max(1+2*l_max,2*n_polar_a);
[legendre_evaluate_ljm___,legendre_evaluate_mlj___,expil__,expi__] = legendre_evaluate_ljm___0(l_max,cos(linspace(0,pi,n_polar_a)),n_azimu_b);
tensor_to_scatter__ = cg_interpolate_n_1(n_order,n_polar_a,n_azimu_b,n_w*n_M,k_p_polar_a__(:),k_p_azimu_b__(:));
scatter_to_tensor__ = ctranspose(tensor_to_scatter__); %<-- this conjugation is not necessary, since the matrix should be real. ;
An__ = @(a_k_Y_) CTF_k_p_.*(tensor_to_scatter__*reshape(cg_evaluate_n_1(l_max,convert_spharm_to_spharm__0(l_max,a_k_Y_),n_polar_a,n_azimu_b,legendre_evaluate_ljm___),[n_polar_a*n_azimu_b,1]));
At__ = @(a_k_X_) convert_spharm__to_spharm_0(l_max,cg_evaluate_t_1(n_polar_a,n_azimu_b,reshape(scatter_to_tensor__*(conj(CTF_k_p_).*a_k_X_),[n_polar_a,n_azimu_b]),l_max,legendre_evaluate_mlj___,expil__,expi__));
AtAn__ = @(a_k_Y_) At__(An__(a_k_Y_));
%[a_k_Y_(1+Y_ij_),~] = pcg(AtAn__,At__(reshape(M_k_p__(1+M_ij_,:),[n_w*n_M,1])));
S_k_p__(1+M_ij_,:) = reshape(An__(a_k_Y_0lsq_(1+Y_ij_)),[n_w_(1+nk_p_r),n_M]);
end;%for nk_p_r=0:n_k_p_r-1;
R_k_p__ = M_k_p__ - S_k_p__;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% R_k_p__: %0.3fs',tmp_t)); end;
RR_ = zeros(n_M,1);
for nM=0:n_M-1;
RR_(1+nM) = real(innerproduct_p_quad(n_k_p_r,[],weight_2d_k_p_r_,n_w_,n_w_sum,R_k_p__(:,1+nM),R_k_p__(:,1+nM)));
end;%for nM=0:n_M-1;
%%%%%%%%;
% Initialize residual clustering. ;
%%%%%%%%;
M_loading__ = randn(n_residual_loading,n_M);
for niteration=0:n_residual_iteration-1;
disp(sprintf(' %% niteration %d/%d',niteration,n_residual_iteration));
%%%%;
% define G_. ;
%%%%;
tmp_t = tic();
G_ = zeros(n_lm_sum,n_residual_loading);
for nloading=0:n_residual_loading-1;
G_(:,1+nloading) = cg_lsq_3(n_order,n_k_p_r,l_max_,n_w_,n_M,R_k_p__*sparse(1:n_M,1:n_M,M_loading__(1+nloading,:),n_M,n_M),CTF_index_,CTF_k_p__,euler_polar_a_,euler_azimu_b_,euler_gamma_z_);
end;%for nloading=0:n_residual_loading-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% G_: %0.3fs',tmp_t)); end;
%%%%;
% orthonormalize G_. ;
%%%%;
tmp_t = tic();
for nloading_0=0:n_residual_loading-1;
G_(:,1+nloading_0) = G_(:,1+nloading_0)/fnorm(G_(:,1+nloading_0));
for nloading_1=nloading_0+1:n_residual_loading-1;
G_(:,1+nloading_1) = G_(:,1+nloading_1) - (ctranspose(G_(:,1+nloading_0))*G_(:,1+nloading_1))*G_(:,1+nloading_0);
end;%for nloading_1=nloading_0+1:n_residual_loading-1;
end;%for nloading_0=0:n_residual_loading-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% orthonormalize G_: %0.3fs',tmp_t)); end;
%%%%;
% update loadings. ;
%%%%;
tmp_t = tic();
EG___ = zeros(n_w_sum,n_M,n_residual_loading);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
Y_ij_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
n_w = n_w_(1+nk_p_r);
M_ij_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
n_polar_a = ceil(n_w/2);
n_azimu_b = max(1+2*l_max,2*n_polar_a);
[legendre_evaluate_ljm___,legendre_evaluate_mlj___,expil__,expi__] = legendre_evaluate_ljm___0(l_max,cos(linspace(0,pi,n_polar_a)),n_azimu_b);
for nM=0:n_M-1;
if (isempty(CTF_index_) | numel(CTF_index_)==1);
CTF_k_p_ = reshape(repmat(CTF_k_p__(1+M_ij_),[1,1]),[n_w*1,1]);
 else;
CTF_k_p_ = reshape(CTF_k_p__(1+M_ij_,1+CTF_index_(1+nM)),[n_w*1,1]);
end;%if (isempty(CTF_index_) | numel(CTF_index_)==1);
[k_p_polar_a__,k_p_azimu_b__] = cg_rhs_1(1,n_w,euler_polar_a_(1+nM),euler_azimu_b_(1+nM),+euler_gamma_z_(1+nM));
tensor_to_scatter__ = cg_interpolate_n_1(n_order,n_polar_a,n_azimu_b,n_w*1,k_p_polar_a__(:),k_p_azimu_b__(:));
scatter_to_tensor__ = ctranspose(tensor_to_scatter__); %<-- this conjugation is not necessary, since the matrix should be real. ;
An__ = @(a_k_Y_) CTF_k_p_.*(tensor_to_scatter__*reshape(cg_evaluate_n_1(l_max,convert_spharm_to_spharm__0(l_max,a_k_Y_),n_polar_a,n_azimu_b,legendre_evaluate_ljm___),[n_polar_a*n_azimu_b,1]));
%At__ = @(a_k_X_) convert_spharm__to_spharm_0(l_max,cg_evaluate_t_1(n_polar_a,n_azimu_b,reshape(scatter_to_tensor__*(conj(CTF_k_p_).*a_k_X_),[n_polar_a,n_azimu_b]),l_max,legendre_evaluate_mlj___,expil__,expi__));
for nloading=0:n_residual_loading-1;
EG___(1+M_ij_,1+nM,1+nloading) = An__(G_(1+Y_ij_,1+nloading));
end;%for nloading=0:n_residual_loading-1;
end;%for nM=0:n_M-1;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% EG__: %0.3fs',tmp_t)); end;
%%%%;
tmp_t = tic();
for nM=0:n_M-1;
for nloading=0:n_residual_loading-1;
M_loading__(1+nloading,1+nM) = real(innerproduct_p_quad(n_k_p_r,[],weight_2d_k_p_r_,n_w_,n_w_sum,R_k_p__(:,1+nM),EG___(:,1+nM,1+nloading)))/RR_(1+nM);
end;%for nloading=0:n_residual_loading-1;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% M_loading__: %0.3fs',tmp_t)); end;
%%%%;
flag_plot=1;
if flag_plot;
subplot(4,8,1+niteration);
scatter(M_loading__(1,:),M_loading__(2,:),15,label_,'filled'); axisnotick; drawnow;
end;%if flag_plot;
%%%%;
end;%for niteration=0:n_residual_iteration-1;

