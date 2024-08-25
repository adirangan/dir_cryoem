%%%%%%%%;
% This is a simple test of eig_ddssnll_lanczos_1 using a surrogate hessian. ;
%%%%%%%%;

flag_verbose = 1;
flag_disp = 1; nf=0;

n_k_p_r = [];
k_p_r_ = [];
k_p_r_max = [];
l_max_ = [];
a_k_Y_quad_yk_ = [];
n_k_all = [];
n_k_all_csum_ = [];
k_p_r_all_ = [];
k_p_azimu_b_all_ = [];
k_p_polar_a_all_ = [];
weight_3d_k_all_ = [];
weight_shell_k_ = [];
weight_3d_k_p_r_ = [];
a_k_p_quad_ = [];
n_w_ = [];
weight_2d_k_p_r_ = [];
weight_2d_wk_ = [];
n_S = [];
S_k_p_q2d_wkS__ = [];
viewing_polar_a_S_ = [];
viewing_azimu_b_S_ = [];
viewing_weight_S_ = [];
n_viewing_polar_a = [];
viewing_polar_a_ = [];
n_viewing_azimu_b_ = [];
n_M = [];
M_k_p_wkM__ = [];
n_CTF = [];
index_nCTF_from_nM_ = [];
CTF_k_p_r_kC__ = [];
CTF_k_p_wkC__ = [];
n_eta = [];
index_neta_from_nM_ = [];
eta_k_p_r_ke__ = [];
eta_k_p_wke__ = [];
euler_polar_a_M_ = [];
euler_azimu_b_M_ = [];
euler_gamma_z_M_ = [];
KAPPA = [];
Ylm_uklma___ = [];
k_p_azimu_b_sub_uka__ = [];
k_p_polar_a_sub_uka__ = [];
l_max_uk_ = [];
index_nu_n_k_per_shell_from_nk_p_r_ = [];
index_k_per_shell_uka__ = [];
V_lmm___ = [];
L_lm__ = [];
d0W_betazeta_mlma____ = [];
d1W_betazeta_mlma____ = [];
d2W_betazeta_mlma____ = [];
U_SmallRotation_Delta_ykabc3__ = [];
v_tilde_ykabci__  = [];
w_tilde_ykabc_  = [];
alph_tilde_i_ = [];
beta_tilde_i_ = []; 

%%%%%%%%;
rng(0);
flag_w = 1;
flag_U = 1;
parameter = struct('type','parameter');
parameter.flag_surrogate = 1;
n_k_p_r = 3;
k_p_r_ = transpose(1:n_k_p_r);
k_p_r_max = max(k_p_r_);
l_max_ = 1+k_p_r_;
n_lm_ = (l_max_+1).^2;
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
%%%%;
n_M = 5;
scaling_volumetric = (4*pi)^2 * sqrt(pi/2);
if flag_w==0;
weight_3d_k_p_r_ = ones(n_k_p_r,1);
weight_2d_k_p_r_ = ones(n_k_p_r,1) * scaling_volumetric ;
weight_3d_riesz_k_p_r_ = weight_2d_k_p_r_;
end;%if flag_w==0;
if flag_w==1;
weight_3d_k_p_r_ = 1+rand(n_k_p_r,1);
weight_2d_k_p_r_ = 1+rand(n_k_p_r,1);
weight_3d_riesz_k_p_r_ = weight_2d_k_p_r_;
end;%if flag_w==1;
%%%%;
Y_l_val_ = zeros(n_lm_sum,1);
Y_m_val_ = zeros(n_lm_sum,1);
Y_k_val_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
tmp_l_val_ = zeros(n_lm_(1+nk_p_r),1);
tmp_m_val_ = zeros(n_lm_(1+nk_p_r),1);
na=0; 
for l_val=0:l_max;
for m_val=-l_val:+l_val;
tmp_l_val_(1+na) = l_val;
tmp_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
Y_l_val_(1+tmp_index_) = tmp_l_val_;
Y_m_val_(1+tmp_index_) = tmp_m_val_;
Y_k_val_(1+tmp_index_) = k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
weight_Y_ = zeros(n_lm_sum,1);
weight_3d_riesz_yk_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
weight_Y_(1+tmp_index_) = weight_3d_k_p_r_(1+nk_p_r);
weight_3d_riesz_yk_(1+tmp_index_) = weight_3d_riesz_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
n_ykabc = n_lm_sum + n_M*3;
weight_3d_riesz_ykabc_ = cat(1,weight_3d_riesz_yk_/scaling_volumetric,ones(3*n_M,1));
numerator_root_weight_3d_riesz_ykabc_ = reshape(sqrt(weight_3d_riesz_ykabc_),[n_ykabc,1]);
denomator_root_weight_3d_riesz_ykabc_ = 1./max(1e-12,reshape(sqrt(weight_3d_riesz_ykabc_),[n_ykabc,1]));
%%%%%%%%;

n_k_all_csum_ = zeros(n_k_p_r+1,1);
a_k_Y_quad_yk_ = zeros(n_lm_sum,1);
a_k_p_quad_ = [0];
n_ykabc = n_lm_sum+n_M*3;
n_3 = 3;
%%%%;
if flag_U==1;
[U_SmallRotation_Delta_ykabc3__,~] = qr(randn(n_ykabc,3) + i*randn(n_ykabc,3));
U_SmallRotation_Delta_ykabc3__ = U_SmallRotation_Delta_ykabc3__(:,1:3); %<-- provide externally. ;
for n3=0:n_3-1;
if (flag_verbose>0); disp(sprintf(' %% n3 %d: uu pre = %0.6f',n3,local_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_SmallRotation_Delta_ykabc3__(:,1+n3),U_SmallRotation_Delta_ykabc3__(:,1+n3)))); end;
U_SmallRotation_Delta_ykabc3__(:,1+n3) = local_normalize_f_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_SmallRotation_Delta_ykabc3__(:,1+n3));
if (flag_verbose>0); disp(sprintf(' %% n3 %d: uu pos = %0.6f',n3,local_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_SmallRotation_Delta_ykabc3__(:,1+n3),U_SmallRotation_Delta_ykabc3__(:,1+n3)))); end;
end;%for n3=0:n_3-1;
U_tilde_SmallRotation_Delta_ykabc3__ = bsxfun(@times,numerator_root_weight_3d_riesz_ykabc_,U_SmallRotation_Delta_ykabc3__);
for n3=0:n_3-1;
if (flag_verbose>0); disp(sprintf(' %% n3 %d: uu tilde = %0.6f',n3,local_weightless_f_bar_dot_g_(n_k_p_r,l_max_,n_M,U_tilde_SmallRotation_Delta_ykabc3__(:,1+n3),U_tilde_SmallRotation_Delta_ykabc3__(:,1+n3)))); end;
U_tilde_SmallRotation_Delta_ykabc3__(:,1+n3) = local_weightless_normalize_f_(n_k_p_r,l_max_,n_M,U_tilde_SmallRotation_Delta_ykabc3__(:,1+n3));
end;%for n3=0:n_3-1;
end;%if flag_U==1;
%%%%;
if flag_U==0;
U_SmallRotation_Delta_ykabc3__ = zeros(n_ykabc,[]); %<-- skip entirely. ;
U_tilde_SmallRotation_Delta_ykabc3__ = zeros(n_ykabc,[]); %<-- skip entirely. ;
end;%if flag_U==0;
%%%%;
tmp_w_ = transpose(linspace(0,pi/2,n_ykabc));
surrogate_S_DHD_ = sort(sin(tmp_w_).^2 + 0.250*cos(2*tmp_w_).^2 + 0.125*sin(4*tmp_w_).^2);
[surrogate_V_DHD__,~] = qr(rand(n_ykabc,n_ykabc) + i*rand(n_ykabc,n_ykabc));
surrogate_DHD__ = surrogate_V_DHD__*diag(surrogate_S_DHD_)*ctranspose(surrogate_V_DHD__);
surrogate_DGD__ = surrogate_DHD__ - U_tilde_SmallRotation_Delta_ykabc3__*ctranspose(U_tilde_SmallRotation_Delta_ykabc3__);
[surrogate_V_DGD__,surrogate_S_DGD__] = eigs(surrogate_DGD__,n_ykabc); surrogate_S_DGD_ = diag(real(surrogate_S_DGD__));
parameter.lanczos_n_iteration_max = n_ykabc;
parameter.surrogate_H__ = diag(denomator_root_weight_3d_riesz_ykabc_)*surrogate_DGD__*diag(numerator_root_weight_3d_riesz_ykabc_);
parameter.flag_verbose = 0;
[ ...
 parameter ...
,U_tilde_SmallRotation_Delta_ykabc3__ ...
,v_tilde_ykabci__  ...
,w_tilde_ykabc_  ...
,alph_tilde_i_ ...
,beta_tilde_i_ ... 
] = ...
eig_ddssnll_lanczos_1( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_quad_yk_ ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,weight_3d_k_p_r_ ...
,a_k_p_quad_ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,S_k_p_q2d_wkS__ ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,n_M ...
,M_k_p_wkM__ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,CTF_k_p_wkC__ ...
,n_eta ...
,index_neta_from_nM_ ...
,eta_k_p_r_ke__ ...
,eta_k_p_wke__ ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,KAPPA ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlma____ ...
,d1W_betazeta_mlma____ ...
,d2W_betazeta_mlma____ ...
,U_tilde_SmallRotation_Delta_ykabc3__ ...
,v_tilde_ykabci__  ...
,w_tilde_ykabc_  ...
,alph_tilde_i_ ...
,beta_tilde_i_ ... 
);
%%%%%%%%;

%%%%%%%%;
n_iteration = numel(alph_tilde_i_);
T_tilde__ = real(spdiags([circshift(beta_tilde_i_,-1),alph_tilde_i_,beta_tilde_i_],[-1,0,+1],n_iteration,n_iteration));
lambda_xi__ = -Inf*ones(n_iteration,n_iteration);
for niteration=0:n_iteration-1;
T_tilde_sub__ = T_tilde__(1:1+niteration,1:1+niteration);
lambda_sub_ = eigs(T_tilde_sub__,[],1+niteration);
lambda_xi__(1:1+niteration,1+niteration) = sort(lambda_sub_,'ascend');
end;%for niteration=0:n_iteration-1;
S_x_ = sort(eigs(T_tilde__,[],n_iteration),'ascend');
S_x_min = min(S_x_);
S_x_max = max(S_x_);
%%%%;
vv_ns4___ = zeros(n_iteration,n_iteration,4);
for niteration=0:n_iteration-1;
T_tilde_sub__ = T_tilde__(1:1+niteration,1:1+niteration);
[TV_tilde_sub__,lambda_sub__] = eigs(T_tilde_sub__,[],1+niteration);
lambda_sub_ = diag(lambda_sub__);
[lambda_srt_,ij_srt_] = sort(lambda_sub_,'ascend');
for index_lambda=0:1+niteration-1;
ij_use = ij_srt_(1+index_lambda);
lambda_use = lambda_srt_(1+index_lambda);
TV_tilde_eig_ = TV_tilde_sub__(:,ij_use);
v_tilde_eig_ykabc_ = v_tilde_ykabci__(:,1:1+niteration)*TV_tilde_eig_;
[v_tilde_eig_dvol_yk_,v_tilde_eig_polar_a_M_use_,v_tilde_eig_azimu_b_M_use_,v_tilde_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M,v_tilde_eig_ykabc_);
[tmp_vv,tmp_vv_dvol,tmp_vv_a,tmp_vv_b,tmp_vv_c] = local_weightless_f_bar_dot_g_(n_k_p_r,l_max_,n_M,v_tilde_eig_ykabc_,v_tilde_eig_ykabc_);
str_vv = sprintf('tmp_vv %0.2f,tmp_vv_dvol %0.2f,tmp_vv_a %0.2f,tmp_vv_b %0.2f,tmp_vv_c %0.2f',tmp_vv,tmp_vv_dvol,tmp_vv_a,tmp_vv_b,tmp_vv_c);
if (flag_verbose>1); disp(sprintf(' %% %s',str_vv)); end;
vv_ns4___(1+niteration,1+index_lambda,:) = [tmp_vv_dvol;tmp_vv_a;tmp_vv_b;tmp_vv_c];
end;%for index_lambda=0:1+niteration-1;
end;%for niteration=0:n_iteration-1;
%%%%%%%%;

%%%%%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,768]);fig81s;
str_title = 'test_eig_ddssnll_lanczos_surrogate';
markersize_use = 8;
linewidth_sml = 0.5;
linewidth_big = 2;
%%%%;
subplot(1,5,[1,2,3]);
hold on;
plot(repmat([0;n_iteration],[1,n_iteration]),repmat(reshape(S_x_,[1,n_iteration]),[2,1]),'-','Color',0.85*[1,1,1],'LineWidth',linewidth_sml);
ni_xi__ = repmat([1:n_iteration],[n_iteration,1]);
tmp_index_ = efind(isfinite(lambda_xi__));
plot(ni_xi__(1+tmp_index_),lambda_xi__(1+tmp_index_),'r.','MarkerSize',markersize_use);
plot((n_iteration+0.5)*ones(n_ykabc,1),surrogate_S_DHD_,'ks','MarkerFaceColor','c');
hold off;
xlabel('iteration'); ylabel('sigma');
xlim([0,1+n_iteration]);
ylim([S_x_min-0.25,S_x_max+0.25]);
title(str_title,'Interpreter','none');
%%%%;
subplot(1,5,[4]);
plot(sort(surrogate_S_DHD_),sort(lambda_xi__(:,end)) - sort(surrogate_S_DHD_),'o','MarkerFaceColor','g')
xlabel('true eigenvalue');
xlim([S_x_min-0.25,S_x_max+0.25]);
ylabel('lamda difference');
grid on;
title('lanczos vs H');
%%%%;
subplot(1,5,[5]);
plot(sort(surrogate_S_DGD_),sort(lambda_xi__(:,end)) - sort(surrogate_S_DGD_),'o','MarkerFaceColor','g')
xlabel('true eigenvalue');
xlim([S_x_min-0.25,S_x_max+0.25]);
ylabel('lamda difference');
grid on;
title('lanczos vs G');
%%%%%%%%;

