%%%%%%%%;
% Sets up a simple volume (in spherical-harmonic-coordinates), ;
% then tests pm_template_2.m. ;
%%%%%%%%;

platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;

flag_verbose = 1;
flag_recalc = 0;
flag_replot = 0;
flag_center = 1;
flag_invert = 0;
tolerance_master = 1e-2;
nf=0;

dir_base = sprintf('/%s/rangan/dir_cryoem/dir_rangan_playroom',string_root);
dir_jpg = sprintf('%s/dir_jpg',dir_base);
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

%%%%%%%%;
% Define spatial grid. ;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = half_diameter_x_c;
n_x_c = 64;
x_c_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
x_c_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
x_c_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
[x_c_0___,x_c_1___,x_c_2___] = ndgrid(x_c_0_,x_c_1_,x_c_2_); n_xxx_c = n_x_c^3;
xxx_c_weight = (2*x_p_r_max/n_x_c)^3;

%%%%%%%%;
% Now set up k-quadrature on sphere. ;
%%%%%%%%;
k_p_r_max = 48.0/(2*pi); k_eq_d = 1.0/(2*pi); str_T_vs_L = 'L';
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_T_vs_L ...
);
%%%%%%%%;

%%%%%%%%;
% Now set up polar-quadrature-weights on disk. ;
%%%%%%%%;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
n_w_max = 2*(l_max_max+1); n_w_0in_ = n_w_max*ones(n_k_p_r,1);
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_0in_ ...
);
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;

%%%%%%%%;
% Random templates. ;
%%%%%%%%;
n_S = 2370;
S_k_p_wkS__ = rand(n_w_sum,n_S) + i*rand(n_w_sum,n_S);
S_k_q_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_);
S_k_q_wkS__(:,1+nS) = S_k_q_wk_;
end;%for nS=0:n_S-1;
S_k_q_wSk___ = permute(reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
%%%%%%%%;
n_M = 24*2; %n_M = 1024;
M_k_p_wkM__ = rand(n_w_sum,n_M) + i*rand(n_w_sum,n_M);
M_k_q_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
M_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_);
M_k_q_wkM__(:,1+nM) = M_k_q_wk_;
end;%for nM=0:n_M-1;
M_k_q_wMk___ = permute(reshape(M_k_q_wkM__,[n_w_max,n_k_p_r,n_M]),[1,3,2]);
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% n_S %.4d n_M %.4d',n_S,n_M)); end;

%%%%%%%%;
% Arbitrary CTF. ;
%%%%%%%%;
n_CTF = 1;
%CTF_alpha_C_ = transpose(linspace(0.5,1.0,n_CTF)); %<-- changes sign. ;
CTF_alpha_C_ = transpose(linspace(0.05,2.4048/k_p_r_max,n_CTF)); %<-- stays positive, since first root of besselj(0,.) is ~2.4048 ;
CTF_k_p_r_kC__ = zeros(n_k_p_r,n_CTF);
for nCTF=0:n_CTF-1;
CTF_k_p_r_kC__(:,1+nCTF) = besselj(0,CTF_alpha_C_(1+nCTF)*k_p_r_);
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
pm_n_UX_rank = n_k_p_r; %<-- keep the same for now. ;
UX_knC___ = zeros(n_k_p_r,pm_n_UX_rank,n_CTF);
X_weight_rC__ = zeros(n_k_p_r,n_CTF);
UX_kn__ = eye(n_k_p_r);
X_weight_r_ = sqrt(weight_2d_k_p_r_);
for nCTF=0:n_CTF-1;
UX_knC___(:,:,1+nCTF) = UX_kn__;
X_weight_rC__(:,1+nCTF) = X_weight_r_;
end;%for nCTF=0:n_CTF-1;
clear nCTF UX_kn__ X_weight_r_ ;
%%%%%%%%;
tmp_t = tic();
n_gamma_z = n_w_max;
gamma_z_ = transpose(linspace(0,2*pi,n_gamma_z+1)); gamma_z_ = gamma_z_(1:n_gamma_z);
index_nCTF_from_nM_ = zeros(n_M,1);
index_nS_from_nM_ = zeros(n_M,1);
index_nw_from_nM_ = zeros(n_M,1);
euler_azimu_b_true_M_ = zeros(n_M,1);
euler_polar_a_true_M_ = zeros(n_M,1);
euler_gamma_z_true_M_ = zeros(n_M,1);
rng(0);
for nM=0:n_M-1;
nCTF = mod(nM,n_CTF);
CTF_k_p_r_k_ = CTF_k_p_r_kC__(:,1+nCTF);
index_nCTF_from_nM_(1+nM) = nCTF;
nS = max(0,min(n_S-1,floor(n_S*rand())));
index_nS_from_nM_(1+nM) = nS;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% M_k_p_wkM__: time %0.6fs',tmp_t)); end;
clear nCTF CTF_k_p_r_k_ nS azimu_b polar_a index_gamma_z gamma_z S_k_p_wk_ T_k_p_wk_ M_k_p_wk_ M_k_q_wk_ ;
%%%%%%%%;
 
%%%%%%%%;
delta_r_max = 0.1; svd_eps = 1e-3; n_delta_v_requested = 1024;
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% ampmh_FTK_1: time %0.6fs',tmp_t)); end;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% FTK.n_svd_l %.3d FTK.n_delta_v %.4d ',FTK.n_svd_l,FTK.n_delta_v)); end;

%%%%%%%%;
% Now calculate correlation X_wSM___ and innerproduct Z_wSM___
%%%%%%%%;
X_wSM___ = zeros(n_w_max,n_S,n_M);
gamma_z_wSM___ = zeros(n_w_max,n_S,n_M);
Z_wSM___ = zeros(n_w_max,n_S,n_M);
UX_M_l2_M_ = zeros(n_M,1);
CTF_UX_S_l2_SC__ = zeros(n_S,n_CTF);
nCTF=0; %for nCTF=0:n_CTF-1;
if (flag_verbose); disp(sprintf(' %% nCTF %d/%d',nCTF,n_CTF)); end;
UX_kn__ = UX_knC___(:,:,1+nCTF);
X_weight_r_ = X_weight_rC__(:,1+nCTF);
CTF_k_p_r_k_ = CTF_k_p_r_kC__(:,1+nCTF);
index_M_sub_ = efind(index_nCTF_from_nM_==nCTF); n_M_sub = numel(index_M_sub_);
tmp_t = tic();
svd_VUXM_sub_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M_sub,M_k_q_wkM__(:,1+index_M_sub_),pm_n_UX_rank,UX_kn__,X_weight_r_);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% tmpmh_VUXM_lwnM____3: time %0.6fs',tmp_t)); end;
tmp_t = tic();
UX_M_sub_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M_sub,pm_n_UX_rank,svd_VUXM_sub_lwnM____);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% ampmh_UX_M_sub_l2_dm__1: time %0.6fs',tmp_t)); end;
tmp_index_d_ = intersect(efind(FTK.delta_x_==0),efind(FTK.delta_y_==0)); %<-- should be zero-displacement. ;
UX_M_sub_l2_M_ = reshape(UX_M_sub_l2_dM__(1+tmp_index_d_(1+0),:),[n_M_sub,1]);
[UX_M_sub_k_q_wnM___,UX_M_sub_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,n_M_sub,svd_VUXM_sub_lwnM____);
%%%%;
tmp_t = tic();
CTF_UX_S_k_q_wnS__ = reshape(permute(reshape(reshape(S_k_q_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_weight_r_.*CTF_k_p_r_k_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
CTF_UX_S_l2_S_ = sum(abs(CTF_UX_S_k_q_wnS__).^2,1)/n_w_max;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% CTF_UX_S_k_q_wnS__: time %0.6fs',tmp_t)); end;
%%%%;
tmp_t = tic();
parameter_ampmh = struct('type','parameter');
[ ...
 parameter_ampmh ...
,X_sub_wSM___ ...
,~ ...
,~ ...
,gamma_z_sub_wSM___ ...
,~ ...
] = ...
ampmh_X_wSM___8( ...
 parameter_ampmh ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_S_ ...
,n_M_sub ...
,svd_VUXM_sub_lwnM____ ...
,UX_M_sub_l2_dM__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% ampmh_X_wSM___8: time %0.6fs',tmp_t)); end;
UX_M_l2_M_(1+index_M_sub_) = UX_M_sub_l2_M_;
CTF_UX_S_l2_SC__(:,1+nCTF) = CTF_UX_S_l2_S_;
Z_sub_wSM___ = bsxfun(@times,bsxfun(@times,X_sub_wSM___,reshape(sqrt(CTF_UX_S_l2_S_),[1,n_S,1])),reshape(sqrt(UX_M_sub_l2_M_),[1,1,n_M_sub]));
X_wSM___(:,:,1+index_M_sub_) = X_sub_wSM___;
gamma_z_wSM___(:,:,1+index_M_sub_) = gamma_z_sub_wSM___;
Z_wSM___(:,:,1+index_M_sub_) = Z_sub_wSM___;
%%%%;
clear UX_kn__ X_weight_r_ CTF_k_p_r_k_ index_M_sub_ svd_VUXM_sub_lwnM____ UX_M_sub_l2_dM__ UX_M_sub_l2_M_ ;
clear UX_M_sub_k_q_wnM___ UX_M_sub_k_p_wnM___ ;
clear CTF_UX_S_k_q_wnS__ CTF_UX_S_l2_S_ ;
clear X_sub_wSM___ gamma_z_sub_wSM___ Z_sub_wSM___ ;
%end;%for nCTF=0:n_CTF-1;
clear nCTF;
%%%%%%%%;





