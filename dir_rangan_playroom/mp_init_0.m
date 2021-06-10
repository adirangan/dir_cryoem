function ...
[ ...
 mp_init_image_polar_a_0est_ ...
,mp_init_image_azimu_b_0est_ ...
,mp_init_image_gamma_z_0est_ ...
,mp_init_image_delta_0est_dM__ ...
,mp_init_image_sigma_gamma_z_0est_ ...
,mp_init_image_sigma_delta_0est_ ...
,mp_init_image_flag_update_0est_ ...
,mp_init_fnorm_d_image_gamma_z_0est_ ...
,mp_init_fnorm_d_image_delta_0est_ ...
] = ...
mp_init_0( ...
 dir_trunk ...
,string_infix ...
,n_M ...
,dat_M_uni_k_p__ ...
,dat_n_UX_rank ...
,dat_nUX_rank ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_uni_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,delta_r_max_use ...
,svd_eps ...
,n_delta_v_requested ...
,n_CTF_rank ...
,CTF_idx_ ...
,CTF_k_p_r__ ...
,UX__ ...
,X_weight_r_ ...
);

verbose=1;

date_diff_threshold = 1.5;

str_mp_init = sprintf('nUX%.3d',dat_nUX_rank);
infix_mp_init = sprintf('mp_init_%s_%s',string_infix,str_mp_init);
% ;
fname_mat_A = sprintf('%s_mat/%s_A.mat',dir_trunk,infix_mp_init);
fname_mat_B = sprintf('%s_mat/%s_B.mat',dir_trunk,infix_mp_init);
fname_mat_C = sprintf('%s_mat/%s_C.mat',dir_trunk,infix_mp_init);
fname_tmp_ABC = sprintf('%s_mat/%s_ABC.tmp',dir_trunk,infix_mp_init);
% ;
flag_exist = ( exist(fname_mat_A,'file') &  exist(fname_mat_B,'file') & exist(fname_mat_C,'file') );
if ( flag_exist);
disp(sprintf(' %% %s found, not creating',fname_mat_A));
disp(sprintf(' %% %s found, not creating',fname_mat_B));
disp(sprintf(' %% %s found, not creating',fname_mat_C));
flag_skip=0;
end;%if ( flag_exist);
if (~flag_exist);
flag_skip=0;
if ( exist(fname_tmp_ABC,'file'));
tmp_date_diff = datenum(clock) - datenum(dir(fname_tmp_ABC).date);
if (tmp_date_diff< date_diff_threshold);
disp(sprintf(' %% %s found, tmp_date_diff = %0.2f = recent, skipping',fname_tmp_ABC,tmp_date_diff));
flag_skip=1;
end;%if (tmp_date_diff< date_diff_threshold);
if (tmp_date_diff>=date_diff_threshold);
disp(sprintf(' %% %s found, tmp_date_diff = %0.2f = stale, deleting',fname_tmp_ABC,tmp_date_diff));
delete(fname_tmp_ABC);
flag_skip=0;
end;%if (tmp_date_diff>=date_diff_threshold);
end;%if ( exist(fname_tmp_ABC,'file'));
end;%if (~flag_exist);

if (~flag_skip);
save(fname_tmp_ABC,'fname_tmp_ABC');

delta_r_max = delta_r_max_use;
delta_r_p = 0.05;
delta_r_s = delta_r_max/sqrt(2*log(1/delta_r_p));
delta_r_N = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2));
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
if (verbose>1); disp(sprintf(' %% p-val %0.4f delta_r_max %0.6f sigma %0.4f N_pixel %0.4f --> FTK.n_svd_l %d, n_delta_v_requested %d',delta_r_p,delta_r_max,delta_r_s,delta_r_N,FTK.n_svd_l,n_delta_v_requested)); end;
% ;
n_w_max = max(n_w_uni_);
n_w_uni_max = max(n_w_uni_);
n_w_uni_sum = sum(n_w_uni_);
n_w_uni_csum_ = cumsum([0;n_w_uni_]);
% ;
pm_n_UX_rank = dat_n_UX_rank;
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_N_pixel = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2)) ; 
% ;
init_viewing_k_eq_d = 1/(2*pi)/sqrt(1);
init_pm_n_UX_rank = 1+dat_nUX_rank;
init_pm_n_k_p_r = init_pm_n_UX_rank;
init_pm_k_p_r_ = ones(init_pm_n_k_p_r,1);
init_pm_k_p_r_max = 1;
init_pm_n_w_ = n_w_max*ones(init_pm_n_k_p_r,1);
init_pm_n_w_max = n_w_max;
init_pm_n_w_sum = sum(init_pm_n_w_);
init_pm_n_w_csum_ = cumsum([0;init_pm_n_w_]);
init_pm_weight_3d_k_p_r_ = ones(init_pm_n_k_p_r,1);
init_pm_weight_2d_k_p_r_ = ones(init_pm_n_k_p_r,1);
flag_MS_vs_SM = 1;
% ;
[UCTF_kc__,SCTF_c__,VCTF_Mc__] = svds(CTF_k_p_r__(:,1+CTF_idx_(1:n_M)),n_CTF_rank);
VSCTF_Mc__ = VCTF_Mc__*SCTF_c__;
% ;
%%%%%%%%;
% Construct M_uni_k_q__ without taking into account the translations. ;
%% Construct M_uni_k_q__ while taking into account the translations. ;
%%%%%%%%;
tmp_t = tic();
dat_M_uni_k_q__ = zeros(n_w_uni_sum,n_M);
for nM=0:n_M-1;
dat_M_uni_k_p_ = ...
transf_p_to_p( ...
 n_k_p_r ...
,k_p_r_ ...
,n_w_uni_ ...
,n_w_uni_sum ...
,dat_M_uni_k_p__(:,1+nM) ...
,+0.0 ...
,+0.0 ...
);
dat_M_uni_k_q__(:,1+nM) = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_uni_ ...
,n_w_uni_sum ...
,dat_M_uni_k_p_ ...
);
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% M_uni_k_q__: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now form svd_VUXM_lwnM____ using these translated images. ;
%%%%%%%%;
tmp_t = tic();
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_uni_,n_M,dat_M_uni_k_q__,init_pm_n_UX_rank,UX__,X_weight_r_);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the non-translated images. ;
%%%%%%%%;
tmp_t = tic();
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_uni_,n_M,init_pm_n_UX_rank,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now, form principal-images (using the zero-displacement). ;
%%%%%%%%;
tmp_t = tic();
[UX_M_k_q_wnM___,UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_uni_,init_pm_n_UX_rank,n_M,svd_VUXM_lwnM____,zeros(n_M,1),zeros(n_M,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_k_q_wnM___: %0.6fs',tmp_t)); end;
UX_M_k_p_wnM__ = reshape(UX_M_k_p_wnM___(:,1:init_pm_n_k_p_r,:),[n_w_uni_max*init_pm_n_k_p_r,n_M]);
UX_M_k_q_wnM__ = reshape(UX_M_k_q_wnM___(:,1:init_pm_n_k_p_r,:),[n_w_uni_max*init_pm_n_k_p_r,n_M]);
% ;
%%%%%%%%;
tmp_t = tic();
SM_k_q_ = svd(UX_M_k_q_wnM__);
n_M_rank = min(find(SM_k_q_/SM_k_q_(1)<1e-3)); if isempty(n_M_rank); n_M_rank = min(size(UX_M_k_q_wnM__)); end;
[UM_k_q__,SM_k_q__,VM_k_q__] = svds(UX_M_k_q_wnM__,n_M_rank);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% SM_k_q_: %0.3fs',tmp_t)); end;
if (verbose>1); disp(sprintf(' %% n_M %d --> n_M_rank %d',n_M,n_M_rank)); end;
%%%%%%%%;
% Now calculate correlations across principal-images. ;
%%%%%%%%;
if ( exist(fname_mat_A,'file')); disp(sprintf(' %% %s found, not creating',fname_mat_A)); end;
if (~exist(fname_mat_A,'file')); disp(sprintf(' %% %s not found, creating',fname_mat_A));
%%%%%%%%;
tmp_index = efind( (FTK.delta_x_==0) & (FTK.delta_y_==0) );
UX_M_l2_M_ = UX_M_l2_dM__(1+tmp_index,:);
n_M_Mbatch = 24;
tmp_t = tic();
[ ...
 X_wMM___ ...
,delta_j_wMM___ ...
,I_value_wMM___ ...
] = ...
ampmh_X_wSM___4( ...
 FTK ...
,n_w_max ...
,init_pm_n_UX_rank ...
,n_M ...
,n_M_rank ...
,n_M_Mbatch ...
,UM_k_q__ ...
,SM_k_q__ ...
,VM_k_q__ ...
,UX_M_l2_M_ ...
,n_M ...
,n_M_Mbatch ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wMM___: %0.3fs',tmp_t)); end;
X_wMM___ = real(X_wMM___);
%%%%%%%%;
% Extract similarity-matrix using optimal gamma_z for each pair. ;
%%%%%%%%;
X_MM__ = zeros(n_M,n_M);
d_gamma_z_MM__ = zeros(n_M,n_M);
d_delta_j_MM__ = zeros(n_M,n_M);
d_delta_dMM___ = zeros(2,n_M,n_M);
for nM0=0:n_M-1;
for nM1=0:n_M-1;
[X_MM,tmp_index] = max(X_wMM___(:,1+nM0,1+nM1)); tmp_index = tmp_index-1;
assert(X_MM==X_wMM___(1+tmp_index,1+nM0,1+nM1));
d_delta_j_MM__(1+nM0,1+nM1) = delta_j_wMM___(1+tmp_index,1+nM0,1+nM1);
d_delta_dMM___(:,1+nM0,1+nM1) = [FTK.delta_x_(1+d_delta_j_MM__(1+nM0,1+nM1)) ; FTK.delta_y_(1+d_delta_j_MM__(1+nM0,1+nM1))];
d_gamma_z_MM__(1+nM0,1+nM1) = 2*pi*tmp_index/n_w_max; %<-- uniformly spaced inplane-angle gamma_z. ;
X_MM__(1+nM0,1+nM1) = X_MM;
end;%for nM1=0:n_M-1;
end;%for nM0=0:n_M-1;
%%%%%%%%;
save(fname_mat_A,'UX_M_l2_M_','d_delta_j_MM__','d_delta_dMM___','d_gamma_z_MM__','X_MM__');
end;%if (~exist(fname_mat_A,'file'));
load(fname_mat_A);

if ( exist(fname_mat_B,'file')); disp(sprintf(' %% %s found, not creating',fname_mat_B)); end;
if (~exist(fname_mat_B,'file')); disp(sprintf(' %% %s not found, creating',fname_mat_B));
%%%%%%%%;
% project images onto principal-modes of similarity-matrix. ;
% Then center on the median. ;
% Then project onto the sphere. ;
% Then push points towards uniform distribution on the sphere. ;
%%%%%%%%;
[UX_MM__,SX_MM__,VX_MM__] = svds(X_MM__,n_M_rank); SX_MM_ = diag(SX_MM__);
UX_MM_median_ = adi_median(transpose(UX_MM__(:,1:3)));
P0_k_c_0_ = UX_MM__(:,1+0) - UX_MM_median_(1+0);
P0_k_c_1_ = UX_MM__(:,1+1) - UX_MM_median_(1+1);
P0_k_c_2_ = UX_MM__(:,1+2) - UX_MM_median_(1+2);
P0_k_c_012_ = sqrt( P0_k_c_0_.^2 + P0_k_c_1_.^2 + P0_k_c_2_.^2 );
P1_k_c_0_ = P0_k_c_0_./P0_k_c_012_;
P1_k_c_1_ = P0_k_c_1_./P0_k_c_012_;
P1_k_c_2_ = P0_k_c_2_./P0_k_c_012_;
P1_k_c_01_ = sqrt( P1_k_c_0_.^2 + P1_k_c_1_.^2 );
P1_azimu_b_ = pi+atan2(P1_k_c_1_,P1_k_c_0_);
P1_polar_a_ = atan2(P1_k_c_01_,P1_k_c_2_);
P2_k_c_0_ = sin(P1_polar_a_).*cos(P1_azimu_b_);
P2_k_c_1_ = sin(P1_polar_a_).*sin(P1_azimu_b_);
P2_k_c_2_ = cos(P1_polar_a_);
P3__ = mp_init_surface_repulsion_wrap_0([P2_k_c_0_,P2_k_c_1_,P2_k_c_2_]);
P3_k_c_0_ = P3__(:,1+0);
P3_k_c_1_ = P3__(:,1+1);
P3_k_c_2_ = P3__(:,1+2);
P3_k_c_01_ = sqrt( P3_k_c_0_.^2 + P3_k_c_1_.^2 );
P3_azimu_b_ = pi+atan2(P3_k_c_1_,P3_k_c_0_);
P3_polar_a_ = atan2(P3_k_c_01_,P3_k_c_2_);
%%%%%%%%;
save(fname_mat_B ...
     ,'UX_MM__','SX_MM_','UX_MM_median_' ...
     ,'P2_k_c_0_','P2_k_c_1_','P2_k_c_2_' ...
     ,'P3__','P3_azimu_b_','P3_polar_a_' ...
     );
end;%if (~exist(fname_mat_B,'file'));
load(fname_mat_B);
mp_init_image_polar_a_0est_ = P3_polar_a_;
mp_init_image_azimu_b_0est_ = P3_azimu_b_;
P3_k_c_0_ = sin(P3_polar_a_).*cos(P3_azimu_b_);
P3_k_c_1_ = sin(P3_polar_a_).*sin(P3_azimu_b_);
P3_k_c_2_ = cos(P3_polar_a_);
P3__ = [ P3_k_c_0_ , P3_k_c_1_ , P3_k_c_2_ ];

%%%%%%%%;
% Check the relationship between delta_ and gamma_z for pairs of nearby images. ;
%%%%%%%%;
n_neighbor = 6;
P3_knn_index__ = knnsearch(P3__,P3__,'K',1+n_neighbor) - 1;
X_MM_knn_avg_ = zeros(n_M,1);
for nM0=0:n_M-1;
nM1_ = P3_knn_index__(1+nM0,1+1+[0:n_neighbor-1]);
X_MM_knn_avg_(1+nM0) = mean(X_MM__(1+nM0,1+nM1_));
end;%for nM0=0:n_M-1;
%%%%%%%%;
% find tightly knit local group. ;
%%%%%%%%;
[~,nM0] = max(X_MM_knn_avg_); nM0 = nM0 - 1;
nM1_ = P3_knn_index__(1+nM0,:);
if (verbose); disp(sprintf(' %% tightly knit local group nM0 %d',nM0)); disp(X_MM__(1+nM0,1+nM1_)); end;
%%%%%%%%;
% check relationship between true delta_ and gamma_z. ;
% We expect that X_MM__(1+nM0,1+nM1) measures: ;
% dot( rotate(M_0_ , +tmp_d_gamma_z) , transf(M_1_ , +tmp_d_delta_) );
% So if M_0_ and M_1_ are both drawn from the same template, ;
% we expect that (ignoring the CTF) : ;
% M_0_ = transf(rotate(S_0_,+tmp_gamma_z_true_0),-tmp_delta_true_0_);
% M_1_ = transf(rotate(S_1_,+tmp_gamma_z_true_1),-tmp_delta_true_1_);
% X_MM = dot( rotate(M_0_ , +tmp_d_gamma_z) , transf(M_1_ , +tmp_d_delta_) );
% Which mean that, if X_MM is close to 1, then. ;
% tmp_delta_true_1_ = rotate(+tmp_d_gamma_z) * tmp_delta_true_0_ + tmp_d_delta_ ;
% tmp_gamma_z_true_1 = tmp_gamma_z_true_0 + tmp_d_gamma_z ;
%%%%%%%%;

if ( exist(fname_mat_C,'file')); disp(sprintf(' %% %s found, not creating',fname_mat_C)); end;
if (~exist(fname_mat_C,'file')); disp(sprintf(' %% %s not found, creating',fname_mat_C));
%%%%%%%%;
% Now run a simple message-passing algorithm to set mp_init_gamma_z_0est_ and mp_init_delta_x_0est_ and mp_init_delta_y_0est_. ;
% Given any clique of 1+n_neighbor nodes, ;
% we are searching for an assignment of mp_init_gamma_z_0est_ and mp_init_delta_x_0est_ and mp_init_delta_y_0est_ ;
% such that the errors: ;
% tmp_error_gamma__(nM0,nM1) = fnorm(periodize(tmp_gamma_z_0est_1 - tmp_gamma_z_true_1,-pi/2,+pi/2));
%%%%%%%%;
mp_init_n_iteration = 512;
mp_init_eta_step = 0.125;
%%%%%%%%;
mp_init_image_flag_update_0est_ = zeros(n_M,1);
mp_init_image_gamma_z_0est_ = zeros(n_M,1);
mp_init_image_delta_0est_dM__ = zeros(2,n_M);
%%%%;
% seed most similar clique with null data. ;
%%%%;
[~,nM0] = max(X_MM_knn_avg_); nM0 = nM0 - 1;
mp_init_image_gamma_z_0est_(1+nM0) = 0.0; 
mp_init_image_delta_0est_dM__(:,1+nM0) = [0.0;0.0];
mp_init_image_flag_update_0est_(1+nM0) = 1;
%%%%;
[ ...
 mp_init_image_gamma_z_0est_ ...
,mp_init_image_delta_0est_dM__ ...
,mp_init_image_sigma_gamma_z_0est_ ...
,mp_init_image_sigma_delta_0est_ ...
,mp_init_image_flag_update_0est_ ...
,mp_init_fnorm_d_image_gamma_z_0est_ ...
,mp_init_fnorm_d_image_delta_0est_ ...
] = ...
mp_init_align_gamma_z_delta_0( ...
 n_M ...
,P3__ ...
,d_gamma_z_MM__ ...
,d_delta_dMM___ ...
,mp_init_image_flag_update_0est_ ...
,mp_init_image_gamma_z_0est_ ...
,mp_init_image_delta_0est_dM__ ...
,mp_init_n_iteration ...
,mp_init_eta_step ...
);
%%%%%%%%;
save(fname_mat_C ...
,'mp_init_image_gamma_z_0est_' ...
,'mp_init_image_delta_0est_dM__' ...
,'mp_init_image_sigma_gamma_z_0est_' ...
,'mp_init_image_sigma_delta_0est_' ...
,'mp_init_image_flag_update_0est_' ...
,'mp_init_fnorm_d_image_gamma_z_0est_' ...
,'mp_init_fnorm_d_image_delta_0est_' ...
     );
end;%if (~exist(fname_mat_C,'file'));

delete(fname_tmp_ABC);
end;%if (~flag_skip);

load(fname_mat_C);
