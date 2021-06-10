function ...
tpmutam_X_true_0( ...
 fname_0 ...
,fname_1 ...
,fname_2 ...
,n_M ...
,M_k_p__ ...
,n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,CTF_idx_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_quad_ ...
);

% updates mat-files with X_true_. ;
% test_principled_marching_updated_translation_alternating_minimization_X_true_.m ;

date_diff_threshold = 0.1;
verbose=1;
flag_recalculate=0;

%%%%%%%%;
fname_mat = sprintf('%s_%s_%s.mat',fname_0,fname_1,fname_2);
if (~exist(fname_mat,'file')); disp(sprintf(' %% Warning. %s not found in tpmutam_X_true_0.',fname_mat)); end;
if ( exist(fname_mat,'file')); disp(sprintf(' %% %s found.',fname_mat)); end;
tmp_ = load(fname_mat);
flag_found = 0;
fname_field = sprintf('X_true_%s_',fname_1);
if ( isfield(tmp_,fname_field)); disp(sprintf(' %% %s found, not creating',fname_field)); end;
if (flag_recalculate | ~isfield(tmp_,fname_field));

flag_skip=0;
fname_mat_upd = sprintf('%s.upd',fname_mat);
if ( exist(fname_mat_upd,'file'));
upd_date_diff = datenum(clock) - datenum(dir(fname_mat_upd).date);
if (upd_date_diff< date_diff_threshold);
disp(sprintf(' %% %s found, upd_date_diff = %0.2f = recent, skipping',fname_mat_upd,upd_date_diff));
flag_skip=1;
end;%if (upd_date_diff< date_diff_threshold);
if (upd_date_diff>=date_diff_threshold);
disp(sprintf(' %% %s found, upd_date_diff = %0.2f = stale, deleting',fname_mat_upd,upd_date_diff));
delete(fname_mat_upd);
flag_skip=0;
end;%if (upd_date_diff>=date_diff_threshold);
end;%if ( exist(fname_mat_upd,'file'));

if (~flag_skip);
disp(sprintf(' %% %s.%s not found, creating',fname_mat,fname_field));
save(fname_mat_upd,'fname_mat');

if strcmp(fname_1,'SM');
mat_X_best_ = tmp_.X_best_SM_;
mat_a_CTF_UX_Y_0lsq__ = tmp_.a_CTF_UX_Y_0lsq_SM__;
mat_euler_polar_a__ = tmp_.euler_polar_a_SM__;
mat_euler_azimu_b__ = tmp_.euler_azimu_b_SM__;
mat_euler_gamma_z__ = tmp_.euler_gamma_z_SM__;
mat_image_delta_x__ = tmp_.image_delta_x_SM__;
mat_image_delta_y__ = tmp_.image_delta_y_SM__;
mat_image_I_value__ = tmp_.image_I_value_SM__;
mat_image_X_value__ = tmp_.image_X_value_SM__;
mat_image_S_index__ = tmp_.image_S_index_SM__;
mat_X_best_time = tmp_.X_best_SM_time;
mat_dat_M_loading__ = tmp_.dat_M_loading_SM__;
mat_dat_M_loading_time = tmp_.dat_M_loading_SM_time;
end;%if strcmp(fname_1,'SM');
if strcmp(fname_1,'MS');
mat_X_best_ = tmp_.X_best_MS_;
mat_a_CTF_UX_Y_0lsq__ = tmp_.a_CTF_UX_Y_0lsq_MS__;
mat_euler_polar_a__ = tmp_.euler_polar_a_MS__;
mat_euler_azimu_b__ = tmp_.euler_azimu_b_MS__;
mat_euler_gamma_z__ = tmp_.euler_gamma_z_MS__;
mat_image_delta_x__ = tmp_.image_delta_x_MS__;
mat_image_delta_y__ = tmp_.image_delta_y_MS__;
mat_image_I_value__ = tmp_.image_I_value_MS__;
mat_image_X_value__ = tmp_.image_X_value_MS__;
mat_image_S_index__ = tmp_.image_S_index_MS__;
mat_X_best_time = tmp_.X_best_MS_time;
mat_dat_M_loading__ = tmp_.dat_M_loading_MS__;
mat_dat_M_loading_time = tmp_.dat_M_loading_MS_time;
end;%if strcmp(fname_1,'MS');
clear tmp_;
%%%%%%%%;

n_iteration = numel(mat_X_best_);
X_true_ = zeros(n_iteration,1);
for niteration=0:n_iteration-1;
tmp_euler_polar_a_ = mat_euler_polar_a__(:,1+niteration);
tmp_euler_azimu_b_ = mat_euler_azimu_b__(:,1+niteration);
tmp_euler_gamma_z_ = mat_euler_gamma_z__(:,1+niteration);
tmp_image_delta_x_ = mat_image_delta_x__(:,1+niteration);
tmp_image_delta_y_ = mat_image_delta_y__(:,1+niteration);
%%%%%%%%;
tmp_t = tic;
c_k_Y_reco_ = ...
cg_lsq_4( ...
 n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p__ ...
,CTF_idx_ ...
,CTF_k_p__ ...
,tmp_euler_polar_a_ ...
,tmp_euler_azimu_b_ ...
,tmp_euler_gamma_z_ ...
,tmp_image_delta_x_ ...
,tmp_image_delta_y_ ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% niteration %.2d/%.2d: M_k_p__ --> c_k_Y_reco_ time %0.2fs',niteration,n_iteration,tmp_t)); end;
[tmp_X_true_orig,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,c_k_Y_reco_);
[tmp_X_true_flip,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,flipY(n_k_p_r,l_max_,c_k_Y_reco_));
if (tmp_X_true_orig>=tmp_X_true_flip); flag_flip = 0; end;
if (tmp_X_true_orig< tmp_X_true_flip); flag_flip = 1; end;
tmp_X_true_reco = max(tmp_X_true_orig,tmp_X_true_flip);
if (verbose); disp(sprintf(' %% niteration %.2d/%.2d: tmp_X_true_reco %0.3f <-- flag_flip %d',niteration,n_iteration,tmp_X_true_reco,flag_flip)); end;
X_true_(1+niteration) = tmp_X_true_reco;
end;%for niteration=0:n_iteration-1;

if strcmp(fname_1,'SM');
X_true_SM_ = X_true_;
save(fname_mat,'-append','X_true_SM_');
end;%if strcmp(fname_1,'SM');
if strcmp(fname_1,'MS');
X_true_MS_ = X_true_;
save(fname_mat,'-append','X_true_MS_');
end;%if strcmp(fname_1,'MS');

delete(fname_mat_upd);
end;%if (~flag_skip);

end;%if (flag_recalculate | ~isfield(tmp_,fname_field));
