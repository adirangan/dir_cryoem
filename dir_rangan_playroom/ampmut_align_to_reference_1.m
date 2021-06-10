function ...
[ ...
 parameter ...
,X_best_ ...
,X_best_flag_flip_ ...
] = ...
ampmut_align_to_reference_1( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_M ...
,M_k_p__ ...
,CTF_index_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_quad_ ...
,n_iteration ...
,euler_polar_a_calc__ ...
,euler_azimu_b_calc__ ...
,euler_gamma_z_calc__ ...
,image_delta_x_calc__ ...
,image_delta_y_calc__ ...
);

verbose=2;
if (verbose); disp(sprintf(' %% [entering ampmut_align_to_reference_1]')); end;

%%%%%%%%;
if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'cg_lsq_n_order')); parameter.cg_lsq_n_order = 5; end; %<-- parameter_bookmark. ;
cg_lsq_n_order = parameter.cg_lsq_n_order;

%%%%%%%%;
X_best_ = zeros(n_iteration,1);
X_best_flag_flip_ = zeros(n_iteration,1);
for niteration=0:n_iteration-1;
%%%%%%%%;
tmp_euler_polar_a_ = euler_polar_a_calc__(:,1+niteration);
tmp_euler_azimu_b_ = euler_azimu_b_calc__(:,1+niteration);
tmp_euler_gamma_z_ = euler_gamma_z_calc__(:,1+niteration);
tmp_image_delta_x_ = image_delta_x_calc__(:,1+niteration);
tmp_image_delta_y_ = image_delta_y_calc__(:,1+niteration);
%%%%%%%%;
tmp_t = tic;
c_k_Y_reco_ = ...
cg_lsq_4( ...
 cg_lsq_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p__ ...
,CTF_index_ ...
,CTF_k_p__ ...
,tmp_euler_polar_a_ ...
,tmp_euler_azimu_b_ ...
,tmp_euler_gamma_z_ ...
,tmp_image_delta_x_ ...
,tmp_image_delta_y_ ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% niteration %.2d/%.2d: M_k_p__ --> c_k_Y_reco_ time %0.2fs',niteration,n_iteration,tmp_t)); end;
%%%%%%%%;
tmp_t = tic;
[ ...
 X_best ...
,flag_flip ...
] = ...
register_spharm_to_spharm_wigner_wrap_1( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,0 ...
,l_max_ ...
,a_k_Y_quad_ ...
,c_k_Y_reco_ ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% niteration %.2d/%.2d: X_best time %0.2fs',niteration,n_iteration,tmp_t)); end;
%%%%%%%%;
if (verbose); disp(sprintf(' %% niteration %.2d/%.2d: X_best %0.4f X_best_flag_flip %d',niteration,n_iteration,X_best,flag_flip)); end;
X_best_(1+niteration) = X_best;
X_best_flag_flip_(1+niteration) = flag_flip;
%%%%%%%%;
end;%for niteration=0:n_iteration-1;

if (verbose); disp(sprintf(' %% [finished ampmut_align_to_reference_1]')); end;
