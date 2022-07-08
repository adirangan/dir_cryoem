%%%%%%%%;
% based on test_viewing_angle_distribution_0.m ;
% trying ring-restriction applied to LetB1. ;
%%%%%%%%;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;
dir_base = sprintf('/%s/rangan/dir_cryoem',string_root);
flag_replot = 0;
flag_replot_vol = 1;
flag_verbose = 1; nf=0;
date_diff_threshold = 0.25;
flag_force_create_mat=0; flag_force_create_tmp=0;
flag_qbp_vs_lsq = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
table_data__ = { ...
%'p28hRPT1_x0' , 'p28hRP' , 0.98 , 'emd_8674.map' , 'T1.star' ; ...
'ISWINCP_x0' , 'ISWINCP' , 1.07 , 'emd_9718.map' , 'ADPBeF.star' ; ...
'trpv1_x0' , 'trpv1' , 1.2156 , 'emd_5778.mrc' , 'tv1_relion_data.star' ; ...
'rib80s_x0' , 'rib80s' , 1.34 , 'emd_2660.mrc' , 'shiny_2sets.star' ; ...
'MlaFEDB_x0' , 'MlaFEDB' , 0.832 , 'emd_22116.map' , 'Empiar_10536_00_to_23.star' ; ...
'LetB1_x0' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1.star' ; ...
'TMEM16F_x0' , 'TMEM16F' , 1.059 , 'emd_20244.map' , 'All_8192.star' ; ...
'LSUbl17dep_x0' , 'LSUbl17dep' , 1.31 , 'emd_8434.map' , 'Parameters_negated.star' ; ...
'ps1_x0' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
};
n_experiment = size(table_data__,1);
ampm__ = cell(n_experiment,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_ = clock;rng(tmp_(end));
%for nexperiment=(randperm(n_experiment)-1);
nexperiment = 4; %<-- LetB1. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
na=0;
fname_prefix = table_data__{1+nexperiment,1+na}; na=na+1;
dir_nopath_data_star = table_data__{1+nexperiment,1+na}; na=na+1;
Pixel_Spacing = table_data__{1+nexperiment,1+na}; na=na+1;
fname_nopath_volume = table_data__{1+nexperiment,1+na}; na=na+1;
fname_nopath_star = table_data__{1+nexperiment,1+na}; na=na+1;
if (flag_verbose); 
disp(sprintf(' %% nexperiment %d/%d: %16s %16s %0.3f %16s %32s' ...
	     ,nexperiment,n_experiment ...
	     ,fname_prefix,dir_nopath_data_star,Pixel_Spacing,fname_nopath_volume,fname_nopath_star ...
	     ));
end;%if (flag_verbose); 
global_parameter = struct('type','parameter');
global_parameter.flag_replot=flag_replot;
global_parameter.flag_replot_vol=flag_replot_vol;
global_parameter.flag_invert=0;
global_parameter.flag_center_image=0;
global_parameter.tolerance_master = 1e-2;
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); global_parameter.flag_invert = 0; end;
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); global_parameter.flag_center_image = 1; end;
if (strcmp(dir_nopath_data_star,'precatalytic_spliceosome')); global_parameter.flag_center_image = 1; end;
dir_pm = sprintf('%s/dir_%s/dir_pm',dir_base,fname_prefix);
if (flag_verbose); disp(sprintf(' %% fname_prefix: %s',fname_prefix)); end;
[~,ampm_] = ampm_fromdisk_3([],dir_pm);
ampm__{1+nexperiment} = ampm_;

tolerance_master = global_parameter.tolerance_master;
flag_replot = global_parameter.flag_replot;
flag_replot_vol = global_parameter.flag_replot_vol;
N_k_p_use__ = ampm_.M_k_p__; if (global_parameter.flag_center_image==1); N_k_p_use__ = ampm_.N_k_p__; end;
M_k_p_wkM__ = N_k_p_use__;
tmp_delta_r_max = 0.030;
tmp_delta_r_upb = 0.250;
n_M = ampm_.n_M;
n_iteration_orig = size(ampm_.corr_a_k_Y_i_,1);
n_iteration_seco = max(1,min(n_iteration_orig-2,6)); %<-- keep at 6?. ;
%%%%%%%%;
% indices. ;
%%%%%%%%;
n_k_all = ampm_.n_k_all;
n_k_all_csum_ = ampm_.n_k_all_csum_;
k_p_r_all_ = ampm_.k_p_r_all_;
k_p_azimu_b_all_ = ampm_.k_p_azimu_b_all_;
k_p_polar_a_all_ = ampm_.k_p_polar_a_all_;
weight_3d_k_all_ = ampm_.weight_3d_k_all_;
weight_shell_k_ = ampm_.weight_shell_k_;
n_k_p_r = ampm_.n_k_p_r;
k_p_r_ = ampm_.k_p_r_;
k_p_r_max = ampm_.k_p_r_max;
weight_2d_k_p_r_ = ampm_.weight_2d_k_p_r_;
weight_3d_k_p_r_ = ampm_.weight_3d_k_p_r_;
l_max_ = ampm_.l_max_;
n_w_ = ampm_.n_w_;
n_w_max = ampm_.n_w_max;
%%%%%%%%;
n_k_all = ampm_.n_k_all;
weight_3d_k_all_ = ampm_.weight_3d_k_all_;
k_c_0_all_ = ampm_.k_c_0_all_;
k_c_1_all_ = ampm_.k_c_1_all_;
k_c_2_all_ = ampm_.k_c_2_all_;
k_p_r_max = ampm_.k_p_r_max;
n_xxx_u = ampm_.n_xxx_u;
xxx_u_weight_ = ampm_.xxx_u_weight_;
x_u_0___ = ampm_.x_u_0___;
x_u_1___ = ampm_.x_u_1___;
x_u_2___ = ampm_.x_u_2___;
x_p_r_max = ampm_.x_p_r_max;
%%%%%%%%;
index_nCTF_from_nM_ = ampm_.index_nCTF_from_nM_;
n_CTF = ampm_.n_CTF;
CTF_k_p_r_kC__ = ampm_.CTF_k_p_r_kC__;
%%%%%%%%;
a_k_Y_true_yk_ = ampm_.a_k_Y_quad_;
euler_polar_a_true_ = ampm_.euler_polar_a_Mi__(:,end);
euler_azimu_b_true_ = ampm_.euler_azimu_b_Mi__(:,end);
euler_gamma_z_true_ = ampm_.euler_gamma_z_Mi__(:,end);
image_delta_x_true_ = ampm_.image_delta_x_acc_Mi__(:,end) + ampm_.image_delta_x_upd_Mi__(:,end);
image_delta_y_true_ = ampm_.image_delta_y_upd_Mi__(:,end) + ampm_.image_delta_y_upd_Mi__(:,end);
%%%%%%%%;
sample_sphere_k_eq_d = 1/(2*pi);
half_diameter_x_c = ampm_.half_diameter_x_c;
n_x_u_pack = ampm_.n_x_u_pack;
%%%%%%%%;
qbp_eps = tolerance_master;

%%%%%%%%;
% sort. ;
%%%%%%%%;
[~,ij_X_best_ampm_srt_] = sort(ampm_.X_best_ampm_ia__(end,:),'descend'); index_X_best_ampm_srt_ = ij_X_best_ampm_srt_-1;

%%%%%%%%;
% Choose the na. ;
%%%%%%%%;
na_ = index_X_best_ampm_srt_([1:9,round(numel(ij_X_best_ampm_srt_)/2)-14:round(numel(ij_X_best_ampm_srt_)/2)+14]); %<-- take top 9 best as well as 29 near the median. ;
n_na = numel(na_);
%%%%%%%%%%%%%%%%;
for nna=0:n_na-1;
%%%%%%%%%%%%%%%%;
na = na_(1+nna);
str_fname_mat = ampm_.str_fname_mat_a_{1+na};
str_fname_mat_pre = str_fname_mat(1:strfind(str_fname_mat,'.mat')-1);

%%%%%%%%;
% Now run ampmut once again, this time using ring. ;
%%%%%%%%;
tmp_str_ = ampm_.str_fname_mat_a_{1+na};
tmp_key_ = 'X_2d_xcor_d0_';
tmp_str_xfix = tmp_str_(strfind(tmp_str_,tmp_key_) + length(tmp_key_):strfind(tmp_str_,'.mat')-1);
XC_fname_pre = sprintf('%s_mat/X_2d_xcor_d0_r%s',ampm_.dir_pm,tmp_str_xfix);
[XC_flag_skip,XC_fname_mat] = open_fname_tmp(XC_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);

XC_fname_align_a_k_Y_pre = sprintf('%s_mat/X_2d_xcor_d0_r%s_align_a_k_Y_',ampm_.dir_pm,tmp_str_xfix);
XC_fname_align_a_k_Y_jpg = sprintf('%s.jpg',XC_fname_align_a_k_Y_pre);
XC_fname_snapshot_pre = sprintf('%s_mat/X_2d_xcor_d0_r%s_snapshot',ampm_.dir_pm,tmp_str_xfix);
XC_fname_snapshot_jpg = sprintf('%s.jpg',XC_fname_snapshot_pre);
XC_fname_compare_image_rank_pre = sprintf('%s_mat/X_2d_xcor_d0_r%s_compare_image_rank',ampm_.dir_pm,tmp_str_xfix);
XC_fname_compare_image_rank_jpg = sprintf('%s.jpg',XC_fname_compare_image_rank_pre);


if ~XC_flag_skip;
%%%%%%%%;
%try;
disp(sprintf(' %% %s not found, creating',XC_fname_mat));
XC_parameter = struct('type','parameter');
XC_parameter.rseed = 0;
XC_parameter.flag_rank_vs_tolerance = 0;
XC_parameter.tolerance_pm = 1e-2;
XC_parameter.delta_r_max = tmp_delta_r_max;
XC_parameter.delta_r_upb = tmp_delta_r_upb;
XC_parameter.dir_pm = ampm_.dir_pm;
XC_parameter.flag_alternate_MS_vs_SM = 1;
XC_parameter.fname_pre = XC_fname_pre;
XC_parameter.flag_euler_polar_a_restrict = 1;
ampmut_wrap_5( ...
 XC_parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_max ...
,l_max_ ...
,index_nCTF_from_nM_ ...
,n_CTF ...
,CTF_k_p_r_kC__ ...
,n_M ...
,M_k_p_wkM__ ...
,ampm_.euler_polar_a_ampm_Ma__(:,1+na) ...
,ampm_.euler_azimu_b_ampm_Ma__(:,1+na) ...
,ampm_.euler_gamma_z_ampm_Ma__(:,1+na) ...
,ampm_.image_delta_x_ampm_Ma__(:,1+na) ...
,ampm_.image_delta_y_ampm_Ma__(:,1+na) ...
,ampm_.a_k_Y_ampm_yka__(:,1+na) ...
);
%catch; disp(sprintf(' %% WARNING: error generating %s',XC_fname_mat)); end;%try;
close_fname_tmp(XC_fname_pre);
end;%if ~XC_flag_skip;

%%%%%%%%;
if ( exist(XC_fname_mat,'file'));
%%%%%%%%;
% First collect statistics. ;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% %s found, aligning',XC_fname_mat)); end;
tmp_XC_ = load(XC_fname_mat);
if (~isfield(tmp_XC_.parameter,'fname_pre')); tmp_XC_.parameter.fname_pre = XC_fname_pre; end; %<-- parameter_bookmark. ;
%%%%%%%%;
tmp_XC_fname_pre = sprintf('%s_align_a_CTF_avg_UX_Y_',tmp_XC_.parameter.fname_pre);
tmp_XC_fname_pre = rootswitch(tmp_XC_fname_pre,string_root,'rangan');
tmp_XC_.parameter.fname_align_a_CTF_avg_UX_Y_pre = tmp_XC_fname_pre;
[tmp_XC_flag_skip,tmp_XC_fname_mat] = open_fname_tmp(tmp_XC_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_XC_flag_skip;
% do nothing. ;
close_fname_tmp(tmp_XC_fname_pre);
end;%if ~tmp_XC_flag_skip;
%%%%%%%%;
tmp_XC_fname_pre = sprintf('%s_align_a_k_Y_',tmp_XC_.parameter.fname_pre);
tmp_XC_fname_pre = rootswitch(tmp_XC_fname_pre,string_root,'rangan');
tmp_XC_.parameter.fname_align_a_k_Y_pre = tmp_XC_fname_pre;
[tmp_XC_flag_skip,tmp_XC_fname_mat] = open_fname_tmp(tmp_XC_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_XC_flag_skip;
%try;
disp(sprintf(' %% %s not found, creating',tmp_XC_fname_mat));
ampmut_align_to_a_k_Y_1( ...
 tmp_XC_.parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_M ...
,M_k_p_wkM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,l_max_ ...
,a_k_Y_true_yk_ ...
,[] ...
,tmp_XC_.euler_polar_a_Mi__ ...
,tmp_XC_.euler_azimu_b_Mi__ ...
,tmp_XC_.euler_gamma_z_Mi__ ...
,tmp_XC_.image_delta_x_acc_Mi__ + tmp_XC_.image_delta_x_upd_Mi__ ...
,tmp_XC_.image_delta_y_acc_Mi__ + tmp_XC_.image_delta_y_upd_Mi__ ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
);
%catch; disp(sprintf(' %% WARNING: error generating %s',tmp_XC_fname_mat)); end;%try;
close_fname_tmp(tmp_XC_fname_pre);
end;%if ~tmp_XC_flag_skip;
%%%%%%%%;
% Now use the final step to reconstruct the molecule. ;
%%%%%%%%;
tmp_XC_euler_polar_a_ = tmp_XC_.euler_polar_a_Mi__(:,end);
tmp_XC_euler_azimu_b_ = tmp_XC_.euler_azimu_b_Mi__(:,end);
tmp_XC_euler_gamma_z_ = tmp_XC_.euler_gamma_z_Mi__(:,end);
tmp_XC_image_delta_x_ = tmp_XC_.image_delta_x_acc_Mi__(:,end) + tmp_XC_.image_delta_x_upd_Mi__(:,end);
tmp_XC_image_delta_y_ = tmp_XC_.image_delta_y_acc_Mi__(:,end) + tmp_XC_.image_delta_y_upd_Mi__(:,end);
fname_XC_k_Y_mat = sprintf('%s_a_k_Y_.mat',tmp_XC_.parameter.fname_pre);
fname_XC_k_Y_mat = rootswitch(fname_XC_k_Y_mat,string_root,'rangan');
if (~exist(fname_XC_k_Y_mat,'file'));
if (flag_verbose); disp(sprintf(' %% %s not found, creating',fname_XC_k_Y_mat)); end;
if flag_qbp_vs_lsq==0;
a_k_Y_reco_ = ...
cg_lsq_6( ...
 cg_lsq_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,tmp_XC_euler_polar_a_ ...
,tmp_XC_euler_azimu_b_ ...
,tmp_XC_euler_gamma_z_ ...
,tmp_XC_image_delta_x_ ...
,tmp_XC_image_delta_y_ ...
);
end;%if flag_qbp_vs_lsq==0;
if flag_qbp_vs_lsq==1;
a_k_Y_reco_ = ...
qbp_6( ...
 qbp_eps ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,tmp_XC_euler_polar_a_ ...
,tmp_XC_euler_azimu_b_ ...
,tmp_XC_euler_gamma_z_ ...
,tmp_XC_image_delta_x_ ...
,tmp_XC_image_delta_y_ ...
);
end;%if flag_qbp_vs_lsq==1;
save(fname_XC_k_Y_mat,'a_k_Y_reco_');
end;%if (~exist(fname_XC_k_Y_mat,'file'));
%tmp_ = load(fname_XC_k_Y_mat); tmp_XC_a_k_Y_reco_ = tmp_.a_k_Y_reco_; clear tmp_; %<-- not used. ;
%%%%%%%%;
end;%if ( exist(XC_fname_mat,'file'));
%%%%%%%%;

%%%%%%%%;
if ( exist(XC_fname_mat,'file') & ~exist(XC_fname_snapshot_jpg) );
[XC_fname_snapshot_flag_skip,XC_fname_snapshot_jpg] = open_fname_tmp(XC_fname_snapshot_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp,'jpg');
if ~XC_fname_snapshot_flag_skip;
%%%%%%%%;
if (flag_verbose) disp(sprintf(' %% %s not found, creating',XC_fname_snapshot_jpg)); end;
tmp_XC_ = load(XC_fname_mat);
flag_found = 1; tmp_XC_align_a_k_Y_ = []; tmp_XC_a_k_Y_ = [];
%%%%;
tmp_XC_fname_pre = sprintf('%s_align_a_k_Y_',tmp_XC_.parameter.fname_pre);
tmp_XC_fname_pre = rootswitch(tmp_XC_fname_pre,string_root,'rangan');
tmp_XC_.parameter.fname_align_a_k_Y_pre = tmp_XC_fname_pre;
tmp_XC_fname_mat = sprintf('%s.mat',tmp_XC_fname_pre);
if (~exist(tmp_XC_fname_mat,'file'));
disp(sprintf(' %% Warning, %s not found in ampmut_wrap_wrap_5',tmp_XC_fname_mat));
flag_found = 0;
end;%if (~exist(tmp_XC_fname_mat,'file'));
if ( exist(tmp_XC_fname_mat,'file'));
tmp_XC_align_a_k_Y_ = load(tmp_XC_fname_mat);
end;%if ( exist(tmp_XC_fname_mat,'file'));
%%%%;
tmp_XC_fname_pre = sprintf('%s_a_k_Y_',tmp_XC_.parameter.fname_pre);
tmp_XC_fname_pre = rootswitch(tmp_XC_fname_pre,string_root,'rangan');
tmp_XC_.parameter.fname_a_k_Y_pre = tmp_XC_fname_pre;
tmp_XC_fname_mat = sprintf('%s.mat',tmp_XC_fname_pre);
if (~exist(tmp_XC_fname_mat,'file'));
disp(sprintf(' %% Warning, %s not found in ampmut_wrap_wrap_5',tmp_XC_fname_mat));
flag_found = 0;
end;%if (~exist(tmp_XC_fname_mat,'file'));
if ( exist(tmp_XC_fname_mat,'file'));
tmp_XC_a_k_Y_ = load(tmp_XC_fname_mat);
end;%if ( exist(tmp_XC_fname_mat,'file'));
%%%%;
if ( flag_found & ~isempty(tmp_XC_align_a_k_Y_) & ~isempty(tmp_XC_a_k_Y_) );
%%%%;
tmp_X_best = tmp_XC_align_a_k_Y_.X_best_(end);
tmp_flag_flip = tmp_XC_align_a_k_Y_.X_best_flag_flip_(end);
tmp_polar_a_best = tmp_XC_align_a_k_Y_.polar_a_best_(end);
tmp_azimu_b_best = tmp_XC_align_a_k_Y_.azimu_b_best_(end);
tmp_gamma_z_best = tmp_XC_align_a_k_Y_.gamma_z_best_(end);
tmp_delta_x_best = tmp_XC_align_a_k_Y_.delta_best__(1+0,end);
tmp_delta_y_best = tmp_XC_align_a_k_Y_.delta_best__(1+1,end);
tmp_a_k_Y_reco_ = tmp_XC_a_k_Y_.a_k_Y_reco_;
%%%%;
[ ... 
 tmp_b_k_Y_reco_ ...
] = ...
spharm_register_and_rotate_2( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_true_yk_ ...
,tmp_a_k_Y_reco_ ...
,0 ...
,tmp_X_best ...
,tmp_flag_flip ...
,tmp_polar_a_best ...
,tmp_azimu_b_best ...
,tmp_gamma_z_best ...
,tmp_delta_x_best ...
,tmp_delta_y_best ...
);
%%%%;
[ ... 
  tmp_b_x_u_reco_ ...
] = ...
convert_spharm_to_x_c_3( ...
 sample_sphere_k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,tmp_b_k_Y_reco_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
%%%%;
figure(5);clf;figbig;
subplot(1,2,1);
isosurface_f_x_u_0(reshape(real(tmp_b_x_u_reco_),n_x_u_pack,n_x_u_pack,n_x_u_pack),[90,95,99]);
title(sprintf('tmp_b_x_u_: corr %0.4f',tmp_X_best),'Interpreter','none');
subplot(1,2,2);
isosurface_f_x_u_0(reshape(real(tmp_b_x_u_reco_),n_x_u_pack,n_x_u_pack,n_x_u_pack),[98.5]);
title(sprintf('tmp_b_x_u_: corr %0.4f',tmp_X_best),'Interpreter','none');
sgtitle(XC_fname_snapshot_jpg,'Interpreter','none');
disp(sprintf(' %% writing %s',XC_fname_snapshot_jpg));
print('-djpeg',XC_fname_snapshot_jpg);
close(gcf);
%%%%%%%%;
end;%if ( flag_found & ~isempty(tmp_XC_align_a_k_Y_) & ~isempty(tmp_XC_a_k_Y_) );
%%%%%%%%;
close_fname_tmp(XC_fname_snapshot_pre);
end;%if ~XC_fname_snapshot_flag_skip;
end;%if ( exist(XC_fname_mat,'file') & ~exist(XC_fname_snapshot_jpg) );
%%%%%%%%;

%%%%%%%%;
if ( exist(XC_fname_mat,'file') & ~exist(XC_fname_compare_image_rank_jpg) );
[XC_fname_compare_image_rank_flag_skip,XC_fname_compare_image_rank_jpg] = open_fname_tmp(XC_fname_compare_image_rank_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp,'jpg');
if ~XC_fname_compare_image_rank_flag_skip;
%%%%%%%%;
if (flag_verbose) disp(sprintf(' %% %s not found, creating',XC_fname_compare_image_rank_jpg)); end;
[ ...
 XC_parameter ...
] = ...
ampmut_compare_image_rank_1( ...
 XC_parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,n_CTF ...
,1 ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,l_max_ ...
,a_k_Y_true_yk_ ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
,[] ...
,XC_fname_mat ...
,XC_fname_compare_image_rank_pre ...
);
%%%%%%%%;
close_fname_tmp(XC_fname_compare_image_rank_pre);
end;%if ~XC_fname_compare_image_rank_flag_skip;
end;%if ( exist(XC_fname_mat,'file') & ~exist(XC_fname_compare_image_rank_jpg) );
%%%%%%%%;

%%%%%%%%;
% Now display a 'before' and 'after' picture. ;
%%%%%%%%;
flag_continue=0;
flag_continue = 1 ...
  & exist(XC_fname_mat,'file') ...
  & exist(sprintf('%s.mat',XC_fname_align_a_k_Y_pre),'file') ...
  ;
if flag_continue;
%%%%%%%%;
flag_disp = 1;
if flag_disp;
%%%%%%%%;
fname_ring_vol_pre = sprintf('%s_vol',XC_fname_pre);
fname_ring_vol_fig_pre = sprintf('%s_FIGD',fname_ring_vol_pre);
fname_ring_vol_fig_jpg = sprintf('%s.jpg',fname_ring_vol_fig_pre);
if (flag_replot_vol | ~exist(fname_ring_vol_fig_jpg,'file'));
%%%%%%%%;
test_viewing_angle_ring_vol_helper_0;
%%%%%%%%;
if strfind(ampm_.dir_pm,'ISWINCP');    percent_threshold_ = 94.50; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'trpv1');      percent_threshold_ = 91.25; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'rib80s');     percent_threshold_ = 86.25; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'MlaFEDB');    percent_threshold_ = 95.00; tmp_nx = 8; end;
if strfind(ampm_.dir_pm,'LetB1');      percent_threshold_ = 91.75; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'TMEM16F');    percent_threshold_ = 94.50; tmp_nx = 16; end;
if strfind(ampm_.dir_pm,'LSUbl17dep'); percent_threshold_ = 86.25; tmp_nx = 15; end;
if strfind(ampm_.dir_pm,'ps1');        percent_threshold_ = 96.00; tmp_nx = 11; end;
%percent_threshold_ = 85:1.25:97.5;
n_percent_threshold = numel(percent_threshold_);
ncrop = 40;
%%%%;
n_x_u_pack = ampm_.n_x_u_pack;
tmp_window_ = zeros(n_x_u_pack,n_x_u_pack,n_x_u_pack);
tmp_index_ = tmp_nx:n_x_u_pack-1-tmp_nx;
tmp_window_(1+tmp_index_,1+tmp_index_,1+tmp_index_)=1;
tmp_index_ = efind(tmp_window_);
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;
fontsize_use = 12;
if (n_percent_threshold> 1); p_row = 2; p_col = n_percent_threshold; end;
if (n_percent_threshold==1); p_row = 1; p_col = 2; end;
for npt=0:n_percent_threshold-1;
tmp_percent_threshold = percent_threshold_(1+npt); nnp=0;
%%%%;
subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
isosurface_f_x_u_0(a_x_u_ampm_alig_(1+tmp_index_),tmp_percent_threshold);
title(sprintf('AMPM %.3f (Z %0.2f)',tmp_percent_threshold,corr_full_reco_vs_crop_ampm_x_(1+ncrop))); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
isosurface_f_x_u_0(a_x_u_ampm_ring_alig_(1+tmp_index_),tmp_percent_threshold);
title(sprintf('AMPM (ring) %.3f (Z %0.2f)',tmp_percent_threshold,corr_full_reco_vs_crop_ampm_ring_x_(1+ncrop))); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
end;%for npt=0:n_percent_threshold-1;
set(gcf,'Position',1+[0,0,512*2,512+64]);
sgtitle(fname_ring_vol_fig_pre,'Interpreter','none');
print('-djpeg',fname_ring_vol_fig_jpg);
close gcf;
%%%%%%%%;
end;%if (flag_replot_vol | ~exist(fname_ring_vol_fig_jpg,'file'));
end;%if flag_disp;
%%%%%%%%;
end;%if flag_continue;

%%%%%%%%%%%%%%%%;
end;%for nna=0:n_na-1;
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%end;%for nexperiment=0:n_experiment-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;



