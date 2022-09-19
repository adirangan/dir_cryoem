clear;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

dir_base = sprintf('/%s/rangan/dir_cryoem',string_root);
dir_jpg = sprintf('%s/dir_jpg',dir_base);
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;
dir_mat = sprintf('%s/dir_mat',dir_base);
if (~exist(dir_mat,'dir')); disp(sprintf(' %% mkdir %s',dir_mat)); mkdir(dir_mat); end;

verbose=1;
flag_recoll = 0;
flag_recalc = 0;
flag_replot = 0;
flag_center_volume = 0;
flag_center_image = 0;
flag_invert = 0;
flag_crop = 0;
tolerance_master = 1e-2;
nf=0;

table_data__ = { ...
 0 , 'p28hRPT1_x0' , 'p28hRPT1' , 0.98 , 'emd_8674.map' , 'T1.star' , +2.5 ; ...
 0 , 'trpv1_x0' , 'trpv1' , 1.2156 , 'emd_5778.mrc' , 'tv1_relion_data.star' , +2.5 ; ...
 1 , 'trpv1c'   , 'trpv1' , 1.2156 , 'emd_5778.mrc' , 'tv1_relion_data.star' , +2.5 ; ...
 0 , 'rib80s_x0' , 'rib80s' , 1.34 , 'emd_2660.mrc' , 'shiny_2sets.star' , +2.5 ; ...
 1 , 'rib80s' , 'rib80s' , 1.34 , 'emd_2660.mrc' , 'shiny_2sets.star' , +2.5 ; ...
 0 , 'MlaFEDB_x0' , 'MlaFEDB' , 0.832 , 'emd_22116.map' , 'Empiar_10536_00_to_23.star' , +2.5 ; ...
 1 , 'MlaFEDB' , 'MlaFEDB' , 0.832 , 'emd_22116.map' , 'Empiar_10536_00_to_23.star' , +2.5 ; ...
 0 , 'ISWINCP_x0' , 'ISWINCP' , 1.07 , 'emd_9718.map' , 'ADPBeF.star' , +2.5 ; ...
 1 , 'ISWINCP' , 'ISWINCP' , 1.07 , 'emd_9718.map' , 'ADPBeF.star' , +2.5 ; ...
 0 , 'TMEM16F_x0' , 'TMEM16F' , 1.059 , 'emd_20244.map' , 'All_8192.star' , +2.5 ; ...
 1 , 'TMEM16F' , 'TMEM16F' , 1.059 , 'emd_20244.map' , 'All_8192.star' , +2.5 ; ...
 0 , 'LSUbl17dep_x0' , 'LSUbl17dep' , 1.31 , 'emd_8434.map' , 'Parameters.star' , +2.5 ; ...
 1 , 'LSUbl17dep' , 'LSUbl17dep' , 1.31 , 'emd_8434.map' , 'Parameters.star' , +2.5 ; ...
 0 , 'ps1_x0' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' , +2.5 ; ...
 1 , 'ps0' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' , +2.5 ; ...
% 1 , 'ps1' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' , +2.5 ; ...
 0 , 'LetB1_x0' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1_350MB.star' , +2.5 ; ...
 1 , 'LetB1' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1_350MB.star' , +2.5 ; ...
};
n_experiment = size(table_data__,1);
collect__ = cell(n_experiment,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nexperiment=0:n_experiment-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
na=0;
flag_load_type = table_data__{1+nexperiment,1+na}; na=na+1;
fname_prefix = table_data__{1+nexperiment,1+na}; na=na+1;
dir_nopath_data_star = table_data__{1+nexperiment,1+na}; na=na+1;
Pixel_Spacing = table_data__{1+nexperiment,1+na}; na=na+1;
fname_nopath_volume = table_data__{1+nexperiment,1+na}; na=na+1;
fname_nopath_star = table_data__{1+nexperiment,1+na}; na=na+1;
sheres_nls_opt_from_quad = table_data__{1+nexperiment,1+na}; na=na+1; %<-- best nls for bayesian-inference reconstruction starting with a_k_Y_quad_. ;
if (verbose);
disp(sprintf(' %% nexperiment %d/%d (load_type %d): %16s %16s %0.3f %16s %32s' ...
	     ,nexperiment,n_experiment,flag_load_type ...
	     ,fname_prefix,dir_nopath_data_star,Pixel_Spacing,fname_nopath_volume,fname_nopath_star ...
	     ));
end;%if (verbose);
%%%%;
flag_invert = 0;
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); flag_invert = 1; end;
flag_center_image = 0;
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); flag_center_image = 1; end;
flag_crop = 1;
%%%%;
fname_prefix_xfix = sprintf('%s',fname_prefix);
dir_pm = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_pm',string_root,fname_prefix_xfix);
if (~exist(sprintf('%s_mat',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_pm)); mkdir(sprintf('%s_mat',dir_pm)); end;
if (~exist(sprintf('%s_jpg',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_pm)); mkdir(sprintf('%s_jpg',dir_pm)); end;
dir_relion = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_relion',string_root,fname_prefix_xfix);
if (~exist(sprintf('%s_mat',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_relion)); mkdir(sprintf('%s_mat',dir_relion)); end;
if (~exist(sprintf('%s_jpg',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_relion)); mkdir(sprintf('%s_jpg',dir_relion)); end;
dir_data_star = sprintf('/%s/rangan/dir_cryoem/dir_%s',string_root,dir_nopath_data_star);
dir_pm_mat = sprintf('%s_mat',dir_pm);
%%%%%%%%;

parameter_sub = struct('type','parameter');
parameter_sub.flag_store_S_k_p__=0;
parameter_sub.flag_store_M_k_p__=0;
parameter_sub.flag_exclude_ampm=1;
[~,ampm_sub_] = ampm_fromdisk_3(parameter_sub,dir_pm);

corr_crop_reco_x_ = [];
if flag_crop;
fname_mat = sprintf('%s/test_pm_collect_corr_crop_reco_.mat',dir_pm_mat);
if ( flag_recoll | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
tmp_ = load(sprintf('%s/a_k_p_quad_.mat',dir_pm_mat),'n_k_p_r','k_p_r_','k_p_r_max','weight_3d_k_p_r_','a_x_u_reco_');
n_k_p_r = tmp_.n_k_p_r;
k_p_r_ = tmp_.k_p_r_;
k_p_r_max = tmp_.k_p_r_max;
weight_3d_k_p_r_ = tmp_.weight_3d_k_p_r_;
a_x_u_reco_ = tmp_.a_x_u_reco_;
clear tmp_;
tmp_ = load(sprintf('%s/a_k_Y_quad_.mat',dir_pm_mat),'l_max_','a_k_Y_quad_');
l_max_ = tmp_.l_max_;
a_k_Y_quad_ = tmp_.a_k_Y_quad_;
clear tmp_;
if flag_load_type==0;
tmp_ = load(sprintf('%s/a_x_u_base_.mat',dir_pm_mat),'half_diameter_x_c','x_u_0_','x_u_1_','x_u_2_','n_x_u_pack','a_x_u_base_');
a_x_u_base_ = tmp_.a_x_u_base_;
end;%if flag_load_type==0;
if flag_load_type==1;
tmp_ = load(sprintf('%s/a_x_u_pack_.mat',dir_pm_mat),'half_diameter_x_c','x_u_0_','x_u_1_','x_u_2_','n_x_u_pack','a_x_u_pack_');
a_x_u_base_ = tmp_.a_x_u_pack_;
end;%if flag_load_type==1;
half_diameter_x_c = tmp_.half_diameter_x_c;
x_u_0_ = tmp_.x_u_0_; x_u_1_ = tmp_.x_u_1_; x_u_2_ = tmp_.x_u_2_; n_x_u_pack = tmp_.n_x_u_pack;
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;
x_u_r___ = sqrt(x_u_0___.^2 + x_u_1___.^2 + x_u_2___.^2);
clear tmp_;
if flag_load_type==0;
if flag_center_image==0; tmp_ = load(sprintf('%s/a_k_Y_reco_from_M__.mat',dir_pm_mat),'corr_a_k_Y_i_','a_k_Y_reco_yki__'); end;
if flag_center_image==1; tmp_ = load(sprintf('%s/a_k_Y_reco_from_N__.mat',dir_pm_mat),'corr_a_k_Y_i_','a_k_Y_reco_yki__'); end;
corr_a_k_Y_i_ = tmp_.corr_a_k_Y_i_;
corr_a_k_Y = corr_a_k_Y_i_(end);
a_k_Y_reco_ = tmp_.a_k_Y_reco_yki__(:,end);
end;%if flag_load_type==0;
if flag_load_type==1;
tmp_ = load(sprintf('%s/c_k_Y_.mat',dir_pm_mat),'c_k_Y_reco_','X_best_reco');
corr_a_k_Y_i_ = tmp_.X_best_reco;
corr_a_k_Y = tmp_.X_best_reco;
a_k_Y_reco_ = tmp_.c_k_Y_reco_;
end;%if flag_load_type==1;
clear tmp_;
%%%%;
corr_crop_reco_x_ = zeros(n_x_u_pack,1);
sample_sphere_k_eq_d = 1/(2*pi);
a_k_Y_true_yk_ = a_k_Y_quad_;
%%%%;
[ ... 
 b_k_Y_reco_ ...
] = ...
spharm_register_and_rotate_2( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_true_yk_ ...
,a_k_Y_reco_ ...
,0 ...
);
%%%%;
[ ... 
  b_x_u_reco_ ...
] = ...
convert_spharm_to_x_c_3( ...
 sample_sphere_k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,b_k_Y_reco_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
%%%%;
n_crop = n_x_u_pack;
for ncrop=0:n_crop-1;
r_crop = 1.0*ncrop/(n_crop-1);
c_x_u_reco_ = b_x_u_reco_.*reshape(x_u_r___<=r_crop,[n_xxx_u,1]);
corr_crop_reco_x_(1+ncrop) = real(corr(c_x_u_reco_(:),a_x_u_reco_(:)));
end;%for ncrop=0:n_crop-1;
if (verbose); disp(sprintf(' %% max(corr_crop_reco_x_): %0.3f',max(corr_crop_reco_x_))); end;
%%%%;
save(fname_mat,'corr_crop_reco_x_');
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
end;%if flag_crop;

corr_crop_xa__ = [];
if flag_crop;
fname_mat = sprintf('%s/test_pm_collect_corr_crop_.mat',dir_pm_mat);
if ( flag_recoll | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
tmp_ = load(sprintf('%s/a_k_p_quad_.mat',dir_pm_mat),'n_k_p_r','k_p_r_','k_p_r_max','weight_3d_k_p_r_','a_x_u_reco_');
n_k_p_r = tmp_.n_k_p_r;
k_p_r_ = tmp_.k_p_r_;
k_p_r_max = tmp_.k_p_r_max;
weight_3d_k_p_r_ = tmp_.weight_3d_k_p_r_;
a_x_u_reco_ = tmp_.a_x_u_reco_;
clear tmp_;
tmp_ = load(sprintf('%s/a_k_Y_quad_.mat',dir_pm_mat),'l_max_','a_k_Y_quad_');
l_max_ = tmp_.l_max_;
a_k_Y_quad_ = tmp_.a_k_Y_quad_;
clear tmp_;
if flag_load_type==0;
tmp_ = load(sprintf('%s/a_x_u_base_.mat',dir_pm_mat),'half_diameter_x_c','x_u_0_','x_u_1_','x_u_2_','n_x_u_pack','a_x_u_base_');
a_x_u_base_ = tmp_.a_x_u_base_;
end;%if flag_load_type==0;
if flag_load_type==1;
tmp_ = load(sprintf('%s/a_x_u_pack_.mat',dir_pm_mat),'half_diameter_x_c','x_u_0_','x_u_1_','x_u_2_','n_x_u_pack','a_x_u_pack_');
a_x_u_base_ = tmp_.a_x_u_pack_;
end;%if flag_load_type==1;if flag_load_type==1;
half_diameter_x_c = tmp_.half_diameter_x_c;
x_u_0_ = tmp_.x_u_0_; x_u_1_ = tmp_.x_u_1_; x_u_2_ = tmp_.x_u_2_; n_x_u_pack = tmp_.n_x_u_pack;
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;
x_u_r___ = sqrt(x_u_0___.^2 + x_u_1___.^2 + x_u_2___.^2);
clear tmp_;
%%%%;
tmp_X_2d_xcor_d0_a1t_align_ = ls(sprintf('%s/X_2d_xcor_d0_a1t*_align_a_k_Y_.mat',dir_pm_mat));
tmp_index_start_ = strfind(tmp_X_2d_xcor_d0_a1t_align_,sprintf('/%s',string_root))-1;
tmp_index_final_ = strfind(tmp_X_2d_xcor_d0_a1t_align_,'.mat')+4-1;
n_X_2d_xcor_d0_a1t_align = min(numel(tmp_index_start_),numel(tmp_index_final_));
corr_crop_xa__ = zeros(n_x_u_pack,n_X_2d_xcor_d0_a1t_align);
a_k_Y_true_yk_ = a_k_Y_quad_;
sample_sphere_k_eq_d = 1/(2*pi);
%%%%;
for na=0:n_X_2d_xcor_d0_a1t_align-1;
tmp_fname_align_a_k_Y_mat = tmp_X_2d_xcor_d0_a1t_align_(1+tmp_index_start_(1+na)+0:1+tmp_index_final_(1+na)-1);
if (tmp_fname_align_a_k_Y_mat(1)==sprintf('\n')); tmp_fname_align_a_k_Y_mat = tmp_fname_align_a_k_Y_mat(2:end); end;
tmp_fname_a_k_Y_mat = tmp_fname_align_a_k_Y_mat;
tmp_ij = strfind(tmp_fname_a_k_Y_mat,'_align_a_k_Y_.mat');
tmp_fname_a_k_Y_mat = tmp_fname_a_k_Y_mat([1:tmp_ij,tmp_ij+7:end]);
%%%%;
tmp_fname_fsc_crop_kx__pre = tmp_fname_align_a_k_Y_mat;
tmp_ij = strfind(tmp_fname_fsc_crop_kx__pre,'_align_a_k_Y_.mat');
tmp_fname_fsc_crop_kx__pre = tmp_fname_fsc_crop_kx__pre([1:tmp_ij-1]);
tmp_fname_fsc_crop_kx__pre = sprintf('%s_fsc_crop_kx__',tmp_fname_fsc_crop_kx__pre);
tmp_fname_fsc_crop_kx__mat = sprintf('%s.mat',tmp_fname_fsc_crop_kx__pre);
flag_corr_crop = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if  exist(tmp_fname_fsc_crop_kx__mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_ = load(tmp_fname_fsc_crop_kx__mat);
corr_crop_ = zeros(64,1);
if ( isfield(tmp_,'corr_full_reco_vs_crop_ampm_x_')); corr_crop_ = tmp_.corr_full_reco_vs_crop_ampm_x_; flag_corr_crop = 1; end;
clear tmp_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if  exist(tmp_fname_fsc_crop_kx__mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ~exist(tmp_fname_fsc_crop_kx__mat,'file') | ~flag_corr_crop;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_fname_align_crop_a_k_Y_pre = tmp_fname_align_a_k_Y_mat;
tmp_ij = strfind(tmp_fname_align_crop_a_k_Y_pre,'_align_a_k_Y_.mat');
tmp_fname_align_crop_a_k_Y_pre = tmp_fname_align_crop_a_k_Y_pre([1:tmp_ij-1]);
tmp_fname_align_crop_a_k_Y_pre = sprintf('%s_align_crop_a_k_Y_',tmp_fname_align_crop_a_k_Y_pre);
[tmp_fname_align_crop_a_k_Y_skip,tmp_fname_align_crop_a_k_Y_mat] = open_fname_tmp(tmp_fname_align_crop_a_k_Y_pre);
if ~tmp_fname_align_crop_a_k_Y_skip;
if (  exist(tmp_fname_align_a_k_Y_mat,'file') &  exist(tmp_fname_a_k_Y_mat,'file') );
tmp_XB_align_a_k_Y_ = load(tmp_fname_align_a_k_Y_mat);
tmp_XB_a_k_Y_ = load(tmp_fname_a_k_Y_mat);
%%%%;
tmp_X_best = tmp_XB_align_a_k_Y_.X_best_(end);
tmp_flag_flip = tmp_XB_align_a_k_Y_.X_best_flag_flip_(end);
tmp_polar_a_best = tmp_XB_align_a_k_Y_.polar_a_best_(end);
tmp_azimu_b_best = tmp_XB_align_a_k_Y_.azimu_b_best_(end);
tmp_gamma_z_best = tmp_XB_align_a_k_Y_.gamma_z_best_(end);
tmp_delta_x_best = tmp_XB_align_a_k_Y_.delta_best__(1+0,end);
tmp_delta_y_best = tmp_XB_align_a_k_Y_.delta_best__(1+1,end);
tmp_a_k_Y_reco_ = tmp_XB_a_k_Y_.a_k_Y_reco_;
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
n_crop = n_x_u_pack;
corr_crop_ = zeros(n_crop,1);
for ncrop=0:n_crop-1;
r_crop = 1.0*ncrop/(n_crop-1);
tmp_c_x_u_reco_ = tmp_b_x_u_reco_.*reshape(x_u_r___<=r_crop,[n_xxx_u,1]);
corr_crop_(1+ncrop) = real(corr(tmp_c_x_u_reco_(:),a_x_u_reco_(:)));
end;%for ncrop=0:n_crop-1;
if (verbose); disp(sprintf(' %% na %d/%d: %s: max(corr_crop_): %0.3f',na,n_X_2d_xcor_d0_a1t_align,tmp_fname_align_crop_a_k_Y_mat,max(corr_crop_))); end;
%%%%;
save(tmp_fname_align_crop_a_k_Y_mat,'corr_crop_');
%%%%;
end;%if (  exist(tmp_fname_align_a_k_Y_mat,'file') &  exist(tmp_fname_a_k_Y_mat,'file') );
close_fname_tmp(tmp_fname_align_crop_a_k_Y_pre);
end;%if ~tmp_fname_align_crop_a_k_Y_skip;
corr_crop_ = zeros(n_x_u_pack,1);
if (  exist(tmp_fname_align_crop_a_k_Y_mat,'file'));
tmp_ = load(tmp_fname_align_crop_a_k_Y_mat); corr_crop_ = tmp_.corr_crop_; clear tmp_;
end;%if (  exist(tmp_fname_align_crop_a_k_Y_mat,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ~exist(tmp_fname_fsc_crop_kx__mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%;
corr_crop_xa__(:,1+na) = corr_crop_;
end;%for na=0:n_X_2d_xcor_d0_a1t_align-1;
%%%%;
save(fname_mat,'corr_crop_xa__');
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
end;%if flag_crop;

fname_mat = sprintf('%s/test_pm_collect.mat',dir_pm_mat);
if ( flag_recoll | ~exist(fname_mat,'file'));
%%%%%%%%;
disp(sprintf(' %% %s not found, creating',fname_mat));
dir_X_relion = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_relion_mat',string_root,dir_nopath_data_star);
fname_X_relion_mat = sprintf('%s/X_relion_.mat',dir_X_relion);
if (~exist(fname_X_relion_mat,'file'));
fname_X_relion_mat = sprintf('%s/job_1024/X_relion_.mat',dir_X_relion);
end;%if (~exist(fname_X_relion_mat,'file'));
if (~exist(fname_X_relion_mat,'file'));
disp(sprintf(' %% %s not found',fname_X_relion_mat));
end;%if (~exist(fname_X_relion_mat,'file'));
X_relion_ = [0];
if ( exist(fname_X_relion_mat,'file'));
tmp_ = load(fname_X_relion_mat);
X_relion_ = tmp_.X_relion_;
clear tmp_;
end;%if ( exist(fname_X_relion_mat,'file'));
X_relion = X_relion_(end);
%%%%;
if flag_load_type==0;
if flag_center_image==0; tmp_ = load(sprintf('%s/a_k_Y_reco_from_M__.mat',dir_pm_mat),'corr_a_k_Y_i_'); end;
if flag_center_image==1; tmp_ = load(sprintf('%s/a_k_Y_reco_from_N__.mat',dir_pm_mat),'corr_a_k_Y_i_'); end;
corr_a_k_Y_i_ = tmp_.corr_a_k_Y_i_;
corr_a_k_Y = corr_a_k_Y_i_(end);
end;%if flag_load_type==0;
if flag_load_type==1;
tmp_ = load(sprintf('%s/c_k_Y_.mat',dir_pm_mat),'X_best_reco');
corr_a_k_Y_i_ = tmp_.X_best_reco;
corr_a_k_Y = tmp_.X_best_reco;
end;%if flag_load_type==1;
clear tmp_;
tmp_ = load(sprintf('%s/M_k_p__.mat',dir_pm_mat),'n_M');
n_M = tmp_.n_M;
clear tmp_;
tmp_ = load(sprintf('%s/a_k_p_quad_.mat',dir_pm_mat),'n_k_p_r','k_p_r_');
n_k_p_r = tmp_.n_k_p_r;
k_p_r_ = tmp_.k_p_r_;
clear tmp_;
tmp_ = load(sprintf('%s/a_k_Y_quad_.mat',dir_pm_mat),'l_max_');
l_max_ = tmp_.l_max_;
clear tmp_;
if flag_load_type==0;
tmp_ = load(sprintf('%s/CTF_k_p_wkC__.mat',dir_pm_mat),'index_nCTF_from_nM_','CTF_avg_k_p_r_k_');
index_nCTF_from_nM_ = tmp_.index_nCTF_from_nM_;
CTF_avg_k_p_r_k_ = tmp_.CTF_avg_k_p_r_k_;
end;%if flag_load_type==0;
if flag_load_type==1;
tmp_ = load(sprintf('%s/CTF_k_p__.mat',dir_pm_mat),'CTF_index_','CTF_avg_k_p_r_');
index_nCTF_from_nM_ = tmp_.CTF_index_;
CTF_avg_k_p_r_k_ = tmp_.CTF_avg_k_p_r_;
end;%if flag_load_type==1;
CTF_avg_k_p_r_kk__ = CTF_avg_k_p_r_k_*transpose(CTF_avg_k_p_r_k_);
clear tmp_;
tmp_ = load(sprintf('%s/S_k_p__.mat',dir_pm_mat),'weight_2d_k_p_r_');
weight_2d_k_p_r_ = tmp_.weight_2d_k_p_r_;
clear tmp_;
%%%%;
tmp_X_2d_xcor_d0_a1t_align_ = ls(sprintf('%s/X_2d_xcor_d0_a1t*_align_a_k_Y_.mat',dir_pm_mat));
tmp_index_start_ = strfind(tmp_X_2d_xcor_d0_a1t_align_,sprintf('/%s',string_root))-1;
tmp_index_final_ = strfind(tmp_X_2d_xcor_d0_a1t_align_,'.mat')+4-1;
n_X_2d_xcor_d0_a1t_align = min(numel(tmp_index_start_),numel(tmp_index_final_));
X_rand_a_ = zeros(n_X_2d_xcor_d0_a1t_align,1);
X_best_a_ = zeros(n_X_2d_xcor_d0_a1t_align,1);
SX_ka__ = zeros(n_k_p_r,n_X_2d_xcor_d0_a1t_align);
delta_sigma_use_a_ = zeros(n_X_2d_xcor_d0_a1t_align,1);
flag_rank_vs_tolerance_a_ = zeros(n_X_2d_xcor_d0_a1t_align,1);
tolerance_pm_a_ = zeros(n_X_2d_xcor_d0_a1t_align,1);
rank_pm_a_ = zeros(n_X_2d_xcor_d0_a1t_align,1);
rseed_a_ = zeros(n_X_2d_xcor_d0_a1t_align,1);
disc_p25_a_ = zeros(n_X_2d_xcor_d0_a1t_align,1);
disc_p50_a_ = zeros(n_X_2d_xcor_d0_a1t_align,1);
disc_p75_a_ = zeros(n_X_2d_xcor_d0_a1t_align,1);
for na=0:n_X_2d_xcor_d0_a1t_align-1;
tmp_fname_align_a_k_Y_mat = tmp_X_2d_xcor_d0_a1t_align_(1+tmp_index_start_(1+na)+0:1+tmp_index_final_(1+na)-1);
if (tmp_fname_align_a_k_Y_mat(1)==sprintf('\n')); tmp_fname_align_a_k_Y_mat = tmp_fname_align_a_k_Y_mat(2:end); end;
tmp_ = load(tmp_fname_align_a_k_Y_mat,'X_best_'); tmp_X_best_ = tmp_.X_best_; clear tmp_;
tmp_X_best = tmp_X_best_(end);
tmp_fname_X_2d_Memp_d1_align_a_k_Y_mat = tmp_fname_align_a_k_Y_mat;
tmp_ij = strfind(tmp_fname_X_2d_Memp_d1_align_a_k_Y_mat,'xcor_d0');
tmp_fname_X_2d_Memp_d1_align_a_k_Y_mat(tmp_ij+[0:6]) = 'Memp_d1';
tmp_ = load(tmp_fname_X_2d_Memp_d1_align_a_k_Y_mat,'X_best_'); tmp_X_rand = tmp_.X_best_(1); clear tmp_;
tmp_fname_a_k_Y_mat = tmp_fname_align_a_k_Y_mat;
tmp_ij = strfind(tmp_fname_a_k_Y_mat,'_align_a_k_Y_.mat');
tmp_fname_a_k_Y_mat = tmp_fname_a_k_Y_mat([1:tmp_ij,tmp_ij+7:end]);
tmp_ = load(tmp_fname_a_k_Y_mat,'a_k_Y_reco_'); tmp_a_k_Y_reco_ = tmp_.a_k_Y_reco_; clear tmp_;
tmp_fname_compare_image_rank_mat = tmp_fname_align_a_k_Y_mat;
tmp_ij = strfind(tmp_fname_compare_image_rank_mat,'_align_a_k_Y_.mat');
tmp_fname_compare_image_rank_mat = tmp_fname_compare_image_rank_mat([1:tmp_ij]);
tmp_fname_compare_image_rank_mat = sprintf('%scompare_image_rank.mat',tmp_fname_compare_image_rank_mat);
tmp_ = load(tmp_fname_compare_image_rank_mat,'tmp_ij_min_threshold_vs_X_TM_avg_');
tmp_ij_min_threshold_vs_X_TM_avg_ = tmp_.tmp_ij_min_threshold_vs_X_TM_avg_;
%tmp_ij_min_threshold_vs_X_TM_avg_Mi__ = tmp_.tmp_ij_min_threshold_vs_X_TM_avg_Mi__;
tmp_disc_p25 = tmp_ij_min_threshold_vs_X_TM_avg_(1+max(0,min(n_M-1,floor(n_M*1/4))));
tmp_disc_p50 = tmp_ij_min_threshold_vs_X_TM_avg_(1+max(0,min(n_M-1,floor(n_M*2/4))));
tmp_disc_p75 = tmp_ij_min_threshold_vs_X_TM_avg_(1+max(0,min(n_M-1,floor(n_M*3/4))));
clear tmp_;
%%%%;
delta_sigma_base = 0.0;
[ ...
 X_2d_xavg_dx_kk__ ...
,X_2d_xavg_dx_weight_r_ ...
] = ...
principled_marching_cost_matrix_6( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,l_max_ ...
,[] ...
,[] ...
,tmp_a_k_Y_reco_ ...
,CTF_avg_k_p_r_kk__ ...
,delta_sigma_base ...
);
[tmp_UX__,tmp_SX__,tmp_VX__] = svd(X_2d_xavg_dx_kk__); tmp_SX_ = diag(tmp_SX__);
%%%%;
tmp_infix = tmp_fname_align_a_k_Y_mat(strfind(tmp_fname_align_a_k_Y_mat,'X_2d_xcor_d0_a1t'):end);
tmp_infix = tmp_infix(numel('X_2d_xcor_d0_a1t')+1:end);
delta_sigma_use = str2num(tmp_infix(1:4)); delta_sigma_use = delta_sigma_use/1000; tmp_infix = tmp_infix(5:end);
flag_rank_vs_tolerance = 0;
if (tmp_infix(1)=='p');
flag_rank_vs_tolerance = 0;
tolerance_pm = str2num(tmp_infix(2:3)); tolerance_pm = 10.^(-tolerance_pm/10);
rank_pm = max(find(tmp_SX_/max(tmp_SX_)> tolerance_pm));
end;%if (tmp_infix(1)=='p');
if (tmp_infix(1)=='n');
flag_rank_vs_tolerance = 1;
rank_pm = str2num(tmp_infix(2:3));
tolerance_pm = sqrt(tmp_SX_(rank_pm+0)*tmp_SX_(rank_pm+1))/max(tmp_SX_);
end;%if (tmp_infix(1)=='n');
tmp_infix = tmp_infix(4:end);
rseed = str2num(tmp_infix(2));
tmp_infix = tmp_infix(3:end);
assert(strcmp(tmp_infix,'_align_a_k_Y_.mat'));
%%%%;
disp(sprintf(' %% %s: flag_rank_vs_tolerance %d tolerance_pm %0.3f rank_pm %.2d rseed %d X_rand %0.3f X_best %0.3f',tmp_fname_align_a_k_Y_mat,flag_rank_vs_tolerance,tolerance_pm,rank_pm,rseed,tmp_X_rand,tmp_X_best));
%%%%;
X_rand_a_(1+na)=tmp_X_rand;
X_best_a_(1+na)=tmp_X_best;
SX_ka__(:,1+na)=tmp_SX_;
delta_sigma_use_a_(1+na)=delta_sigma_use;
flag_rank_vs_tolerance_a_(1+na)=flag_rank_vs_tolerance;
tolerance_pm_a_(1+na)=tolerance_pm;
rank_pm_a_(1+na)=rank_pm;
rseed_a_(1+na)=rseed;
disc_p25_a_(1+na)=tmp_disc_p25;
disc_p50_a_(1+na)=tmp_disc_p50;
disc_p75_a_(1+na)=tmp_disc_p75;
%%%%;
end;%for na=0:n_X_2d_xcor_d0_a1t_align-1;
%%%%;
save(fname_mat ...
,'fname_prefix' ...
,'dir_nopath_data_star' ...
,'Pixel_Spacing' ...
,'fname_nopath_volume' ...
,'fname_nopath_star' ...
,'flag_invert' ...
,'flag_center_image' ...
,'dir_pm' ...
,'dir_relion' ...
,'dir_data_star' ...
,'dir_pm_mat' ...
,'dir_X_relion' ...
,'fname_X_relion_mat' ...
,'X_relion_' ...
,'X_relion' ...
,'corr_a_k_Y_i_' ...
,'corr_a_k_Y' ...
,'n_M' ...
,'n_k_p_r' ...
,'k_p_r_' ...
,'l_max_' ...
,'index_nCTF_from_nM_' ...
,'CTF_avg_k_p_r_k_' ...
,'weight_2d_k_p_r_' ...
,'X_rand_a_' ...
,'X_best_a_' ...
,'SX_ka__' ...
,'delta_sigma_use_a_' ...
,'flag_rank_vs_tolerance_a_' ...
,'tolerance_pm_a_' ...
,'rank_pm_a_' ...
,'rseed_a_' ...
,'disc_p25_a_' ...
,'disc_p50_a_' ...
,'disc_p75_a_' ...
);
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
collect__{1+nexperiment} = load(fname_mat);
end;%if ( exist(fname_mat,'file'));

fsc_crop_reco_stab_alig_from_x__ = [];
if flag_center_image==1; tmp_str = 'N'; end; if flag_center_image==0; tmp_str = 'M'; end;
str_fsc_crop_reco_stab_alig_from_x__ = sprintf('%s/fsc_crop_reco_stab_alig_from_%s__.mat',dir_pm_mat,tmp_str);
if ( exist(str_fsc_crop_reco_stab_alig_from_x__,'file'));
fsc_crop_reco_stab_alig_from_x__ = load(str_fsc_crop_reco_stab_alig_from_x__);
end;%if ( exist(str_fsc_crop_reco_stab_alig_from_x__,'file'));

rseed_ = 0:2; n_rseed = numel(rseed_);
sigma_sheres_ = exp([-3.0:0.5:-1.0]); n_sigma_sheres = numel(sigma_sheres_);
fsc_sher_rs__ = cell(n_rseed,n_sigma_sheres);
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
for nsigma_sheres=0:n_sigma_sheres-1;
fsc_sher_rs__{1+nrseed,1+nsigma_sheres} = struct('type','fsc');
sigma_sheres = sigma_sheres_(1+nsigma_sheres);
str_nls = sprintf('nls%.2d',round(10*-log(sigma_sheres)));
str_rseed = sprintf('r%d',rseed);
tmp_str_infix = sprintf('X_2d_Memp_d1_%st0100p20%s',str_nls,str_rseed);
tmp_fname_mat = sprintf('%s/%s_fsc.mat',dir_pm_mat,tmp_str_infix);
if  exist(tmp_fname_mat,'file');
if (verbose>0); disp(sprintf(' %% %s found',tmp_fname_mat)); end;
tmp_ = load(tmp_fname_mat);
fsc_sher_rs__{1+nrseed,1+nsigma_sheres}.corr_full_reco_vs_crop_reco_stab_x_ = tmp_.corr_full_reco_vs_crop_reco_stab_x_;
clear tmp_;
end;%if  exist(tmp_fname_mat,'file');
end;%for nsigma_sheres=0:n_sigma_sheres-1;
end;%for nrseed=0:n_rseed-1;

collect__{1+nexperiment}.sheres_nls_opt_from_quad = sheres_nls_opt_from_quad;
collect__{1+nexperiment}.corr_crop_xa__ = [];
collect__{1+nexperiment}.corr_crop_reco_x_ = [];
if flag_crop; collect__{1+nexperiment}.corr_crop_xa__ = corr_crop_xa__; end;
if flag_crop; collect__{1+nexperiment}.corr_crop_reco_x_ = corr_crop_reco_x_; end;
collect__{1+nexperiment}.fsc_crop_reco_stab_alig_from_x__ = fsc_crop_reco_stab_alig_from_x__;
collect__{1+nexperiment}.fsc_sher_rs__ = fsc_sher_rs__;
collect__{1+nexperiment}.ampm_sub_ = ampm_sub_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nexperiment=0:n_experiment-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% accumulate across molecules. ;
%%%%%%%%;
dir_nopath_data_star_ = table_data__(:,1+2);
u_dir_nopath_data_star_ = unique(dir_nopath_data_star_);
u_dir_nopath_data_star_ = cat( 1 , setdiff(u_dir_nopath_data_star_,'p28hRPT1') , 'p28hRPT1' );
n_u_experiment = numel(u_dir_nopath_data_star_);
index_nexperiment_from_nu__ = cell(n_u_experiment,1);
for nu_experiment=0:n_u_experiment-1;
dir_nopath_data_star = u_dir_nopath_data_star_{1+nu_experiment};
index_nexperiment_from_nu__{1+nu_experiment} = efind(strcmp(dir_nopath_data_star_,dir_nopath_data_star));
end;%for nu_experiment=0:n_u_experiment-1;

ylim_ = [-0.125,1.00];
fontsize_use = 16;
linewidth_use = 0.25;
tmp_horiz_gap = 0.125;
facecolor_reli_ = [0.00,0.00,0.85];
facecolor_sher_ = [0.00,0.85,0.85];
facecolor_ampm_ = [0.85,0.00,0.85];
facecolor_corr_OAM_unst_ = [0.15,0.15,0.15];
facecolor_corr_OAM_stab_ = [0.35,0.85,0.00];
facecolor_corr_OBI_stab_ = [0.00,0.50,0.50];

error('stopping for now');

%%%%%%%%%%%%%%%%;
for nu_experiment=0:n_u_experiment-1-1;%for nu_experiment=0:n_u_experiment-1;
index_nexperiment_from_nu_ = index_nexperiment_from_nu__{1+nu_experiment};
%%%%%%%%%%%%%%%%;

u_dir_nopath_data_star = u_dir_nopath_data_star_{1+nu_experiment};
flag_rank_vs_tolerance_use = 1-strcmp(u_dir_nopath_data_star,'p28hRPT1');
%%%%%%%%;
tmp_collect_corr_OAM_unst = -Inf;
tmp_collect_rand = -Inf;
tmp_collect_reli = -Inf;
tmp_collect_corr_OAM_stab = -Inf;
tmp_collect_corr_OBI_stab = -Inf;
na_ampm_tot=0;
tmp_X_ampm_a_ = [];
tmp_tolerance_pm_a_ = [];
tmp_delta_sigma_use_a_ = [];
tmp_flag_rank_vs_tolerance_a_ = [];
for nl=0:numel(index_nexperiment_from_nu_)-1;
nexperiment = index_nexperiment_from_nu_(1+nl);
%%%%%%%%;
collect_ = collect__{1+nexperiment};
tmp_X_ = collect_.X_best_a_;
if (~isempty(collect_.corr_crop_xa__));
tmp_X_ = collect_.corr_crop_xa__;
[~,tmp_index_crop] = max(mean(tmp_X_,2)); tmp_index_crop = tmp_index_crop-1;
if (verbose); disp(sprintf(' %% %s: crop at index %d',collect_.dir_nopath_data_star,tmp_index_crop)); end;
tmp_X_ = tmp_X_(1+tmp_index_crop,:);
end;%if (~isempty(collect_.corr_crop_xa__));
%%%%%%%%;
n_a = numel(collect_.tolerance_pm_a_);
tmp_index_ = na_ampm_tot + [0:n_a-1];
tmp_X_ampm_a_(1 + tmp_index_) = tmp_X_;
tmp_tolerance_pm_a_(1 + tmp_index_) = collect_.tolerance_pm_a_;
tmp_delta_sigma_use_a_(1 + tmp_index_) = collect_.delta_sigma_use_a_;
tmp_flag_rank_vs_tolerance_a_(1 + tmp_index_) = collect_.flag_rank_vs_tolerance_a_;
%%%%;
if ( isempty(collect_.corr_crop_reco_x_)); tmp_collect_corr_OAM_unst = max(tmp_collect_corr_OAM_unst,real(collect_.corr_a_k_Y)); end;
if (~isempty(collect_.corr_crop_reco_x_)); tmp_collect_corr_OAM_unst = max(tmp_collect_corr_OAM_unst,real(collect_.corr_crop_reco_x_(1+tmp_index_crop))); end;
tmp_collect_rand = max(tmp_collect_rand,mean(real(collect_.X_rand_a_)));
tmp_collect_reli = max(tmp_collect_reli,collect_.X_relion);
%%%%%%%%;
tmp_c = -Inf;
if (~isempty(collect_.fsc_crop_reco_stab_alig_from_x__));
tmp_c = collect_.fsc_crop_reco_stab_alig_from_x__.corr_full_reco_vs_crop_reco_stab_x_(1+tmp_index_crop);
end;%if (~isempty(collect_.fsc_crop_reco_stab_alig_from_x__));
tmp_collect_corr_OAM_stab = max(tmp_collect_corr_OAM_stab,tmp_c);
%%%%%%%%;
tmp_c = -Inf;
if ~isempty(collect_.ampm_sub_);
if ~isempty(collect_.ampm_sub_.sher_from_quad_);
if ~isempty(collect_.ampm_sub_.sher_from_quad_.corr_full_reco_vs_crop_reco_stab_x_);
tmp_c = collect_.ampm_sub_.sher_from_quad_.corr_full_reco_vs_crop_reco_stab_x_(1+tmp_index_crop);
end;%if ~isempty(collect_.ampm_sub_.sher_from_quad_.corr_full_reco_vs_crop_reco_stab_x_);
end;%if ~isempty(collect_.ampm_sub_.sher_from_quad_);
end;%if ~isempty(collect_.ampm_sub_);
tmp_collect_corr_OBI_stab = max(tmp_collect_corr_OBI_stab,tmp_c);
%%%%%%%%;
na_ampm_tot = na_ampm_tot + n_a;
end;%for nl=0:numel(index_nexperiment_from_nu_)-1;
n_a_ampm_tot = na_ampm_tot;
%%%%%%%%;
na_sher_tot=0;
tmp_X_sher_a_=[];
tmp_nls_sher_a_=[];
tmp_nls_quad_a_=[];
for nl=0:numel(index_nexperiment_from_nu_)-1;
nexperiment = index_nexperiment_from_nu_(1+nl);
collect_ = collect__{1+nexperiment};
fsc_sher_rs__ = collect_.fsc_sher_rs__;
for nsigma_sheres=0:n_sigma_sheres-1;
sigma_sheres = sigma_sheres_(1+nsigma_sheres);
lsigma_sheres = log(sigma_sheres);
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
if isfield(fsc_sher_rs__{1+nrseed,1+nsigma_sheres},'corr_full_reco_vs_crop_reco_stab_x_');
tmp_X = fsc_sher_rs__{1+nrseed,1+nsigma_sheres}.corr_full_reco_vs_crop_reco_stab_x_(1+tmp_index_crop);
tmp_X_sher_a_(1+na_sher_tot) = tmp_X;
tmp_nls_sher_a_(1+na_sher_tot) = -lsigma_sheres;
tmp_nls_quad_a_(1+na_sher_tot) = collect__{1+nexperiment}.sheres_nls_opt_from_quad;
na_sher_tot = na_sher_tot+1;
end;%if isfield(fsc_sher_rs__{1+nrseed,1+nsigma_sheres},'corr_full_reco_vs_crop_reco_stab_x_');
end;%for nrseed=0:n_rseed-1;
end;%for nsigma_sheres=0:n_sigma_sheres-1;
end;%for nl=0:numel(index_nexperiment_from_nu_)-1;
n_a_sher_tot = na_sher_tot;
%%%%%%%%;
tmp_index_ampm_ = efind( ...
 tmp_X_ampm_a_~=0 ...
& tmp_tolerance_pm_a_<=10^-1.75 ...
& tmp_delta_sigma_use_a_> 0.00 ...
& tmp_flag_rank_vs_tolerance_a_==flag_rank_vs_tolerance_use ...
);
n_a_ampm_use = numel(tmp_index_ampm_);
tmp_index_sher_ = efind( ...
 tmp_X_sher_a_~=0 ...
& abs(tmp_nls_sher_a_-tmp_nls_quad_a_)<=1e-3 ...
);
n_a_sher_use = numel(tmp_index_sher_);
%%%%%%%%;
if (verbose);
disp(sprintf(' %% nu_experiment %d/%d: %s',nu_experiment,n_u_experiment,u_dir_nopath_data_star));
disp(sprintf(' %% %% tmp_collect_corr_OAM_unst %0.6f',tmp_collect_corr_OAM_unst));
disp(sprintf(' %% %% tmp_collect_rand %0.6f',tmp_collect_rand));
disp(sprintf(' %% %% tmp_collect_reli %0.6f',tmp_collect_reli));
disp(sprintf(' %% %% tmp_collect_corr_OAM_stab %0.6f',tmp_collect_corr_OAM_stab));
disp(sprintf(' %% %% tmp_collect_corr_OBI_stab %0.6f',tmp_collect_corr_OBI_stab));
disp(sprintf(' %% %% n_a_ampm_tot %d <-- use %d',n_a_ampm_tot,n_a_ampm_use));
disp(sprintf(' %% %% n_a_sher_tot %d <-- use %d',n_a_sher_tot,n_a_sher_use));
end;%if (verbose);

tmp_X_ampm_use_srt_ = sort(tmp_X_ampm_a_(1+tmp_index_ampm_),'ascend');
tmp_X_sher_use_srt_ = sort(tmp_X_sher_a_(1+tmp_index_sher_),'ascend');

fname_fig_pre = sprintf('%s/test_pm_collect_%s_FIGJ',dir_jpg,u_dir_nopath_data_star);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
figure(1+nf);nf=nf+1;clf;figsml;
nh=0;nx=0;
hold on;
xtick_pre = 0;
xtick_min_ = [];
xtick_max_ = [];
str_xtick_ = {};
%%%%;
tmp_width = 1; tmp_horiz = xtick_pre; tmp_verti = tmp_collect_reli;
patch(tmp_horiz + tmp_width*[0;0;1;1;0],tmp_verti*[0;1;1;0;0],0,'LineWidth',linewidth_use,'EdgeColor','k','FaceColor',facecolor_reli_);
str_xtick_{1+nx} = 'RDN';
xtick_min_(1+nx) = xtick_pre;
xtick_max_(1+nx) = xtick_pre + tmp_width;
xtick_pre = xtick_max_(1+nx) + tmp_horiz_gap;
nx=nx+1;nh = nh+1;
%%%%;
tmp_width = 1.5*1/max(1,n_a_sher_use);
for na=0:n_a_sher_use-1;
tmp_horiz = xtick_pre + na*tmp_width; tmp_verti = tmp_X_sher_use_srt_(1+na);
patch(tmp_horiz + tmp_width*[0;0;1;1;0],tmp_verti*[0;1;1;0;0],0,'LineWidth',linewidth_use,'EdgeColor','k','FaceColor',facecolor_sher_);
end;%for na=0:n_a_sher_use-1;
str_xtick_{1+nx} = 'BI';
xtick_min_(1+nx) = xtick_pre;
xtick_max_(1+nx) = xtick_pre + (1+na)*tmp_width;
xtick_pre = xtick_max_(1+nx) + tmp_horiz_gap;
nx=nx+1;nh = nh+tmp_width*n_a_sher_use;
%%%%;
tmp_width = 1; tmp_horiz = xtick_pre; tmp_verti = tmp_collect_corr_OBI_stab;
patch(tmp_horiz + tmp_width*[0;0;1;1;0],tmp_verti*[0;1;1;0;0],0,'LineWidth',linewidth_use,'EdgeColor','k','FaceColor',facecolor_corr_OBI_stab_);
str_xtick_{1+nx} = 'OBIS';
xtick_min_(1+nx) = xtick_pre;
xtick_max_(1+nx) = xtick_pre + tmp_width;
xtick_pre = xtick_max_(1+nx) + tmp_horiz_gap;
nx=nx+1;nh = nh+1;
%%%%;
tmp_width = 6*1/max(1,n_a_ampm_use);
for na=0:n_a_ampm_use-1;
tmp_horiz = xtick_pre + na*tmp_width; tmp_verti = tmp_X_ampm_use_srt_(1+na);
patch(tmp_horiz + tmp_width*[0;0;1;1;0],tmp_verti*[0;1;1;0;0],0,'LineWidth',linewidth_use,'EdgeColor','k','FaceColor',facecolor_ampm_);
end;%for na=0:n_a_ampm_use-1;
str_xtick_{1+nx} = 'AMPM';
xtick_min_(1+nx) = xtick_pre;
xtick_max_(1+nx) = xtick_pre + (1+na)*tmp_width;
xtick_pre = xtick_max_(1+nx) + tmp_horiz_gap;
nx=nx+1;nh = nh+tmp_width*n_a_ampm_use;
%%%%;
tmp_width = 1; tmp_horiz = xtick_pre; tmp_verti = tmp_collect_corr_OAM_stab;
patch(tmp_horiz + tmp_width*[0;0;1;1;0],tmp_verti*[0;1;1;0;0],0,'LineWidth',linewidth_use,'EdgeColor','k','FaceColor',facecolor_corr_OAM_stab_);
str_xtick_{1+nx} = 'OAMS';
xtick_min_(1+nx) = xtick_pre;
xtick_max_(1+nx) = xtick_pre + tmp_width;
xtick_pre = xtick_max_(1+nx) + tmp_horiz_gap;
nx=nx+1;nh = nh+1;
%%%%;
tmp_width = 1; tmp_horiz = xtick_pre; tmp_verti = tmp_collect_corr_OAM_unst;
patch(tmp_horiz + tmp_width*[0;0;1;1;0],tmp_verti*[0;1;1;0;0],0,'LineWidth',linewidth_use,'EdgeColor','k','FaceColor',facecolor_corr_OAM_unst_);
str_xtick_{1+nx} = 'OAM';
xtick_min_(1+nx) = xtick_pre;
xtick_max_(1+nx) = xtick_pre + tmp_width;
xtick_pre = xtick_max_(1+nx) + tmp_horiz_gap;
nx=nx+1;nh = nh+1;
%%%%;
hold off;
xtick_mid_ = (xtick_min_ + xtick_max_)/2;
ylabel('correlation (masked)'); %xlabel('trial');
xlim_ = [0-tmp_horiz_gap,xtick_pre+tmp_horiz_gap]; xlim(xlim_);
ylim(ylim_); set(gca,'YTick',-0.1:+0.1:1.0);
title(sprintf('%s',u_dir_nopath_data_star),'Interpreter','none');
set(gca,'TickLength',[0,0]);
set(gca,'XTickLabel',str_xtick_,'XTick',xtick_mid_); xtickangle(90);
set(gca,'FontSize',fontsize_use);
sgtitle(fname_fig_pre,'Interpreter','none');
grid on; grid minor;
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
sgtitle('');
tmp_dir = sprintf('/%s/rangan/dir_cryoem/dir_ampm_manuscript/dir_ampm_fig_collect',string_root);
fname_fig_jpg_strip = sprintf('%s/test_pm_collect_%s_FIGJ_strip.jpg',tmp_dir,u_dir_nopath_data_star);
disp(sprintf(' %% writing %s',fname_fig_jpg_strip));
print('-djpeg',sprintf('%s',fname_fig_jpg_strip));
%close(gcf);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));

%%%%%%%%%%%%%%%%;
end;%for nu_experiment=0:n_u_experiment-1;
%%%%%%%%%%%%%%%%;

disp('returning');return;

%%%%%%%%;
% Now calculate crosscorrelations between reconstructed volumes. ;
%%%%%%%%;
%%%%%%%%%%%%%%%%;
for nu_experiment=0:n_u_experiment-1-1;%for nu_experiment=0:n_u_experiment-1;
index_nexperiment_from_nu_ = index_nexperiment_from_nu__{1+nu_experiment};
%%%%%%%%%%%%%%%%;
u_dir_nopath_data_star = u_dir_nopath_data_star_{1+nu_experiment};
fname_mat_pre = sprintf('%s/test_pm_collect_%s_ampm_alig_crosscorrelation',dir_mat,u_dir_nopath_data_star);
fname_mat = sprintf('%s.mat',fname_mat_pre);
%%%%%%%%%%%%;
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%%%%%;
flag_rank_vs_tolerance_use = 1-strcmp(u_dir_nopath_data_star,'p28hRPT1');
n_experiment_local = numel(index_nexperiment_from_nu_);
disp(sprintf(' %% nu_experiment %d/%d u_dir_nopath_data_star %s flag_rank_vs_tolerance_use %d n_experiment_local %d',nu_experiment,n_u_experiment,u_dir_nopath_data_star,flag_rank_vs_tolerance_use,n_experiment_local));
%%%%;
na_ampm_tot=0;
tmp_X_ampm_a_ = [];
tmp_tolerance_pm_a_ = [];
tmp_delta_sigma_use_a_ = [];
tmp_flag_rank_vs_tolerance_a_ = [];
ampm_local__ = cell(n_experiment_local,1);
%%%%%%%%;
for nl=0:n_experiment_local-1;
nexperiment = index_nexperiment_from_nu_(1+nl);
%%%%;
collect_ = collect__{1+nexperiment};
tmp_X_ = collect_.X_best_a_;
if (~isempty(collect_.corr_crop_xa__));
tmp_X_ = collect_.corr_crop_xa__;
[~,tmp_index_crop] = max(mean(tmp_X_,2)); tmp_index_crop = tmp_index_crop-1;
if (verbose); disp(sprintf(' %% %s: crop at index %d',collect_.dir_nopath_data_star,tmp_index_crop)); end;
tmp_X_ = tmp_X_(1+tmp_index_crop,:);
end;%if (~isempty(collect_.corr_crop_xa__));
%%%%;
n_a = numel(collect_.tolerance_pm_a_);
tmp_index_ = na_ampm_tot + [0:n_a-1];
tmp_X_ampm_a_(1 + tmp_index_) = tmp_X_;
tmp_tolerance_pm_a_(1 + tmp_index_) = collect_.tolerance_pm_a_;
tmp_delta_sigma_use_a_(1 + tmp_index_) = collect_.delta_sigma_use_a_;
tmp_flag_rank_vs_tolerance_a_(1 + tmp_index_) = collect_.flag_rank_vs_tolerance_a_;
%%%%;
na_ampm_tot = na_ampm_tot + n_a;
%%%%;
[~,ampm_local__{1+nl}] = ampm_fromdisk_3(struct('type','parameter','flag_store_S_k_p__',0,'flag_store_M_k_p__',0),collect_.dir_pm);
%%%%;
end;%for nl=0:n_experiment_local-1;
n_a_ampm_tot = na_ampm_tot;
%%%%%%%%;
tmp_index_ampm_ = efind( ...
 tmp_X_ampm_a_~=0 ...
& tmp_tolerance_pm_a_<=10^-1.75 ...
& tmp_delta_sigma_use_a_> 0.00 ...
& tmp_flag_rank_vs_tolerance_a_==flag_rank_vs_tolerance_use ...
);
n_a_ampm_use = numel(tmp_index_ampm_);
%%%%%%%%;
if (verbose);
disp(sprintf(' %% nu_experiment %d/%d: %s',nu_experiment,n_u_experiment,u_dir_nopath_data_star));
disp(sprintf(' %% %% n_a_ampm_tot %d <-- use %d',n_a_ampm_tot,n_a_ampm_use));
end;%if (verbose);
%%%%%%%%;
tmp_a_k_Y_ampm_alig_yka__ = zeros(ampm_local__{1}.n_lm_sum,0);
na=0;
for nl=0:n_experiment_local-1;
nexperiment = index_nexperiment_from_nu_(1+nl);
tmp_ampm_local_ = ampm_local__{1+nl};
tmp_n_a = numel(tmp_ampm_local_.str_fname_mat_a_);
for tmp_na=0:tmp_n_a-1;
if (mod(tmp_na,8)==0); disp(sprintf(' %% tmp_na %d/%d',tmp_na,tmp_n_a)); end;
[ ... 
tmp_a_k_Y_ampm_alig_yka__(:,1+na) ...
] = ...
spharm_register_and_rotate_2( ...
 tmp_ampm_local_.n_k_p_r ...
,tmp_ampm_local_.k_p_r_ ...
,tmp_ampm_local_.k_p_r_max ...
,tmp_ampm_local_.weight_3d_k_p_r_ ...
,tmp_ampm_local_.l_max_ ...
,tmp_ampm_local_.a_k_Y_quad_ ...
,tmp_ampm_local_.a_k_Y_ampm_yka__(:,1+tmp_na) ...
,0 ...
,tmp_ampm_local_.X_best_ampm_ia__(end,1+tmp_na) ...
,tmp_ampm_local_.X_best_flag_flip_ampm_ia__(end,1+tmp_na) ...
,tmp_ampm_local_.polar_a_best_ampm_ia__(end,1+tmp_na) ...
,tmp_ampm_local_.azimu_b_best_ampm_ia__(end,1+tmp_na) ...
,tmp_ampm_local_.gamma_z_best_ampm_ia__(end,1+tmp_na) ...
,tmp_ampm_local_.delta_best_ampm_dia__(1+0,end,1+tmp_na) ...
,tmp_ampm_local_.delta_best_ampm_dia__(1+1,end,1+tmp_na) ...
,tmp_ampm_local_.delta_best_ampm_dia__(1+2,end,1+tmp_na) ...
);
na=na+1;
end;%for tmp_na=0:n_a-1;
end;%for nl=0:n_experiment_local-1;
%%%%%%%%;
[tmp_X_ampm_use_srt_,tmp_ij_X_ampm_use_srt_] = sort(tmp_X_ampm_a_(1+tmp_index_ampm_),'ascend');
tmp_index_use_ = tmp_index_ampm_(tmp_ij_X_ampm_use_srt_); n_index_use = numel(tmp_index_use_);
tmp_a_k_Y_ampm_alig_use_yka__ = tmp_a_k_Y_ampm_alig_yka__(:,1+tmp_index_use_);
%%%%%%%%;
% Now perform stripped down calculation of innerproducts. ;
%%%%%%%%;
n_k_p_r = tmp_ampm_local_.n_k_p_r;
k_p_r_ = tmp_ampm_local_.k_p_r_;
weight_3d_k_p_r_ = tmp_ampm_local_.weight_3d_k_p_r_;
n_w_max = tmp_ampm_local_.n_w_max;
n_m_max = tmp_ampm_local_.n_m_max;
l_max_ = tmp_ampm_local_.l_max_; l_max_max = max(l_max_); n_l_max = l_max_max + 1;
%%%%%%%%;
% calculate d_. ;
%%%%%%%%;
n_beta = n_w_max; beta_ = transpose(linspace(-pi,+pi,n_beta+1)); beta_ = beta_(1:end-1);
t_0in = tic;
d_mmlb____ = zeros(n_m_max,n_m_max,n_l_max,n_beta);
for nbeta=0:n_beta-1;
if (verbose); if (mod(nbeta,16)==0); disp(sprintf(' %% nbeta %d/%d',nbeta,n_beta)); end; end;
beta = beta_(1+nbeta);
W_ = wignerd_b(l_max_max,-beta);
n_W_ = zeros(1,1+l_max_max); for (l_val=0:l_max_max); n_W_(1+l_val) = numel(W_{1+l_val}); end;
d_mml___ = zeros(n_m_max,n_m_max,1+l_max_max);
for l_val=0:l_max_max;
d_mml___(1+l_max_max + [-l_val:+l_val],1+l_max_max + [-l_val:+l_val],1+l_val) = W_{1+l_val};
end;%for l_val=0:l_max_max;
d_mmlb____(:,:,:,1+nbeta) = d_mml___;
end;%for nbeta=0:n_beta-1;
t_out = toc(t_0in); if (verbose>1); disp(sprintf(' %% calculate d_: t %0.6f',t_out)); end;
%%%%%%%%;
t_0in = tic;
tmp_a_k_Y_ampm_alig_use_mkla____ = zeros(n_m_max,n_k_p_r,n_l_max,n_index_use);
for nindex_use=0:n_index_use-1;
tmp_a_k_Y_mlk_ = tmp_a_k_Y_ampm_alig_use_yka__(:,1+nindex_use);
tmp_a_k_Y_mlk_ = spharm_normalize_2(n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,tmp_a_k_Y_mlk_);
tmp_a_k_Y_mlk___ = convert_spharm_to_spharm___0(n_k_p_r,l_max_,tmp_a_k_Y_mlk_);
tmp_a_k_Y_mkl___ = permute(tmp_a_k_Y_mlk___,[1,3,2]);
tmp_a_k_Y_ampm_alig_use_mkla____(:,:,:,1+nindex_use) = tmp_a_k_Y_mkl___;
clear tmp_a_k_Y_mlk_ tmp_a_k_Y_mlk___ tmp_a_k_Y_mkl___;
end;%for nindex_use=0:n_index_use-1;
t_out = toc(t_0in); if (verbose>1); disp(sprintf(' %% calculate tmp_a_k_Y_ampm_alig_use_mkla_____: t %0.6f',t_out)); end;
%%%%%%%%;
t_0in = tic;
n_sum = n_index_use*(n_index_use-1)/2; nsum=0;
tmp_XX_use__ = zeros(n_index_use,n_index_use);
for nindex_use_0=0:n_index_use-1;
tmp_a_k_Y_0_yk_ = tmp_a_k_Y_ampm_alig_use_yka__(:,1+nindex_use_0);
tmp_a_k_Y_0_mkl___ = tmp_a_k_Y_ampm_alig_use_mkla____(:,:,:,1+nindex_use_0);
for nindex_use_1=nindex_use_0+0:n_index_use-1;
tmp_a_k_Y_1_yk_ = tmp_a_k_Y_ampm_alig_use_yka__(:,1+nindex_use_1);
tmp_a_k_Y_1_mkl___ = tmp_a_k_Y_ampm_alig_use_mkla____(:,:,:,1+nindex_use_1);
if (mod(nsum,16)==0); if (verbose); disp(sprintf(' %% nsum %.6d/%.6d (%.3d,%.3d)',nsum,n_sum,nindex_use_0,nindex_use_1)); end; end;
[ ...
 tmp_X0_mmb___ ...
] = ...
register_spharm_to_spharm_single_beta_3_stripped_0( ...
 verbose ...
,n_m_max ...
,n_k_p_r ...
,n_l_max ...
,weight_3d_k_p_r_ ...
,tmp_a_k_Y_0_mkl___ ...
,tmp_a_k_Y_1_mkl___ ...
,[] ...
,[] ...
,n_beta ...
,d_mmlb____ ...
);
tmp_X0_mmb___ = real(tmp_X0_mmb___);
tmp_XX_use__(1+nindex_use_0,1+nindex_use_1) = max(tmp_X0_mmb___,[],'all');
tmp_XX_use__(1+nindex_use_1,1+nindex_use_0) = tmp_XX_use__(1+nindex_use_0,1+nindex_use_1);
nsum = nsum+1;
end;%for nindex_use_1=nindex_use_0+0:n_index_use-1;
end;%for nindex_use_0=0:n_index_use-1;
t_out = toc(t_0in); if (verbose>1); disp(sprintf(' %% calculate tmp_X0_mmb___: t %0.6f',t_out)); end;
%%%%%%%%;
tmp_YY__ = real(corr(tmp_a_k_Y_ampm_alig_yka__));
tmp_YY_use__ = tmp_YY__(1+tmp_index_use_,1+tmp_index_use_);
%%%%%%%%;
save(fname_mat ...
     ,'n_a_ampm_tot' ...
     ,'flag_rank_vs_tolerance_use' ...
     ,'tmp_index_crop' ...
     ,'tmp_X_ampm_a_' ...
     ,'tmp_tolerance_pm_a_' ...
     ,'tmp_delta_sigma_use_a_' ...
     ,'tmp_flag_rank_vs_tolerance_a_' ...
     ,'tmp_index_ampm_' ...
     ,'tmp_X_ampm_use_srt_' ...
     ,'tmp_ij_X_ampm_use_srt_' ...
     ,'tmp_index_use_' ...
     ,'tmp_XX_use__' ...
     ,'tmp_YY__','tmp_YY_use__' ...
     );
%%%%%%%%%%%%;
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%%%%%;
if  exist(fname_mat,'file');
%%%%%%%%%%%%;
load(fname_mat);
%%%%%%%%;
fname_fig_pre = sprintf('%s/test_pm_collect_%s_ampm_alig_crosscorrelation_FIGK',dir_jpg,u_dir_nopath_data_star);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot>1 | ~exist(fname_fig_jpg,'file'));
figure(1+nf);nf=nf+1;clf;figsml;figbeach();
fontsize_use = 12;
%tmp_YY_use__ = tmp_YY__(1+tmp_index_use_,1+tmp_index_use_);
tmp_index_offdiagonal_ = efind(1-eye(numel(tmp_index_use_)));
tmp_YY_med = median(tmp_YY_use__(1+tmp_index_offdiagonal_));
tmp_YY_avg = mean(tmp_YY_use__(1+tmp_index_offdiagonal_));
imagesc(tmp_YY_use__,[0,1]);
axisnotick; axis image; tmp_c_ = colorbar; set(tmp_c_,'Ticks',[0,0.5,1]);
title(sprintf('%s (%0.2f) [%0.2f]',u_dir_nopath_data_star,tmp_YY_avg,tmp_YY_med),'Interpreter','none');
sgtitle(fname_fig_pre,'Interpreter','none');
set(gca,'FontSize',fontsize_use);
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (flag_replot>1 | ~exist(fname_fig_jpg,'file'));
%%%%%%%%;
fname_fig_pre = sprintf('%s/test_pm_collect_%s_ampm_alig_crosscorrelation_FIGL',dir_jpg,u_dir_nopath_data_star);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
figure(1+nf);nf=nf+1;clf;figsml;figbeach();
fontsize_use = 24;
%tmp_XX_use__ = tmp_XX__(1+tmp_index_use_,1+tmp_index_use_);
tmp_index_0ondiagonal_ = efind(0+eye(numel(tmp_index_use_)));
tmp_XX_use__(1+tmp_index_0ondiagonal_) = 1;
tmp_index_offdiagonal_ = efind(1-eye(numel(tmp_index_use_)));
tmp_XX_med = median(tmp_XX_use__(1+tmp_index_offdiagonal_));
tmp_XX_avg = mean(tmp_XX_use__(1+tmp_index_offdiagonal_));
imagesc(tmp_XX_use__,[0,1]);
axisnotick; axis image; tmp_c_ = colorbar; set(tmp_c_,'Ticks',[0,0.5,1],'TickLength',[0]);
title(sprintf('%s (%0.2f) [%0.2f]',u_dir_nopath_data_star,tmp_XX_avg,tmp_XX_med),'Interpreter','none');
sgtitle(fname_fig_pre,'Interpreter','none');
set(gca,'FontSize',fontsize_use);
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
sgtitle('');
tmp_str = u_dir_nopath_data_star; if strcmp(tmp_str,'precatalytic_spliceosome'); tmp_str = 'ps1'; end;
title(sprintf('%s (%0.2f) [%0.2f]',tmp_str,tmp_XX_avg,tmp_XX_med),'Interpreter','none');
tmp_dir = sprintf('/%s/rangan/dir_cryoem/dir_ampm_manuscript/dir_ampm_fig_bootstrap',string_root);
fname_fig_jpg_strip = sprintf('%s/test_pm_collect_%s_ampm_alig_crosscorrelation_FIGL_strip.jpg',tmp_dir,u_dir_nopath_data_star);
disp(sprintf(' %% writing %s',fname_fig_jpg_strip));
print('-djpeg',sprintf('%s',fname_fig_jpg_strip));
close(gcf);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%%%%%%%%%;
end;%if  exist(fname_mat,'file');
%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%;
end;%for nu_experiment=0:n_u_experiment-1-1;%for nu_experiment=0:n_u_experiment-1;
%%%%%%%%%%%%%%%%;


%%%%%%%%;
% Now cluster the LSUbl17dep results. ;
% assume nu_experiment==1. ;
%%%%%%%%;
nu_experiment = 1
index_nexperiment_from_nu_ = index_nexperiment_from_nu__{1+nu_experiment};
%%%%%%%%%%%%%%%%;
u_dir_nopath_data_star = u_dir_nopath_data_star_{1+nu_experiment};
fname_mat_pre = sprintf('%s/test_pm_collect_%s_ampm_alig_crosscorrelation',dir_mat,u_dir_nopath_data_star);
fname_mat = sprintf('%s.mat',fname_mat_pre);
flag_rank_vs_tolerance_use = 1-strcmp(u_dir_nopath_data_star,'p28hRPT1');
n_experiment_local = numel(index_nexperiment_from_nu_);
disp(sprintf(' %% nu_experiment %d/%d u_dir_nopath_data_star %s flag_rank_vs_tolerance_use %d n_experiment_local %d',nu_experiment,n_u_experiment,u_dir_nopath_data_star,flag_rank_vs_tolerance_use,n_experiment_local));
%%%%;
na_ampm_tot=0;
tmp_X_ampm_a_ = [];
tmp_tolerance_pm_a_ = [];
tmp_delta_sigma_use_a_ = [];
tmp_flag_rank_vs_tolerance_a_ = [];
ampm_local__ = cell(n_experiment_local,1);
%%%%;
for nl=0:n_experiment_local-1;
nexperiment = index_nexperiment_from_nu_(1+nl);
%%%%;
collect_ = collect__{1+nexperiment};
tmp_X_ = collect_.X_best_a_;
if (~isempty(collect_.corr_crop_xa__));
tmp_X_ = collect_.corr_crop_xa__;
[~,tmp_index_crop] = max(mean(tmp_X_,2)); tmp_index_crop = tmp_index_crop-1;
if (verbose); disp(sprintf(' %% %s: crop at index %d',collect_.dir_nopath_data_star,tmp_index_crop)); end;
tmp_X_ = tmp_X_(1+tmp_index_crop,:);
end;%if (~isempty(collect_.corr_crop_xa__));
%%%%;
n_a = numel(collect_.tolerance_pm_a_);
tmp_index_ = na_ampm_tot + [0:n_a-1];
tmp_X_ampm_a_(1 + tmp_index_) = tmp_X_;
tmp_tolerance_pm_a_(1 + tmp_index_) = collect_.tolerance_pm_a_;
tmp_delta_sigma_use_a_(1 + tmp_index_) = collect_.delta_sigma_use_a_;
tmp_flag_rank_vs_tolerance_a_(1 + tmp_index_) = collect_.flag_rank_vs_tolerance_a_;
%%%%;
na_ampm_tot = na_ampm_tot + n_a;
%%%%;
[~,ampm_local__{1+nl}] = ampm_fromdisk_3(struct('type','parameter','flag_store_S_k_p__',0,'flag_store_M_k_p__',0),collect_.dir_pm);
%%%%;
end;%for nl=0:n_experiment_local-1;
n_a_ampm_tot = na_ampm_tot;
%%%%%%%%;
tmp_index_ampm_ = efind( ...
 tmp_X_ampm_a_~=0 ...
& tmp_tolerance_pm_a_<=10^-1.75 ...
& tmp_delta_sigma_use_a_> 0.00 ...
& tmp_flag_rank_vs_tolerance_a_==flag_rank_vs_tolerance_use ...
);
n_a_ampm_use = numel(tmp_index_ampm_);
%%%%%%%%;
if (verbose);
disp(sprintf(' %% nu_experiment %d/%d: %s',nu_experiment,n_u_experiment,u_dir_nopath_data_star));
disp(sprintf(' %% %% n_a_ampm_tot %d <-- use %d',n_a_ampm_tot,n_a_ampm_use));
end;%if (verbose);
%%%%%%%%;
load(fname_mat);
%%%%;
tmp_ampm_local_ = ampm_local__{1+0};
dir_pm = tmp_ampm_local_.dir_pm;
n_k_p_r = tmp_ampm_local_.n_k_p_r;
k_p_r_ = tmp_ampm_local_.k_p_r_;
k_p_r_max = tmp_ampm_local_.k_p_r_max;
weight_3d_k_p_r_ = tmp_ampm_local_.weight_3d_k_p_r_;
weight_3d_k_all_ = tmp_ampm_local_.weight_3d_k_all_;
weight_shell_k_ = tmp_ampm_local_.weight_shell_k_;
n_w_max = tmp_ampm_local_.n_w_max;
n_m_max = tmp_ampm_local_.n_m_max;
l_max_ = tmp_ampm_local_.l_max_; l_max_max = max(l_max_); n_l_max = l_max_max + 1;
n_k_all = tmp_ampm_local_.n_k_all;
n_k_all_csum_ = tmp_ampm_local_.n_k_all_csum_;
k_p_r_all_ = tmp_ampm_local_.k_p_r_all_;
k_p_azimu_b_all_ = tmp_ampm_local_.k_p_azimu_b_all_;
k_p_polar_a_all_ = tmp_ampm_local_.k_p_polar_a_all_;
k_c_0_all_ = tmp_ampm_local_.k_c_0_all_;
k_c_1_all_ = tmp_ampm_local_.k_c_1_all_;
k_c_2_all_ = tmp_ampm_local_.k_c_2_all_;
n_xxx_u = tmp_ampm_local_.n_xxx_u;
xxx_u_weight_ = tmp_ampm_local_.xxx_u_weight_;
x_u_0___ = tmp_ampm_local_.x_u_0___;
x_u_1___ = tmp_ampm_local_.x_u_1___;
x_u_2___ = tmp_ampm_local_.x_u_2___;
x_p_r_max = tmp_ampm_local_.x_p_r_max;
n_x_u_pack = tmp_ampm_local_.n_x_u_pack;
if ~exist('Ylm_klma___','var');
[ ...
 Ylm_klma___ ...
] = ...
get_Ylm_wrap_0( ...
 verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,n_k_p_r ...
,l_max_ ...
);
end;%if ~exist('Ylm_klma___','var');
%%%%%%%%;
tmp_a_k_Y_ampm_alig_yka__ = zeros(ampm_local__{1}.n_lm_sum,0);
na=0;
for nl=0:n_experiment_local-1;
nexperiment = index_nexperiment_from_nu_(1+nl);
tmp_ampm_local_ = ampm_local__{1+nl};
tmp_n_a = numel(tmp_ampm_local_.str_fname_mat_a_);
for tmp_na=0:tmp_n_a-1;
if (mod(tmp_na,8)==0); disp(sprintf(' %% tmp_na %d/%d',tmp_na,tmp_n_a)); end;
[ ... 
tmp_a_k_Y_ampm_alig_yka__(:,1+na) ...
] = ...
spharm_register_and_rotate_2( ...
 tmp_ampm_local_.n_k_p_r ...
,tmp_ampm_local_.k_p_r_ ...
,tmp_ampm_local_.k_p_r_max ...
,tmp_ampm_local_.weight_3d_k_p_r_ ...
,tmp_ampm_local_.l_max_ ...
,tmp_ampm_local_.a_k_Y_quad_ ...
,tmp_ampm_local_.a_k_Y_ampm_yka__(:,1+tmp_na) ...
,0 ...
,tmp_ampm_local_.X_best_ampm_ia__(end,1+tmp_na) ...
,tmp_ampm_local_.X_best_flag_flip_ampm_ia__(end,1+tmp_na) ...
,tmp_ampm_local_.polar_a_best_ampm_ia__(end,1+tmp_na) ...
,tmp_ampm_local_.azimu_b_best_ampm_ia__(end,1+tmp_na) ...
,tmp_ampm_local_.gamma_z_best_ampm_ia__(end,1+tmp_na) ...
,tmp_ampm_local_.delta_best_ampm_dia__(1+0,end,1+tmp_na) ...
,tmp_ampm_local_.delta_best_ampm_dia__(1+1,end,1+tmp_na) ...
,tmp_ampm_local_.delta_best_ampm_dia__(1+2,end,1+tmp_na) ...
);
na=na+1;
end;%for tmp_na=0:n_a-1;
end;%for nl=0:n_experiment_local-1;
%%%%%%%%;
[tmp_X_ampm_use_srt_,tmp_ij_X_ampm_use_srt_] = sort(tmp_X_ampm_a_(1+tmp_index_ampm_),'ascend');
tmp_index_use_ = tmp_index_ampm_(tmp_ij_X_ampm_use_srt_); n_index_use = numel(tmp_index_use_);
tmp_a_k_Y_ampm_alig_use_yka__ = tmp_a_k_Y_ampm_alig_yka__(:,1+tmp_index_use_);
%%%%%%%%;
[tmp_U__,tmp_S__,tmp_V__] = svds(tmp_XX_use__,2);
[~,tmp_ij_] = sort(tmp_U__(:,2),'ascend');
imagesc(tmp_XX_use__(tmp_ij_,tmp_ij_),[0.75,1]);figbeach();
tmp_a_k_Y_ampm_alig_0_yk_ = mean(tmp_a_k_Y_ampm_alig_use_yka__(:,tmp_ij_( 1:20)),2);
tmp_a_k_Y_ampm_alig_1_yk_ = mean(tmp_a_k_Y_ampm_alig_use_yka__(:,tmp_ij_(21:48)),2);
%%%%%%%%;
[ ...
 tmp_a_k_p_ampm_alig_0_ ...
,Ylm_klma___ ...
] = ...
convert_spharm_to_k_p_3( ...
 verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,tmp_a_k_Y_ampm_alig_0_yk_ ...
,Ylm_klma___ ...
);
[ ...
 tmp_a_x_c_ampm_alig_0_ ...
] = ...
convert_k_p_to_x_c_1( ...
 verbose ...
,n_k_all ...
,weight_3d_k_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,k_p_r_max ...
,n_xxx_u ...
,xxx_u_weight_ ...
,x_u_0___ ...
,x_u_1___ ...
,x_u_2___ ...
,x_p_r_max ...
,tmp_a_k_p_ampm_alig_0_ ...
);
%%%%;
[ ...
 tmp_a_k_p_ampm_alig_1_ ...
,Ylm_klma___ ...
] = ...
convert_spharm_to_k_p_3( ...
 verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,tmp_a_k_Y_ampm_alig_1_yk_ ...
,Ylm_klma___ ...
);
[ ...
 tmp_a_x_c_ampm_alig_1_ ...
] = ...
convert_k_p_to_x_c_1( ...
 verbose ...
,n_k_all ...
,weight_3d_k_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,k_p_r_max ...
,n_xxx_u ...
,xxx_u_weight_ ...
,x_u_0___ ...
,x_u_1___ ...
,x_u_2___ ...
,x_p_r_max ...
,tmp_a_k_p_ampm_alig_1_ ...
);
%%%%;

fname_fig_pre = sprintf('/%s/rangan/dir_cryoem/dir_jpg/LSUbl17dep_ampm_vol_cluster_FIGM',string_root);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
if strfind(tmp_ampm_local_.dir_pm,'ISWINCP');    percent_threshold_ = 94.50; tmp_nx = 14; end;
if strfind(tmp_ampm_local_.dir_pm,'trpv1');      percent_threshold_ = 91.25; tmp_nx = 14; end;
if strfind(tmp_ampm_local_.dir_pm,'rib80s');     percent_threshold_ = 86.25; tmp_nx = 14; end;
if strfind(tmp_ampm_local_.dir_pm,'MlaFEDB');    percent_threshold_ = 95.00; tmp_nx = 8; end;
if strfind(tmp_ampm_local_.dir_pm,'LetB1');      percent_threshold_ = 91.75; tmp_nx = 14; end;
if strfind(tmp_ampm_local_.dir_pm,'TMEM16F');    percent_threshold_ = 94.50; tmp_nx = 16; end;
if strfind(tmp_ampm_local_.dir_pm,'LSUbl17dep'); percent_threshold_ = 86.25; tmp_nx = 15; end;
if strfind(tmp_ampm_local_.dir_pm,'ps1');        percent_threshold_ = 96.00; tmp_nx = 11; end;
tmp_window_ = zeros(n_x_u_pack,n_x_u_pack,n_x_u_pack);
tmp_index_ = tmp_nx:n_x_u_pack-1-tmp_nx;
tmp_window_(1+tmp_index_,1+tmp_index_,1+tmp_index_)=1;
tmp_index_ = efind(tmp_window_);
figure(1+nf);nf=nf+1;clf;figmed;
fontsize_use = 24;
subplot(1,2,1); isosurface_f_x_u_0(tmp_a_x_c_ampm_alig_0_(1+tmp_index_),percent_threshold_(1+0)); axis off; title('Class  ''Low'''); set(gca,'FontSize',fontsize_use); view([-150,-50]);
subplot(1,2,2); isosurface_f_x_u_0(tmp_a_x_c_ampm_alig_1_(1+tmp_index_),percent_threshold_(1+0)); axis off; title('Class ''High'''); set(gca,'FontSize',fontsize_use); view([-150,-22]);
set(gcf,'Position',1+[0,0,512+256,512]);
sgtitle(fname_fig_jpg,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
sgtitle('');
tmp_dir = sprintf('/%s/rangan/dir_cryoem/dir_ampm_manuscript/dir_ampm_fig_bootstrap',string_root);
fname_fig_jpg_strip = sprintf('%s/LSUbl17dep_ampm_vol_cluster_FIGM_strip.jpg',tmp_dir);
disp(sprintf(' %% writing %s',fname_fig_jpg_strip));
print('-djpeg',sprintf('%s',fname_fig_jpg_strip));
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));







