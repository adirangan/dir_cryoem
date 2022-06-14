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

verbose=1;
flag_recoll = 0;
flag_replot = 1;
flag_center_volume = 0;
flag_center_image = 0;
flag_invert = 0;
flag_crop = 0;
tolerance_master = 1e-2;
nf=0;

table_data__ = { ...
 0 , 'p28hRPT1_x0' , 'p28hRPT1' , 0.98 , 'emd_8674.map' , 'T1.star' ; ...
 0 , 'trpv1_x0' , 'trpv1' , 1.2156 , 'emd_5778.mrc' , 'tv1_relion_data.star' ; ...
 1 , 'trpv1c'   , 'trpv1' , 1.2156 , 'emd_5778.mrc' , 'tv1_relion_data.star' ; ...
 0 , 'rib80s_x0' , 'rib80s' , 1.34 , 'emd_2660.mrc' , 'shiny_2sets.star' ; ...
 1 , 'rib80s' , 'rib80s' , 1.34 , 'emd_2660.mrc' , 'shiny_2sets.star' ; ...
 0 , 'MlaFEDB_x0' , 'MlaFEDB' , 0.832 , 'emd_22116.map' , 'Empiar_10536_00_to_23.star' ; ...
 1 , 'MlaFEDB' , 'MlaFEDB' , 0.832 , 'emd_22116.map' , 'Empiar_10536_00_to_23.star' ; ...
 0 , 'ISWINCP_x0' , 'ISWINCP' , 1.07 , 'emd_9718.map' , 'ADPBeF.star' ; ...
 1 , 'ISWINCP' , 'ISWINCP' , 1.07 , 'emd_9718.map' , 'ADPBeF.star' ; ...
 0 , 'TMEM16F_x0' , 'TMEM16F' , 1.059 , 'emd_20244.map' , 'All_8192.star' ; ...
 1 , 'TMEM16F' , 'TMEM16F' , 1.059 , 'emd_20244.map' , 'All_8192.star' ; ...
 0 , 'LSUbl17dep_x0' , 'LSUbl17dep' , 1.31 , 'emd_8434.map' , 'Parameters.star' ; ...
 1 , 'LSUbl17dep' , 'LSUbl17dep' , 1.31 , 'emd_8434.map' , 'Parameters.star' ; ...
 0 , 'ps1_x0' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
 1 , 'ps0' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
% 1 , 'ps1' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
 0 , 'LetB1_x0' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1_350MB.star' ; ...
 1 , 'LetB1' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1_350MB.star' ; ...
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

collect__{1+nexperiment}.corr_crop_xa__ = [];
collect__{1+nexperiment}.corr_crop_reco_x_ = [];
if flag_crop; collect__{1+nexperiment}.corr_crop_xa__ = corr_crop_xa__; end;
if flag_crop; collect__{1+nexperiment}.corr_crop_reco_x_ = corr_crop_reco_x_; end;
collect__{1+nexperiment}.fsc_crop_reco_stab_alig_from_x__ = fsc_crop_reco_stab_alig_from_x__;
collect__{1+nexperiment}.fsc_sher_rs__ = fsc_sher_rs__;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nexperiment=0:n_experiment-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

fname_fig = sprintf('%s/test_pm_collect_plus_X_best_FIGA',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
tlim_ = [1,3];%tlim_ = [1.25,2.5];
dlim_ = [0,0.15];
figure(1+nf);nf=nf+1;clf;figbig; fig80s;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
markersize_use = 8;
str_symbol_ = {'o','^','s','p','h'};
p_row = 3; p_col = ceil(n_experiment/p_row); np=0;
for nexperiment=0:n_experiment-1;
collect_ = collect__{1+nexperiment};
tmp_X_ = collect_.X_best_a_;
if (~isempty(collect_.corr_crop_xa__));
tmp_X_ = collect_.corr_crop_xa__;
[~,tmp_index_crop] = max(mean(tmp_X_,2)); tmp_index_crop = tmp_index_crop-1;
if (verbose); disp(sprintf(' %% %s: crop at index %d',collect_.dir_nopath_data_star,tmp_index_crop)); end;
tmp_X_ = tmp_X_(1+tmp_index_crop,:);
end;%if (~isempty(collect_.corr_crop_xa__));
subplot(p_row,p_col,1+np);np=np+1;
n_a = numel(collect_.tolerance_pm_a_);
[~,tmp_index_] = sort(collect_.tolerance_pm_a_,'ascend'); tmp_index_ = tmp_index_-1;
hold on;
for na=0:n_a-1;
tmp_delta_sigma_use = collect_.delta_sigma_use_a_(1+tmp_index_(1+na));
tmp_flag_rank_vs_tolerance = collect_.flag_rank_vs_tolerance_a_(1+tmp_index_(1+na));
if (1 | tmp_delta_sigma_use> 0);
if (1 | tmp_flag_rank_vs_tolerance==1);
str_symbol = str_symbol_{1+tmp_flag_rank_vs_tolerance};
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*(tmp_delta_sigma_use-min(dlim_))/diff(dlim_))));
tmp_X = tmp_X_(1+tmp_index_(1+na));
if (tmp_X~=0);
plot(-log10(collect_.tolerance_pm_a_(1+tmp_index_(1+na))),tmp_X,str_symbol,'MarkerFaceColor',c_80s__(1+nc_80s,:),'MarkerSize',markersize_use);
end;%if (tmp_X~=0);
end;%if (tmp_flag_rank_vs_tolerance==1);
end;%if (tmp_delta_sigma_use> 0);
end;%for na=0:n_a-1;
if ( isempty(collect_.corr_crop_reco_x_)); plot(tlim_,real(collect_.corr_a_k_Y)*[1,1],'k-','LineWidth',3); end;
if (~isempty(collect_.corr_crop_reco_x_)); plot(tlim_,real(collect_.corr_crop_reco_x_(1+tmp_index_crop))*[1,1],'k-','LineWidth',3); end;
plot(tlim_,mean(real(collect_.X_rand_a_))*[1,1],'k:','LineWidth',3);
plot(tlim_,collect_.X_relion*[1,1],'b-','LineWidth',2);
hold off;
xlabel('tolerance');
ylabel('correlation');
xlim(tlim_);
ylim([-0.2,+1.0]);
grid on;
title(collect_.dir_nopath_data_star,'Interpreter','none');
end;%for nexperiment=0:n_experiment-1;
sgtitle(fname_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/test_pm_collect_plus_X_best_crop_FIGB',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 3;
n_crop = 64;
for nexperiment=0:n_experiment-1;
collect_ = collect__{1+nexperiment};
subplot(p_row,ceil(n_experiment/p_row),1+nexperiment);
hold on;
plot(1:n_crop,collect_.X_relion*ones(n_crop,1),'k:');
plot(1:n_crop,collect_.corr_crop_reco_x_,'ko-');
plot(1:n_crop,mean(collect_.corr_crop_xa__(:,:),2),'m.');
tmp_index_ = efind( (collect_.delta_sigma_use_a_> 0) & (collect_.flag_rank_vs_tolerance_a_==0) );
plot(1:n_crop,mean(collect_.corr_crop_xa__(:,1+tmp_index_),2),'bx');
tmp_index_ = efind( (collect_.delta_sigma_use_a_> 0) & (collect_.flag_rank_vs_tolerance_a_==1) );
plot(1:n_crop,mean(collect_.corr_crop_xa__(:,1+tmp_index_),2),'rx');
hold off;
ylim([-0.2,1]); ylabel('correlation');
xlim([1,n_crop]); xlabel('crop');
grid on;
title(collect_.dir_nopath_data_star,'Interpreter','none');
end;%for nexperiment=0:n_experiment-1;
sgtitle(fname_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/test_pm_collect_plus_X_disc_FIGC',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
n_M = 1024;
tlim_ = [1,3];%tlim_ = [1.25,2.5];
dlim_ = [0,0.15];
figure(1+nf);nf=nf+1;clf;figbig; fig80s;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
markersize_use = 8;
str_symbol_ = {'o','^','s','p','h'};
p_row = 3; p_col = ceil(n_experiment/p_row); np=0;
for nexperiment=0:n_experiment-1;
collect_ = collect__{1+nexperiment};
subplot(p_row,p_col,1+np);np=np+1;
n_a = numel(collect_.tolerance_pm_a_);
[~,tmp_index_] = sort(collect_.tolerance_pm_a_,'ascend'); tmp_index_ = tmp_index_-1;
hold on;
for na=0:n_a-1;
tmp_delta_sigma_use = collect_.delta_sigma_use_a_(1+tmp_index_(1+na));
tmp_flag_rank_vs_tolerance = collect_.flag_rank_vs_tolerance_a_(1+tmp_index_(1+na));
if (1 | tmp_delta_sigma_use> 0);
if (1 | tmp_flag_rank_vs_tolerance==1);
str_symbol = str_symbol_{1+tmp_flag_rank_vs_tolerance};
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*(tmp_delta_sigma_use-min(dlim_))/diff(dlim_))));
tmp_X_disc = collect_.disc_p25_a_(1+tmp_index_(1+na))/n_M;
plot(-log10(collect_.tolerance_pm_a_(1+tmp_index_(1+na))),tmp_X_disc,str_symbol,'MarkerFaceColor',c_80s__(1+nc_80s,:),'MarkerSize',markersize_use/4);
tmp_X_disc = collect_.disc_p50_a_(1+tmp_index_(1+na))/n_M;
plot(-log10(collect_.tolerance_pm_a_(1+tmp_index_(1+na))),tmp_X_disc,str_symbol,'MarkerFaceColor',c_80s__(1+nc_80s,:),'MarkerSize',markersize_use/2);
tmp_X_disc = collect_.disc_p75_a_(1+tmp_index_(1+na))/n_M;
plot(-log10(collect_.tolerance_pm_a_(1+tmp_index_(1+na))),tmp_X_disc,str_symbol,'MarkerFaceColor',c_80s__(1+nc_80s,:),'MarkerSize',markersize_use/1);
end;%if (tmp_flag_rank_vs_tolerance==1);
end;%if (tmp_delta_sigma_use> 0);
end;%for na=0:n_a-1;
hold off;
xlabel('tolerance');
ylabel('average rank');
xlim(tlim_);
ylim([0.35,+1.0]);
grid on;
title(collect_.dir_nopath_data_star,'Interpreter','none');
end;%for nexperiment=0:n_experiment-1;
sgtitle(fname_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));

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
%%%%%%%%%%%%%%%%%%%%%%%;
fname_fig = sprintf('%s/test_pm_collect_plus_X_best_FIGD',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
tlim_ = [1,3];%tlim_ = [1.25,2.5];
lslim_ = [-3,-1]; 
dlim_ = [0,0.15];
figure(1+nf);nf=nf+1;clf;figbig; fig80s;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
markersize_use = 8;
str_symbol_ = {'o','^','s','p','h'};  
p_row = 3; p_col = ceil(2*n_u_experiment/p_row); np=0;
%%%%%%%%%%%%%%%%;
for nu_experiment=0:n_u_experiment-1;
index_nexperiment_from_nu_ = index_nexperiment_from_nu__{1+nu_experiment};
%%%%%%%%%%%%%%%%;
subplot(p_row,p_col,1+np);np=np+1;
tmp_collect_corr = -Inf;
tmp_collect_rand = -Inf;
tmp_collect_reli = -Inf;
tmp_collect_corr_reco_stab = -Inf;
hold on;
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
[~,tmp_index_] = sort(collect_.tolerance_pm_a_,'ascend'); tmp_index_ = tmp_index_-1;
for na=0:n_a-1;
tmp_delta_sigma_use = collect_.delta_sigma_use_a_(1+tmp_index_(1+na));
tmp_flag_rank_vs_tolerance = collect_.flag_rank_vs_tolerance_a_(1+tmp_index_(1+na));
if (0 | tmp_delta_sigma_use> 0);
if ( ...
     ( strcmp(collect_.dir_nopath_data_star,'p28hRPT1') & tmp_flag_rank_vs_tolerance==0) ...
     | ...
     (~strcmp(collect_.dir_nopath_data_star,'p28hRPT1') & tmp_flag_rank_vs_tolerance==1) ...
     );
str_symbol = str_symbol_{1+tmp_flag_rank_vs_tolerance};
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*(tmp_delta_sigma_use-min(dlim_))/diff(dlim_))));
tmp_X = tmp_X_(1+tmp_index_(1+na));
if (tmp_X~=0);
plot(-log10(collect_.tolerance_pm_a_(1+tmp_index_(1+na))),tmp_X,str_symbol,'MarkerFaceColor',c_80s__(1+nc_80s,:),'MarkerSize',markersize_use);
end;%if (tmp_X~=0);
end;%if (tmp_flag_rank_vs_tolerance==1);
end;%if (tmp_delta_sigma_use> 0);
end;%for na=0:n_a-1;
%%%%;
if ( isempty(collect_.corr_crop_reco_x_)); tmp_collect_corr = max(tmp_collect_corr,real(collect_.corr_a_k_Y)); end;
if (~isempty(collect_.corr_crop_reco_x_)); tmp_collect_corr = max(tmp_collect_corr,real(collect_.corr_crop_reco_x_(1+tmp_index_crop))); end;
tmp_collect_rand = max(tmp_collect_rand,mean(real(collect_.X_rand_a_)));
tmp_collect_reli = max(tmp_collect_reli,collect_.X_relion);
%%%%%%%%;
tmp_c = -Inf;
if (~isempty(collect_.fsc_crop_reco_stab_alig_from_x__));
tmp_c = collect_.fsc_crop_reco_stab_alig_from_x__.corr_full_reco_vs_crop_reco_stab_x_(1+tmp_index_crop);
end;%if (~isempty(collect_.fsc_crop_reco_stab_alig_from_x__));
tmp_collect_corr_reco_stab = max(tmp_collect_corr_reco_stab,tmp_c);
%%%%%%%%;
end;%for nl=0:numel(index_nexperiment_from_nu_)-1;
%%%%;
plot(tlim_,tmp_collect_corr*[1,1],'k-','LineWidth',3);
%plot(tlim_,tmp_collect_rand*[1,1],'k:','LineWidth',3);
plot(tlim_,tmp_collect_reli*[1,1],'b-','LineWidth',2);
plot(tlim_,tmp_collect_corr_reco_stab*[1,1],'c-','LineWidth',2);
hold off;
%%%%;
title(collect_.dir_nopath_data_star,'Interpreter','none');
xlabel('tolerance');
ylabel('correlation');
xlim(tlim_);
ylim([-0.2,+1.0]);
grid on;
%%%%%%%%%%%%%%%%;
subplot(p_row,p_col,1+np);np=np+1;
tmp_collect_corr = -Inf;
tmp_collect_rand = -Inf;
tmp_collect_reli = -Inf;
tmp_collect_corr_reco_stab = -Inf;
hold on;
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
str_symbol = 'x';
fsc_sher_rs__ = collect_.fsc_sher_rs__;
for nsigma_sheres=0:n_sigma_sheres-1;
sigma_sheres = sigma_sheres_(1+nsigma_sheres);
lsigma_sheres = log(sigma_sheres);
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
if isfield(fsc_sher_rs__{1+nrseed,1+nsigma_sheres},'corr_full_reco_vs_crop_reco_stab_x_');
tmp_X = fsc_sher_rs__{1+nrseed,1+nsigma_sheres}.corr_full_reco_vs_crop_reco_stab_x_(1+tmp_index_crop);
if (tmp_X~=0);
plot(lsigma_sheres,tmp_X,str_symbol,'MarkerFaceColor','g','MarkerSize',markersize_use);
end;%if (tmp_X~=0);
end;%if isfield(fsc_sher_rs__{1+nrseed,1+nsigma_sheres},'corr_full_reco_vs_crop_reco_stab_x_');
end;%for nrseed=0:n_rseed-1;
end;%for nsigma_sheres=0:n_sigma_sheres-1;
%%%%;
if ( isempty(collect_.corr_crop_reco_x_)); tmp_collect_corr = max(tmp_collect_corr,real(collect_.corr_a_k_Y)); end;
if (~isempty(collect_.corr_crop_reco_x_)); tmp_collect_corr = max(tmp_collect_corr,real(collect_.corr_crop_reco_x_(1+tmp_index_crop))); end;
tmp_collect_rand = max(tmp_collect_rand,mean(real(collect_.X_rand_a_)));
tmp_collect_reli = max(tmp_collect_reli,collect_.X_relion);
%%%%%%%%;
tmp_c = -Inf;
if (~isempty(collect_.fsc_crop_reco_stab_alig_from_x__));
tmp_c = collect_.fsc_crop_reco_stab_alig_from_x__.corr_full_reco_vs_crop_reco_stab_x_(1+tmp_index_crop);
end;%if (~isempty(collect_.fsc_crop_reco_stab_alig_from_x__));
tmp_collect_corr_reco_stab = max(tmp_collect_corr_reco_stab,tmp_c);
%%%%%%%%;
end;%for nl=0:numel(index_nexperiment_from_nu_)-1;
%%%%;
plot(lslim_,tmp_collect_corr*[1,1],'k-','LineWidth',3);
%plot(lslim_,tmp_collect_rand*[1,1],'k:','LineWidth',3);
plot(lslim_,tmp_collect_reli*[1,1],'b-','LineWidth',2);
plot(lslim_,tmp_collect_corr_reco_stab*[1,1],'c-','LineWidth',2);
hold off;
%%%%;
title(collect_.dir_nopath_data_star,'Interpreter','none');
xlabel('$\log(\sigma)$','Interpreter','latex');
ylabel('correlation');
xlim(lslim_);
ylim([-0.2,+1.0]);
grid on;
%%%%%%%%%%%%%%%%;
end;%for nu_experiment=0:n_u_experiment-1;
%%%%%%%%%%%%%%%%;
sgtitle(fname_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% plot residual discriminability. ;
% ignore the _x0 runs. ;
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
%%%%%%%%%%%%%%%%%%%%%%%;
fname_fig = sprintf('%s/test_pm_collect_plus_X_disc_FIGE',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
tlim_ = [1,3];%tlim_ = [1.25,2.5];
dlim_ = [0,0.15];
figure(1+nf);nf=nf+1;clf;figbig; fig80s;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
markersize_use = 8;
str_symbol_ = {'o','^','s','p','h'};  
p_row = 3; p_col = ceil(n_u_experiment/p_row); np=0;
%%%%%%%%;
for nu_experiment=0:n_u_experiment-1;
index_nexperiment_from_nu_ = index_nexperiment_from_nu__{1+nu_experiment};
subplot(p_row,p_col,1+np);np=np+1;
hold on;
for nl=0:numel(index_nexperiment_from_nu_)-1;
nexperiment = index_nexperiment_from_nu_(1+nl);
collect_ = collect__{1+nexperiment};
%%%%%%%%;
n_a = numel(collect_.tolerance_pm_a_);
[~,tmp_index_] = sort(collect_.tolerance_pm_a_,'ascend'); tmp_index_ = tmp_index_-1;
for na=0:n_a-1;
tmp_delta_sigma_use = collect_.delta_sigma_use_a_(1+tmp_index_(1+na));
tmp_flag_rank_vs_tolerance = collect_.flag_rank_vs_tolerance_a_(1+tmp_index_(1+na));
if ( ...
     strcmp(collect_.dir_nopath_data_star,'p28hRPT1') ...
     | ...
     (~strcmp(collect_.dir_nopath_data_star,'p28hRPT1') ...
      & ...
      isempty(strfind(collect_.fname_prefix,'_x0')) ...
      ) ...
     );
if (0 | tmp_delta_sigma_use> 0);
if ( ...
     ( strcmp(collect_.dir_nopath_data_star,'p28hRPT1') & tmp_flag_rank_vs_tolerance==0) ...
     | ...
     (~strcmp(collect_.dir_nopath_data_star,'p28hRPT1') & tmp_flag_rank_vs_tolerance==1) ...
     );
str_symbol = str_symbol_{1+tmp_flag_rank_vs_tolerance};
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*(tmp_delta_sigma_use-min(dlim_))/diff(dlim_))));
tmp_X_disc = collect_.disc_p25_a_(1+tmp_index_(1+na))/n_M;
plot(-log10(collect_.tolerance_pm_a_(1+tmp_index_(1+na))),tmp_X_disc,str_symbol,'MarkerFaceColor',c_80s__(1+nc_80s,:),'MarkerSize',markersize_use/4);
tmp_X_disc = collect_.disc_p50_a_(1+tmp_index_(1+na))/n_M;
plot(-log10(collect_.tolerance_pm_a_(1+tmp_index_(1+na))),tmp_X_disc,str_symbol,'MarkerFaceColor',c_80s__(1+nc_80s,:),'MarkerSize',markersize_use/2);
tmp_X_disc = collect_.disc_p75_a_(1+tmp_index_(1+na))/n_M;
plot(-log10(collect_.tolerance_pm_a_(1+tmp_index_(1+na))),tmp_X_disc,str_symbol,'MarkerFaceColor',c_80s__(1+nc_80s,:),'MarkerSize',markersize_use/1);
end;%if (tmp_flag_rank_vs_tolerance==1);
end;%if (tmp_delta_sigma_use> 0);
end;%if retain ;
end;%for na=0:n_a-1;
%%%%%%%%;
end;%for nl=0:numel(index_nexperiment_from_nu_)-1;
%%%%;
plot(tlim_,0.5*ones(size(tlim_)),'k:','LineWidth',3);
hold off;
xlabel('tolerance');
ylabel('average rank');
xlim(tlim_);
ylim([0.35,+1.0]);
grid on;
title(collect_.dir_nopath_data_star,'Interpreter','none');
end;%for nu_experiment=0:n_u_experiment-1;
%%%%%%%%;
sgtitle(fname_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));


disp('returning');return;


