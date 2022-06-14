function ...
[ ...
 global_parameter ...
] = ...
test_pm_23_clean_20220307( ...
 global_parameter ...
,fname_prefix ...
,dir_nopath_data_star ...
,Pixel_Spacing ...
,fname_nopath_volume ...
,fname_nopath_star ...
);

%{

  dir_base = '/data/rangan/dir_cryoem';
  for fname_prefix_ = {'ISWINCP_x0','p28hRPT1_x0','trpv1_x0','rib80s_x0','MlaFEDB_x0','LetB1_x0','TMEM16F_x0','LSUbl17dep_x0','ps1_x0','LetB1_x0'};
  fname_prefix = fname_prefix_{1};
  dir_pm = sprintf('%s/dir_%s/dir_pm',dir_base,fname_prefix);
  %%%%;
  fname_mat = sprintf('%s_mat/fsc_crop_reco_stab_alig_from_M__.mat',dir_pm);
  if (exist(fname_mat,'file'));
  tmp_ = load(fname_mat);
  ncrop = 40;
  tmp_c0 = tmp_.corr_reco_vs_reco_stab;
  tmp_c1 = tmp_.corr_full_reco_vs_crop_reco_stab_x_(1+ncrop);
  disp(sprintf(' %% %s: \n %% %% %% corr_reco_vs_reco_stab %0.2f corr_full_reco_vs_crop_reco_stab_40 %0.2f',fname_mat,tmp_c0,tmp_c1));
  end;%if (exist(fname_mat,'file'));
  %%%%;
  fname_mat = sprintf('%s_mat/fsc_crop_reco_stab_alig_from_N__.mat',dir_pm);
  if (exist(fname_mat,'file'));
  tmp_ = load(fname_mat);
  ncrop = 40;
  tmp_c0 = tmp_.corr_reco_vs_reco_stab;
  tmp_c1 = tmp_.corr_full_reco_vs_crop_reco_stab_x_(1+ncrop);
  disp(sprintf(' %% %s: \n %% %% %% corr_reco_vs_reco_stab %0.2f corr_full_reco_vs_crop_reco_stab_40 %0.2f',fname_mat,tmp_c0,tmp_c1));
  end;%if (exist(fname_mat,'file'));
  %%%%;
  end;% for fname_prefix_;
  
  %}


if (nargin<1);
table_data__ = { ...
'ISWINCP_x0' , 'ISWINCP' , 1.07 , 'emd_9718.map' , 'ADPBeF.star' ; ...
%'p28hRPT1_x0' , 'p28hRP' , 0.98 , 'emd_8674.map' , 'T1.star' ; ...
'trpv1_x0' , 'trpv1' , 1.2156 , 'emd_5778.mrc' , 'tv1_relion_data.star' ; ...
'rib80s_x0' , 'rib80s' , 1.34 , 'emd_2660.mrc' , 'shiny_2sets.star' ; ...
'MlaFEDB_x0' , 'MlaFEDB' , 0.832 , 'emd_22116.map' , 'Empiar_10536_00_to_23.star' ; ...
'LetB1_x0' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1.star' ; ...
'TMEM16F_x0' , 'TMEM16F' , 1.059 , 'emd_20244.map' , 'All_8192.star' ; ...
'LSUbl17dep_x0' , 'LSUbl17dep' , 1.31 , 'emd_8434.map' , 'Parameters_negated.star' ; ...
'ps1_x0' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
'LetB1_x0' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1_350MB.star' ; ...
};
n_experiment = size(table_data__,1);
%%%%%%%%;
tmp_ = clock;rng(tmp_(end));
for nexperiment=(randperm(n_experiment)-1);
na=0;
fname_prefix = table_data__{1+nexperiment,1+na}; na=na+1;
dir_nopath_data_star = table_data__{1+nexperiment,1+na}; na=na+1;
Pixel_Spacing = table_data__{1+nexperiment,1+na}; na=na+1;
fname_nopath_volume = table_data__{1+nexperiment,1+na}; na=na+1;
fname_nopath_star = table_data__{1+nexperiment,1+na}; na=na+1;
disp(sprintf(' %% nexperiment %d/%d: %16s %16s %0.3f %16s %32s' ...
	     ,nexperiment,n_experiment ...
	     ,fname_prefix,dir_nopath_data_star,Pixel_Spacing,fname_nopath_volume,fname_nopath_star ...
	     ));
global_parameter = struct('type','parameter');
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); global_parameter.flag_invert = 0; end;
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); global_parameter.flag_center_image = 1; end;
if (strcmp(dir_nopath_data_star,'precatalytic_spliceosome')); global_parameter.flag_center_image = 1; end;
global_parameter = ...
test_pm_23_clean_20220307( ...
 global_parameter ...
,fname_prefix ...
,dir_nopath_data_star ...
,Pixel_Spacing ...
,fname_nopath_volume ...
,fname_nopath_star ...
);
end;%for nexperiment=0:n_experiment-1;
%%%%%%%%;
disp('returning');return;
end;%if (nargin<1);

% try: ;
% global_parameter=[];fname_prefix='trpv1_x0';dir_nopath_data_star='trpv1';Pixel_Spacing=1.2156;fname_nopath_volume='emd_5778.mrc';fname_nopath_star='tv1_relion_data.star';
% global_parameter=[];fname_prefix='rib80s_x0';dir_nopath_data_star='rib80s';Pixel_Spacing=1.34;fname_nopath_volume='emd_2660.mrc';fname_nopath_star='shiny_2sets.star';
% global_parameter=[];fname_prefix='LetB1_x0';dir_nopath_data_star='LetB1';Pixel_Spacing=1.31;fname_nopath_volume='emd_20993.map';fname_nopath_star='job_569_model_1_350MB.star';

verbose=1;
if (verbose); disp(sprintf(' %% [entering test_pm_23_clean_20220307]')); end;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

if isempty(global_parameter); global_parameter = struct('type','parameter'); end;
if (~isfield(global_parameter,'flag_recalc')); global_parameter.flag_recalc = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_replot')); global_parameter.flag_replot = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_center_volume')); global_parameter.flag_center_volume = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_center_image')); global_parameter.flag_center_image = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_invert')); global_parameter.flag_invert = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'tolerance_master')); global_parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
flag_recalc = global_parameter.flag_recalc;
flag_replot = global_parameter.flag_replot;
flag_center_volume = global_parameter.flag_center_volume;
flag_center_image = global_parameter.flag_center_image;
flag_invert = global_parameter.flag_invert;
tolerance_master = global_parameter.tolerance_master;
nf=0;

%%%%%%%%;
fname_prefix_xfix = sprintf('%s',fname_prefix);
dir_pm = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_pm',string_root,fname_prefix_xfix);
if (~exist(sprintf('%s_mat',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_pm)); mkdir(sprintf('%s_mat',dir_pm)); end;
if (~exist(sprintf('%s_jpg',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_pm)); mkdir(sprintf('%s_jpg',dir_pm)); end;
dir_relion = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_relion',string_root,fname_prefix_xfix);
if (~exist(sprintf('%s_mat',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_relion)); mkdir(sprintf('%s_mat',dir_relion)); end;
if (~exist(sprintf('%s_jpg',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_relion)); mkdir(sprintf('%s_jpg',dir_relion)); end;
string_rusty_root = 'mnt/home';
dir_relion_bin = sprintf('/%s/rangan/relion/build/bin',string_rusty_root);
dir_data_star = sprintf('/%s/rangan/dir_cryoem/dir_%s',string_root,dir_nopath_data_star);
%%%%%%%%;
% all classes and subclasses. ;
%%%%%%%%;
fname_nopath_volume_ = ...
{ ...
 fname_nopath_volume ... %<-- class 0. ;
};
n_volume = numel(fname_nopath_volume_);
flag_het = 0; if (n_volume> 1); flag_het = 1; end;

flag_clean = 1;
flag_clean_really = 0;

dir_pm_mat = sprintf('%s_mat',dir_pm);
fname_mat_ = ...
  { ...
'a_k_Y_0lsq_reco_.mat' ...
'a_k_Y_0qbp_reco_.mat' ...
'a_k_Y_reco_from_M__.mat' ...
'a_k_Y_reco_from_N__.mat' ...
'a_k_Y_reco_stab_alig_from_M__.mat' ...
'a_k_Y_reco_stab_alig_from_N__.mat' ...
'a_k_Y_reco_stab_from_M__.mat' ...
'a_k_Y_reco_stab_from_N__.mat' ...
'fsc_crop_reco_stab_alig_from_M__.mat' ...
'fsc_crop_reco_stab_alig_from_N__.mat' ...
'X_TM_.mat' ...
'test_pm_collect_corr_crop_.mat' ...
'test_pm_collect_corr_crop_reco_.mat' ...
'test_pm_collect_fsc_crop_ampm_.mat' ...
'test_pm_collect_fsc_crop_reco_.mat' ...
'test_pm_collect.mat' ...
  };
for nl=0:numel(fname_mat_)-1;
fname_mat = sprintf('%s/%s',dir_pm_mat,fname_mat_{1+nl});
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found',fname_mat));
if flag_clean;
disp(sprintf(' %% %s found, deleting',fname_mat));
if flag_clean_really; delete(fname_mat); end;%if flag_clean_really;
end;%if flag_clean;
end;%if ( exist(fname_mat,'file'));
end;%for nl=0:numel(fname_mat_)-1;

delta_sigma_use = 0.0;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
tolerance_pm_ = 0.1.^[1:0.5:3];
n_tolerance_pm = numel(tolerance_pm_);
delta_r_max_factor_ = [0.00,0.125,0.25,0.50];
n_delta_r_max_factor = numel(delta_r_max_factor_);
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
if ( exist(sprintf('%s_mat/delta_r_max_legacy.txt',dir_pm)));
delta_r_max_legacy = textread(sprintf('%s_mat/delta_r_max_legacy.txt',dir_pm));
disp(sprintf(' %% loading delta_r_max_legacy: %0.16f',delta_r_max_legacy));
delta_r_max_use = delta_r_max_factor * delta_r_max_legacy;
end;%if ( exist(sprintf('%s_mat/delta_r_max_legacy.txt',dir_pm)));
for flag_alternate_MS_vs_SM = [1];%for flag_alternate_MS_vs_SM = [0:1];
for ntolerance_pm=0:n_tolerance_pm-1;
tolerance_pm = tolerance_pm_(1+ntolerance_pm);
rank_pm = 12 + 2*ntolerance_pm;
str_strategy = ''; if (flag_alternate_MS_vs_SM); str_strategy = sprintf('%sa1',str_strategy); end;
str_delta_r_max = sprintf('t%.4d',floor(1000*delta_r_max_use));
str_tolerance_pm = sprintf('p%.2d',floor(10*-log10(tolerance_pm)));
str_rank_pm = sprintf('n%.2d',rank_pm);
str_rseed = sprintf('r%d',dat_rseed);
for flag_rank_vs_tolerance=0:1;
if (flag_rank_vs_tolerance==0); str_xfix = sprintf('%s%s%s%s',str_strategy,str_delta_r_max,str_tolerance_pm,str_rseed); end;
if (flag_rank_vs_tolerance==1); str_xfix = sprintf('%s%s%s%s',str_strategy,str_delta_r_max,str_rank_pm,str_rseed); end;
for flag_xcor_vs_Memp=0:1;
if flag_xcor_vs_Memp==0; str_infix = 'X_2d_Memp_d1_'; end;
if flag_xcor_vs_Memp==1; str_infix = 'X_2d_xcor_d0_'; end;
XA_fname_mat = sprintf('%s_mat/%s%s.mat',dir_pm,str_infix,str_xfix);
XA_fname_align_a_k_Y_pre = sprintf('%s_mat/%s%s_align_a_k_Y_',dir_pm,str_infix,str_xfix);
XA_fname_align_a_k_Y_mat = sprintf('%s.mat',XA_fname_align_a_k_Y_pre);
XA_fname_align_a_k_Y_jpg = sprintf('%s.jpg',XA_fname_align_a_k_Y_pre);
XA_fname_snapshot_pre = sprintf('%s_mat/%s%s_snapshot',dir_pm,str_infix,str_xfix);
XA_fname_snapshot_mat = sprintf('%s.mat',XA_fname_snapshot_pre);
XA_fname_snapshot_jpg = sprintf('%s.jpg',XA_fname_snapshot_pre);
XA_fname_compare_image_rank_pre = sprintf('%s_mat/%s%s_compare_image_rank',dir_pm,str_infix,str_xfix);
XA_fname_compare_image_rank_mat = sprintf('%s.mat',XA_fname_compare_image_rank_pre);
XA_fname_compare_image_rank_jpg = sprintf('%s.jpg',XA_fname_compare_image_rank_pre);
XA_fname_align_crop_a_k_Y_mat = sprintf('%s_mat/%s%s_align_crop_a_k_Y_.mat',dir_pm,str_infix,str_xfix);
XA_fname_fsc_crop_kx__mat = sprintf('%s_mat/%s%s_fsc_crop_kx__.mat',dir_pm,str_infix,str_xfix);
if (~exist(XA_fname_mat,'file')); disp(sprintf(' %% %s not found',XA_fname_mat)); end;%if (~exist(XA_fname_mat,'file'));
if ( exist(XA_fname_mat,'file'));
if (verbose>2);
if (~exist(XA_fname_align_a_k_Y_mat,'file')); disp(sprintf(' %% %s not found',XA_fname_align_a_k_Y_mat)); end;
if (~exist(XA_fname_align_a_k_Y_jpg,'file')); disp(sprintf(' %% %s not found',XA_fname_align_a_k_Y_jpg)); end;
if (~exist(XA_fname_snapshot_mat,'file')); disp(sprintf(' %% %s not found',XA_fname_snapshot_mat)); end;
if (~exist(XA_fname_snapshot_jpg,'file')); disp(sprintf(' %% %s not found',XA_fname_snapshot_jpg)); end;
if (~exist(XA_fname_compare_image_rank_mat,'file')); disp(sprintf(' %% %s not found',XA_fname_compare_image_rank_mat)); end;
if (~exist(XA_fname_compare_image_rank_jpg,'file')); disp(sprintf(' %% %s not found',XA_fname_compare_image_rank_jpg)); end;
if (~exist(XA_fname_align_crop_a_k_Y_mat,'file')); disp(sprintf(' %% %s not found',XA_fname_align_crop_a_k_Y_mat)); end;
if (~exist(XA_fname_fsc_crop_kx__mat,'file')); disp(sprintf(' %% %s not found',XA_fname_fsc_crop_kx__mat)); end;
end;%if (verbose>2);
%%%%;
if flag_clean;
if ( exist(XA_fname_align_a_k_Y_mat,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_align_a_k_Y_mat));
if flag_clean_really; delete(XA_fname_align_a_k_Y_mat); end; end;
if ( exist(XA_fname_align_a_k_Y_jpg,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_align_a_k_Y_jpg));
if flag_clean_really; delete(XA_fname_align_a_k_Y_jpg); end; end;
if ( exist(XA_fname_snapshot_mat,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_snapshot_mat));
if flag_clean_really; delete(XA_fname_snapshot_mat); end; end;
if ( exist(XA_fname_snapshot_jpg,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_snapshot_jpg));
if flag_clean_really; delete(XA_fname_snapshot_jpg); end; end;
if ( exist(XA_fname_compare_image_rank_mat,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_compare_image_rank_mat));
if flag_clean_really; delete(XA_fname_compare_image_rank_mat); end; end;
if ( exist(XA_fname_compare_image_rank_jpg,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_compare_image_rank_jpg));
if flag_clean_really; delete(XA_fname_compare_image_rank_jpg); end; end;
if ( exist(XA_fname_align_crop_a_k_Y_mat,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_align_crop_a_k_Y_mat));
if flag_clean_really; delete(XA_fname_align_crop_a_k_Y_mat); end; end;
if ( exist(XA_fname_fsc_crop_kx__mat,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_fsc_crop_kx__mat));
if flag_clean_really; delete(XA_fname_fsc_crop_kx__mat); end; end;
end;%if flag_clean;
%%%%;
end;%if ( exist(XA_fname_mat,'file'));
end;%for flag_xcor_vs_Memp=0:1;
end;%for flag_rank_vs_tolerance=0:1;
end;%for flag_alternate_MS_vs_SM = [1];
end;%for ntolerance_pm=0:n_tolerance_pm-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for ndat_rseed=0:n_dat_rseed-1;

verbose=1;
if (verbose); disp(sprintf(' %% [finished test_pm_23_clean_20220307]')); end;

disp('returning'); return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
