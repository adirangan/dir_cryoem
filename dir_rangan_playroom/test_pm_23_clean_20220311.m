function ...
[ ...
 global_parameter ...
] = ...
test_pm_23_clean_20220311( ...
 global_parameter ...
,fname_prefix ...
,dir_nopath_data_star ...
,Pixel_Spacing ...
,fname_nopath_volume ...
,fname_nopath_star ...
);

%{

scp -p dir_ISWINCP_x0/dir_pm_mat/delta_r_max_legacy.txt dir_ISWINCP_nM0000to0255/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_ISWINCP_x0/dir_pm_mat/delta_r_max_legacy.txt dir_ISWINCP_nM0000to0511/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_ISWINCP_x0/dir_pm_mat/delta_r_max_legacy.txt dir_ISWINCP_nM0000to0767/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_trpv1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_trpv1_nM0000to0255/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_trpv1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_trpv1_nM0000to0511/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_trpv1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_trpv1_nM0000to0767/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_rib80s_x0/dir_pm_mat/delta_r_max_legacy.txt dir_rib80s_nM0000to0255/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_rib80s_x0/dir_pm_mat/delta_r_max_legacy.txt dir_rib80s_nM0000to0511/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_rib80s_x0/dir_pm_mat/delta_r_max_legacy.txt dir_rib80s_nM0000to0767/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_MlaFEDB_x0/dir_pm_mat/delta_r_max_legacy.txt dir_MlaFEDB_nM0000to0255/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_MlaFEDB_x0/dir_pm_mat/delta_r_max_legacy.txt dir_MlaFEDB_nM0000to0511/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_MlaFEDB_x0/dir_pm_mat/delta_r_max_legacy.txt dir_MlaFEDB_nM0000to0767/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_LetB1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_LetB1_nM0000to0255/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_LetB1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_LetB1_nM0000to0511/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_LetB1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_LetB1_nM0000to0767/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_TMEM16F_x0/dir_pm_mat/delta_r_max_legacy.txt dir_TMEM16F_nM0000to0255/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_TMEM16F_x0/dir_pm_mat/delta_r_max_legacy.txt dir_TMEM16F_nM0000to0511/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_TMEM16F_x0/dir_pm_mat/delta_r_max_legacy.txt dir_TMEM16F_nM0000to0767/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_LSUbl17dep_x0/dir_pm_mat/delta_r_max_legacy.txt dir_LSUbl17dep_nM0000to0255/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_LSUbl17dep_x0/dir_pm_mat/delta_r_max_legacy.txt dir_LSUbl17dep_nM0000to0511/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_LSUbl17dep_x0/dir_pm_mat/delta_r_max_legacy.txt dir_LSUbl17dep_nM0000to0767/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_ps1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_ps1_nM0000to0255/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_ps1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_ps1_nM0000to0511/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_ps1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_ps1_nM0000to0767/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_ps1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_ps1_x1/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_ps1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_ps1_x2/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_ps1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_ps1_x3/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_ps1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_ps1_x4/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_ps1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_ps1_x5/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_ps1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_ps1_x6/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_ps1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_ps1_x7/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_ps1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_ps1_x0to7/dir_pm_mat/delta_r_max_legacy.txt
scp -p dir_ps1_x0/dir_pm_mat/delta_r_max_legacy.txt dir_ps1_x0to7_combine/dir_pm_mat/delta_r_max_legacy.txt
  
  %}

if (nargin<1);
table_data__ = { ...
'ISWINCP_nM0000to0255' , 'ISWINCP' , 1.07 , 'emd_9718.map' , 'ADPBeF.star' ; ...
'trpv1_nM0000to0255' , 'trpv1' , 1.2156 , 'emd_5778.mrc' , 'tv1_relion_data.star' ; ...
'rib80s_nM0000to0255' , 'rib80s' , 1.34 , 'emd_2660.mrc' , 'shiny_2sets.star' ; ...
'MlaFEDB_nM0000to0255' , 'MlaFEDB' , 0.832 , 'emd_22116.map' , 'Empiar_10536_00_to_23.star' ; ...
'LetB1_nM0000to0255' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1.star' ; ...
'TMEM16F_nM0000to0255' , 'TMEM16F' , 1.059 , 'emd_20244.map' , 'All_8192.star' ; ...
'LSUbl17dep_nM0000to0255' , 'LSUbl17dep' , 1.31 , 'emd_8434.map' , 'Parameters_negated.star' ; ...
'ps1_nM0000to0255' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
'ISWINCP_nM0000to0511' , 'ISWINCP' , 1.07 , 'emd_9718.map' , 'ADPBeF.star' ; ...
'trpv1_nM0000to0511' , 'trpv1' , 1.2156 , 'emd_5778.mrc' , 'tv1_relion_data.star' ; ...
'rib80s_nM0000to0511' , 'rib80s' , 1.34 , 'emd_2660.mrc' , 'shiny_2sets.star' ; ...
'MlaFEDB_nM0000to0511' , 'MlaFEDB' , 0.832 , 'emd_22116.map' , 'Empiar_10536_00_to_23.star' ; ...
'LetB1_nM0000to0511' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1.star' ; ...
'TMEM16F_nM0000to0511' , 'TMEM16F' , 1.059 , 'emd_20244.map' , 'All_8192.star' ; ...
'LSUbl17dep_nM0000to0511' , 'LSUbl17dep' , 1.31 , 'emd_8434.map' , 'Parameters_negated.star' ; ...
'ps1_nM0000to0511' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
'ISWINCP_nM0000to0767' , 'ISWINCP' , 1.07 , 'emd_9718.map' , 'ADPBeF.star' ; ...
'trpv1_nM0000to0767' , 'trpv1' , 1.2156 , 'emd_5778.mrc' , 'tv1_relion_data.star' ; ...
'rib80s_nM0000to0767' , 'rib80s' , 1.34 , 'emd_2660.mrc' , 'shiny_2sets.star' ; ...
'MlaFEDB_nM0000to0767' , 'MlaFEDB' , 0.832 , 'emd_22116.map' , 'Empiar_10536_00_to_23.star' ; ...
'LetB1_nM0000to0767' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1.star' ; ...
'TMEM16F_nM0000to0767' , 'TMEM16F' , 1.059 , 'emd_20244.map' , 'All_8192.star' ; ...
'LSUbl17dep_nM0000to0767' , 'LSUbl17dep' , 1.31 , 'emd_8434.map' , 'Parameters_negated.star' ; ...
'ps1_nM0000to0767' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
'ps1_x1' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
'ps1_x2' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
'ps1_x3' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
'ps1_x4' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
'ps1_x5' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
'ps1_x6' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
'ps1_x7' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
'ps1_x0to7' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
'ps1_x0to7_combine' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
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
test_pm_23_clean_20220311( ...
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

verbose=1;
if (verbose); disp(sprintf(' %% [entering test_pm_23_clean_20220311]')); end;

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

flag_clean = 1;
flag_clean_really = 0;

delta_sigma_use = 0.05;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
tolerance_pm_ = 0.1.^[2:0.5:3];
n_tolerance_pm = numel(tolerance_pm_);
rank_pm_ = [16,18];
n_rank_pm = numel(rank_pm_);
%delta_r_max_factor_ = [0.00,0.125,0.25,0.50];
delta_r_max_factor_ = 1.00;
n_delta_r_max_factor = numel(delta_r_max_factor_);
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
%{
if ( exist(sprintf('%s_mat/delta_r_max_legacy.txt',dir_pm)));
delta_r_max_legacy = textread(sprintf('%s_mat/delta_r_max_legacy.txt',dir_pm));
disp(sprintf(' %% loading delta_r_max_legacy: %0.16f',delta_r_max_legacy));
delta_r_max_use = delta_r_max_factor * delta_r_max_legacy;
end;%if ( exist(sprintf('%s_mat/delta_r_max_legacy.txt',dir_pm)));
 %}
for flag_alternate_MS_vs_SM = [1];%for flag_alternate_MS_vs_SM = [0:1];
for flag_rank_vs_tolerance=[0:1];
%%%%%%%%;
if flag_rank_vs_tolerance==0;
for ntolerance_pm=0:n_tolerance_pm-1;
tolerance_pm = tolerance_pm_(1+ntolerance_pm);
str_strategy = ''; if (flag_alternate_MS_vs_SM); str_strategy = sprintf('%sa1',str_strategy); end;
str_delta_r_max = sprintf('t%.4d',floor(1000*delta_r_max_use));
str_tolerance_pm = sprintf('p%.2d',floor(10*-log10(tolerance_pm)));
str_rseed = sprintf('r%d',dat_rseed);
if (flag_rank_vs_tolerance==0); str_xfix = sprintf('%s%s%s%s',str_strategy,str_delta_r_max,str_tolerance_pm,str_rseed); end;
if (flag_rank_vs_tolerance==1); str_xfix = sprintf('%s%s%s%s',str_strategy,str_delta_r_max,str_rank_pm,str_rseed); end;
test_pm_23_clean_helper_0(verbose,flag_clean,flag_clean_really,dir_pm,str_xfix);
end;%for ntolerance_pm=0:n_tolerance_pm-1;
end;%if flag_rank_vs_tolerance==0;
%%%%%%%%;
if flag_rank_vs_tolerance==1;
for nrank_pm=0:n_rank_pm-1;
rank_pm = rank_pm_(1+nrank_pm);
str_strategy = ''; if (flag_alternate_MS_vs_SM); str_strategy = sprintf('%sa1',str_strategy); end;
str_delta_r_max = sprintf('t%.4d',floor(1000*delta_r_max_use));
str_rank_pm = sprintf('n%.2d',rank_pm);
str_rseed = sprintf('r%d',dat_rseed);
if (flag_rank_vs_tolerance==0); str_xfix = sprintf('%s%s%s%s',str_strategy,str_delta_r_max,str_tolerance_pm,str_rseed); end;
if (flag_rank_vs_tolerance==1); str_xfix = sprintf('%s%s%s%s',str_strategy,str_delta_r_max,str_rank_pm,str_rseed); end;
test_pm_23_clean_helper_0(verbose,flag_clean,flag_clean_really,dir_pm,str_xfix);
end;%for nrank_pm=0:n_rank_pm-1;
end;%if flag_rank_vs_tolerance==1;
%%%%%%%%;
end;%for flag_rank_vs_tolerance=[0:1];
end;%for flag_alternate_MS_vs_SM = [1];
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for ndat_rseed=0:n_dat_rseed-1;

verbose=1;
if (verbose); disp(sprintf(' %% [finished test_pm_23_clean_20220311]')); end;

disp('returning'); return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
