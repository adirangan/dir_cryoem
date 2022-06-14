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
% 0 , 'p28hRPT1_x0' , 'p28hRPT1' , 0.98 , 'emd_8674.map' , 'T1.star' ; ...
 0 , 'trpv1_x0' , 'trpv1' , 1.2156 , 'emd_5778.mrc' , 'tv1_relion_data.star' ; ...
% 1 , 'trpv1c'   , 'trpv1' , 1.2156 , 'emd_5778.mrc' , 'tv1_relion_data.star' ; ...
 0 , 'rib80s_x0' , 'rib80s' , 1.34 , 'emd_2660.mrc' , 'shiny_2sets.star' ; ...
% 1 , 'rib80s' , 'rib80s' , 1.34 , 'emd_2660.mrc' , 'shiny_2sets.star' ; ...
 0 , 'MlaFEDB_x0' , 'MlaFEDB' , 0.832 , 'emd_22116.map' , 'Empiar_10536_00_to_23.star' ; ...
% 1 , 'MlaFEDB' , 'MlaFEDB' , 0.832 , 'emd_22116.map' , 'Empiar_10536_00_to_23.star' ; ...
 0 , 'ISWINCP_x0' , 'ISWINCP' , 1.07 , 'emd_9718.map' , 'ADPBeF.star' ; ...
% 1 , 'ISWINCP' , 'ISWINCP' , 1.07 , 'emd_9718.map' , 'ADPBeF.star' ; ...
 0 , 'TMEM16F_x0' , 'TMEM16F' , 1.059 , 'emd_20244.map' , 'All_8192.star' ; ...
% 1 , 'TMEM16F' , 'TMEM16F' , 1.059 , 'emd_20244.map' , 'All_8192.star' ; ...
 0 , 'LSUbl17dep_x0' , 'LSUbl17dep' , 1.31 , 'emd_8434.map' , 'Parameters.star' ; ...
% 1 , 'LSUbl17dep' , 'LSUbl17dep' , 1.31 , 'emd_8434.map' , 'Parameters.star' ; ...
 0 , 'ps1_x0' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
% 1 , 'ps0' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
% 1 , 'ps1' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
 0 , 'LetB1_x0' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1_350MB.star' ; ...
% 1 , 'LetB1' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1_350MB.star' ; ...
};
n_experiment = size(table_data__,1);
collect__ = cell(n_experiment,1);
collect_aX__ = cell(n_experiment,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_ = clock;rng(tmp_(end));
for nexperiment=(randperm(n_experiment)-1);
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

test_pm_collect_excerpt_fsc_x_0( ...
 struct('Pixel_Spacing',Pixel_Spacing) ...
,dir_pm ...
);

test_pm_collect_excerpt_fsc_xa__0( ...
 struct('Pixel_Spacing',Pixel_Spacing) ...
,dir_pm ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nexperiment=0:n_experiment-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Compare fsc from optimal-reconstruction with those from experiments. ;
%%%%%%%%;
fname_fig = sprintf('%s/test_pm_collect_plus_fsc_FIGI',dir_jpg);

disp('returning');return;


