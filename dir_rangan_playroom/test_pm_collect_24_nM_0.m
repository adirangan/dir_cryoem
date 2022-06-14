%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;
dir_base = sprintf('/%s/rangan/dir_cryoem',string_root);
flag_verbose = 1; nf=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
table_data__ = { ...
'ISWINCP' , 'ISWINCP' , 1.07 , 'emd_9718.map' , 'ADPBeF.star' ; ...
%'p28hRPT1' , 'p28hRP' , 0.98 , 'emd_8674.map' , 'T1.star' ; ...
'trpv1' , 'trpv1' , 1.2156 , 'emd_5778.mrc' , 'tv1_relion_data.star' ; ...
'rib80s' , 'rib80s' , 1.34 , 'emd_2660.mrc' , 'shiny_2sets.star' ; ...
'MlaFEDB' , 'MlaFEDB' , 0.832 , 'emd_22116.map' , 'Empiar_10536_00_to_23.star' ; ...
'LetB1' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1.star' ; ...
'TMEM16F' , 'TMEM16F' , 1.059 , 'emd_20244.map' , 'All_8192.star' ; ...
'LSUbl17dep' , 'LSUbl17dep' , 1.31 , 'emd_8434.map' , 'Parameters_negated.star' ; ...
'ps1' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
};
n_experiment = size(table_data__,1);
ampm___ = cell(n_experiment,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nexperiment=0:n_experiment-1;
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
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); global_parameter.flag_invert = 0; end;
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); global_parameter.flag_center_image = 1; end;
if (strcmp(dir_nopath_data_star,'precatalytic_spliceosome')); global_parameter.flag_center_image = 1; end;
n_M_ = [256,512,768,1024]; n_n_M = numel(n_M_);
ampm__ = cell(n_n_M,1);
for nn_M=0:n_n_M-1;
n_M = n_M_(1+nn_M);
global_parameter.nM_start = n_M*0;
global_parameter.nM_final = global_parameter.nM_start + n_M*1 - 1;
if (n_M==1024);
fname_prefix_plus_nM=sprintf('%s_x0',fname_prefix);
end;%if (n_M==1024);
if (n_M<1024);
fname_prefix_plus_nM=sprintf('%s_nM%.4dto%.4d',fname_prefix,global_parameter.nM_start,global_parameter.nM_final);
end;%if (n_M<1024);
dir_pm = sprintf('%s/dir_%s/dir_pm',dir_base,fname_prefix_plus_nM);
if (flag_verbose); disp(sprintf(' %% fname_prefix_plus_nM: %s',fname_prefix_plus_nM)); end;
[~,ampm_] = ampm_fromdisk_1(struct('type','parameter','flag_store_S_k_p__',0,'flag_store_M_k_p__',0),dir_pm);
ampm__{1+nn_M} = ampm_;
end;%for nn_M=0:n_n_M-1;
ampm___{1+nexperiment} = ampm__;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nexperiment=0:n_experiment-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

figure(1+nf);nf=nf+1;figbig;clf;
p_row = 2; p_col = ceil(n_experiment/p_row); np=0;
for nexperiment=0:n_experiment-1;
fname_prefix = table_data__{1+nexperiment,1+0};;
ampm__ = ampm___{1+nexperiment};
subplot(p_row,p_col,1+np);np=np+1;
str_symbol_ = {'o','^','s','p','h'}; n_str_symbol = numel(str_symbol_);
markersize_use = 8;
hold on;
%%%%;
for nn_M=0:n_n_M-1;
n_M = n_M_(1+nn_M);
ampm_ = ampm__{1+nn_M};
for na=0:ampm_.n_str_filter-1;
t_sval = ampm_.t_sval_a_(1+na);
p_vs_n = ampm_.p_vs_n_a_(1+na);
n_sval = ampm_.n_sval_a_(1+na);
p_sval = ampm_.p_sval_a_(1+na);
r_sval = ampm_.r_sval_a_(1+na);
X_best_ampm = ampm_.X_best_ampm_ia__(end,1+na);
str_symbol = str_symbol_{1+(p_vs_n==0)};
tmp_c = 'm'; if (p_vs_n==1); tmp_c = 'c'; end;
dx = 0.125-0.25*p_vs_n;
flag_include = ...
 (t_sval> 0) ...
  & ( ((p_vs_n==0)&&(n_sval>=14)) | ((p_vs_n==1)&&(p_sval>=20)) ) ...
  & ( X_best_ampm~=0 ) ...
  ;
if (flag_include); plot(nn_M+dx,X_best_ampm,str_symbol,'MarkerEdgeColor','k','MarkerFaceColor',tmp_c,'MarkerSize',markersize_use); end;
end;%for na=0:ampm_.n_str_filter-1;
end;%for nn_M=0:n_n_M-1;
%%%%;
hold off;
set(gca,'Xtick',[0:n_n_M-1],'XTickLabel',n_M_); xtickangle(90);
set(gca,'YTick',0:0.1:1.0);
xlabel('n_M','Interpreter','none');
ylabel('correlation (unmasked)');
xlim([0-0.5,n_n_M-1+0.5]);
ylim([-0.125,+1.125]);
grid on;
set(gca,'FontSize',12);
title(sprintf('%s',fname_prefix),'Interpreter','none');
drawnow();
end;%for nexperiment=0:n_experiment-1;
%%%%%%%%;
sgtitle(sprintf('test_pm_collect_24_nM_0'),'Interpreter','none');
fname_fig = sprintf('%s/dir_jpg/test_pm_collect_24_nM_FIGA',dir_base);
if (flag_verbose); disp(sprintf(' %% writing %s',fname_fig)); end;
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
