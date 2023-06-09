%%%%%%%%;
% Derived from test_viewing_angle_distribution_0.m ;
% Also tries to plot figures for relion, etc. ;
%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;
dir_base = sprintf('/%s/rangan/dir_cryoem',string_root);
flag_recalc = 0;
flag_replot = 0;
flag_verbose = 1; nf=0;

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
for nexperiment=(randperm(n_experiment)-1);
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
global_parameter.flag_recalc=flag_recalc;
global_parameter.flag_replot=flag_replot;
global_parameter.flag_invert=0;
global_parameter.flag_center_image=0;
global_parameter.tolerance_master = 1e-2;
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); global_parameter.flag_invert = 0; end;
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); global_parameter.flag_center_image = 1; end;
if (strcmp(dir_nopath_data_star,'precatalytic_spliceosome')); global_parameter.flag_center_image = 1; end;
dir_pm = sprintf('%s/dir_%s/dir_pm',dir_base,fname_prefix);
if (flag_verbose); disp(sprintf(' %% fname_prefix: %s',fname_prefix)); end;
[~,ampm_] = ampm_fromdisk_3(struct('type','parameter','flag_exclude_ampm',1),dir_pm);
ampm__{1+nexperiment} = ampm_;

tolerance_master = global_parameter.tolerance_master;
flag_recalc = global_parameter.flag_recalc;
flag_replot = global_parameter.flag_replot;
N_k_p_use__ = ampm_.M_k_p__; if (global_parameter.flag_center_image==1); N_k_p_use__ = ampm_.N_k_p__; end;
tmp_delta_r_max = 0.030;
tmp_delta_r_upb = 0.250;
n_S = ampm_.n_S; n_M = ampm_.n_M; 
n_iteration_seco = 6;

str_exte_type = 'relion';
%%%%%%%%%%%%%%%%;
if strcmp(str_exte_type,'relion'); 
dir_exte = sprintf('%s/dir_%s/dir_relion_mat/job_1024',dir_base,dir_nopath_data_star); 
fname_exte_mrc = sprintf('%s/run_it300_class001.mrc',dir_exte); 
end;%if strcmp(str_exte_type,'relion');
fname_exte_pre = sprintf('%s/test_viewing_angle_distribution_%s',dir_exte,str_exte_type);
[flag_exte_skip,fname_exte_mat] = open_fname_tmp(fname_exte_pre);
if flag_recalc | ~flag_exte_skip;
%%%%%%%%%%%%%%%%;
% set default indices. ;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
n_x_u_pack = 64;
%%%%%%%%;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;
xxx_u_weight_ = (2*x_p_r_max/n_x_u_pack)^3;
%%%%%%%%;
% load stored indices. ;
%%%%%%%%;
n_w_ = ampm_.n_w_;
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
weight_3d_k_p_r_ = ampm_.weight_3d_k_p_r_;
l_max_ = ampm_.l_max_;
n_lm_ = ampm_.n_lm_;
n_lm_max = ampm_.n_lm_max;
n_lm_sum = ampm_.n_lm_sum;
n_lm_csum_ = ampm_.n_lm_csum_;
l_max_max = ampm_.l_max_max;
m_max_ = ampm_.m_max_;
n_m_max = ampm_.n_m_max;
Y_l_val_ = ampm_.Y_l_val_;
Y_m_val_ = ampm_.Y_m_val_;
Y_k_val_ = ampm_.Y_k_val_;
weight_Y_ = ampm_.weight_Y_;
weight_2d_k_p_r_ = ampm_.weight_2d_k_p_r_;
weight_2d_k_all_ = ampm_.weight_2d_k_all_;
%%%%%%%%;
n_k_all = ampm_.n_k_all;
weight_3d_k_all_ = ampm_.weight_3d_k_all_;
k_c_0_all_ = ampm_.k_c_0_all_;
k_c_1_all_ = ampm_.k_c_1_all_;
k_c_2_all_ = ampm_.k_c_2_all_;
k_p_r_max = ampm_.k_p_r_max;
n_x_u_pack = ampm_.n_x_u_pack;
n_xxx_u = ampm_.n_xxx_u;
xxx_u_weight_ = ampm_.xxx_u_weight_;
x_u_0_ = ampm_.x_u_0_;
x_u_1_ = ampm_.x_u_1_;
x_u_2_ = ampm_.x_u_2_;
x_u_0___ = ampm_.x_u_0___;
x_u_1___ = ampm_.x_u_1___;
x_u_2___ = ampm_.x_u_2___;
x_p_r_max = ampm_.x_p_r_max;
%%%%%%%%;
% Now cluster the CTF based on tolerance_cluster. ;
%%%%%%%%;
parameter_CTF = struct('type','parameter');
parameter_CTF.tolerance_master = tolerance_master;
[ ...
 parameter_CTF ...
,ampm_.index_ncluster_from_nCTF_ ...
] = ...
knn_cluster_CTF_k_p_r_kC__0( ...
 parameter_CTF ...
,ampm_.n_k_p_r ...
,ampm_.k_p_r_ ...
,ampm_.weight_2d_k_p_r_ ...
,ampm_.n_CTF ...
,ampm_.CTF_k_p_r_kC__ ...
);
%%%%%%%%;
ampm_.n_cluster = 1+max(ampm_.index_ncluster_from_nCTF_);
%%%%%%%%;
% load volume. ;
%%%%%%%%;
a_x_u_exte_ = cast(ReadMRC(fname_exte_mrc),'double');
if ndims(a_x_u_exte_)< 3; n_x_u_exte = round(numel(a_x_u_exte_).^(1/3)); a_x_u_exte_ = reshape(a_x_u_exte_,[n_x_u_exte,n_x_u_exte,n_x_u_exte]); end;
n_x_u_exte = size(a_x_u_exte_,1);
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
n_pack_from_exte = n_x_u_exte/n_x_u_pack;
tmp_pack_row_ij_ = zeros(n_x_u_pack,1);
tmp_pack_col_ij_ = zeros(n_x_u_pack,1);
tmp_pack_val_ij_ = zeros(n_x_u_pack,1);
na=0;
for nx_u=0:n_x_u_exte-1;
tmp_pack_row_ij_(1+na) = 1+nx_u;
tmp_pack_col_ij_(1+na) = 1+floor(nx_u/n_pack_from_exte);
tmp_pack_val_ij_(1+na) = 1/n_pack_from_exte;
na=na+1;
end;%for nx_u=0:n_x_u_exte-1;
x_u_pack_from_exte_ = sparse(tmp_pack_row_ij_,tmp_pack_col_ij_,tmp_pack_val_ij_,n_x_u_exte,n_x_u_pack);
a_x_u_pack_from_exte_ = reshape(a_x_u_exte_,[n_x_u_exte*n_x_u_exte,n_x_u_exte])*x_u_pack_from_exte_;
a_x_u_pack_from_exte_ = reshape(permute(reshape(a_x_u_pack_from_exte_,[n_x_u_exte,n_x_u_exte,n_x_u_pack]),[3,1,2]),[n_x_u_exte*n_x_u_pack,n_x_u_exte])*x_u_pack_from_exte_;
a_x_u_pack_from_exte_ = reshape(permute(reshape(a_x_u_pack_from_exte_,[n_x_u_pack,n_x_u_exte,n_x_u_pack]),[3,1,2]),[n_x_u_pack*n_x_u_pack,n_x_u_exte])*x_u_pack_from_exte_;
a_x_u_pack_from_exte_ = permute(reshape(a_x_u_pack_from_exte_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),[3,1,2]);
%%%%%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_exte_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_x_u_pack_from_exte_(:).*xxx_u_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_p_exte_ time %0.2fs',tmp_t));
%%%%%%%%;
tmp_t = tic;
[a_k_Y_exte_] = convert_k_p_to_spharm_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_exte_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_exte_ time %0.2fs',tmp_t));
%%%%%%%%;
tmp_t = tic;
[ ... 
 a_k_Y_exte_alig_ ...
,X_best_exte_alig ...
] = ...
spharm_register_and_rotate_2( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,ampm_.a_k_Y_quad_ ...
,a_k_Y_exte_ ...
,0 ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_exte_alig_ time %0.2fs',tmp_t));
%%%%%%%%;
if ~exist('Ylm_klma___','var'); Ylm_klma___ = []; end;
tmp_t = tic();
[ ...
 a_k_p_exte_alig_  ...
,Ylm_klma___ ...
] =  ...
convert_spharm_to_k_p_3( ...
 0*flag_verbose ...
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
,a_k_Y_exte_alig_ ...
,Ylm_klma___ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% convert_spharm_to_k_p_3: %0.6fs',tmp_t));
%%%%%%%%;
tmp_t = tic(); 
[ ...
 a_x_u_exte_alig_ ...
] = ...
convert_k_p_to_x_c_1( ...
 0*flag_verbose ...
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
,a_k_p_exte_alig_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% convert_k_p_to_x_c_1: %0.6fs',tmp_t));
%%%%%%%%;
% Now align M_k_p__ to a_k_Y_exte_. ;
%%%%%%%%;
tmp_t = tic();
parameter_align = struct('type','parameter');
parameter_align.fname_align_a_k_Y_pre = ''; %<-- empty for now. ;
parameter_align.delta_r_max = tmp_delta_r_max;
parameter_align.n_delta_v_requested = 24;
parameter_align.delta_r_upb = tmp_delta_r_upb;
parameter_align.n_iteration = n_iteration_seco;
[ ...
 parameter_align ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,X_exte_alig_SMi___ ...
,~ ...
,~ ...
,~ ...
,~ ...
] = ...
pm_align_M_k_p_to_a_k_Y_3( ...
 parameter_align ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,N_k_p_use__ ...
,ampm_.n_cluster ...
,ampm_.index_ncluster_from_nCTF_ ...
,ampm_.n_CTF ...
,ampm_.index_nCTF_from_nM_ ...
,ampm_.CTF_k_p_r_kC__ ...
,l_max_ ...
,a_k_Y_exte_alig_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_align_M_k_p_to_a_k_Y_3: %0.2fs',tmp_t));
%%%%%%%%;
% load true X_SM__. ;
%%%%%%%%;
tmp_str_filter_ = ls(sprintf('%s/*viewing_angle_distribution.mat',ampm_.dir_pm_mat));
tmp_index_start_ = strfind(tmp_str_filter_,sprintf('/%s',string_root))-1;
tmp_index_final_ = strfind(tmp_str_filter_,'.mat')+4-1;
tmp_str_filter_ = tmp_str_filter_(1+tmp_index_start_(1):1+tmp_index_final_(1)-1);
tmp_true_ = load(tmp_str_filter_,'X_SM__');
if (flag_verbose> 0); disp(sprintf(' %% X_SM__ <-- %s',tmp_str_filter_)); end;
X_true_SM__ = tmp_true_.X_SM__;
X_exte_alig_SM__ = X_exte_alig_SMi___(:,:,end);
%%%%;
P_true_SM__ = X_true_SM__;
for nM=0:n_M-1;
tmp_S_ = P_true_SM__(:,1+nM);
[~,tmp_ij_] = sort(tmp_S_); [~,tmp_ij_] = sort(tmp_ij_);
P_true_SM__(:,1+nM) = tmp_ij_;
end;%for nM=0:n_M-1;
%%%%;
P_exte_alig_SM__ = X_exte_alig_SM__;
for nM=0:n_M-1;
tmp_S_ = P_exte_alig_SM__(:,1+nM);
[~,tmp_ij_] = sort(tmp_S_); [~,tmp_ij_] = sort(tmp_ij_);
P_exte_alig_SM__(:,1+nM) = tmp_ij_;
end;%for nM=0:n_M-1;
%%%%;
Q_true_SM__ = X_true_SM__;
for nS=0:n_S-1;
tmp_S_ = Q_true_SM__(1+nS,:);
[~,tmp_ij_] = sort(tmp_S_); [~,tmp_ij_] = sort(tmp_ij_);
Q_true_SM__(1+nS,:) = tmp_ij_;
end;%for nS=0:n_S-1;
%%%%;
Q_exte_alig_SM__ = X_exte_alig_SM__;
for nS=0:n_S-1;
tmp_S_ = Q_exte_alig_SM__(1+nS,:);
[~,tmp_ij_] = sort(tmp_S_); [~,tmp_ij_] = sort(tmp_ij_);
Q_exte_alig_SM__(1+nS,:) = tmp_ij_;
end;%for nS=0:n_S-1;
%%%%;
CX_S_ = zeros(n_S,1);
CP_S_ = zeros(n_S,1);
CQ_S_ = zeros(n_S,1);
for nS=0:n_S-1;
CX_S_(1+nS) = corr(transpose(X_true_SM__(1+nS,:)),transpose(X_exte_alig_SM__(1+nS,:)));
CP_S_(1+nS) = corr(transpose(P_true_SM__(1+nS,:)),transpose(P_exte_alig_SM__(1+nS,:)));
CQ_S_(1+nS) = corr(transpose(Q_true_SM__(1+nS,:)),transpose(Q_exte_alig_SM__(1+nS,:)));
end;%for nS=0:n_S-1;
CX_M_ = zeros(n_M,1);
CP_M_ = zeros(n_M,1);
CQ_M_ = zeros(n_M,1);
for nM=0:n_M-1;
CX_M_(1+nM) = corr((X_true_SM__(:,1+nM)),(X_exte_alig_SM__(:,1+nM)));
CP_M_(1+nM) = corr((P_true_SM__(:,1+nM)),(P_exte_alig_SM__(:,1+nM)));
CQ_M_(1+nM) = corr((Q_true_SM__(:,1+nM)),(Q_exte_alig_SM__(:,1+nM)));
end;%for nM=0:n_M-1;
%%%%%%%%;
save(fname_exte_mat ...
     ,'str_exte_type','' ...
     ,'a_k_Y_exte_alig_','X_best_exte_alig','a_k_p_exte_alig_','a_x_u_exte_alig_' ...
     ,'X_true_SM__','X_exte_alig_SM__' ...
     ,'P_true_SM__','P_exte_alig_SM__' ...
     ,'Q_true_SM__','Q_exte_alig_SM__' ...
     ,'CX_S_','CP_S_','CQ_S_' ...
     ,'CX_M_','CP_M_','CQ_M_' ...
     ,'-v7.3' ...
     );
close_fname_tmp(fname_exte_pre);
end;%if flag_recalc | ~flag_exte_skip;
%%%%%%%%%%%%%%%%;
if  exist(fname_exte_mat,'file');
%%%%%%%%%%%%%%%%;
tmp_exte_ = load(fname_exte_mat);
if strfind(ampm_.dir_pm,'p28hRPT1');   percent_threshold_ = 95.25; tmp_nx = 5; end;
if strfind(ampm_.dir_pm,'ISWINCP');    percent_threshold_ = 94.50; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'trpv1');      percent_threshold_ = 91.25; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'rib80s');     percent_threshold_ = 86.25; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'MlaFEDB');    percent_threshold_ = 95.00; tmp_nx = 8; end;
if strfind(ampm_.dir_pm,'LetB1');      percent_threshold_ = 91.75; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'TMEM16F');    percent_threshold_ = 94.50; tmp_nx = 16; end;
if strfind(ampm_.dir_pm,'LSUbl17dep'); percent_threshold_ = 86.25; tmp_nx = 15; end;
if strfind(ampm_.dir_pm,'ps1');        percent_threshold_ = 96.00; tmp_nx = 11; end;
n_x_u_pack = ampm_.n_x_u_pack;
tmp_window_ = zeros(n_x_u_pack,n_x_u_pack,n_x_u_pack);
tmp_index_ = tmp_nx:n_x_u_pack-1-tmp_nx;
tmp_window_(1+tmp_index_,1+tmp_index_,1+tmp_index_)=1;
tmp_index_ = efind(tmp_window_);
%%%%%%%%;
fname_exte_fig_pre = sprintf('%s/%s_%s_viewing_angle_distribution_FIGABCD',dir_exte,dir_nopath_data_star,tmp_exte_.str_exte_type);
fname_exte_fig_jpg = sprintf('%s.jpg',fname_exte_fig_pre);
if (flag_replot>=1) | ~exist(fname_exte_fig_jpg,'file');
if (flag_verbose); disp(sprintf(' %% %s not found, creating',fname_exte_fig_jpg)); end;
n_S = ampm_.n_S; n_M = ampm_.n_M; 
figure(1+nf);nf=nf+1;clf;figmed;
fontsize_use = 16;
tmp_percent_threshold = percent_threshold_(1);
n_h = 128; dh = 1/max(1,n_h); n_z = 1.0/max(1,n_h*n_h + numel(tmp_exte_.X_true_SM__))/(dh*dh); %<-- now integrates to 1. ;
hlim_ = [-4,+4];
%%%%;
subplot_{1} = subplot(1,4,1);
isosurface_f_x_u_0(tmp_exte_.a_x_u_exte_alig_(1+tmp_index_),tmp_percent_threshold);
%title(sprintf('%s',tmp_exte_.str_exte_type)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
xlabel('');ylabel('');zlabel(sprintf('%s',tmp_exte_.str_exte_type)); set(gca,'FontSize',fontsize_use);
%%%%;
subplot_{2} = subplot(1,4,2);
imagesc(log(n_z) + log(1+hist2d_0(tmp_exte_.X_true_SM__(:),tmp_exte_.X_exte_alig_SM__(:),n_h,n_h,[0,+1],[0,+1])),hlim_);
set(gca,'Ydir','normal'); tmp_c_= colorbar; set(tmp_c_,'Ticks',[-4,+4]);
xlabel('corr (true)','Interpreter','none');
ylabel(sprintf('corr (%s)',tmp_exte_.str_exte_type),'Interpreter','none');
title('Correlation','Interpreter','none');
axis image; set(gca,'XTick',[1,n_h],'XTickLabel',[0,1],'YTick',[1,n_h],'YTickLabel',[0,1]);
set(gca,'FontSize',fontsize_use);
%%%%;
subplot_{3} = subplot(1,4,3);
imagesc(log(n_z) + log(1+hist2d_0(tmp_exte_.P_true_SM__(:),tmp_exte_.P_exte_alig_SM__(:),n_h,n_h,[1,n_S],[1,n_S])),hlim_);
set(gca,'Ydir','normal'); tmp_c_= colorbar; set(tmp_c_,'Ticks',[-4,+4]);
xlabel('rank (true)','Interpreter','none');
ylabel(sprintf('rank (%s)',tmp_exte_.str_exte_type),'Interpreter','none');
title('template rank','Interpreter','none');
axis image; set(gca,'XTick',[1,n_h],'XTickLabel',[0,1],'YTick',[1,n_h],'YTickLabel',[0,1]);
set(gca,'FontSize',fontsize_use);
%%%%;
subplot_{4} = subplot(1,4,4);
imagesc(log(n_z) + log(1+hist2d_0(tmp_exte_.Q_true_SM__(:),tmp_exte_.Q_exte_alig_SM__(:),n_h,n_h,[1,n_M],[1,n_M])),hlim_);
set(gca,'Ydir','normal'); tmp_c_= colorbar; set(tmp_c_,'Ticks',[-4,+4]);
xlabel('rank (true)','Interpreter','none');
ylabel(sprintf('rank (%s)',tmp_exte_.str_exte_type),'Interpreter','none');
title('image rank','Interpreter','none');
axis image; set(gca,'XTick',[1,n_h],'XTickLabel',[0,1],'YTick',[1,n_h],'YTickLabel',[0,1]);
set(gca,'FontSize',fontsize_use);
%%%%;
colormap(subplot_{2},colormap_80s);
colormap(subplot_{3},colormap_80s);
colormap(subplot_{4},colormap_80s);
set(gcf,'Position',1+[0,0,512*4,512]);
%%%%;
sgtitle(fname_exte_fig_pre,'Interpreter','none');
print('-djpeg',fname_exte_fig_jpg);
sgtitle('');
tmp_dir = sprintf('/%s/rangan/dir_cryoem/dir_ampm_manuscript/dir_ampm_exte_fig',string_root);
if ~exist(tmp_dir,'dir'); disp(sprintf(' %% mkdir %s',tmp_dir)); mkdir(tmp_dir); end;
fname_fig_jpg_strip = sprintf('%s/%s_%s_viewing_angle_distribution_FIGABCD.jpg',tmp_dir,dir_nopath_data_star,tmp_exte_.str_exte_type);
disp(sprintf(' %% writing %s',fname_fig_jpg_strip));
print('-djpeg',sprintf('%s',fname_fig_jpg_strip));
close gcf;
end;%if (flag_replot>=1) | ~exist(sprintf('%s.jpg',fname_exte_fig_jpg),'file');
%%%%%%%%;
%%%%%%%%%%%%%%%%;
end;%if  exist(fname_exte_mat,'file');
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nexperiment=0:n_experiment-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;





