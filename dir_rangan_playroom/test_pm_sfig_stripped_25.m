%%%%%%%%;
% Stripped down version of test_pm_sfig_24.m ;
% Restricted to trpv1 (for now). ;
% Involves translations. ;
%%%%%%%%;

str_thisfunction = 'test_pm_sfig_stripped_25';

global_parameter=[];
fname_prefix='trpv1_x0';
dir_nopath_data_star='trpv1';
Pixel_Spacing=1.2156;
fname_nopath_volume='emd_5778.mrc';
fname_nopath_star='tv1_relion_data.star';

flag_verbose=1;
if (flag_verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

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

flag_exist = ...
    exist(sprintf('%s/%s',dir_data_star,fname_nopath_volume),'file') ...
  & exist(sprintf('%s/%s',dir_data_star,fname_nopath_star),'file') ...
  ;

if (~flag_exist);
disp(sprintf(' %% Warning, %s and %s not found',fname_nopath_volume,fname_nopath_star));
end;%if (~flag_exist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ( flag_exist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
disp(sprintf(' %% %s and %s found, processing...',fname_nopath_volume,fname_nopath_star));

[~,ampm] = ampm_fromdisk_3(struct('type','parameter','flag_exclude_ampm',1),dir_pm);
tolerance_master = 0.001;
rank_pm_ = [1:16,ampm.n_k_p_r-1]; n_rank_pm = numel(rank_pm_);
%%%%%%%%%%%%%%%%;
for delta_r_max = [ 0.010 , 0.025 , 0.050 , 0.075 , 0.100 ];
%%%%%%%%%%%%%%%%;
str_delta_r_max = sprintf('d%.4d',round(1000*delta_r_max));

%%%%%%%%;
fname_pre = sprintf('%s_mat/pm_fig_X_2d_xcor_d0_%s_FIGL_',dir_pm,str_delta_r_max);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre);
%%%%%%%%;
if ~flag_skip;
%%%%%%%%;
X_rwSM____ = cell(n_rank_pm,1);
for nrank_pm=0:n_rank_pm-1;
rank_pm = rank_pm_(1+nrank_pm);
parameter = struct('type','parameter');
parameter.tolerance_master = 0.001;
parameter.delta_r_max = delta_r_max;
parameter.flag_xcor_vs_Memp = 1; %<-- use volumetric principal-modes. ;
parameter.flag_rank_vs_tolerance = 1; %<-- use rank_pm;
parameter.rank_pm = rank_pm;
parameter.qbp_sheres_stop = 1; %<-- stop after innerproduct calculation. ;
parameter.flag_verbose = 1;
M_k_p_wkM__ = ampm.M_k_p__;
image_delta_x_prev_M_ = ampm.image_center_delta_x_c_0_;
image_delta_y_prev_M_ = ampm.image_center_delta_x_c_1_;
a_k_Y_base_yk_ = ampm.a_k_Y_quad_;
index_ncluster_from_nCTF_ = [];
%%%%;
tmp_t = tic();
[ ...
 parameter ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,X_rwSM____{1+nrank_pm} ...
] = ...
qbp_sheres_0( ...
 parameter ...
,ampm.n_k_p_r ...
,ampm.k_p_r_ ...
,ampm.k_p_r_max ...
,ampm.weight_3d_k_p_r_ ...
,ampm.weight_2d_k_p_r_ ...
,ampm.n_w_max ...
,ampm.l_max_ ...
,ampm.n_M ...
,M_k_p_wkM__ ...
,image_delta_x_prev_M_ ...
,image_delta_y_prev_M_ ...
,ampm.n_CTF ...
,ampm.index_nCTF_from_nM_ ...
,ampm.CTF_k_p_r_kC__ ...
,index_ncluster_from_nCTF_ ...
,a_k_Y_base_yk_ ...
);
tmp_t = toc(tmp_t);
if (flag_verbose); disp(sprintf(' %% qbp_sheres_0 [stop 1]; rank_pm %d: %0.6fs',rank_pm,tmp_t)); end;
end;%for nrank_pm=0:n_rank_pm-1;
%%%%%%%%;
n_S = ampm.n_S; n_M = ampm.n_M; n_w_max = ampm.n_w_max;
l2er_X_X_r_ = zeros(n_rank_pm,1);
corr_X_X_r_ = zeros(n_rank_pm,1);
prct_X_X_r_ = zeros(n_rank_pm,1);
X_wSMe___ = X_rwSM____{end};
for nrank_pm=0:n_rank_pm-1;
rank_pm = rank_pm_(1+nrank_pm);
X_wSM0___ = X_rwSM____{1+nrank_pm}; 
l2er_X_X = fnorm(X_wSMe___ - X_wSM0___) / fnorm(X_wSMe___);
l2er_X_X_r_(1+nrank_pm) = l2er_X_X;
corr_X_X = corr(X_wSMe___(:),X_wSM0___(:));
corr_X_X_r_(1+nrank_pm) = corr_X_X;
[~,tmp_nw_SM_] = max(X_wSM0___,[],1); tmp_nw_SM_ = tmp_nw_SM_-1;
tmp_X_SM___ = zeros(1,n_S,n_M);
for nM=0:n_M-1; for nS=0:n_S-1;
tmp_nw = tmp_nw_SM_(1,1+nS,1+nM);
tmp_X = X_wSMe___(1+tmp_nw,1+nS,1+nM);
tmp_X_SM___(1,1+nS,1+nM) = tmp_X;
end;end;%for nM=0:n_M-1; for nS=0:n_S-1;
tmp_X_SM___ = repmat(tmp_X_SM___,[n_w_max,1,1]);
tmp_X_SM___ = X_wSMe___ > tmp_X_SM___;
tmp_p__ = sum(tmp_X_SM___,1)/n_w_max;
tmp_p = mean(tmp_p__,'all');
prct_X_X_r_(1+nrank_pm) = tmp_p;
end;%for nrank_pm=0:n_rank_pm-1;
%%%%%%%%;
save(fname_mat ...
     ,'tolerance_master' ...
     ,'rank_pm_' ...
     ,'rank_pm','nrank_pm' ...
     ,'delta_r_max','str_delta_r_max' ...
     ,'l2er_X_X_r_','corr_X_X_r_','prct_X_X_r_' ...
     );
close_fname_tmp(fname_pre);
%%%%%%%%;
end;%if ~flag_skip;
%%%%%%%%;

if  exist(fname_mat,'file');
load(fname_mat);
%%%%%%%%;
fname_fig = sprintf('%s_jpg/pm_fig_X_2d_xcor_d0_%s_FIGL_',dir_pm,str_delta_r_max);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;figsml;
tmp_ij_ = 1:numel(rank_pm_)-1;
rank_pm_use_ = rank_pm_(tmp_ij_);
%%%%;
subplot(1,1,1);
hold on;
%plot(rank_pm_use_,l2er_X_X_r_(tmp_ij_),'ko-','MarkerFaceColor','k');
plot(rank_pm_use_,0+corr_X_X_r_(tmp_ij_),'ko-','MarkerFaceColor','r');
plot(rank_pm_use_,1-prct_X_X_r_(tmp_ij_),'ko-','MarkerFaceColor','c');
hold off;
xlim([0,rank_pm_use_(end)+1]); set(gca,'XTick',0:1:rank_pm_use_(end)); xlabel('rank $H$','Interpreter','latex');
ylim([0.5,1]); ylabel('value'); set(gca,'YTick',0.5:0.05:1.0); grid on;
ylabel('correlation'); grid on;
%legend({'frob','corr','prct'});
legend({'corr','prct'});
title('correlation');
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%if  exist(fname_mat,'file');

%%%%%%%%%%%%%%%%;
end;%for delta_r_max = [ 0.010 , 0.025 , 0.050 , 0.075 , 0.100 ];
%%%%%%%%%%%%%%%%;

fname_fig = sprintf('%s_jpg/pm_fig_X_2d_xcor_d0_dxxxx_FIGL_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;figsml;
linewidth_use = 2; markersize_use = 8; fontsize_use = 12;
tmp_ij_ = 1:numel(rank_pm_)-1;
rank_pm_use_ = rank_pm_(tmp_ij_);
delta_r_max_ = [ 0.010 , 0.025 , 0.050 , 0.075 , 0.100 ]; n_delta_r_max = numel(delta_r_max_);
%tmp_symbol_ = {'o-','^-','s-','p-','h-'};
tmp_symbol_ = {'o-','o-','o-','o-','o-'};
tmp_intens_ = 0.90.^[0:n_delta_r_max-1];
tmp_legend_ = ...
{ ...
   '\delta = 0.010 corr' ...
   '\delta = 0.025 corr' ...
   '\delta = 0.050 corr' ...
   '\delta = 0.075 corr' ...
   '\delta = 0.100 corr' ...
   '\delta = 0.010 prct' ...
   '\delta = 0.025 prct' ...
   '\delta = 0.050 prct' ...
   '\delta = 0.075 prct' ...
   '\delta = 0.100 prct' ...
};
tmp_d_ = cell(n_delta_r_max,1);
for ndelta_r_max=0:n_delta_r_max-1;
delta_r_max = delta_r_max_(1+ndelta_r_max);
str_delta_r_max = sprintf('d%.4d',round(1000*delta_r_max));
tmp_fname_mat = sprintf('%s_mat/pm_fig_X_2d_xcor_d0_%s_FIGL_.mat',dir_pm,str_delta_r_max);
if  exist(tmp_fname_mat,'file'); tmp_d_{1+ndelta_r_max} = load(tmp_fname_mat); end;
end;%for ndelta_r_max=0:n_delta_r_max-1;
%%%%;
subplot(1,1,1);
hold on;
for ndelta_r_max=0:n_delta_r_max-1;
tmp_symbol = tmp_symbol_{1+ndelta_r_max}; tmp_intens = tmp_intens_(1+ndelta_r_max);
plot(rank_pm_use_,0+tmp_d_{1+ndelta_r_max}.corr_X_X_r_(tmp_ij_),tmp_symbol,'MarkerEdgeColor','k','MarkerFaceColor',[1,0,0]*tmp_intens,'Color',[1,0,0]*tmp_intens,'LineWidth',linewidth_use,'MarkerSize',markersize_use);
end;%for ndelta_r_max=0:n_delta_r_max-1;
for ndelta_r_max=0:n_delta_r_max-1;
tmp_symbol = tmp_symbol_{1+ndelta_r_max}; tmp_intens = tmp_intens_(1+ndelta_r_max);
plot(rank_pm_use_,1-tmp_d_{1+ndelta_r_max}.prct_X_X_r_(tmp_ij_),tmp_symbol,'MarkerEdgeColor','k','MarkerFaceColor',[0,1,1]*tmp_intens,'Color',[0,1,1]*tmp_intens,'LineWidth',linewidth_use,'MarkerSize',markersize_use);
end;%for ndelta_r_max=0:n_delta_r_max-1;
hold off;
%xlim([0,rank_pm_use_(end)+1]); set(gca,'XTick',0:1:rank_pm_use_(end)); xlabel('rank $H$','Interpreter','latex');
xlim([0,rank_pm_use_(8)+1]); set(gca,'XTick',0:1:rank_pm_use_(8)); xlabel('rank $H$','Interpreter','latex');
ylim([0.5,1]); ylabel('value'); set(gca,'YTick',0.5:0.05:1.0); grid on;
ylabel('correlation'); grid on;
%legend({'frob','corr','prct'});
legend(tmp_legend_,'location','SouthEast');
set(gca,'FontSize',fontsize_use);
title('correlation');
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ( flag_exist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

disp('returning'); return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
