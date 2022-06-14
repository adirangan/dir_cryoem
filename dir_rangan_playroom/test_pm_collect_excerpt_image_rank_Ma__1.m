function ...
[ ...
 parameter ...
,str_fname_align_a_k_Y_mat_a_ ...
,image_X_value_Ma__ ...
,R_k_p_l2_Ma__ ...
,ij_X_value_nsort_from_nM_Ma__ ...
,ij_R_k_p_l2_nsort_from_nM_Ma__ ...
,ij_min_nsort_from_nM_Ma__ ...
,ij_max_nsort_from_nM_Ma__ ...
,X_value_rank_avg_from_nt_Ma__ ...
,R_k_p_l2_rank_avg_from_nt_Ma__ ...
,min_rank_avg_from_nt_Ma__ ...
,max_rank_avg_from_nt_Ma__ ...
,X_TM_ ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
,corr_a_k_Y_true ...
,image_X_value_true_ ...
] = ...
test_pm_collect_excerpt_image_rank_Ma__1( ...
 parameter ...
,dir_pm ...
);

str_thisfunction = 'test_pm_collect_excerpt_image_rank_Ma__1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if nargin<1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

dir_base = '/data/rangan/dir_cryoem';
parameter = struct('type','parameter');
parameter.flag_center_image = 1;
parameter.flag_recalc = 0;
n_octile = 8; na = 0;
for noctile=0:n_octile-1;
dir_pm = sprintf('%s/dir_ps1_x%d/dir_pm',dir_base,noctile);
[ ...
 parameter ...
,tmp_str_fname_align_a_k_Y_mat_a_ ...
,tmp_image_X_value_Ma__ ...
,tmp_R_k_p_l2_Ma__ ...
,tmp_ij_X_value_nsort_from_nM_Ma__ ...
,tmp_ij_R_k_p_l2_nsort_from_nM_Ma__ ...
,tmp_ij_min_nsort_from_nM_Ma__ ...
,tmp_ij_max_nsort_from_nM_Ma__ ...
,tmp_X_value_rank_avg_from_nt_Ma__ ...
,tmp_R_k_p_l2_rank_avg_from_nt_Ma__ ...
,tmp_min_rank_avg_from_nt_Ma__ ...
,tmp_max_rank_avg_from_nt_Ma__ ...
,tmp_X_TM_ ...
,tmp_euler_polar_a_true_ ...
,tmp_euler_azimu_b_true_ ...
,tmp_euler_gamma_z_true_ ...
,tmp_image_delta_x_true_ ...
,tmp_image_delta_y_true_ ...
,tmp_corr_a_k_Y_true ...
,tmp_image_X_value_true_ ...
] = ...
test_pm_collect_excerpt_image_rank_Ma__1( ...
 parameter ...
,dir_pm ...
);
%%%%;
tmp_f_p = @(sA) ~isempty(strfind(sA,'X_2d_xcor_d0_a1t0122p25r0'));
tmp_index_ = efind(cellfun(tmp_f_p,tmp_str_fname_align_a_k_Y_mat_a_));
n_a_local = numel(tmp_index_);
n_a_local_(1+noctile) = n_a_local;
image_X_value_Moa___(:,1+noctile,1 + [0:n_a_local-1]) = tmp_image_X_value_Ma__(:,1+tmp_index_);
R_k_p_l2_Moa___(:,1+noctile,1 + [0:n_a_local-1]) = tmp_R_k_p_l2_Ma__(:,1+tmp_index_);
image_X_value_true_Mo__(:,1+noctile) = tmp_image_X_value_true_;
euler_polar_a_true_Mo__(:,1+noctile) = tmp_euler_polar_a_true_;
euler_azimu_b_true_Mo__(:,1+noctile) = tmp_euler_azimu_b_true_;
euler_gamma_z_true_Mo__(:,1+noctile) = tmp_euler_gamma_z_true_;
na = na + n_a_local;
end;%for noctile=0:n_octile-1;
n_a = na;
%%%%;
image_X_value_true_Moa___ = repmat(image_X_value_true_Mo__,[1,1,n_a_local]);
euler_polar_a_true_Moa___ = repmat(euler_polar_a_true_Mo__,[1,1,n_a_local]);
euler_azimu_b_true_Moa___ = repmat(euler_azimu_b_true_Mo__,[1,1,n_a_local]);
euler_gamma_z_true_Moa___ = repmat(euler_gamma_z_true_Mo__,[1,1,n_a_local]);
%%%%;
n_M = size(tmp_X_TM_,1); n_M_sum = n_M*n_octile;
image_X_value_ = reshape(mean(image_X_value_Moa___,3),[n_M*n_octile,1]);
image_X_value_true_ = reshape(mean(image_X_value_true_Moa___,3),[n_M*n_octile,1]);
euler_polar_a_true_ = reshape(mean(euler_polar_a_true_Moa___,3),[n_M*n_octile,1]);
euler_azimu_b_true_ = reshape(mean(euler_azimu_b_true_Moa___,3),[n_M*n_octile,1]);
euler_gamma_z_true_ = reshape(mean(euler_gamma_z_true_Moa___,3),[n_M*n_octile,1]);
%%%%;
%[~,ij_image_X_value_nM_from_nsort_] = sort(image_X_value_,'ascend');
[~,ij_image_X_value_nM_from_nsort_] = sort(image_X_value_true_,'ascend');
[~,ij_image_X_value_nsort_from_nM_] = sort(ij_image_X_value_nM_from_nsort_,'ascend');
ij_from_noctile__ = cell(n_octile,1);
for noctile=0:n_octile-1;
ij_from_noctile__{1+noctile} = ij_image_X_value_nM_from_nsort_(1+[0:n_M-1] + floor((n_M_sum-n_M)*noctile/(n_octile-1)));
end;%for noctile=0:n_octile-1;
flag_disp=0;
if flag_disp;
figure(1);clf;figsml;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
hold on;
for noctile=0:n_octile-1;
ij_from_noctile_ = ij_from_noctile__{1+noctile};
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*noctile/n_octile)));
plot(noctile*n_M+[0:n_M-1],image_X_value_(ij_from_noctile_),'.','Color',c_80s__(1+nc_80s,:));
end;%for noctile=0:n_octile-1;
hold off;
xlim([0,n_M_sum]);
xlabel('nMsum','Interpreter','none');
end;%if flag_disp;
%%%%%%%%;
[ ...
 tmp_n_all ...
,tmp_azimu_b_all_ ...
,tmp_polar_a_all_ ...
,tmp_weight_all_ ...
,tmp_k_c_0_all_ ...
,tmp_k_c_1_all_ ...
,tmp_k_c_2_all_ ...
,tmp_n_polar_a ...
,tmp_polar_a_ ...
,tmp_n_azimu_b_ ...
] = ...
sample_shell_5( ...
 1.0 ...
,1/(2*pi/2) ...
) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
nf=0;

flag_replot = 1;
fname_fig = sprintf('%s/dir_ps1_x0to7_combine/dir_pm_jpg/test_pm_collect_excerpt_image_rank_true_1',dir_base);
if ( flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file') );
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
colormap(1-colormap('gray'));
p_row = 2; p_col = 5; np=0;
markersize_sml = 4;
n_h_sml = 32;
n_h_big = 16;
%%%%;
xlim_ = [0,0.7];
subplot(p_row,p_col,1+np); np=np+1;
plot(image_X_value_true_,image_X_value_,'k.','MarkerSize',markersize_sml);
xlim(xlim_); xlabel('X_value_true','Interpreter','none'); set(gca,'Xtick',linspace(0,0.7,8));
ylim(xlim_); ylabel('X_value_ampm','Interpreter','none'); set(gca,'Ytick',linspace(0,0.7,8));
axis square; grid on;
%%%%;
subplot(p_row,p_col,1+np); np=np+1;
tmp_x_ = image_X_value_true_;
tmp_y_ = image_X_value_;
imagesc(log2(1+hist2d_0(tmp_x_,tmp_y_,n_h_sml,n_h_sml,xlim_,xlim_)));
set(gca,'ydir','normal');
axis image; axisnotick;
xlabel('X_value_true','Interpreter','none');
ylabel('X_value_ampm','Interpreter','none');
%%%%;
for noctile=0:n_octile-1;
subplot(p_row,p_col,1+np); np=np+1;
ij_from_noctile_ = ij_from_noctile__{1+noctile};
tmp_a_ = euler_polar_a_true_(ij_from_noctile_);
tmp_b_ = euler_azimu_b_true_(ij_from_noctile_);
%imagesc(log2(1+hist2d_0(tmp_b_,tmp_a_,2*n_h_big,1*n_h_big,[0,2*pi],[0,1*pi])));
imagesc(hist2d_0(tmp_b_,tmp_a_,2*n_h_big,1*n_h_big,[0,2*pi],[0,1*pi]));
%{
hist2d_polar_a_azimu_b_0( ...
 tmp_polar_a_all_ ...
,tmp_azimu_b_all_ ...
,tmp_weight_all_ ...
,tmp_a_ ...
,tmp_b_ ...
,[0,0.55*1024] ...
,colormap_beach ...
,1 ...
,0 ...
);
 %}
xlabel('azimu_b'); ylabel('polar_a');
axis image; axisnotick;
title(sprintf('1+noctile %d/%d',1+noctile,n_octile));
end;%for noctile=0:n_octile-1;
%%%%;
sgtitle(sprintf('%s',fname_fig),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if ( flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file') );
if ( exist(sprintf('%s.jpg',fname_fig),'file') );
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file') );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
disp('returning'); return;
end;%if nargin<1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); dir_pm=[]; end; na=na+1;
if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_recalc'); parameter.flag_recalc=0; end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
if ~isfield(parameter,'flag_center_image'); parameter.flag_center_image=0; end;
flag_recalc = parameter.flag_recalc;
flag_verbose = parameter.flag_verbose;
flag_center_image = parameter.flag_center_image;
if isempty(dir_pm); dir_pm = pwd; end;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); string_root = 'home'; end;
if (strcmp(platform,'eval1')); string_root = 'home'; end;
if (strcmp(platform,'rusty')); string_root = 'mnt/home'; end;
%%%%%%%%;
dir_base = sprintf('/%s/rangan/dir_cryoem',string_root);

dir_pm_mat = sprintf('%s_mat',dir_pm);

%%%%%%%%;
tmp_fname = sprintf('%s/M_k_p__.mat',dir_pm_mat);
if ~exist(tmp_fname,'file'); disp(sprintf(' %% Warning, %s not found in %s',tmp_fname,str_thisfunction)); return; end;
tmp_ = load(tmp_fname,'n_M');
n_M = tmp_.n_M;
clear tmp_;
%%%%%%%%;
tmp_fname = sprintf('%s/X_TM_.mat',dir_pm_mat);
if ~exist(tmp_fname,'file'); disp(sprintf(' %% Warning, %s not found in %s',tmp_fname,str_thisfunction)); return; end;
tmp_ = load(tmp_fname);
X_TM_ = tmp_.X_TM_;
[~,ij_X_TM_nM_from_nsort_] = sort(X_TM_,'ascend'); [~,ij_X_TM_nsort_from_nM_] = sort(ij_X_TM_nM_from_nsort_,'ascend');
clear tmp_;
%%%%%%%%;
flag_N_vs_M = flag_center_image;
if flag_N_vs_M==0; tmp_str = 'M'; end;%if flag_N_vs_M==0; 
if flag_N_vs_M==1; tmp_str = 'N'; end;%if flag_N_vs_M==1;
tmp_fname = sprintf('%s/a_k_Y_reco_from_%s__.mat',dir_pm_mat,tmp_str);
if ~exist(tmp_fname,'file'); disp(sprintf(' %% Warning, %s not found in %s',tmp_fname,str_thisfunction)); return; end;
tmp_ = load(tmp_fname);
euler_polar_a_true_ = tmp_.euler_polar_a_Mi__(:,end);
euler_azimu_b_true_ = tmp_.euler_azimu_b_Mi__(:,end);
euler_gamma_z_true_ = tmp_.euler_gamma_z_Mi__(:,end);
image_delta_x_true_ = tmp_.image_delta_x_acc_Mi__(:,end) + tmp_.image_delta_x_upd_Mi__(:,end);
image_delta_y_true_ = tmp_.image_delta_y_acc_Mi__(:,end) + tmp_.image_delta_y_upd_Mi__(:,end);
corr_a_k_Y_true = tmp_.corr_a_k_Y_i_(end);
image_X_value_true_ = tmp_.image_X_value_Mi__(:,end-1);
clear tmp_;
%%%%%%%%;

%%%%%%%%;
% find estimated image_rank. ;
%%%%%%%%;
fname_pre = sprintf('%s/test_pm_collect_compare_image_rank_ampm_',dir_pm_mat);
[fname_skip,fname_mat] = open_fname_tmp(fname_pre);
if ( flag_recalc | ~fname_skip );
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
tmp_X_2d_xcor_d0_a1t_align_ = ls(sprintf('%s/X_2d_xcor_d0_a1t*_align_a_k_Y_.mat',dir_pm_mat));
tmp_index_start_ = strfind(tmp_X_2d_xcor_d0_a1t_align_,sprintf('/%s',string_root))-1;
tmp_index_final_ = strfind(tmp_X_2d_xcor_d0_a1t_align_,'.mat')+4-1;
n_X_2d_xcor_d0_a1t_align = min(numel(tmp_index_start_),numel(tmp_index_final_));
%%%%;
image_X_value_Ma__ = zeros(n_M,n_X_2d_xcor_d0_a1t_align);
R_k_p_l2_Ma__ = zeros(n_M,n_X_2d_xcor_d0_a1t_align);
ij_X_value_nsort_from_nM_Ma__ = zeros(n_M,n_X_2d_xcor_d0_a1t_align);
ij_R_k_p_l2_nsort_from_nM_Ma__ = zeros(n_M,n_X_2d_xcor_d0_a1t_align);
ij_min_nsort_from_nM_Ma__ = zeros(n_M,n_X_2d_xcor_d0_a1t_align);
ij_max_nsort_from_nM_Ma__ = zeros(n_M,n_X_2d_xcor_d0_a1t_align);
X_value_rank_avg_from_nt_Ma__ = zeros(n_M,n_X_2d_xcor_d0_a1t_align);
R_k_p_l2_rank_avg_from_nt_Ma__ = zeros(n_M,n_X_2d_xcor_d0_a1t_align);
min_rank_avg_from_nt_Ma__ = zeros(n_M,n_X_2d_xcor_d0_a1t_align);
max_rank_avg_from_nt_Ma__ = zeros(n_M,n_X_2d_xcor_d0_a1t_align);
%%%%%%%%;
for na=0:n_X_2d_xcor_d0_a1t_align-1;
tmp_fname_align_a_k_Y_mat = tmp_X_2d_xcor_d0_a1t_align_(1+tmp_index_start_(1+na)+0:1+tmp_index_final_(1+na)-1);
str_fname_align_a_k_Y_mat_a_{1+na} = tmp_fname_align_a_k_Y_mat;
if (tmp_fname_align_a_k_Y_mat(1)==sprintf('\n')); tmp_fname_align_a_k_Y_mat = tmp_fname_align_a_k_Y_mat(2:end); end;
tmp_fname_a_k_Y_mat = tmp_fname_align_a_k_Y_mat;
tmp_ij = strfind(tmp_fname_a_k_Y_mat,'_align_a_k_Y_.mat');
tmp_fname_a_k_Y_mat = tmp_fname_a_k_Y_mat([1:tmp_ij,tmp_ij+7:end]);
tmp_fname_compare_image_rank_pre = tmp_fname_align_a_k_Y_mat;
tmp_ij = strfind(tmp_fname_compare_image_rank_pre,'_align_a_k_Y_.mat');
tmp_fname_compare_image_rank_pre = tmp_fname_compare_image_rank_pre([1:tmp_ij-1]);
tmp_fname_compare_image_rank_pre = sprintf('%s_compare_image_rank',tmp_fname_compare_image_rank_pre);
tmp_fname_compare_image_rank_mat = sprintf('%s.mat',tmp_fname_compare_image_rank_pre);
tmp_image_X_value_ = zeros(n_M,1);
tmp_R_k_p_l2_ = zeros(n_M,1);
tmp_ij_X_value_nsort_from_nM_ = zeros(n_M,1);
tmp_ij_R_k_p_l2_nsort_from_nM_ = zeros(n_M,1);
tmp_ij_min_nsort_from_nM_ = zeros(n_M,1);
tmp_ij_max_nsort_from_nM_ = zeros(n_M,1);
tmp_X_value_rank_avg_from_nt_ = zeros(n_M,1);
tmp_R_k_p_l2_rank_avg_from_nt_ = zeros(n_M,1);
tmp_min_rank_avg_from_nt_ = zeros(n_M,1);
tmp_max_rank_avg_from_nt_ = zeros(n_M,1);
%%%%%%%%;
if (  exist(tmp_fname_align_a_k_Y_mat,'file') &  exist(tmp_fname_a_k_Y_mat,'file') &  exist(tmp_fname_compare_image_rank_mat,'file') );
tmp_ = load(tmp_fname_compare_image_rank_mat);
assert(tmp_.n_M==n_M);
tmp_error = fnorm(tmp_.X_TM_-X_TM_)/max(1e-12,fnorm(X_TM_));
if (tmp_error>1e-2); disp(sprintf(' %% Warning, X_TM_ vs tmp_.X_TM_: %0.16f in %s',tmp_error,tmp_fname_align_a_k_Y_mat)); end;
tmp_image_X_value_ = tmp_.tmp_image_X_value_;
tmp_R_k_p_l2_ = tmp_.tmp_R_k_p_l2_;
[~,tmp_ij_X_value_nM_from_nsort_] = sort(tmp_image_X_value_,'ascend'); [~,tmp_ij_X_value_nsort_from_nM_] = sort(tmp_ij_X_value_nM_from_nsort_,'ascend');
[~,tmp_ij_R_k_p_l2_nM_from_nsort_] = sort(tmp_R_k_p_l2_,'descend'); [~,tmp_ij_R_k_p_l2_nsort_from_nM_] = sort(tmp_ij_R_k_p_l2_nM_from_nsort_,'ascend');
tmp_ij_min_ = min(tmp_ij_X_value_nsort_from_nM_,tmp_ij_R_k_p_l2_nsort_from_nM_);
[~,tmp_ij_min_nM_from_nsort_] = sort(tmp_ij_min_,'ascend'); [~,tmp_ij_min_nsort_from_nM_] = sort(tmp_ij_min_nM_from_nsort_,'ascend');
tmp_ij_max_ = max(tmp_ij_X_value_nsort_from_nM_,tmp_ij_R_k_p_l2_nsort_from_nM_);
[~,tmp_ij_max_nM_from_nsort_] = sort(tmp_ij_max_,'ascend'); [~,tmp_ij_max_nsort_from_nM_] = sort(tmp_ij_max_nM_from_nsort_,'ascend');
tmp_X_value_rank_avg_from_nt_ = nsort_from_nv_threshold_0(1,tmp_ij_X_value_nsort_from_nM_,ij_X_TM_nsort_from_nM_);
tmp_R_k_p_l2_rank_avg_from_nt_ = nsort_from_nv_threshold_0(1,tmp_ij_R_k_p_l2_nsort_from_nM_,ij_X_TM_nsort_from_nM_);
tmp_min_rank_avg_from_nt_ = nsort_from_nv_threshold_0(1,tmp_ij_min_nsort_from_nM_,ij_X_TM_nsort_from_nM_);
tmp_max_rank_avg_from_nt_ = nsort_from_nv_threshold_0(1,tmp_ij_max_nsort_from_nM_,ij_X_TM_nsort_from_nM_);
end;%if (  exist(tmp_fname_align_a_k_Y_mat,'file') &  exist(tmp_fname_a_k_Y_mat,'file') &  exist(tmp_fname_compare_image_rank_mat,'file') );
%%%%%%%%;
image_X_value_Ma__(:,1+na) = tmp_image_X_value_;
R_k_p_l2_Ma__(:,1+na) = tmp_R_k_p_l2_;
ij_X_value_nsort_from_nM_Ma__(:,1+na) = tmp_ij_X_value_nsort_from_nM_;
ij_R_k_p_l2_nsort_from_nM_Ma__(:,1+na) = tmp_ij_R_k_p_l2_nsort_from_nM_;
ij_min_nsort_from_nM_Ma__(:,1+na) = tmp_ij_min_nsort_from_nM_;
ij_max_nsort_from_nM_Ma__(:,1+na) = tmp_ij_max_nsort_from_nM_;
X_value_rank_avg_from_nt_Ma__(:,1+na) = tmp_X_value_rank_avg_from_nt_;
R_k_p_l2_rank_avg_from_nt_Ma__(:,1+na) = tmp_R_k_p_l2_rank_avg_from_nt_;
min_rank_avg_from_nt_Ma__(:,1+na) = tmp_min_rank_avg_from_nt_;
max_rank_avg_from_nt_Ma__(:,1+na) = tmp_max_rank_avg_from_nt_;
%%%%%%%%;
end;%for na=0:n_X_2d_xcor_d0_a1t_align-1;
%%%%;
save( ...
fname_mat ...
,'str_fname_align_a_k_Y_mat_a_' ...
,'image_X_value_Ma__' ...
,'R_k_p_l2_Ma__' ...
,'ij_X_value_nsort_from_nM_Ma__' ...
,'ij_R_k_p_l2_nsort_from_nM_Ma__' ...
,'ij_min_nsort_from_nM_Ma__' ...
,'ij_max_nsort_from_nM_Ma__' ...
,'X_value_rank_avg_from_nt_Ma__' ...
,'R_k_p_l2_rank_avg_from_nt_Ma__' ...
,'min_rank_avg_from_nt_Ma__' ...
,'max_rank_avg_from_nt_Ma__' ...
,'X_TM_' ...
,'euler_polar_a_true_' ...
,'euler_azimu_b_true_' ...
,'euler_gamma_z_true_' ...
,'image_delta_x_true_' ...
,'image_delta_y_true_' ...
,'corr_a_k_Y_true' ...
,'image_X_value_true_' ...
);
%%%%%%%%;
close_fname_tmp(fname_pre);
end;%if ( flag_recalc | ~fname_skip );
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
