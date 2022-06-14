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
test_pm_collect_excerpt_image_rank_Ma__0( ...
 parameter ...
,dir_pm ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if nargin<1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_center_image = 1;
parameter.flag_recalc = 0;
dir_pm = '/data/rangan/dir_cryoem/dir_ps1_x0/dir_pm';
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
test_pm_collect_excerpt_image_rank_Ma__0( ...
 parameter ...
,dir_pm ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
nf=0;
n_a = numel(str_fname_align_a_k_Y_mat_a_);
n_M = numel(X_TM_);
%tmp_f_n = @(sA) ~isempty(strfind(sA,'n16')) | ~isempty(strfind(sA,'n18')) | ~isempty(strfind(sA,'n20')) ;
%tmp_f_n = @(sA) ~isempty(strfind(sA,'n12')) | ~isempty(strfind(sA,'n14')) ;
%tmp_index_ = efind(cellfun(tmp_f_n,str_fname_align_a_k_Y_mat_a_));
tmp_f_p = @(sA)  ~isempty(strfind(sA,'p20')) | ~isempty(strfind(sA,'p25')) | ~isempty(strfind(sA,'p30')) ;
%tmp_f_p = @(sA) ~isempty(strfind(sA,'p10')) | ~isempty(strfind(sA,'p15')) ;
tmp_index_ = efind(cellfun(tmp_f_p,str_fname_align_a_k_Y_mat_a_));
disp(sprintf(' %% using %d/%d',numel(tmp_index_),n_a));
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
markersize_use = 16;
subplot(1,2,1); hold on;
for nl=0:numel(tmp_index_)-1;
na = tmp_index_(1+nl);
plot3(euler_polar_a_true_,euler_azimu_b_true_,X_value_rank_avg_from_nt_Ma__(:,1+na),'.','MarkerEdgeColor','k','MarkerSize',markersize_use);
end;%for nl=0:numel(tmp_index_)-1;
axis vis3d;
xlabel('polar_a','Interpreter','none');
ylabel('azimu_b','Interpreter','none');
zlabel('X_value_rank_avg','Interpreter','none');
subplot(1,2,2); hold on;
plot3(euler_polar_a_true_,euler_azimu_b_true_,mean(X_value_rank_avg_from_nt_Ma__(:,1+tmp_index_),2),'.','MarkerEdgeColor','k','MarkerSize',markersize_use);
axis vis3d;
xlabel('polar_a','Interpreter','none');
ylabel('azimu_b','Interpreter','none');
zlabel('X_value_rank_avg','Interpreter','none');
return;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 2; p_col = 3; np=0;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot(1:n_M,X_value_rank_avg_from_nt_Ma__(:,1+tmp_index_)/n_M,'k-','LineWidth',1);
plot(1:n_M,mean(X_value_rank_avg_from_nt_Ma__(:,1+tmp_index_)/n_M,2),'b-','LineWidth',3);
hold off;
xlim([1,n_M]);ylim([0.5,1.0]);
xlabel('X threshold'); ylabel('true rank avg');
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot(1:n_M,R_k_p_l2_rank_avg_from_nt_Ma__(:,1+tmp_index_)/n_M,'k-','LineWidth',1);
plot(1:n_M,mean(R_k_p_l2_rank_avg_from_nt_Ma__(:,1+tmp_index_)/n_M,2),'g-','LineWidth',3);
hold off;
xlim([1,n_M]);ylim([0.5,1.0]);
xlabel('R threshold'); ylabel('true rank avg');
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot(1:n_M,min_rank_avg_from_nt_Ma__(:,1+tmp_index_)/n_M,'k-','LineWidth',1);
plot(1:n_M,mean(min_rank_avg_from_nt_Ma__(:,1+tmp_index_)/n_M,2),'r-','LineWidth',3);
hold off;
xlim([1,n_M]);ylim([0.5,1.0]);
xlabel('min threshold'); ylabel('true rank avg');
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot(1:n_M,max_rank_avg_from_nt_Ma__(:,1+tmp_index_)/n_M,'k-','LineWidth',1);
plot(1:n_M,mean(max_rank_avg_from_nt_Ma__(:,1+tmp_index_)/n_M,2),'c-','LineWidth',3);
hold off;
xlim([1,n_M]);ylim([0.5,1.0]);
xlabel('max threshold'); ylabel('true rank avg');
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot(1:n_M,mean(X_value_rank_avg_from_nt_Ma__(:,1+tmp_index_)/n_M,2),'b-','LineWidth',3);
plot(1:n_M,mean(R_k_p_l2_rank_avg_from_nt_Ma__(:,1+tmp_index_)/n_M,2),'g-','LineWidth',3);
plot(1:n_M,mean(min_rank_avg_from_nt_Ma__(:,1+tmp_index_)/n_M,2),'r-','LineWidth',3);
plot(1:n_M,mean(max_rank_avg_from_nt_Ma__(:,1+tmp_index_)/n_M,2),'c-','LineWidth',3);
hold off;
legend({'X','R','min','max'},'Location','NorthWest');
xlim([1,n_M]);ylim([0.5,1.0]);
xlabel('threshold'); ylabel('true rank avg');

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
if ~exist(tmp_fname,'file'); disp(sprintf(' %% Warning, %s not found in test_pm_collect_excerpt_image_rank_Ma__0',tmp_fname)); return; end;
tmp_ = load(tmp_fname,'n_M');
n_M = tmp_.n_M;
clear tmp_;
%%%%%%%%;
tmp_fname = sprintf('%s/X_TM_.mat',dir_pm_mat);
if ~exist(tmp_fname,'file'); disp(sprintf(' %% Warning, %s not found in test_pm_collect_excerpt_image_rank_Ma__0',tmp_fname)); return; end;
tmp_ = load(tmp_fname);
X_TM_ = tmp_.X_TM_;
[~,ij_X_TM_nM_from_nsort_] = sort(X_TM_,'ascend'); [~,ij_X_TM_nsort_from_nM_] = sort(ij_X_TM_nM_from_nsort_,'ascend');
clear tmp_;
%%%%%%%%;
flag_N_vs_M = flag_center_image;
if flag_N_vs_M==0; tmp_str = 'M'; end;%if flag_N_vs_M==0; 
if flag_N_vs_M==1; tmp_str = 'N'; end;%if flag_N_vs_M==1;
tmp_fname = sprintf('%s/a_k_Y_reco_from_%s__.mat',dir_pm_mat,tmp_str);
if ~exist(tmp_fname,'file'); disp(sprintf(' %% Warning, %s not found in test_pm_collect_excerpt_image_rank_Ma__0',tmp_fname)); return; end;
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
