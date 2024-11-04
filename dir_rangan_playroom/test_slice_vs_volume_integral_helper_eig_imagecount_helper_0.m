%%%%%%%%;
% Note, may have to reset tolerance_pm = 1e-3;
%%%%%%%%;
tolerance_pm = 1e-3;
%%%%%%%%;
% recalculate idealized principal-modes. ;
%%%%%%%%;
n_3 = 3;
if ~exist('KAPPA','var'); KAPPA=[]; end;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___=[]; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__=[]; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__=[]; end;
if ~exist('l_max_uk_','var'); l_max_uk_=[]; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_=[]; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__=[]; end;
if ~exist('V_lmm___','var'); V_lmm___=[]; end;
if ~exist('L_lm__','var'); L_lm__=[]; end;
if ~exist('d0W_betazeta_mlma____','var'); d0W_betazeta_mlma____=[]; end;
if ~exist('d1W_betazeta_mlma____','var'); d1W_betazeta_mlma____=[]; end;
if ~exist('d2W_betazeta_mlma____','var'); d2W_betazeta_mlma____=[]; end;

%%%%%%%%;
% If necessary, calculate the idealized principal-modes for unit CTF. ;
%%%%%%%%;
if ~exist('X_2d_x1_d0_kk__','var');
[X_2d_x1_d0_kk__,X_2d_x1_d0_weight_r_] = principled_marching_cost_matrix_6(n_k_p_r,k_p_r_,weight_2d_k_p_r_,l_max_,[],[],a_k_Y_quad_yk_);
end;%if ~exist('X_2d_x1_d0_kk__','var');
%%%%%%%%;
% Now determine principal-modes. ;
%%%%%%%%;
if ~exist('tolerance_pm','var'); tolerance_pm = 1e-3; end;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
X_kk__ = X_2d_x1_d0_kk__;
[tmp_UX__,tmp_SX__,tmp_VX__] = svds(X_kk__,n_UX_rank); tmp_SX_ = diag(tmp_SX__);
pm_n_UX_rank = max(find(tmp_SX_/max(tmp_SX_)> tolerance_pm));
UX_kn__ = zeros(n_k_p_r,n_UX_rank); SX_k_ = zeros(n_UX_rank,1);
UX_kn__(:,:) = tmp_UX__(:,1+[0:n_UX_rank-1]);
SX_k_(:) = tmp_SX_(1+[0:n_UX_rank-1]);
nlt = -log10(tolerance_pm);
str_tolerance_pm = sprintf('nlt%.2dpm%d',10*nlt,pm_n_UX_rank);
if (flag_verbose>0); disp(sprintf(' %% tolerance_pm %0.6f: pm_n_UX_rank %d/%d --> %s',tolerance_pm,pm_n_UX_rank,n_UX_rank,str_tolerance_pm)); end;
%%%%%%%%;

%%%%%%%%;
% First collect/collate data. ;
%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
str_dir_mat = sprintf('%s_mat',dir_ssnll);
if ~exist(str_dir_mat,'dir'); disp(sprintf(' %% mkdir %s',str_dir_mat)); mkdir(str_dir_mat); end;
str_dir_jpg = sprintf('%s_jpg',dir_ssnll);
if ~exist(str_dir_jpg,'dir'); disp(sprintf(' %% mkdir %s',str_dir_jpg)); mkdir(str_dir_jpg); end;
str_dir_jpg_stripped = sprintf('%s_jpg_stripped',dir_ssnll);
if ~exist(str_dir_jpg_stripped,'dir'); disp(sprintf(' %% mkdir %s',str_dir_jpg_stripped)); mkdir(str_dir_jpg_stripped); end;
lsigma_ = [NaN,-3:0.5:+3]; n_lsigma = numel(lsigma_);
lanczos_n_iteration_max = 32;
ddssnll_mid_q2d_si__ = zeros(n_lsigma,lanczos_n_iteration_max);
ddssnll_dif_q2d_si__ = zeros(n_lsigma,lanczos_n_iteration_max);
ddssnll_lsq_q2d_si__ = zeros(n_lsigma,lanczos_n_iteration_max);
nfound = 0;
for nlsigma=0:n_lsigma-1;
lsigma = lsigma_(1+nlsigma);
if ~isfinite(lsigma); str_infix = sprintf('p_empirical'); end;
if lsigma< -1e-12; str_infix = sprintf('lsigma_n%.3d',fix(100*abs(lsigma))); end;
if lsigma>=-1e-12; str_infix = sprintf('lsigma_p%.3d',fix(100*abs(lsigma))); end;
str_fname_nopath_prefix = sprintf('eig_i1_from_synth_%s_%s',str_tolerance_pm,str_infix);
str_dir_sub_mat = sprintf('%s/dir_%s',str_dir_mat,str_fname_nopath_prefix);
disp(sprintf(' %% %.2d/%.2d: %s',nlsigma,n_lsigma,str_fname_nopath_prefix));
for index_lambda=0:lanczos_n_iteration_max-1;
str_fname_nopath_sub_prefix = sprintf('%s_i%.3dl%.3d',str_fname_nopath_prefix,lanczos_n_iteration_max-1,index_lambda);
tmp_fname_sub_mat = sprintf('%s/%s.mat',str_dir_sub_mat,str_fname_nopath_sub_prefix);
if ~exist(tmp_fname_sub_mat,'file');
disp(sprintf(' %% Warning, %s not found',tmp_fname_sub_mat));
end;%if ~exist(tmp_fname_sub_mat,'file');
if  exist(tmp_fname_sub_mat,'file');
tmp_ = load(tmp_fname_sub_mat,'tmp_ddssnll_mid_q2d','tmp_ddssnll_dif_q2d','tmp_ddssnll_lsq_q2d');
disp(sprintf(' %% %% %s: mid %+16.4f dif %+16.4f lsq %+16.4f',str_fname_nopath_sub_prefix,tmp_.tmp_ddssnll_mid_q2d,tmp_.tmp_ddssnll_dif_q2d,tmp_.tmp_ddssnll_lsq_q2d));
ddssnll_mid_q2d_si__(1+nlsigma,1+index_lambda) = tmp_.tmp_ddssnll_mid_q2d;
ddssnll_dif_q2d_si__(1+nlsigma,1+index_lambda) = tmp_.tmp_ddssnll_dif_q2d;
ddssnll_lsq_q2d_si__(1+nlsigma,1+index_lambda) = tmp_.tmp_ddssnll_lsq_q2d;
nfound  = nfound+1;
clear tmp_;
end;%if  exist(tmp_fname_sub_mat,'file');
end;%for index_lambda=0:lanczos_n_iteration_max-1;
disp(sprintf(' %% '));
end;% for nlsigma=0:n_lsigma-1;
if (flag_verbose>0); disp(sprintf(' %% found %d/%d',nfound,n_lsigma*lanczos_n_iteration_max)); end;
%%%%%%%%;

%%%%%%%%;
% Now visualize. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed; fig81s;
p_row = 1; p_col = 3; np=0;
subplot(p_row,p_col,1+np);np=np+1; 
imagesc(real(ddssnll_mid_q2d_si__)); axisnotick; xlabel('index_lambda','Interpreter','none'); ylabel('index_sigma','Interpreter','none'); title('mid');
subplot(p_row,p_col,1+np);np=np+1; 
imagesc(real(ddssnll_dif_q2d_si__)); axisnotick; xlabel('index_lambda','Interpreter','none'); ylabel('index_sigma','Interpreter','none'); title('dif');
subplot(p_row,p_col,1+np);np=np+1; 
imagesc(real(ddssnll_lsq_q2d_si__)); axisnotick; xlabel('index_lambda','Interpreter','none'); ylabel('index_sigma','Interpreter','none'); title('lsq');
set(gcf,'Position',1+[0,0,1024*2,512]);
fname_fig_pre = sprintf('%s/eig_i1_from_synth_%s_lsigma_xxxx',str_dir_jpg,str_tolerance_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/eig_i1_from_synth_%s_lsigma_xxxx',str_dir_jpg_stripped,str_tolerance_pm);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
print('-depsc',fname_fig_stripped_eps);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
%%%%%%%%;

figure(1+nf);nf=nf+1;clf;figmed; figbeach;
p_row = 1; p_col = 3; np=0;
llim_ = [-6,+6];
subplot(p_row,p_col,1+np);np=np+1; 
imagesc(log(max(1e-12,real(ddssnll_mid_q2d_si__))),llim_); axisnotick; 
xlabel('index_lambda','Interpreter','none'); ylabel('index_sigma','Interpreter','none'); title('log(mid)');
tmp_c_ = colorbar; set(tmp_c_,'Ticks',[-6,-3,0,+3,+6],'TickLength',0);
subplot(p_row,p_col,1+np);np=np+1; 
imagesc(log(max(1e-12,real(ddssnll_dif_q2d_si__))),llim_); axisnotick; 
xlabel('index_lambda','Interpreter','none'); ylabel('index_sigma','Interpreter','none'); title('log(dif)');
tmp_c_ = colorbar; set(tmp_c_,'Ticks',[-6,-3,0,+3,+6],'TickLength',0);
subplot(p_row,p_col,1+np);np=np+1; 
imagesc(log(max(1e-12,real(ddssnll_lsq_q2d_si__))),llim_); axisnotick; 
xlabel('index_lambda','Interpreter','none'); ylabel('index_sigma','Interpreter','none'); title('log(lsq)');
tmp_c_ = colorbar; set(tmp_c_,'Ticks',[-6,-3,0,+3,+6],'TickLength',0);
set(gcf,'Position',1+[0,0,1024*2,512]);
fname_fig_pre = sprintf('%s/eig_i1_from_synth_%s_lsigma_xxxx_log',str_dir_jpg,str_tolerance_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/eig_i1_from_synth_%s_lsigma_xxxx_log',str_dir_jpg_stripped,str_tolerance_pm);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
print('-depsc',fname_fig_stripped_eps);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
%%%%%%%%;

%%%%%%%%;
% reconstruct tmp_a_x_u_reco_. ;
%%%%%%%%;
a_k_Y_reco_yk_ = zeros(n_lm_sum,1);
for pm_nk_p_r=0:pm_n_k_p_r-1;
for nk_p_r=0:n_k_p_r-1;
if (flag_verbose>1); disp(sprintf(' %% adding pm_nk_p_r %d/%d nk_p_r %d/%d',pm_nk_p_r,pm_n_k_p_r,nk_p_r,n_k_p_r)); end;
tmp_l_max = l_max_(1+nk_p_r);
pm_tmp_l_max = pm_l_max_(1+pm_nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
pm_tmp_n_lm = (pm_tmp_l_max+1).^2;
pm_tmp_index_ = pm_n_lm_csum_(1+pm_nk_p_r) + (0:pm_tmp_n_lm-1);
a_k_Y_reco_yk_(1+tmp_index_) = a_k_Y_reco_yk_(1+tmp_index_) + UX_kn__(1+nk_p_r,1+pm_nk_p_r)/max(1e-12,X_weight_r_(1+nk_p_r))*pm_a_k_Y_quad_yk__(1:tmp_n_lm,1+pm_nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for pm_nk_p_r=0:pm_n_k_p_r-1;
%%%%;
tmp_a_k_p_reco_ = zeros(n_k_all,1);
tmp_t = tic;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
 tmp_a_k_p_reco_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_spharm_to_k_p_4( ...
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
,a_k_Y_reco_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% tmp_a_k_p_reco_ time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
tmp_a_x_u_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,tmp_a_k_p_reco_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: tmp_a_x_u_reco_ time %0.2fs',tmp_t));
%%%%%%%%;

lambda_cut = 3.5;
%%%%%%%%;
% Run through each single example and try to visualize the alignment perturbation. ;
if strcmp(dir_nopath_data_star,'rib80s');
val_zoom_use = 1.50;
prct_use = 98.00;
npick_ = { ...
  ,[1,4] ...
  ,[2,1] ...
  ,[3,0] ...
  ,[4,2] ...
  ,[5,1] ...
  ,[6,2] ...
  ,[7,0] ...
  ,[10,0] ...
};
end;%if strcmp(dir_nopath_data_star,'rib80s');
if strcmp(dir_nopath_data_star,'trpv1');
val_zoom_use = 2.00;
prct_use = 98.75;
npick_ = { ...
  [0,4] ...
  ,[1,4] ...
  ,[2,2] ...
  ,[3,1] ...
  ,[3,2] ...
  ,[4,0] ...
  ,[4,1] ...
  ,[6,3] ...
  ,[7,0] ...
  ,[7,2] ...
  ,[11,0] ...
};
end;%if strcmp(dir_nopath_data_star,'TRPV1');
if strcmp(dir_nopath_data_star,'ISWINCP');
val_zoom_use = 1.85;
prct_use = 99.25;
npick_ = { ...
  ,[1,4] ...
  ,[2,1] ...
  ,[3,0] ...
  ,[4,0] ...
  ,[4,1] ...
  ,[4,2] ...
  ,[5,0] ...
  ,[6,2] ...
  ,[7,0] ...
  ,[7,1] ...
  ,[12,0] ...
};
end;%if strcmp(dir_nopath_data_star,'ISWINCP');
if strcmp(dir_nopath_data_star,'MlaFEDB');
val_zoom_use = 1.75;
prct_use = 99.00;
npick_ = { ...
  ,[1,4] ...
  ,[2,1] ...
  ,[3,0] ...
  ,[4,0] ...
  ,[5,0] ...
  ,[6,0] ...
  ,[7,2] ...
  ,[12,0] ...
};
end;%if strcmp(dir_nopath_data_star,'MlaFEDB');
n_pick = numel(npick_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for npick=0:n_pick-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
npick = npick_{1+npick};
nlsigma = npick(1+0); index_lambda = npick(1+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%for nlsigma=0:n_lsigma-1; for index_lambda=0:lanczos_n_iteration_max-1;
%nlsigma= 4; index_lambda=2; %<-- also try: nlsigma=5; index_lambda=1; or nlsigma=4; index_lambda=0;
%nlsigma=11; index_lambda=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
lambda_mid = real(ddssnll_mid_q2d_si__(1+nlsigma,1+index_lambda));
lambda_dif = real(ddssnll_dif_q2d_si__(1+nlsigma,1+index_lambda));
lambda_lsq = real(ddssnll_lsq_q2d_si__(1+nlsigma,1+index_lambda));
lambda_mid_min = min(real(ddssnll_mid_q2d_si__(1+nlsigma,:)));
lambda_dif_min = min(real(ddssnll_dif_q2d_si__(1+nlsigma,:)));
lambda_lsq_min = min(real(ddssnll_lsq_q2d_si__(1+nlsigma,:)));
flag_ismin = (lambda_dif==lambda_dif_min) | (lambda_lsq==lambda_lsq_min);
disp(sprintf(' %% nlsigma %d index_lambda %d lambda_mid +%0.6f lambda_dif +%0.6f lambda_lsq +%0.6f',nlsigma,index_lambda,lambda_mid,lambda_dif,lambda_lsq));
lsigma = lsigma_(1+nlsigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%if (flag_ismin | (lambda_dif<=lambda_cut & lambda_lsq<=lambda_cut));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ~isfinite(lsigma); str_infix = sprintf('p_empirical'); end;
if lsigma< -1e-12; str_infix = sprintf('lsigma_n%.3d',fix(100*abs(lsigma))); end;
if lsigma>=-1e-12; str_infix = sprintf('lsigma_p%.3d',fix(100*abs(lsigma))); end;
str_fname_nopath_prefix = sprintf('eig_i1_from_synth_%s_%s',str_tolerance_pm,str_infix);
str_dir_sub_mat = sprintf('%s/dir_%s',str_dir_mat,str_fname_nopath_prefix);
str_dir_sub_jpg = sprintf('%s/dir_%s',str_dir_jpg,str_fname_nopath_prefix);
disp(sprintf(' %% %.2d/%.2d: %s',nlsigma,n_lsigma,str_fname_nopath_prefix));
str_fname_nopath_sub_prefix = sprintf('%s_i%.3dl%.3d',str_fname_nopath_prefix,lanczos_n_iteration_max-1,index_lambda);
tmp_fname_sub_mat = sprintf('%s/%s.mat',str_dir_sub_mat,str_fname_nopath_sub_prefix);
if ~exist(tmp_fname_sub_mat,'file');
disp(sprintf(' %% Warning, %s not found',tmp_fname_sub_mat));
end;%if ~exist(tmp_fname_sub_mat,'file');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if  exist(tmp_fname_sub_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_ = load(tmp_fname_sub_mat);
%%%%;
from_pm_UX_kn__ = UX_kn__;
pm_n_k_p_r = tmp_.pm_n_k_p_r;
pm_k_p_r_ = tmp_.pm_k_p_r_;
pm_k_p_r_max = tmp_.pm_k_p_r_max;
pm_l_max_ = tmp_.pm_l_max_;
pm_weight_3d_k_p_r_ = tmp_.pm_weight_3d_k_p_r_;
n_M_use = tmp_.n_M_use;
weight_imagecount_M_use_ = tmp_.weight_imagecount_M_use_;
euler_polar_a_M_use_ = tmp_.euler_polar_a_M_use_;
euler_azimu_b_M_use_ = tmp_.euler_azimu_b_M_use_;
euler_gamma_z_M_use_ = tmp_.euler_gamma_z_M_use_;
n_M_imp = tmp_.n_M_imp;
weight_imagecount_M_imp_ = tmp_.weight_imagecount_M_imp_;
scaling_volumetric = tmp_.scaling_volumetric;
pm_weight_3d_riesz_k_p_r_ = tmp_.pm_weight_3d_riesz_k_p_r_;
pm_weight_3d_riesz_weight_imagecount_ykabc_ = tmp_.pm_weight_3d_riesz_weight_imagecount_ykabc_;
pm_v_tilde_eig_ykabc_ = tmp_.pm_v_tilde_eig_ykabc_;
pm_v_eig_ykabc_ = tmp_.pm_v_eig_ykabc_;
%%%%;

%%%%%%%%;
% visualize pm_v_eig_dvol_yk_;
%%%%%%%%;
[pm_v_tilde_eig_dvol_yk_,pm_v_tilde_eig_polar_a_M_use_,pm_v_tilde_eig_azimu_b_M_use_,pm_v_tilde_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(pm_n_k_p_r,pm_l_max_,n_M_imp,pm_v_tilde_eig_ykabc_);
[pm_v_eig_dvol_yk_,pm_v_eig_polar_a_M_use_,pm_v_eig_azimu_b_M_use_,pm_v_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(pm_n_k_p_r,pm_l_max_,n_M_imp,pm_v_eig_ykabc_);
pm_v_k_Y_use_yk_ = pm_v_eig_dvol_yk_;
pm_v_k_Y_use_yk__ = local_yk__from_yk_(pm_n_k_p_r,pm_l_max_,pm_v_k_Y_use_yk_);
v_k_Y_reco_yk_ = zeros(n_lm_sum,1);
pm_n_UX_rank = pm_n_k_p_r;
for nUX_rank=0:pm_n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_l_max = l_max_(1+nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
v_k_Y_reco_yk_(1+tmp_index_) = v_k_Y_reco_yk_(1+tmp_index_) + from_pm_UX_kn__(1+nk_p_r,1+nUX_rank)/max(1e-12,X_weight_r_(1+nk_p_r))*pm_v_k_Y_use_yk__(1:tmp_n_lm,1+nUX_rank);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:pm_n_UX_rank-1;
%%%%;
tmp_t = tic;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
v_k_p_reco_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_spharm_to_k_p_4( ...
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
,v_k_Y_reco_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% v_k_p_reco_ time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
v_x_u_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,v_k_p_reco_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: v_x_u_reco_ time %0.2fs',tmp_t));

flag_disp=0;
if flag_disp;
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGF',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
vmax = max(abs(v_x_u_reco_),[],'all');
amax = max(abs(tmp_a_x_u_reco_),[],'all');
dvol_ = 3e-1*[-1:+1]; p_col = numel(dvol_);
prct_ = [98.5,99.0]; p_row = numel(prct_);
np=0;
for prow=0:p_row-1;
prct = prct_(1+prow);
pcol=0;
for dvol=0:numel(dvol_)-1;
subplot(p_row,p_col,1+np);np=np+1;cla;
dvol = dvol_(1+dvol);
tmp_tmp_a_x_u_reco_ = tmp_a_x_u_reco_/max(1e-12,amax) + dvol*v_x_u_reco_/max(1e-12,vmax);
isosurface_f_x_u_1(struct('percent_threshold_',[prct]),tmp_tmp_a_x_u_reco_);
title(sprintf('prct %5.2f dvol %0.2f',prct,dvol));
pcol=pcol+1;
end;%for dvol=0:numel(dvol_)-1;
end;%for prow=0:p_row-1;
set(gcf,'Position',1+[0,0,1024*2,2*(768-128)]);
str_sgtitle = sprintf('%s index_lambda %d lambda %+0.2f[%+0.2f] = exp(%+0.2f)[exp(%+0.2f)]',str_fname_nopath_prefix,index_lambda,lambda_dif,lambda_lsq,log(abs(lambda_dif)),log(abs(lambda_lsq)));
sgtitle(str_sgtitle,'Interpreter','none');
%sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
%print('-depsc',fname_fig_eps);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%;
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGG',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
prct_ = [98.5,99.0]; p_row = numel(prct_);
p_col = 2; np=0;
for prow=0:p_row-1;
prct = prct_(1+prow);
subplot(p_row,p_col,1+np);np=np+1;cla;
isosurface_f_x_u_1(struct('percent_threshold_',[100-prct,prct],'c_use__',colormap_pm),v_x_u_reco_);
title(sprintf('prct [%5.2f,%5.2f]:',[100-prct,prct]));
subplot(p_row,p_col,1+np);np=np+1;cla;
isosurface_f_x_u_2(struct('percent_threshold_',prct),tmp_a_x_u_reco_,v_x_u_reco_);
title(sprintf('prct [%5.2f]:',[prct]));
end;%for prow=0:p_row-1;
set(gcf,'Position',1+[0,0,1024*1.5,1024*1.5]);
str_sgtitle = sprintf('%s index_lambda %d lambda %+0.2f[%+0.2f] = exp(%+0.2f)[exp(%+0.2f)]',str_fname_nopath_prefix,index_lambda,lambda_dif,lambda_lsq,log(abs(lambda_dif)),log(abs(lambda_lsq)));
sgtitle(str_sgtitle,'Interpreter','none');
%sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
%print('-depsc',fname_fig_eps);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%;
end;%if flag_disp;
 
flag_disp=0;
if flag_disp;
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGH',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
figure(1+nf);nf=nf+1;clf;figmed;
tmp_dtau_euler_fnorm = sum(bsxfun(@times,[tmp_.tmp_dtau_euler_polar_a_M_use_,tmp_.tmp_dtau_euler_azimu_b_M_use_,tmp_.tmp_dtau_euler_gamma_z_M_use_].^2,tmp_.weight_imagecount_M_use_),'all');
tmp_dtau_euler_nrm_polar_a_M_ = tmp_.tmp_dtau_euler_polar_a_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_azimu_b_M_ = tmp_.tmp_dtau_euler_azimu_b_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_gamma_z_M_ = tmp_.tmp_dtau_euler_gamma_z_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
factor_amplify = 2.0;
figbig;
subplot(1,1,1);
hold on;
sphere_compass__0( ...
struct('type','sphere_compass','compass_r_base',1.0/(2*pi)/2*sqrt(0.5),'flag_2d_vs_3d',1,'flag_normalize',1) ...
,tmp_.n_M_use ...
,tmp_.euler_polar_a_M_use_ ...
,tmp_.euler_azimu_b_M_use_ ...
,tmp_dtau_euler_nrm_polar_a_M_ ...
,tmp_dtau_euler_nrm_azimu_b_M_ ...
,tmp_dtau_euler_nrm_gamma_z_M_ ...
);
hold off;
xlim([0,2*pi]); xlabel('azimu_b','Interpreter','none');
ylim([0,1*pi]); ylabel('polar_a','Interpreter','none');
axisnotick;
title('dtau_euler_','Interpreter','none');
str_sgtitle = sprintf('%s index_lambda %d lambda %+0.2f[%+0.2f] = exp(%+0.2f)[exp(%+0.2f)]',str_fname_nopath_prefix,index_lambda,lambda_dif,lambda_lsq,log(abs(lambda_dif)),log(abs(lambda_lsq)));
sgtitle(str_sgtitle,'Interpreter','none');
%sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
%print('-depsc',fname_fig_eps);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%;
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGI',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
figure(1+nf);nf=nf+1;clf;figmed;
tmp_dtau_euler_fnorm = sum(bsxfun(@times,[tmp_.tmp_dtau_euler_polar_a_M_use_,tmp_.tmp_dtau_euler_azimu_b_M_use_,tmp_.tmp_dtau_euler_gamma_z_M_use_].^2,tmp_.weight_imagecount_M_use_),'all');
tmp_dtau_euler_nrm_polar_a_M_ = tmp_.tmp_dtau_euler_polar_a_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_azimu_b_M_ = tmp_.tmp_dtau_euler_azimu_b_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_gamma_z_M_ = tmp_.tmp_dtau_euler_gamma_z_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
factor_amplify = 2.0;
%%%%;
figbig;
for np=0:1;
subplot(1,2,1+np);
hold on;
plot_sphere_grid_0(struct('flag_solid',1));
sphere_compass__0( ...
struct('type','sphere_compass','compass_r_base',1.0/(2*pi)/4,'flag_2d_vs_3d',0,'flag_normalize',1) ...
,tmp_.n_M_use ...
,tmp_.euler_polar_a_M_use_ ...
,tmp_.euler_azimu_b_M_use_ ...
,tmp_dtau_euler_nrm_polar_a_M_ ...
,tmp_dtau_euler_nrm_azimu_b_M_ ...
,tmp_dtau_euler_nrm_gamma_z_M_ ...
);
hold off;
xlabel('x0'); ylabel('x1'); zlabel('x2'); axis equal; axis vis3d; axisnotick3d;
if np==0; view(0,+90); title('dtau_euler_ view(0,+90)','Interpreter','none'); end;
if np==1; view(0,-90); title('dtau_euler_ view(0,-90)','Interpreter','none'); end;
end;%for np=0:1;
%%%%;
sgtitle(str_sgtitle,'Interpreter','none');
%sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
%print('-depsc',fname_fig_eps);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%;
end;%if flag_disp;

%%%%%%%%;
% estimate magnitude of perturbation. ;
%%%%%%%%;
[~,~,tmp_v_std] = spharm_normalize_2(pm_n_k_p_r,pm_k_p_r_,pm_weight_3d_k_p_r_,pm_l_max_,pm_v_eig_dvol_yk_);
[~,~,tmp_a_std] = spharm_normalize_2(pm_n_k_p_r,pm_k_p_r_,pm_weight_3d_k_p_r_,pm_l_max_,pm_a_k_Y_quad_yk_);
if (flag_verbose>0); disp(sprintf(' %% tmp_v_std/tmp_a_std %0.6f',tmp_v_std/tmp_a_std)); end;
[pm_a_k_Y_quad_lr] = sqrt(local_f_bar_dot_g_(pm_n_k_p_r,pm_weight_3d_riesz_k_p_r_,pm_l_max_,0,pm_a_k_Y_quad_yk_,pm_a_k_Y_quad_yk_));
[pm_v_tilde_eig_dvol_lr] = sqrt(local_weightless_f_bar_dot_g_(pm_n_k_p_r,pm_l_max_,0,pm_v_tilde_eig_dvol_yk_,pm_v_tilde_eig_dvol_yk_));
[pm_v_eig_dvol_lr] = sqrt(local_f_bar_dot_g_(pm_n_k_p_r,pm_weight_3d_riesz_k_p_r_,pm_l_max_,0,pm_v_eig_dvol_yk_,pm_v_eig_dvol_yk_));
if (flag_verbose>0); disp(sprintf(' %% pm_v_eig_dvol_lr/pm_a_k_Y_quad_lr %0.6f',pm_v_eig_dvol_lr/pm_a_k_Y_quad_lr)); end;
% So we can imagine multiplying pm_v_eig_dvol_yk_ by (tmp_a_std/tmp_v_std) to produce something with equivalent norm to pm_a_k_Y_quad_yk_. ;
% This then multiplies the value of H by (tmp_a_std/tmp_v_std)^2 --> lambda_mid*(tmp_a_std/tmp_v_std)^2 (in image-driven units). ;
% Note that sum(weight_imagecount_M_use_) --> 1.0, ;
% so our likelihood takes the form of:
% ssnll(a+delta*dvol) - ssnll(a) = N_image * (0.5 * H_scaled * delta^2) ; %<-- here a and dvol have comparable l2-norm. ;
% nll(a+delta*dvol) - nll(a) = sigma^-2 * N_image * (0.5 * H_scaled * delta^2) ; %<-- delta could be thought of as dimensionless. ;
%%%%%%%%;
tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_/max(1e-12,fnorm(tmp_a_x_u_reco_));
v_x_u_reco_fnrm_ = v_x_u_reco_/max(1e-12,fnorm(v_x_u_reco_));
tmp_H_scale = min(lambda_dif,lambda_lsq)*(tmp_a_std/tmp_v_std)^2;
tmp_N_image = 1024;
tmp_inv_sigma_squared = inv(RR_bar/n_x_u_pack^2);
tmp_log20 = 3.0;
tmp_dvol_mag = real(sqrt(tmp_log20 / (tmp_inv_sigma_squared * 0.5 * tmp_H_scale) / tmp_N_image ));
if tmp_dvol_mag>=0.125;
tmp_N_image = ceil(real(tmp_log20 / (tmp_inv_sigma_squared * 0.5 * tmp_H_scale) / 0.125^2 ));
tmp_N_image = max(tmp_N_image,1000*ceil(tmp_N_image/1000));
tmp_dvol_mag = real(sqrt(tmp_log20 / (tmp_inv_sigma_squared * 0.5 * tmp_H_scale) / tmp_N_image ));
end;%if tmp_dvol_mag>=0.125;
if (flag_verbose>0); disp(sprintf(' %% tmp_dvol_mag %+0.6f tmp_N_image %d',tmp_dvol_mag,tmp_N_image)); end;
%%%%%%%%;
 
flag_disp=0;
if flag_disp;
%%%%%%%%;
% Make frames for perturbation movie. ;
%%%%%%%%;
dvol_ = tmp_dvol_mag*[-1:0.5:+1]; n_dvol = numel(dvol_);
prct_ = [prct_use]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
vlim_g_ = prctile(real(v_x_u_reco_(:)),[0.25,99.75]);
val_zoom = val_zoom_use;
for ndvol=0:n_dvol-1;
dvol = dvol_(1+ndvol);
figure(1+nf);nf=nf+1;clf;figsml;
fontsize_use = 16;
tmp_tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_fnrm_ + dvol*v_x_u_reco_fnrm_;
isosurface_f_x_u_2( ...
 struct('vval_',[vval],'vlim_g_',vlim_g_) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),tmp_tmp_a_x_u_reco_fnrm_) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),v_x_u_reco_) ...
);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('prct %5.2f dvol %+.04f',prct,dvol),'Interpreter','latex');
set(gcf,'Position',1+[0,0,1024*2,2*(768-128)]);
set(gca,'FontSize',fontsize_use);
str_sgtitle = sprintf('%s index_lambda %d lambda_dif %0.4f tmp_dvol_mag %0.6f tmp_N_image %d',str_fname_nopath_prefix,index_lambda,lambda_dif,tmp_dvol_mag,tmp_N_image);
sgtitle(str_sgtitle,'Interpreter','none');
str_subplot = sprintf('p%.4d_dvol%d',100*prct,ndvol);
fname_fig_pre = sprintf('%s/%s_FIGK_%s',str_dir_sub_jpg,str_fname_nopath_sub_prefix,str_subplot);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGK_%s_stripped',str_dir_jpg_stripped,str_fname_nopath_sub_prefix,str_subplot);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
end;%for ndvol=0:n_dvol-1;
%%%%%%%%;
end;%if flag_disp;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
%%%%;
subplot(1,2,[1]);
flim_ = 4.0*[0,2.0/(4*pi)];
flag_2d_vs_3d=0;
lsigma_dist = lsigma_dist_(1+nlsigma);
imagesc_polar_a_azimu_b_0( ...
 euler_polar_a_M_use_ ... 
,euler_azimu_b_M_use_ ... 
,factor_imagecount_Ms__(:,1+nlsigma) ...
,flim_ ... 
,colormap_beach ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
title('$\mu(\tau)$','Interpreter','latex');
axisnotick3d; axis equal; axis vis3d;
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(1,2,[2]);
prct_ = [prct_use]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
val_zoom = val_zoom_use;
isosurface_f_x_u_1( ...
 struct('vval_',[vval]) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),tmp_a_x_u_reco_fnrm_) ...
);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('$\\hat{F}$'),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
%%%%;
set(gcf,'Position',1+[0,0,1024*1.0,512]);
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGL',str_dir_sub_jpg,str_fname_nopath_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGL_stripped',str_dir_jpg_stripped,str_fname_nopath_prefix);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
disp(sprintf(' %% writing %s',fname_fig_stripped_pre));
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
%%%%%%%%;  
end;%if flag_disp;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
%%%%;
tmp_dtau_euler_fnorm = sum(bsxfun(@times,[tmp_.tmp_dtau_euler_polar_a_M_use_,tmp_.tmp_dtau_euler_azimu_b_M_use_,tmp_.tmp_dtau_euler_gamma_z_M_use_].^2,tmp_.weight_imagecount_M_use_),'all');
tmp_dtau_euler_nrm_polar_a_M_ = tmp_.tmp_dtau_euler_polar_a_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_azimu_b_M_ = tmp_.tmp_dtau_euler_azimu_b_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_gamma_z_M_ = tmp_.tmp_dtau_euler_gamma_z_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
factor_amplify = 2.0;
%%%%;
for np=0;%for np=0:1;
subplot(1,2,[1]);
hold on;
plot_sphere_grid_0(struct('flag_solid',1));
sphere_compass__0( ...
struct('type','sphere_compass','compass_r_base',1.0/(2*pi)/4,'flag_2d_vs_3d',0,'flag_normalize',1) ...
,tmp_.n_M_use ...
,tmp_.euler_polar_a_M_use_ ...
,tmp_.euler_azimu_b_M_use_ ...
,tmp_dtau_euler_nrm_polar_a_M_ ...
,tmp_dtau_euler_nrm_azimu_b_M_ ...
,tmp_dtau_euler_nrm_gamma_z_M_ ...
);
hold off;
xlabel(''); ylabel(''); zlabel(''); axis vis3d; axisnotick3d; axis([-1,+1,-1,+1,-1,+1]);
if np==0; view(0,+90); title('$\Delta\tau$: azimuthal $+\pi/2$','Interpreter','latex'); end;
if np==1; view(0,-90); title('$\Delta\tau$: azimuthal $-\pi/2$','Interpreter','latex'); end;
set(gca,'FontSize',fontsize_use);
end;%for np=0:1;
%%%%;
subplot(1,2,[2]);
prct = 98.5;
isosurface_f_x_u_1(struct('percent_threshold_',[100-prct,prct],'c_use__',colormap_pm),v_x_u_reco_);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('$\\Delta \\hat{F}$'),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
%%%%;
set(gcf,'Position',1+[0,0,1024*1.0,512]);
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGM',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGM_stripped',str_dir_jpg_stripped,str_fname_nopath_sub_prefix);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
disp(sprintf(' %% writing %s',fname_fig_stripped_pre));
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
%%%%%%%%;  
end;%if flag_disp;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
%%%%;
tmp_dtau_euler_fnorm = sum(bsxfun(@times,[tmp_.tmp_dtau_euler_polar_a_M_use_,tmp_.tmp_dtau_euler_azimu_b_M_use_,tmp_.tmp_dtau_euler_gamma_z_M_use_].^2,tmp_.weight_imagecount_M_use_),'all');
tmp_dtau_euler_nrm_polar_a_M_ = tmp_.tmp_dtau_euler_polar_a_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_azimu_b_M_ = tmp_.tmp_dtau_euler_azimu_b_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_gamma_z_M_ = tmp_.tmp_dtau_euler_gamma_z_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
factor_amplify = 2.0;
%%%%;
subplot(1,3,[1,2]);
hold on;
sphere_compass__0( ...
struct('type','sphere_compass','compass_r_base',1.0/(2*pi)/2*sqrt(0.5),'flag_2d_vs_3d',1,'flag_normalize',1) ...
,tmp_.n_M_use ...
,tmp_.euler_polar_a_M_use_ ...
,tmp_.euler_azimu_b_M_use_ ...
,tmp_dtau_euler_nrm_polar_a_M_ ...
,tmp_dtau_euler_nrm_azimu_b_M_ ...
,tmp_dtau_euler_nrm_gamma_z_M_ ...
);
hold off;
xlim([0,2*pi]); xlabel('azimuthal','Interpreter','latex');
set(gca,'XTick',0:pi/4:2*pi,'XTickLabel',[]);
ylim([0,1*pi]); ylabel('polar','Interpreter','latex');
set(gca,'YTick',0:pi/4:1*pi,'YTickLabel',[]);
grid on;
title('$\Delta\tau$: equatorial','Interpreter','latex');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(1,3,[3]);
prct = 98.5;
isosurface_f_x_u_1(struct('percent_threshold_',[100-prct,prct],'c_use__',colormap_pm),v_x_u_reco_);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('$\\Delta \\hat{F}$'),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
%%%%;
set(gcf,'Position',1+[0,0,1024*1.75,512]);
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGN',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGN_stripped',str_dir_jpg_stripped,str_fname_nopath_sub_prefix);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
disp(sprintf(' %% writing %s',fname_fig_stripped_pre));
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;  
close(gcf);
end;%if flag_disp;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
%%%%;
dvol_ = tmp_dvol_mag*[-1:2.0:+1]; n_dvol = numel(dvol_);
prct_ = [prct_use]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
vlim_g_ = prctile(real(v_x_u_reco_(:)),[0.25,99.75]);
val_zoom = val_zoom_use;
for ndvol=0:n_dvol-1;
%if ndvol==0; subplot(2,5,[ 2, 3, 7, 8]); end;
%if ndvol==1; subplot(2,5,[ 4, 5, 9,10]); end;
if ndvol==0; subplot(1,2,[1]); end;
if ndvol==1; subplot(1,2,[2]); end;
dvol = dvol_(1+ndvol);
tmp_tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_fnrm_ + dvol*v_x_u_reco_fnrm_;
isosurface_f_x_u_2( ...
 struct('vval_',[vval],'vlim_g_',vlim_g_) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),tmp_tmp_a_x_u_reco_fnrm_) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),v_x_u_reco_) ...
);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('$\\Delta \\hat{F}$: $%+0.3f$, $J$%: $%d$',dvol,tmp_N_image),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
end;%for ndvol=0:n_dvol-1;
%%%%;
set(gcf,'Position',1+[0,0,1024*1.0,512]);
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGO',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGO_stripped',str_dir_jpg_stripped,str_fname_nopath_sub_prefix);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
disp(sprintf(' %% writing %s',fname_fig_stripped_pre));
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;  
close(gcf);
end;%if flag_disp;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
interp_k = 1;
%%%%;
dvol_ = tmp_dvol_mag*[-1:2.0:+1]; n_dvol = numel(dvol_);
tmp_tmp_a_x_u_reco_mini_ = min(tmp_a_x_u_reco_fnrm_ + min(dvol_)*v_x_u_reco_fnrm_,tmp_a_x_u_reco_fnrm_ + max(dvol_)*v_x_u_reco_fnrm_);
prct_ = [prct_use]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
vlim_g_ = prctile(real(v_x_u_reco_(:)),[0.25,99.75]);
val_zoom = val_zoom_use;
for ndvol=0:n_dvol-1;
%if ndvol==0; subplot(2,5,[ 2, 3, 7, 8]); end;
%if ndvol==1; subplot(2,5,[ 4, 5, 9,10]); end;
if ndvol==0; subplot(1,2,[1]); c_use__ = flipud(colormap_pm); end;
if ndvol==1; subplot(1,2,[2]); c_use__ =        colormap_pm ; end;
dvol = dvol_(1+ndvol);
tmp_tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_fnrm_ + dvol*v_x_u_reco_fnrm_;
hold on;
isosurface_f_x_u_2( ...
 struct('vval_',[vval],'c_use__',0.95*[1,1,1],'flag_projection',0,'flag_collapse',0,'flag_boxgrid',1,'v_alpha',1.000) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom*1.00),interp3(reshape(tmp_tmp_a_x_u_reco_mini_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline')) ...
);
isosurface_f_x_u_2( ...
 struct('vval_',[vval],'c_use__',     c_use__,'flag_projection',0,'flag_collapse',0,'flag_boxgrid',1,'v_alpha',0.500) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom*0.99),interp3(reshape(tmp_tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline')) ...
);
hold off;
xlabel(''); ylabel(''); zlabel('');
title(sprintf('$\\Delta \\hat{F}$: $%+0.3f$, $J$: $%d$',dvol,tmp_N_image),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
end;%for ndvol=0:n_dvol-1;
%%%%;
set(gcf,'Position',1+[0,0,1024*1.0,512]);
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGP',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGP_stripped',str_dir_jpg_stripped,str_fname_nopath_sub_prefix);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
disp(sprintf(' %% writing %s',fname_fig_stripped_pre));
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;  
close(gcf);
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
interp_k = 1;
dvol_ = tmp_dvol_mag*[-1:2.0:+1]; n_dvol = numel(dvol_);
prct_ = [prct_use]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
vlim_g_ = [0,2.5e-3];
val_zoom = val_zoom_use;
dvol = max(dvol_);
tmp_tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_fnrm_ + dvol*v_x_u_reco_fnrm_;
figure(1+nf);nf=nf+1;clf;figsml;
[~,d_v_] = ...
isosurface_f_x_u_3( ...
 struct('vval_',[vval],'vlim_g_',vlim_g_,'flag_plot',0) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
);
close(gcf);
vlim_g_ = [0,prctile(d_v_,90)];
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
%%%%;
val_zoom = val_zoom_use;
for ndvol=0:n_dvol-1;
%if ndvol==0; subplot(2,5,[ 2, 3, 7, 8]); end;
%if ndvol==1; subplot(2,5,[ 4, 5, 9,10]); end;
if ndvol==0; subplot(1,2,[1]); end;
if ndvol==1; subplot(1,2,[2]); end;
dvol = dvol_(1+ndvol);
tmp_tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_fnrm_ + dvol*v_x_u_reco_fnrm_;
isosurface_f_x_u_3( ...
 struct('vval_',[vval],'vlim_g_',vlim_g_) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('$\\Delta \\hat{F}$: $%+0.3f$, $J$%: $%d$',dvol,tmp_N_image),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
end;%for ndvol=0:n_dvol-1;
%%%%;
set(gcf,'Position',1+[0,0,1024*1.0,512]);
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGQ',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGQ_stripped',str_dir_jpg_stripped,str_fname_nopath_sub_prefix);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
disp(sprintf(' %% writing %s',fname_fig_stripped_pre));
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;  
close(gcf);
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% Make frames for perturbation movie. ;
%%%%%%%%;
interp_k = 1;
dvol_ = tmp_dvol_mag*[-1:0.5:+1]; n_dvol = numel(dvol_);
prct_ = [prct_use]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
vlim_g_ = [0,2.5e-3];
val_zoom = val_zoom_use;
dvol = max(dvol_);
tmp_tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_fnrm_ + dvol*v_x_u_reco_fnrm_;
figure(1+nf);nf=nf+1;clf;figsml;
[~,d_v_] = ...
isosurface_f_x_u_3( ...
 struct('vval_',[vval],'vlim_g_',vlim_g_,'flag_plot',0) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
);
close(gcf);
vlim_g_ = [0,prctile(d_v_,90)];
%%%%%%%%;
for ndvol=0:n_dvol-1;
dvol = dvol_(1+ndvol);
figure(1+nf);nf=nf+1;clf;figsml;
fontsize_use = 16;
tmp_tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_fnrm_ + dvol*v_x_u_reco_fnrm_;
isosurface_f_x_u_3( ...
 struct('vval_',[vval],'vlim_g_',vlim_g_) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('prct %5.2f dvol %+.04f',prct,dvol),'Interpreter','latex');
set(gcf,'Position',1+[0,0,1024*2,2*(768-128)]);
set(gca,'FontSize',fontsize_use);
str_sgtitle = sprintf('%s index_lambda %d lambda_dif %0.4f tmp_dvol_mag %0.6f tmp_N_image %d',str_fname_nopath_prefix,index_lambda,lambda_dif,tmp_dvol_mag,tmp_N_image);
sgtitle(str_sgtitle,'Interpreter','none');
str_subplot = sprintf('p%.4d_dvol%d',100*prct,ndvol);
fname_fig_pre = sprintf('%s/%s_FIGK2_%s',str_dir_sub_jpg,str_fname_nopath_sub_prefix,str_subplot);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGK2_%s_stripped',str_dir_jpg_stripped,str_fname_nopath_sub_prefix,str_subplot);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
end;%for ndvol=0:n_dvol-1;
%%%%%%%%;
end;%if flag_disp;

%%%%;
clear tmp_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if  exist(tmp_fname_sub_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%end;%if (lambda_dif<=lambda_cut & lambda_lsq<=lambda_cut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%end;end;%for nlsigma=0:n_lsigma-1; for index_lambda=0:lanczos_n_iteration_max-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for npick=0:n_pick-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%{

if flag_disp;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
a_k_Y_reco_yk_ = zeros(n_lm_sum,1);
for pm_nk_p_r=0:pm_n_k_p_r-1;
for nk_p_r=0:n_k_p_r-1;
if (flag_verbose>1); disp(sprintf(' %% adding pm_nk_p_r %d/%d nk_p_r %d/%d',pm_nk_p_r,pm_n_k_p_r,nk_p_r,n_k_p_r)); end;
tmp_l_max = l_max_(1+nk_p_r);
pm_tmp_l_max = pm_l_max_(1+pm_nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
pm_tmp_n_lm = (pm_tmp_l_max+1).^2;
pm_tmp_index_ = pm_n_lm_csum_(1+pm_nk_p_r) + (0:pm_tmp_n_lm-1);
a_k_Y_reco_yk_(1+tmp_index_) = a_k_Y_reco_yk_(1+tmp_index_) + UX_kn__(1+nk_p_r,1+pm_nk_p_r)/max(1e-12,X_weight_r_(1+nk_p_r))*pm_a_k_Y_quad_yk__(1:tmp_n_lm,1+pm_nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for pm_nk_p_r=0:pm_n_k_p_r-1;
%%%%;
tmp_a_k_p_reco_ = zeros(n_k_all,1);
tmp_t = tic;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
 tmp_a_k_p_reco_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_spharm_to_k_p_4( ...
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
,a_k_Y_reco_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% tmp_a_k_p_reco_ time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
tmp_a_x_u_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,tmp_a_k_p_reco_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: tmp_a_x_u_reco_ time %0.2fs',tmp_t));
%%%%;
isosurface_f_x_u_1([],tmp_a_x_u_reco_); title(sprintf('pm_nk_p_r %d/%d',pm_nk_p_r,pm_n_k_p_r),'Interpreter','none');
drawnow();
%%%%%%%%;
end;%if flag_disp;

 %}
