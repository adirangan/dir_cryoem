%{
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
%}

%%%%%%%%;
tolerance_pm = 0; pm_n_UX_rank = n_k_p_r;
n_UX_rank = n_k_p_r; %<-- full reconstruction. ;
X_kk__ = eye(n_UX_rank);
UX_kn__ = eye(n_UX_rank);
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
%lsigma_ = [NaN,-3:0.5:+3]; n_lsigma = numel(lsigma_);
lsigma_ = [NaN]; n_lsigma = numel(lsigma_);
n_ctf_select = 3;
%lanczos_n_iteration_max = 32;
ddssnll_mid_q2d_ssi___ = zeros(n_ctf_select,n_lsigma,lanczos_n_iteration_max);
ddssnll_dif_q2d_ssi___ = zeros(n_ctf_select,n_lsigma,lanczos_n_iteration_max);
ddssnll_lsq_q2d_ssi___ = zeros(n_ctf_select,n_lsigma,lanczos_n_iteration_max);
nfound = 0;
for nctf_select=0:n_ctf_select-1;
n_CTF_select = n_CTF_select_s_(1+nctf_select);
CTF_select_k_p_r_kC__ = CTF_select_k_p_r_kCs___(:,:,1+nctf_select);
CTF_select_k_p_wkC__ = CTF_select_k_p_wkCs___(:,:,1+nctf_select);
pm_n_CTF = n_CTF_select;
pm_index_nCTF_from_nM_use_ = reshape(ones(n_M_one,1)*[0:pm_n_CTF-1],[n_M_use,1]);
pm_CTF_k_p_r_kC__ = CTF_select_k_p_r_kC__;
pm_CTF_k_p_wkC__ = CTF_select_k_p_wkC__;
pm_M_use_k_p_wkM__ = pm_M_tri_k_p_wkM__.*pm_CTF_k_p_wkC__(:,1+pm_index_nCTF_from_nM_use_);
DefocusV_CTF_select_C_ = DefocusV_CTF_select_Cs__(:,1+nctf_select);
DefocusV_CTF_select_avg = mean(DefocusV_CTF_select_C_);
str_ctf_select = sprintf('DefV%d',round(DefocusV_CTF_select_avg));
for nlsigma=0:n_lsigma-1;
lsigma = lsigma_(1+nlsigma);
if ~isfinite(lsigma); str_infix = sprintf('p_empirical'); end;
if lsigma< -1e-12; str_infix = sprintf('lsigma_n%.3d',fix(100*abs(lsigma))); end;
if lsigma>=-1e-12; str_infix = sprintf('lsigma_p%.3d',fix(100*abs(lsigma))); end;
if flag_implicit_dtau==0; str_fname_nopath_prefix = sprintf('eig_from_synth_%s_%s_%s',str_tolerance_pm,str_infix,str_ctf_select); end;
if flag_implicit_dtau==1; str_fname_nopath_prefix = sprintf('eig_i1_from_synth_%s_%s_%s',str_tolerance_pm,str_infix,str_ctf_select); end;
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
ddssnll_mid_q2d_ssi___(1+nctf_select,1+nlsigma,1+index_lambda) = tmp_.tmp_ddssnll_mid_q2d;
ddssnll_dif_q2d_ssi___(1+nctf_select,1+nlsigma,1+index_lambda) = tmp_.tmp_ddssnll_dif_q2d;
ddssnll_lsq_q2d_ssi___(1+nctf_select,1+nlsigma,1+index_lambda) = tmp_.tmp_ddssnll_lsq_q2d;
nfound  = nfound+1;
clear tmp_;
end;%if  exist(tmp_fname_sub_mat,'file');
end;%for index_lambda=0:lanczos_n_iteration_max-1;
disp(sprintf(' %% '));
end;% for nlsigma=0:n_lsigma-1;
end;% for nctf_select=0:n_ctf_select-1;
if (flag_verbose>0); disp(sprintf(' %% found %d/%d',nfound,n_ctf_select*n_lsigma*lanczos_n_iteration_max)); end;
%%%%%%%%;

ddssnll_mid_q2d_si__ = reshape(ddssnll_mid_q2d_ssi___,[n_ctf_select*n_lsigma,lanczos_n_iteration_max]);
ddssnll_dif_q2d_si__ = reshape(ddssnll_dif_q2d_ssi___,[n_ctf_select*n_lsigma,lanczos_n_iteration_max]);
ddssnll_lsq_q2d_si__ = reshape(ddssnll_lsq_q2d_ssi___,[n_ctf_select*n_lsigma,lanczos_n_iteration_max]);

%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
subplot(1,1,1);
b_ = bar(1:7,transpose(ddssnll_lsq_q2d_si__(:,1:7)),0.85);
set(b_(1+0),'EdgeColor','k','FaceColor',0.35*[1,1,1]);
set(b_(1+1),'EdgeColor','k','FaceColor',0.65*[1,1,1]);
set(b_(1+2),'EdgeColor','k','FaceColor',0.85*[1,1,1]);
set(gca,'XTick',1:7); set(gca,'YTick',0:0.5:1.5,'YTickLabel',{'','',''});
set(gca,'TickLength',[0,0]); grid on;
xlabel('mode number (defocus 10k, 17k, 23k)'); ylabel('rayleigh-quotient');
fname_fig_pre = sprintf('%s/ctf_select_FIGC',str_dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/ctf_select_FIGC',str_dir_jpg_stripped);
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
% Now visualize. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed; fig81s;
p_row = 1; p_col = 3; np=0;
subplot(p_row,p_col,1+np);np=np+1; 
imagesc(real(ddssnll_mid_q2d_si__)); axisnotick; xlabel('index_lambda','Interpreter','none'); ylabel('index_select','Interpreter','none'); title('mid');
subplot(p_row,p_col,1+np);np=np+1; 
imagesc(real(ddssnll_dif_q2d_si__)); axisnotick; xlabel('index_lambda','Interpreter','none'); ylabel('index_select','Interpreter','none'); title('dif');
subplot(p_row,p_col,1+np);np=np+1; 
imagesc(real(ddssnll_lsq_q2d_si__)); axisnotick; xlabel('index_lambda','Interpreter','none'); ylabel('index_select','Interpreter','none'); title('lsq');
set(gcf,'Position',1+[0,0,1024*2,512]);
fname_fig_pre = sprintf('%s/eig_i1_from_synth_%s_lsigma_xxxx_DefV_xxxx',str_dir_jpg,str_tolerance_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/eig_i1_from_synth_%s_lsigma_xxxx_DefV_xxxx',str_dir_jpg_stripped,str_tolerance_pm);
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
xlabel('index_lambda','Interpreter','none'); ylabel('index_select','Interpreter','none'); title('log(mid)');
tmp_c_ = colorbar; set(tmp_c_,'Ticks',[-6,-3,0,+3,+6],'TickLength',0);
subplot(p_row,p_col,1+np);np=np+1; 
imagesc(log(max(1e-12,real(ddssnll_dif_q2d_si__))),llim_); axisnotick; 
xlabel('index_lambda','Interpreter','none'); ylabel('index_select','Interpreter','none'); title('log(dif)');
tmp_c_ = colorbar; set(tmp_c_,'Ticks',[-6,-3,0,+3,+6],'TickLength',0);
subplot(p_row,p_col,1+np);np=np+1; 
imagesc(log(max(1e-12,real(ddssnll_lsq_q2d_si__))),llim_); axisnotick; 
xlabel('index_lambda','Interpreter','none'); ylabel('index_select','Interpreter','none'); title('log(lsq)');
tmp_c_ = colorbar; set(tmp_c_,'Ticks',[-6,-3,0,+3,+6],'TickLength',0);
set(gcf,'Position',1+[0,0,1024*2,512]);
fname_fig_pre = sprintf('%s/eig_i1_from_synth_%s_lsigma_xxxx_DefV_xxxx_log',str_dir_jpg,str_tolerance_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/eig_i1_from_synth_%s_lsigma_xxxx_DefV_xxxx_log',str_dir_jpg_stripped,str_tolerance_pm);
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
% reconstruct tmp_a_x_u_reco_frompm_. ;
%%%%%%%%;
a_k_Y_reco_frompm_yk_ = zeros(n_lm_sum,1);
for pm_nk_p_r=0:pm_n_k_p_r-1;
for nk_p_r=0:n_k_p_r-1;
if (flag_verbose>1); disp(sprintf(' %% adding pm_nk_p_r %d/%d nk_p_r %d/%d',pm_nk_p_r,pm_n_k_p_r,nk_p_r,n_k_p_r)); end;
tmp_l_max = l_max_(1+nk_p_r);
pm_tmp_l_max = pm_l_max_(1+pm_nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
pm_tmp_n_lm = (pm_tmp_l_max+1).^2;
pm_tmp_index_ = pm_n_lm_csum_(1+pm_nk_p_r) + (0:pm_tmp_n_lm-1);
a_k_Y_reco_frompm_yk_(1+tmp_index_) = a_k_Y_reco_frompm_yk_(1+tmp_index_) + UX_kn__(1+nk_p_r,1+pm_nk_p_r)/max(1e-12,X_weight_r_(1+nk_p_r))*pm_a_k_Y_reco_empi_yk__(1:tmp_n_lm,1+pm_nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for pm_nk_p_r=0:pm_n_k_p_r-1;
%%%%;
tmp_a_k_p_reco_frompm_ = zeros(n_k_all,1);
tmp_t = tic;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
tmp_a_k_p_reco_frompm_ ...
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
,a_k_Y_reco_frompm_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% tmp_a_k_p_reco_frompm_ time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
tmp_a_x_u_reco_frompm_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,tmp_a_k_p_reco_frompm_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: tmp_a_x_u_reco_frompm_ time %0.2fs',tmp_t));
%%%%%%%%;

lambda_cut = 3.5;
%%%%%%%%;
% Run through each single example and try to visualize the alignment perturbation. ;
%%%%%%%%;
if strcmp(dir_nopath_data_star,'rib80s');
val_zoom_use = 1.50;
prct_use = 98.00;
vlim_g_max = 2.5e-3;
npick_ = { ...
};
end;%if strcmp(dir_nopath_data_star,'rib80s');
if strcmp(dir_nopath_data_star,'trpv1');
val_zoom_use = 2.00;
prct_use = 98.75;
vlim_g_max = 2.5e-3;
npick_ = { ...
,[0,0,0],[0,0,1],[0,0,2],[0,0,3],[0,0,4] ...
,[1,0,0],[1,0,1],[1,0,2],[1,0,3] ...
,[2,0,0],[2,0,1],[2,0,2],[2,0,3] ...
};
end;%if strcmp(dir_nopath_data_star,'trpv1');
if strcmp(dir_nopath_data_star,'ISWINCP');
val_zoom_use = 1.85;
prct_use = 99.25;
vlim_g_max = 2.5e-3;
npick_ = { ...
};
end;%if strcmp(dir_nopath_data_star,'ISWINCP');
if strcmp(dir_nopath_data_star,'MlaFEDB');
val_zoom_use = 1.75;
prct_use = 99.00;
vlim_g_max = 2.5e-3;
npick_ = { ...
};
end;%if strcmp(dir_nopath_data_star,'MlaFEDB');
n_pick = numel(npick_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for npick=0:n_pick-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_npick_ = npick_{1+npick};
nctf_select = tmp_npick_(1+0);
nlsigma = tmp_npick_(1+1);
index_lambda = tmp_npick_(1+2);
DefocusV_CTF_select_C_ = DefocusV_CTF_select_Cs__(:,1+nctf_select);
DefocusV_CTF_select_avg = mean(DefocusV_CTF_select_C_);
str_ctf_select = sprintf('DefV%d',round(DefocusV_CTF_select_avg));

%%%%%%%%;
lambda_mid = real(ddssnll_mid_q2d_ssi___(1+nctf_select,1+nlsigma,1+index_lambda));
lambda_dif = real(ddssnll_dif_q2d_ssi___(1+nctf_select,1+nlsigma,1+index_lambda));
lambda_lsq = real(ddssnll_lsq_q2d_ssi___(1+nctf_select,1+nlsigma,1+index_lambda));
lambda_mid_min = min(real(ddssnll_mid_q2d_ssi___(1+nctf_select,1+nlsigma,:)));
lambda_dif_min = min(real(ddssnll_dif_q2d_ssi___(1+nctf_select,1+nlsigma,:)));
lambda_lsq_min = min(real(ddssnll_lsq_q2d_ssi___(1+nctf_select,1+nlsigma,:)));
flag_ismin = (lambda_dif==lambda_dif_min) | (lambda_lsq==lambda_lsq_min);
disp(sprintf(' %% nctf_select %d nlsigma %d index_lambda %d lambda_mid +%0.6f lambda_dif +%0.6f lambda_lsq +%0.6f',nctf_select,nlsigma,index_lambda,lambda_mid,lambda_dif,lambda_lsq));
lsigma = lsigma_(1+nlsigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%if (flag_ismin | (lambda_dif<=lambda_cut & lambda_lsq<=lambda_cut));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ~isfinite(lsigma); str_infix = sprintf('p_empirical'); end;
if lsigma< -1e-12; str_infix = sprintf('lsigma_n%.3d',fix(100*abs(lsigma))); end;
if lsigma>=-1e-12; str_infix = sprintf('lsigma_p%.3d',fix(100*abs(lsigma))); end;
if flag_implicit_dtau==0; str_fname_nopath_prefix = sprintf('eig_from_synth_%s_%s_%s',str_tolerance_pm,str_infix,str_ctf_select); end;
if flag_implicit_dtau==1; str_fname_nopath_prefix = sprintf('eig_i1_from_synth_%s_%s_%s',str_tolerance_pm,str_infix,str_ctf_select); end;
str_dir_sub_mat = sprintf('%s/dir_%s',str_dir_mat,str_fname_nopath_prefix);
str_dir_sub_jpg = sprintf('%s/dir_%s',str_dir_jpg,str_fname_nopath_prefix);
disp(sprintf(' %% %.2d/%.2d: %s',nlsigma,n_lsigma,str_fname_nopath_prefix));
str_fname_nopath_sub_prefix = sprintf('%s_i%.3dl%.3d',str_fname_nopath_prefix,lanczos_n_iteration_max-1,index_lambda);
tmp_fname_sub_mat = sprintf('%s/%s.mat',str_dir_sub_mat,str_fname_nopath_sub_prefix);
if ~exist(tmp_fname_sub_mat,'file');
disp(sprintf(' %% Warning, %s not found',tmp_fname_sub_mat));
end;%if ~exist(tmp_fname_sub_mat,'file');

%pm_a_k_Y_use_yk__ = pm_a_k_Y_quad_yk__;
pm_a_k_Y_use_yk__ = pm_a_k_Y_reco_empi_yk__;
parameter_figure = parameter;
parameter_figure.flag_verbose = flag_verbose;
parameter_figure.flag_disp = flag_disp;
parameter_figure.flag_replot = flag_replot;
parameter_figure.dir_ssnll = dir_ssnll;
parameter_figure.str_fname_nopath_prefix = str_fname_nopath_prefix;
parameter_figure.str_fname_nopath_sub_prefix = str_fname_nopath_sub_prefix;
parameter_figure.tmp_fname_sub_mat = tmp_fname_sub_mat;
parameter_figure.RR_bar = RR_bar;
parameter_figure.index_lambda = index_lambda;
parameter_figure.lambda_mid = lambda_mid;
parameter_figure.lambda_dif = lambda_dif;
parameter_figure.lambda_lsq = lambda_lsq;
parameter_figure.val_zoom_use = val_zoom_use;
parameter_figure.prct_use = prct_use;
parameter_figure.vlim_g_max = vlim_g_max;
[ ...
 parameter_figure ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
test_slice_vs_volume_integral_helper_eig_helper_figure_5(...
 parameter_figure ...
,n_x_u_pack ...
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
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
,k_p_r_max ...
,X_weight_r_ ...
,pm_a_k_Y_use_yk__ ...
,UX_kn__...
,n_viewing_S_use ...
,viewing_polar_a_S_use_ ...
,viewing_azimu_b_S_use_ ...
,viewing_weight_S_use_ ...
,euler_polar_a_tavg_ ...
,euler_azimu_b_tavg_ ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%end;%if (lambda_dif<=lambda_cut & lambda_lsq<=lambda_cut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for npick=0:n_pick-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
