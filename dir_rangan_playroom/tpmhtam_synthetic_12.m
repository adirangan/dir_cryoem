function tpmhtam_synthetic_12(...
 syn_infix...
,syn_rseed...
,syn_n_M...
,syn_snr...
,syn_n_UX_rank...
,syn_n_iteration...
,syn_n_order...
,dir_trunk...
,n_x_u...
,diameter_x_c...
,n_k_p_r...
,k_p_r_...
,k_p_r_max...
,weight_3d_k_p_r_...
,viewing_k_eq_d...
,template_k_eq_d...
,delta_r_max...
,svd_eps...
,n_delta_v_requested...
,n_w_max_0in...
,CTF_uni_avg_k_p_...
,l_max_...
,n_molecule...
,molecule_density_...
,a_k_Y_quad__...
,UX__...
,X_weight_r_...
,a_UX_Y_quad_mavg__...
,f_rand_0in...
,flag_plot...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% test_principled_marching_helper_trpv1_alternating_minimization_synthetic_12. ;
% (similar to tpmhtam_synthetic_9, except enforces uniform n_w_). ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

date_diff_threshold = 1.0;

if (nargin<29) f_rand_0in = 0.05; end;
if (nargin<30) flag_plot = 1; end;

rng(0);

n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
molecule_density_ = molecule_density_/sum(molecule_density_);
verbose=1;

fname_mat = sprintf('%s_mat/tpmhtamsux_%s_M_k_p___n%.3ds%.4dr%.3d.mat',dir_trunk,syn_infix,syn_n_M,floor(1000*syn_snr),syn_n_UX_rank);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
% Generating templates from molecules. ;
%%%%%%%%;
tmp_t = tic();
%[~,n_w_,weight_2d_k_p_r_,weight_2d_k_all_,n_viewing_all,viewing_azimu_b_all_,viewing_polar_a_all_,n_viewing_polar_a,viewing_polar_a_,n_viewing_azimu_b_,template_k_c_0__,template_k_c_1__,template_k_c_2__] = get_template_0(verbose,n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,l_max_,[],viewing_k_eq_d,template_k_eq_d,n_w_uni_);
n_w_uni_ = n_w_max_0in*ones(n_k_p_r,1);
[~,n_w_,weight_2d_k_p_r_,weight_2d_k_all_,n_viewing_all,viewing_azimu_b_all_,viewing_polar_a_all_,n_viewing_polar_a,viewing_polar_a_,n_viewing_azimu_b_,template_k_c_0__,template_k_c_1__,template_k_c_2__] = get_template_0(verbose,n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,l_max_,[],viewing_k_eq_d,0,n_w_uni_);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum = cumsum([0;n_w_]);
S_k_p___ = cell(n_molecule,1);
for nmolecule=0:n_molecule-1;
[S_k_p___{1+nmolecule}] = get_template_0(verbose,n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,l_max_,a_k_Y_quad__(:,1+nmolecule),viewing_k_eq_d,0,n_w_uni_);
end;%for nmolecule=0:n_molecule-1;
n_S = n_viewing_all; %<-- number of templates stored in S_k_p__. ;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% S_k_p___: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Use templates S_k_p___ to estimate typical signal strength (per unit area in k_p). ;
%%%%%%%%;
tmp_t = tic();
template_area_element_ = weight_2d_k_all_*(2*pi)^2; %<-- now sum(template_area_element_)==pi*k_p_r_max^2;
S_k_p_l1_avg_ = zeros(n_molecule,1);
for nmolecule=0:n_molecule-1;
S_k_p_l1_avg_(1+nmolecule) = mean(transpose(template_area_element_)*abs(S_k_p___{1+nmolecule}.*CTF_uni_avg_k_p_))/(pi*k_p_r_max^2);
end;%for nmolecule=0:n_molecule-1;
S_k_p_l1_avg = sum(molecule_density_.*S_k_p_l1_avg_);
if (syn_snr<=0); syn_sigma = 0; end;
if (syn_snr> 0); syn_sigma = S_k_p_l1_avg/syn_snr; end;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% syn_sigma %0.6f: %0.3fs',syn_sigma,tmp_t)); end;
%%%%%%%%
% Set up translation distribution. ;
%%%%%%%%;
% Note that a 2d isotropic gaussian with std delta_sigma has the property ;
% that a fraction (1-exp(-R^2/(2*delta_sigma^2))) is contained within radius R. ;
% For example, try out: ;
% tmp_delta_sigma = 1.45; tmp_R = 3.3; tmp_A_ = randn(1024*64,2)*tmp_delta_sigma; tmp_r_ = sqrt(sum(tmp_A_.^2,2)); tmp_0in = numel(find(tmp_r_<tmp_R))/numel(tmp_r_); tmp_1in = 1-exp(-tmp_R^2/(2*tmp_delta_sigma^2)); disp(sprintf(' %% I_monte: %0.6f I_form: %0.6f',tmp_0in,tmp_1in));
% Thus, to ensure that 95-percent of the true-displacements lie within delta_r_max, ;
% we require that exp(-delta_r_max^2/(2*delta_sigma^2)) = 0.05 ;
% or that delta_sigma = sqrt(delta_r_max^2/log(20^2)). ;
% Roughly speaking, delta_sigma should approximately equal delta_r_max/2.5. ;
% Conversely: if we know delta_sigma, we can conclude that delta_r_max = delta_sigma*sqrt(log(20^2)). ; 
%%%%%%%%;
delta_r_p = 0.05;
delta_r_s = delta_r_max/sqrt(2*log(1/delta_r_p));
delta_r_N = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2));
true_viewing_delta_x_ = delta_r_s*randn(syn_n_M,1);
true_viewing_delta_y_ = delta_r_s*randn(syn_n_M,1);
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
disp(sprintf(' %% p-val %0.4f delta_r_max %0.6f sigma %0.4f N_pixel %0.4f',delta_r_p,delta_r_max,delta_r_s,delta_r_N));
%%%%%%%%;
% Generate synthetic images from templates S_k_p__. ;
%%%%%%%%;
tmp_t = tic();
true_viewing_polar_a_ = zeros(syn_n_M,1);
true_viewing_azimu_b_ = zeros(syn_n_M,1);
true_viewing_gamma_z_ = 2*pi*rand(syn_n_M,1);
syn_M_molecule_type_ = zeros(syn_n_M,1);
syn_M_k_p__ = zeros(n_w_sum,syn_n_M);
syn_UX_M_k_p___ = zeros(n_w_max,syn_n_M,syn_n_UX_rank);
syn_UX_M_k_q___ = zeros(n_w_max,syn_n_M,syn_n_UX_rank);
syn_S_permutation_ = randperm(n_S)-1;
nS=0;
for nM=0:syn_n_M-1;
tmp_d = nM/syn_n_M;
syn_M_molecule_type = length(find(cumsum(molecule_density_)<tmp_d));
syn_M_molecule_type_(1+nM) = syn_M_molecule_type;
if (mod(nM,100)==0); disp(sprintf(' %% nM %d/%d type %d',nM,syn_n_M,syn_M_molecule_type)); end;
tmp_M_k_p_ = S_k_p___{1+syn_M_molecule_type_(1+nM)}(:,1+syn_S_permutation_(1+nS)).*CTF_uni_avg_k_p_; %<-- extract random template. ;
true_viewing_polar_a_(1+nM) = viewing_polar_a_all_(1+syn_S_permutation_(1+nS));
true_viewing_azimu_b_(1+nM) = viewing_azimu_b_all_(1+syn_S_permutation_(1+nS));
tmp_M_k_p_ = rotate_p2p_fx(n_k_p_r,n_w_,n_w_sum,tmp_M_k_p_,+true_viewing_gamma_z_(1+nM));
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_p_,-true_viewing_delta_x_(1+nM),-true_viewing_delta_y_(1+nM));
tmp_M_k_p_ = tmp_M_k_p_ + syn_sigma*randn(n_w_sum,1)./sqrt( template_area_element_ * n_w_sum / (pi*k_p_r_max^2) );
syn_M_k_p__(:,1+nM) = tmp_M_k_p_;
tmp_M_x_c_ = interp_k_p_to_x_c_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_M_k_p_);
tmp_M_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1)) ;
tmp_M_k_p__ = reshape(tmp_M_k_p_,[n_w_max,n_k_p_r]);
for nUX_rank=0:syn_n_UX_rank-1;
tmp_UX_M_k_p_ = zeros(n_w_max,1);
for nk_p_r=0:n_k_p_r-1;
tmp_UX_M_k_p_ = tmp_UX_M_k_p_ + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_M_k_p__(:,1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
syn_UX_M_k_p___(:,1+nM,1+nUX_rank) = tmp_UX_M_k_p_;
syn_UX_M_k_q___(:,1+nM,1+nUX_rank) = interp_p_to_q(1,n_w_max,n_w_max,tmp_UX_M_k_p_);
end;%for nUX_rank=0:syn_n_UX_rank-1;
nS=nS+1; if (nS>=n_S); nS=0; end;
end;%for nM=0:syn_n_M-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% syn_M_k_p__: %0.3fs',tmp_t)); end;
%%%%%%%%;
% prepare precomputation for ampm. ;
%%%%%%%%;
pm_n_UX_rank = syn_n_UX_rank;
tmp_t = tic();
syn_M_k_q__ = zeros(n_w_sum,syn_n_M);
for nM=0:syn_n_M-1;
syn_M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,syn_M_k_p__(:,1+nM));
end;%for nM=0:syn_n_M-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% syn_M_k_q__: %0.3fs',tmp_t)); end;
tmp_t = tic();
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,syn_n_M,syn_M_k_q__,pm_n_UX_rank,UX__,X_weight_r_);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
tmp_t = tic();
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,syn_n_M,pm_n_UX_rank,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
%%%%%%%%;
clear S_k_p___;
for nmolecule=0:n_molecule-1;
tmp_d = length(find(syn_M_molecule_type_==nmolecule));
disp(sprintf(' %% molecule_type: %d <-- %d/%d = %0.2f',nmolecule,tmp_d,syn_n_M,tmp_d/syn_n_M));
end;%for nmolecule=0:n_molecule-1;
%%%%%%%%;
% use true euler-angles and displacements to solve for 'true' model (across all shells). ;
%%%%%%%%;
syn_a_k_Y_0lsq_ = cg_lsq_4(syn_n_order,n_k_p_r,k_p_r_,l_max_,n_w_,syn_n_M,syn_M_k_p__,1,CTF_uni_avg_k_p_,true_viewing_polar_a_,true_viewing_azimu_b_,true_viewing_gamma_z_,true_viewing_delta_x_,true_viewing_delta_y_);
%%%%%%%%;
% Compare current model (across all shells) to a_k_Y_quad_mavg_ (across all shells). ;
%%%%%%%%;
a_k_Y_quad_mavg_ = zeros(n_lm_sum,1);
for nmolecule=0:n_molecule-1;
a_k_Y_quad_mavg_ = a_k_Y_quad_mavg_ + molecule_density_(1+nmolecule)*a_k_Y_quad__(:,1+nmolecule);
end;%for nmolecule=0:n_molecule-1;
[syn_a_k_Y_0lsq_X_best,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_mavg_,syn_a_k_Y_0lsq_);
disp(sprintf(' %% [across all shells] a_k_Y_quad_mavg_ vs syn_a_k_Y_0lsq_: correlation %+0.6f',syn_a_k_Y_0lsq_X_best));
%%%%%%%%;
save(fname_mat ...
     ,'n_lm_','n_lm_max','n_lm_sum','n_lm_csum_','l_max_max','m_max_','n_m_max' ...
     ,'n_w_','weight_2d_k_p_r_','weight_2d_k_all_','n_viewing_all','viewing_azimu_b_all_','viewing_polar_a_all_','n_viewing_polar_a','viewing_polar_a_','n_viewing_azimu_b_','template_k_c_0__','template_k_c_1__','template_k_c_2__' ...
     ,'n_w_max','n_w_sum','n_w_csum' ...
     ,'n_S','S_k_p_l1_avg','syn_sigma' ...
     ,'delta_r_max','delta_r_p','delta_r_s','delta_r_N','svd_eps','n_delta_v_requested','FTK' ...
     ,'true_viewing_delta_x_','true_viewing_delta_y_' ...
     ,'true_viewing_polar_a_' ...
     ,'true_viewing_azimu_b_' ...
     ,'true_viewing_gamma_z_' ...
     ,'molecule_density_' ...
     ,'syn_M_molecule_type_' ...
     ,'syn_M_k_p__' ...
     ,'syn_UX_M_k_p___' ...
     ,'syn_UX_M_k_q___' ...
     ,'syn_S_permutation_' ...
     ,'svd_VUXM_lwnM____','UX_M_l2_dM__' ...
     ,'syn_a_k_Y_0lsq_' ...
     ,'a_k_Y_quad_mavg_' ...
     ,'syn_a_k_Y_0lsq_X_best' ...
     ,'-v7.3' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
fname_fig = sprintf('%s_mat/dir_jpg/tpmhtamsux_%s_M_k_p___n%.3ds%.4dr%.3d',dir_trunk,syn_infix,syn_n_M,floor(1000*syn_snr),syn_n_UX_rank);
if (flag_plot & ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
colormap(colormap_beach());
nUX_rank=1-1; clim_ = 1.5*std(real(syn_UX_M_k_p___(:,:,1+nUX_rank)),1,'all')*[-1,+1];
nUX_rank=min(syn_n_UX_rank-1, 1-1); subplot(1,2,1); imagesc(real(squeeze(syn_UX_M_k_p___(:,:,1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
if syn_n_UX_rank>=2;
nUX_rank=min(syn_n_UX_rank-1, 2-1); subplot(1,2,2); imagesc(real(squeeze(syn_UX_M_k_p___(:,:,1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
end;%if syn_n_UX_rank>=2;
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_mat/dir_jpg/tpmhtamsux_%s_M_k_p___A_n%.3ds%.4dr%.3d',dir_trunk,syn_infix,syn_n_M,floor(1000*syn_snr),syn_n_UX_rank);
if (flag_plot & ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1);
for nUX_rank=0:min(3,syn_n_UX_rank)-1;
subplot(1,min(3,syn_n_UX_rank),1+nUX_rank); hold on;
for nM=0:32:syn_n_M-1;
nc = max(0,min(n_c-1,floor(n_c*nM/syn_n_M)));
plot(2*pi*(0:n_w_max-1)/n_w_max,real(syn_UX_M_k_p___(:,1+nM,1+nUX_rank)),'.','Color',c_(1+nc,:)); 
end;%for nM=0:syn_n_M-1;
xlim([0,2*pi]);
%ylim(4e-8*[-1,+1]); 
xlabel('gamma_z','Interpreter','none');
ylabel('real(syn_UX_M_k_p)','Interpreter','none');
title(sprintf('rank %d',nUX_rank));
end; %for nUX_rank=0:min(15,syn_n_UX_rank)-1;
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

rng(syn_rseed);

fname_0 = sprintf('%s_mat/tpmhtamsux_%s_n%.3ds%.4d',dir_trunk,syn_infix,syn_n_M,floor(1000*syn_snr));

tmp_t = tic();
syn_M_k_q__ = zeros(n_w_sum,syn_n_M);
for nM=0:syn_n_M-1;
syn_M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,syn_M_k_p__(:,1+nM));
end;%for nM=0:syn_n_M-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% syn_M_k_q__: %0.3fs',tmp_t)); end;

pm_n_UX_rank = syn_n_UX_rank;
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_);
pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_N_pixel = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2)) ; 

%%%%%%%%;
% Now set up alternating minimization for 'MS-phase' of successive combinations of principled-image-rings. ;
%%%%%%%%;
for nUX_rank=0:syn_n_UX_rank-1;
fname_2 = sprintf('nUX%.3drng%.3d',nUX_rank,syn_rseed);
MS_fname_pre = sprintf('%s_MS_%s',fname_0,fname_2);
MS_fname_mat = sprintf('%s.mat',MS_fname_pre);
%MS_fname_tmp = sprintf('%s.tmp',MS_fname_pre);
if (~exist(MS_fname_mat,'file'));
flag_skip=0;
MS_fname_tmp = sprintf('%s.tmp',MS_fname_pre);
if ( exist(MS_fname_tmp,'file'));
tmp_date_diff = datenum(clock) - datenum(dir(MS_fname_tmp).date);
if (tmp_date_diff< date_diff_threshold);
disp(sprintf(' %% %s found, tmp_date_diff = %0.2f = recent, skipping',MS_fname_tmp,tmp_date_diff));
flag_skip=1;
end;%if (tmp_date_diff< date_diff_threshold);
if (tmp_date_diff>=date_diff_threshold);
disp(sprintf(' %% %s found, tmp_date_diff = %0.2f = stale, deleting',MS_fname_tmp,tmp_date_diff));
delete(MS_fname_tmp);
flag_skip=0;
end;%if (tmp_date_diff>=date_diff_threshold);
end;%if ( exist(MS_fname_tmp,'file'));
if (~flag_skip);
disp(sprintf(' %% %s not found, creating',MS_fname_pre));
save(MS_fname_tmp,'MS_fname_mat');
%%%%%%%%;
tmp_rseed=syn_rseed;tmp_n_iteration=syn_n_iteration;tmp_n_iteration_register=1;
tmp_viewing_k_eq_d = 1/(2*pi)/sqrt(1);tmp_n_order = syn_n_order;
tmp_pm_n_UX_rank = 1+nUX_rank;
tmp_pm_n_k_p_r = tmp_pm_n_UX_rank;
tmp_pm_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_k_p_r_max = 1;
tmp_pm_n_w_ = n_w_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_w_max = n_w_max;
tmp_pm_n_w_sum = sum(tmp_pm_n_w_);
tmp_pm_n_w_csum_ = cumsum([0;tmp_pm_n_w_]);
tmp_pm_l_max_ = l_max_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_lm_ = (1+tmp_pm_l_max_).^2; tmp_pm_n_lm_sum = sum(tmp_pm_n_lm_);
tmp_pm_weight_3d_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_weight_2d_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
flag_MS_vs_SM = 1;
f_rand = 0;
a_UX_Y_true_ = reshape(a_UX_Y_quad_mavg__(:,1:tmp_pm_n_UX_rank),[n_lm_max*tmp_pm_n_UX_rank,1]); %<-- use as ground truth for this set of principled-image-rings. ;
CTF_index_ = 1;
CTF_k_p__ = ones(n_w_sum,1);
tmp_t = tic();
[X_best_MS_...
,a_UX_Y_0lsq_MS__...
,euler_polar_a_MS__...
,euler_azimu_b_MS__...
,euler_gamma_z_MS__...
,image_delta_x_MS__...
,image_delta_y_MS__...
,image_I_value_MS__...
,image_X_value_MS__...
,image_S_index_MS__...
] = ...
ampm_4(...
 tmp_rseed...
,tmp_n_iteration...
,tmp_n_iteration_register...
,tmp_viewing_k_eq_d...
,tmp_n_order...
,n_k_p_r...
,k_p_r_...
,k_p_r_max...
,weight_3d_k_p_r_...
,weight_2d_k_p_r_...
,delta_r_max...
,svd_eps...
,n_delta_v_requested...
,FTK...
,n_w_...
,tmp_pm_n_UX_rank...
,UX__(:,1:tmp_pm_n_UX_rank)...
,X_weight_r_...
,syn_n_M...
,syn_M_k_p__...
,syn_M_k_q__...
,CTF_index_...
,CTF_k_p__...
,svd_VUXM_lwnM____(:,:,1:tmp_pm_n_UX_rank,:)...
,[]...
,l_max_...
,a_UX_Y_true_...
,[]...
,[]...
,[]...
,[]...
,[]...
,[]...
,flag_MS_vs_SM...
,f_rand...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% X_best_MS_: %0.3fs',tmp_t)); end;
X_best_MS_time = tmp_t;
tmp_t = tic();
tmp_n_residual_loading = 3; tmp_n_residual_iteration = 32;
[UX_M_k_q_wnM___,UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,tmp_pm_n_UX_rank,syn_n_M,svd_VUXM_lwnM____(:,:,1:tmp_pm_n_UX_rank,:),image_delta_x_MS__(:,end),image_delta_y_MS__(:,end));
syn_M_loading_MS__ = get_loading_0(tmp_n_residual_loading,tmp_n_residual_iteration,syn_n_order,tmp_pm_n_k_p_r,tmp_pm_weight_3d_k_p_r_,tmp_pm_l_max_,tmp_pm_n_w_,syn_n_M,reshape(UX_M_k_p_wnM___,[tmp_pm_n_w_sum,syn_n_M]),euler_polar_a_MS__(:,end),euler_azimu_b_MS__(:,end),euler_gamma_z_MS__(:,end));
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% syn_M_loading_MS__: %0.3fs',tmp_t)); end;
syn_M_loading_MS_time = tmp_t;
save(MS_fname_mat ...
     ,'X_best_MS_','X_best_MS_time' ...
     ,'a_UX_Y_0lsq_MS__' ...
     ,'euler_polar_a_MS__','euler_azimu_b_MS__','euler_gamma_z_MS__' ...
     ,'image_delta_x_MS__','image_delta_y_MS__','image_X_value_MS__','image_S_index_MS__','image_I_value_MS__' ...
     ,'syn_M_loading_MS__','syn_M_molecule_type_','syn_M_loading_MS_time' ...
     );
%%%%%%%%;
delete(MS_fname_tmp);
end;%if (~flag_skip);
end;%if (~exist(MS_fname_mat,'file'));
if ( exist(MS_fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',MS_fname_pre));
end;%if ( exist(MS_fname_mat,'file'));
end;%for nUX_rank=0:syn_n_UX_rank-1;

%%%%%%%%;
% Now set up alternating minimization for 'SM-phase' of successive combinations of principled-image-rings. ;
%%%%%%%%;
for nUX_rank=0:syn_n_UX_rank-1;
fname_2 = sprintf('nUX%.3drng%.3d',nUX_rank,syn_rseed);
MS_fname_pre = sprintf('%s_MS_%s',fname_0,fname_2);
MS_fname_mat = sprintf('%s.mat',MS_fname_pre);
%MS_fname_tmp = sprintf('%s.tmp',MS_fname_pre);
tmp_rseed=syn_rseed;tmp_n_iteration=syn_n_iteration;tmp_n_iteration_register=1;
tmp_viewing_k_eq_d = 1/(2*pi)/sqrt(1);tmp_n_order = syn_n_order;
tmp_pm_n_UX_rank = 1+nUX_rank;
tmp_pm_n_k_p_r = tmp_pm_n_UX_rank;
tmp_pm_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_k_p_r_max = 1;
tmp_pm_n_w_ = n_w_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_w_max = n_w_max;
tmp_pm_n_w_sum = sum(tmp_pm_n_w_);
tmp_pm_n_w_csum_ = cumsum([0;tmp_pm_n_w_]);
tmp_pm_l_max_ = l_max_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_lm_ = (1+tmp_pm_l_max_).^2; tmp_pm_n_lm_sum = sum(tmp_pm_n_lm_);
tmp_pm_weight_3d_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_weight_2d_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
if (~exist(MS_fname_mat,'file'));
disp(sprintf(' %% %s not found, skipping',MS_fname_pre));
end;%if (~exist(MS_fname_mat,'file'));
if ( exist(MS_fname_mat,'file'));
disp(sprintf(' %% %s found, loading',MS_fname_pre));
MS_tmp_ = load(MS_fname_mat);
%%%%%%%%;
SM_fname_pre = sprintf('%s_SM_%s',fname_0,fname_2);
SM_fname_mat = sprintf('%s.mat',SM_fname_pre);
%SM_fname_tmp = sprintf('%s.tmp',SM_fname_pre);
if (~exist(SM_fname_mat,'file'));
flag_skip=0;
SM_fname_tmp = sprintf('%s.tmp',SM_fname_pre);
if ( exist(SM_fname_tmp,'file'));
tmp_date_diff = datenum(clock) - datenum(dir(SM_fname_tmp).date);
if (tmp_date_diff< date_diff_threshold);
disp(sprintf(' %% %s found, tmp_date_diff = %0.2f = recent, skipping',SM_fname_tmp,tmp_date_diff));
flag_skip=1;
end;%if (tmp_date_diff< date_diff_threshold);
if (tmp_date_diff>=date_diff_threshold);
disp(sprintf(' %% %s found, tmp_date_diff = %0.2f = stale, deleting',SM_fname_tmp,tmp_date_diff));
delete(SM_fname_tmp);
flag_skip=0;
end;%if (tmp_date_diff>=date_diff_threshold);
end;%if ( exist(SM_fname_tmp,'file'));
if (~flag_skip);
disp(sprintf(' %% %s not found, creating',SM_fname_pre));
save(SM_fname_tmp,'SM_fname_mat');
%%%%%%%%;
tmp_euler_polar_a_MS_ = MS_tmp_.euler_polar_a_MS__(:,end);
tmp_euler_azimu_b_MS_ = MS_tmp_.euler_azimu_b_MS__(:,end);
tmp_euler_gamma_z_MS_ = MS_tmp_.euler_gamma_z_MS__(:,end);
tmp_image_delta_x_MS_ = MS_tmp_.image_delta_x_MS__(:,end);
tmp_image_delta_y_MS_ = MS_tmp_.image_delta_y_MS__(:,end);
tmp_image_I_value_MS_ = MS_tmp_.image_I_value_MS__(:,end);
flag_MS_vs_SM = 0; f_rand = f_rand_0in;
%%%%%%%%;
tmp_rseed=syn_rseed;tmp_n_iteration=syn_n_iteration;tmp_n_iteration_register=1;
tmp_viewing_k_eq_d = 1/(2*pi)/sqrt(1);tmp_n_order = syn_n_order;
tmp_pm_n_UX_rank = 1+nUX_rank;
flag_MS_vs_SM = 0;
a_UX_Y_true_ = reshape(a_UX_Y_quad_mavg__(:,1:tmp_pm_n_UX_rank),[n_lm_max*tmp_pm_n_UX_rank,1]); %<-- use as ground truth for this set of principled-image-rings. ;
CTF_index_ = 1;
CTF_k_p__ = ones(n_w_sum,1);
tmp_t = tic();
[X_best_SM_...
,a_UX_Y_0lsq_SM__...
,euler_polar_a_SM__...
,euler_azimu_b_SM__...
,euler_gamma_z_SM__...
,image_delta_x_SM__...
,image_delta_y_SM__...
,image_I_value_SM__...
,image_X_value_SM__...
,image_S_index_SM__...
] = ...
ampm_4(...
 tmp_rseed...
,tmp_n_iteration...
,tmp_n_iteration_register...
,tmp_viewing_k_eq_d...
,tmp_n_order...
,n_k_p_r...
,k_p_r_...
,k_p_r_max...
,weight_3d_k_p_r_...
,weight_2d_k_p_r_...
,delta_r_max...
,svd_eps...
,n_delta_v_requested...
,FTK...
,n_w_...
,tmp_pm_n_UX_rank...
,UX__(:,1:tmp_pm_n_UX_rank)...
,X_weight_r_...
,syn_n_M...
,syn_M_k_p__...
,syn_M_k_q__...
,CTF_index_...
,CTF_k_p__...
,svd_VUXM_lwnM____(:,:,1:tmp_pm_n_UX_rank,:)...
,[]...
,l_max_...
,a_UX_Y_true_...
,tmp_euler_polar_a_MS_...
,tmp_euler_azimu_b_MS_...
,tmp_euler_gamma_z_MS_...
,tmp_image_delta_x_MS_...
,tmp_image_delta_y_MS_...
,tmp_image_I_value_MS_...
,flag_MS_vs_SM...
,f_rand...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% X_best_SM_: %0.3fs',tmp_t)); end;
X_best_SM_time = tmp_t;
tmp_t = tic();
tmp_n_residual_loading = 3; tmp_n_residual_iteration = 32;
[UX_M_k_q_wnM___,UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,tmp_pm_n_UX_rank,syn_n_M,svd_VUXM_lwnM____(:,:,1:tmp_pm_n_UX_rank,:),image_delta_x_SM__(:,end),image_delta_y_SM__(:,end));
syn_M_loading_SM__ = get_loading_0(tmp_n_residual_loading,tmp_n_residual_iteration,syn_n_order,tmp_pm_n_k_p_r,tmp_pm_weight_3d_k_p_r_,tmp_pm_l_max_,tmp_pm_n_w_,syn_n_M,reshape(UX_M_k_p_wnM___,[tmp_pm_n_w_sum,syn_n_M]),euler_polar_a_SM__(:,end),euler_azimu_b_SM__(:,end),euler_gamma_z_SM__(:,end));
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% syn_M_loading_SM__: %0.3fs',tmp_t)); end;
syn_M_loading_SM_time = tmp_t;
save(SM_fname_mat ...
     ,'X_best_SM_','X_best_SM_time' ...
     ,'a_UX_Y_0lsq_SM__' ...
     ,'euler_polar_a_SM__','euler_azimu_b_SM__','euler_gamma_z_SM__' ...
     ,'image_delta_x_SM__','image_delta_y_SM__','image_X_value_SM__','image_S_index_SM__','image_I_value_SM__' ...
     ,'syn_M_loading_SM__','syn_M_molecule_type_','syn_M_loading_SM_time' ...
     );
%%%%%%%%;
delete(SM_fname_tmp);
end;%if (~flag_skip);
end;%if (~exist(SM_fname_mat,'file'));
if ( exist(SM_fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',SM_fname_pre));
end;%if ( exist(SM_fname_mat,'file'));
%%%%%%%%;
end;%if ( exist(MS_fname_mat,'file'));
end;%for nUX_rank=0:syn_n_UX_rank-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;





