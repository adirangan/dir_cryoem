function tpmham_synthetic_0(syn_rseed,syn_n_M,syn_snr,syn_n_UX_rank,syn_n_iteration,syn_n_order,dir_trunk,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,k_p_r_max,weight_k_p_r_,k_eq_d,S_k_p__,CTF_avg_k_p_,CTF_avg_k_p_r_,n_w_,weight_2d_k_p_r_,weight_2d_k_all_,n_viewing_all,viewing_azimu_b_all_,viewing_polar_a_all_,l_max_,a_k_Y_quad_,UX_,a_UX_Y_quad__);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% test_principled_marching_helper_alternating_minimization_synthetic_0. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

rng(0);

l_max_max = max(l_max_);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum = cumsum([0;n_w_]);

fname_mat = sprintf('%s_mat/tpmham_synthetic_UX_M_k_p___n%.3ds%.3dr%.3d.mat',dir_trunk,syn_n_M,floor(100*syn_snr),syn_n_UX_rank);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
n_S = n_viewing_all; %<-- number of templates stored in S_k_p__. ;
%%%%%%%%;
% Use templates S_k_p__ to estimate typical signal strength (per unit area in k_p). ;
%%%%%%%%;
template_area_element_ = weight_2d_k_all_*(2*pi)^2; %<-- now sum(template_area_element_)==pi*k_p_r_max^2;
S_k_p_l1_avg = mean(transpose(template_area_element_)*abs(S_k_p__.*CTF_avg_k_p_))/(pi*k_p_r_max^2);
if (syn_snr<=0); syn_sigma = 0; end;
if (syn_snr> 0); syn_sigma = S_k_p_l1_avg/syn_snr; end;
%%%%%%%%;
% First generate synthetic images from templates S_k_p__. ;
%%%%%%%%;
true_viewing_polar_a_ = zeros(syn_n_M,1);
true_viewing_azimu_b_ = zeros(syn_n_M,1);
true_viewing_gamma_z_ = 2*pi*rand(syn_n_M,1); 
syn_M_k_p__ = zeros(n_w_sum,syn_n_M);
syn_UX_M_k_p___ = zeros(n_w_max,syn_n_M,syn_n_UX_rank);
syn_UX_M_k_q___ = zeros(n_w_max,syn_n_M,syn_n_UX_rank);
syn_S_permutation_ = randperm(n_S)-1;
nS=0;
for nM=0:syn_n_M-1;
if (mod(nM,100)==0); disp(sprintf(' %% nM %d/%d',nM,syn_n_M)); end;
tmp_M_k_p_ = S_k_p__(:,1+syn_S_permutation_(1+nS)).*CTF_avg_k_p_; %<-- extract random template. ;
%tmp_M_k_p_ = S_k_p__(:,1+syn_S_permutation_(1+nS)); %<-- extract random template. ;
true_viewing_polar_a_(1+nM) = viewing_polar_a_all_(1+syn_S_permutation_(1+nS));
true_viewing_azimu_b_(1+nM) = viewing_azimu_b_all_(1+syn_S_permutation_(1+nS));
tmp_M_k_p_ = rotate_p2p_fx(n_k_p_r,n_w_,n_w_sum,tmp_M_k_p_,+true_viewing_gamma_z_(1+nM));
tmp_M_k_p_ = tmp_M_k_p_ + syn_sigma*randn(n_w_sum,1)./sqrt( template_area_element_ * n_w_sum / (pi*k_p_r_max^2) );
syn_M_k_p__(:,1+nM) = tmp_M_k_p_;
tmp_M_x_c_ = interp_k_p_to_x_c_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_M_k_p_);
tmp_M_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1)) ;
tmp_M_k_p__ = reshape(tmp_M_k_p_,[n_w_max,n_k_p_r]);
for nUX_rank=0:syn_n_UX_rank-1;
tmp_UX_M_k_p_ = zeros(n_w_max,1);
for nk_p_r=0:n_k_p_r-1;
tmp_UX_M_k_p_ = tmp_UX_M_k_p_ + UX_(1+nk_p_r,1+nUX_rank)*tmp_M_k_p__(:,1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
syn_UX_M_k_p___(:,1+nM,1+nUX_rank) = tmp_UX_M_k_p_;
syn_UX_M_k_q___(:,1+nM,1+nUX_rank) = interp_p_to_q(1,n_w_max,n_w_max,tmp_UX_M_k_p_);
end;%for nUX_rank=0:syn_n_UX_rank-1;
nS=nS+1; if (nS>=n_S); nS=0; end;
end;%for nM=0:tmp_n_M-1;
%%%%%%%%;
% use true euler-angles to solve for 'true' model (across all shells). ;
%%%%%%%%;
syn_a_k_Y_0lsq_ = cg_lsq_2(syn_n_order,n_k_p_r,l_max_,n_w_,syn_n_M,syn_M_k_p__,CTF_avg_k_p_r_,true_viewing_polar_a_,true_viewing_azimu_b_,true_viewing_gamma_z_);
%%%%%%%%;
% Compare current model (across all shells) to a_k_Y_quad_ (across all shells). ;
%%%%%%%%;
[syn_a_k_Y_0lsq_X_best,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_k_p_r_,0,l_max_,a_k_Y_quad_,syn_a_k_Y_0lsq_);
disp(sprintf(' %% [across all shells] a_k_Y_quad_ vs syn_a_k_Y_0lsq_: correlation %+0.6f',syn_a_k_Y_0lsq_X_best));
%%%%%%%%;
save(fname_mat ...
     ,'n_S','S_k_p_l1_avg','syn_sigma' ...
     ,'true_viewing_polar_a_' ...
     ,'true_viewing_azimu_b_' ...
     ,'true_viewing_gamma_z_' ...
     ,'syn_M_k_p__' ...
     ,'syn_UX_M_k_p___' ...
     ,'syn_UX_M_k_q___' ...
     ,'syn_S_permutation_' ...
     ,'syn_a_k_Y_0lsq_' ...
     ,'syn_a_k_Y_0lsq_X_best' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
fname_fig = sprintf('%s_mat/dir_jpg/tpmham_synthetic_UX_M_k_p___n%.3ds%.3dr%.3d',dir_trunk,syn_n_M,floor(100*syn_snr),syn_n_UX_rank);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
colormap(colormap_beach());
nUX_rank=1-1; clim_ = 1.5*std(real(syn_UX_M_k_p___(:,:,1+nUX_rank)),1,'all')*[-1,+1];
nUX_rank= 1-1; subplot(2,3,1); imagesc(real(squeeze(syn_UX_M_k_p___(:,:,1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank= 2-1; subplot(2,3,2); imagesc(real(squeeze(syn_UX_M_k_p___(:,:,1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank= 3-1; subplot(2,3,3); imagesc(real(squeeze(syn_UX_M_k_p___(:,:,1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank= 4-1; subplot(2,3,4); imagesc(real(squeeze(syn_UX_M_k_p___(:,:,1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank= 5-1; subplot(2,3,5); imagesc(real(squeeze(syn_UX_M_k_p___(:,:,1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=syn_n_UX_rank-1; subplot(2,3,6); imagesc(real(squeeze(syn_UX_M_k_p___(:,:,1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_mat/dir_jpg/tpmham_synthetic_UX_M_k_p___A_n%.3ds%.3dr%.3d',dir_trunk,syn_n_M,floor(100*syn_snr),syn_n_UX_rank);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1);
for nUX_rank=0:min(15,syn_n_UX_rank)-1;
subplot(3,5,1+nUX_rank); hold on;
for nM=0:32:syn_n_M-1;
nc = max(0,min(n_c-1,floor(n_c*nM/syn_n_M)));
plot(2*pi*(0:n_w_max-1)/n_w_max,real(syn_UX_M_k_p___(:,1+nM,1+nUX_rank)),'.','Color',c_(1+nc,:)); 
end;%for nM=0:syn_n_M-1;
xlim([0,2*pi]);
ylim(4e-8*[-1,+1]); 
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

%%%%%%%%;
% Now set up alternating minimization for each single principled-image-ring. ;
%%%%%%%%;
syn_X_best_singlepc_ = zeros(syn_n_iteration,1);
syn_X_best_allshell_ = zeros(syn_n_iteration,1);
for nUX_rank=0:syn_n_UX_rank-1;
fname_pre = sprintf('%s_mat/tpmham_synthetic_UX_Y_n%.3ds%.3dr%.3d_rng%.3dnUX%.3d',dir_trunk,syn_n_M,floor(100*syn_snr),syn_n_UX_rank,syn_rseed,nUX_rank);
fname_mat = sprintf('%s.mat',fname_pre);
fname_tmp = sprintf('%s.tmp',fname_pre);
if (~exist(fname_mat,'file') & ~exist(fname_tmp,'file'));
disp(sprintf(' %% %s not found, creating',fname_pre));
save(fname_tmp,'fname_mat');
%%%%%%%%;
tmp_verbose=0;
tmp_n_k_p_r = 1;
tmp_a_k_Y_quad_ = a_UX_Y_quad__(:,1+nUX_rank); %<-- use as ground truth for this principled-image-ring. ;
tmp_k_p_r = 1; tmp_n_w = n_w_max; tmp_n_w_max = tmp_n_w; tmp_n_w_sum = tmp_n_w; 
tmp_weight_k_p_r = 1; tmp_weight_2d_k_p_r = 1;
tmp_l_max = l_max_max;
tmp_n_polar_a = max(15,1+2*tmp_l_max); tmp_n_azimu_b = 1+2*tmp_n_polar_a;
[legendre_evaluate_ljm___,legendre_evaluate_mlj___,expil__,expi__] = legendre_evaluate_ljm___0(tmp_l_max,cos(linspace(0,pi,tmp_n_polar_a)),tmp_n_azimu_b);
tmp_M_k_p__ = squeeze(syn_UX_M_k_p___(:,:,1+nUX_rank));
tmp_M_k_q__ = squeeze(syn_UX_M_k_q___(:,:,1+nUX_rank));
tmp_MM_ = zeros(syn_n_M,1);
for nM=0:syn_n_M-1;
tmp_MM_(1+nM) = innerproduct_p_quad(tmp_n_k_p_r,tmp_k_p_r,tmp_weight_2d_k_p_r/(2*pi),tmp_n_w,tmp_n_w_sum,tmp_M_k_p__(:,1+nM),tmp_M_k_p__(:,1+nM));
end;%for nM=0:syn_n_M-1;
%%%%%%%%;
% initialize current euler-angles randomly. ;
%%%%%%%%;
euler_polar_a_ = pi*rand(syn_n_M,1); euler_azimu_b_ = 2*pi*rand(syn_n_M,1); euler_gamma_z_ = 2*pi*rand(syn_n_M,1);
for niteration=0:syn_n_iteration-1;
%%%%%%%%;
% use current euler-angles to solve for current model (on single shell). ;
%%%%%%%%;
[tmp_k_p_polar_a__,tmp_k_p_azimu_b__] = cg_rhs_1(syn_n_M,tmp_n_w,euler_polar_a_,euler_azimu_b_,+euler_gamma_z_);
tensor_to_scatter__ = cg_interpolate_n_1(syn_n_order,tmp_n_polar_a,tmp_n_azimu_b,tmp_n_w*syn_n_M,tmp_k_p_polar_a__(:),tmp_k_p_azimu_b__(:));
scatter_to_tensor__ = transpose(tensor_to_scatter__);
tmp_An__ = @(a_k_Y_) tensor_to_scatter__*reshape(cg_evaluate_n_1(tmp_l_max,convert_spharm_to_spharm__0(tmp_l_max,a_k_Y_),tmp_n_polar_a,tmp_n_azimu_b,legendre_evaluate_ljm___),[tmp_n_polar_a*tmp_n_azimu_b,1]);
tmp_At__ = @(a_k_X_) convert_spharm__to_spharm_0(tmp_l_max,cg_evaluate_t_1(tmp_n_polar_a,tmp_n_azimu_b,reshape(scatter_to_tensor__*a_k_X_,[tmp_n_polar_a,tmp_n_azimu_b]),tmp_l_max,legendre_evaluate_mlj___,expil__,expi__));
tmp_AtAn__ = @(a_k_Y_) tmp_At__(tmp_An__(a_k_Y_));
[tmp_a_k_Y_0lsq_,~] = pcg(tmp_AtAn__,tmp_At__(tmp_M_k_p__(:)));
%%%%%%%%;
% Compare current model (on single shell) to tmp_a_k_Y_quad_ (also on single shell). ;
%%%%%%%%;
[X_best,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(tmp_n_k_p_r,tmp_k_p_r,k_p_r_max,tmp_weight_k_p_r,0,tmp_l_max,tmp_a_k_Y_quad_,tmp_a_k_Y_0lsq_);
disp(sprintf(' %% niteration %d/%d [single shell]: tmp_a_k_Y_quad_ vs tmp_a_k_Y_lsq0_: correlation %+0.6f',niteration,syn_n_iteration,X_best));
syn_X_best_singlepc_(1+niteration) = X_best;
if (mod(niteration,8)==7); %<-- this takes too long. ;
%%%%%%%%;
% use current euler-angles to solve for current model (across all shells). ;
%%%%%%%%;
tmp_b_k_Y_0lsq_ = cg_lsq_2(syn_n_order,n_k_p_r,l_max_,n_w_,syn_n_M,syn_M_k_p__,CTF_avg_k_p_r_,euler_polar_a_,euler_azimu_b_,euler_gamma_z_);
%%%%%%%%;
% Compare current model (across all shells) to syn_a_k_Y_0lsq_ (across all shells). ;
%%%%%%%%;
[X_best,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_k_p_r_,0,l_max_,syn_a_k_Y_0lsq_,tmp_b_k_Y_0lsq_);
disp(sprintf(' %% niteration %d/%d [across all shells] syn_a_k_Y_0lsq_ vs tmp_b_k_Y_0lsq_: correlation %+0.6f',niteration,syn_n_iteration,X_best));
syn_X_best_allshell_(1+niteration) = X_best;
end;%if (mod(niteration,8)==7);
%%%%%%%%;
% use current model (on single shell) to generate current templates. ;
%%%%%%%%;
tmp_viewing_k_eq_d = k_eq_d*4; %<-- make this slightly different from viewing_k_eq_d to test out code. ;
[tmp_S_k_p__,~,~,~,tmp_n_viewing_all,tmp_viewing_azimu_b_all_,tmp_viewing_polar_a_all_,~,~,~,~,~,~] = get_template_0(tmp_verbose,tmp_n_k_p_r,tmp_k_p_r,k_p_r_max,tmp_weight_k_p_r,tmp_l_max,tmp_a_k_Y_0lsq_,tmp_viewing_k_eq_d,-1,tmp_n_w);
tmp_n_S = tmp_n_viewing_all;
tmp_SS_ = zeros(tmp_n_S,1);
for nS=0:tmp_n_S-1;
tmp_SS_(1+nS) = innerproduct_p_quad(tmp_n_k_p_r,tmp_k_p_r,tmp_weight_2d_k_p_r/(2*pi),tmp_n_w,tmp_n_w_sum,tmp_S_k_p__(:,1+nS),tmp_S_k_p__(:,1+nS));
end;%for nS=0:tmp_n_S-1;
for nS=0:tmp_n_S-1; tmp_S_k_q__(:,1+nS) = interp_p_to_q(tmp_n_k_p_r,tmp_n_w,tmp_n_w_sum,tmp_S_k_p__(:,1+nS)); end;%for nS=0:tmp_n_S-1; 
%%%%%%%%;
% Use current templates to calculate current innerproducts/correlations. ;
%%%%%%%%;
tmp_X___ = zeros(tmp_n_w_max,syn_n_M,tmp_n_S);
for nS=0:tmp_n_S-1;
for nM=0:syn_n_M-1;
tmp_X___(:,1+nM,1+nS) = ifft(innerproduct_q_k_stretch_quad_0(tmp_n_k_p_r,tmp_k_p_r,tmp_weight_2d_k_p_r/(2*pi),tmp_n_w,tmp_n_w_sum,tmp_S_k_q__(:,1+nS),tmp_M_k_q__(:,1+nM)))*tmp_n_w_max/sqrt(tmp_MM_(1+nM))/sqrt(tmp_SS_(1+nS)) ; %<-- multiplication by n_w_max not needed in fortran fftw_plan_back. ;
end;%for nM=0:syn_n_M-1;
end;%for nS=0:tmp_n_S-1;
%%%%%%%%;
% Use current correlations to udate current euler-angles. ;
%%%%%%%%;
flag_M_used_ = zeros(syn_n_M,1);
tmp_permutation_ = randperm(tmp_n_S)-1;
nS=0;
while (sum(flag_M_used_)<syn_n_M);
index_M_unused_ = find(flag_M_used_==0)-1;
[~,index_wM_best] = max(real(tmp_X___(:,1+index_M_unused_,1+tmp_permutation(1+nS))),[],'all','linear'); index_wM_best = index_wM_best-1;
[nw_best,index_M_best] = ind2sub([tmp_n_w_max,numel(index_M_unused_)],1+index_wM_best); 
nw_best = nw_best-1; index_M_best = index_M_best-1;
nM_best = index_M_unused_(1+index_M_best);
flag_M_used_(1+nM_best)=1;
euler_polar_a_(1+nM_best) = tmp_viewing_polar_a_all_(1+tmp_permutation_(1+nS));
euler_azimu_b_(1+nM_best) = tmp_viewing_azimu_b_all_(1+tmp_permutation_(1+nS));
euler_gamma_z_(1+nM_best) = 2*pi*nw_best/tmp_n_w_max;
nS = nS+1; if (nS>=tmp_n_S); nS=0; end;
end;%while (sum(flag_M_used_)<syn_n_M);
%%%%%%%%;
% Now return to beginning of loop. ;
%%%%%%%%;
end;%for niteration=0:syn_n_iteration-1;
%%%%%%%%;
save(fname_mat ...
     ,'syn_X_best_singlepc_' ...
     ,'syn_X_best_allshell_' ...
     );
%%%%%%%%;
delete(fname_tmp);
end;%if (~exist(fname_mat,'file') & ~exist(fname_tmp,'file'));
if ( exist(fname_mat,'file') | exist(fname_tmp,'file'));
disp(sprintf(' %% %s found, not creating',fname_pre));
end;%if ( exist(fname_mat,'file') | exist(fname_tmp,'file'));
end;%for nUX_rank=0:syn_n_UX_rank-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;





