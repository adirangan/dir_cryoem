function tpmhham_synthetic_5(...
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
,weight_k_p_r_...
,viewing_k_eq_d...
,template_k_eq_d...
,CTF_avg_k_p_...
,CTF_avg_k_p_r_...
,l_max_...
,n_molecule...
,molecule_density_...
,a_k_Y_quad__...
,UX_...
,X_weight_r_...
,a_UX_Y_quad_mavg__...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% test_principled_marching_heterogeneous_helper_alternating_minimization_synthetic_5. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

rng(0);

n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
molecule_density_ = molecule_density_/sum(molecule_density_);

fname_mat = sprintf('%s_mat/tpmhham_synthetic_5_UX_%s_M_k_p___n%.3ds%.4dr%.3d.mat',dir_trunk,syn_infix,syn_n_M,floor(1000*syn_snr),syn_n_UX_rank);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
verbose=1;
[~,n_w_,weight_2d_k_p_r_,weight_2d_k_all_,n_viewing_all,viewing_azimu_b_all_,viewing_polar_a_all_,n_viewing_polar_a,viewing_polar_a_,n_viewing_azimu_b_,template_k_c_0__,template_k_c_1__,template_k_c_2__] = get_template_0(verbose,n_k_p_r,k_p_r_,k_p_r_max,weight_k_p_r_,l_max_,[],viewing_k_eq_d,template_k_eq_d);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum = cumsum([0;n_w_]);
S_k_p___ = cell(n_molecule,1);
for nmolecule=0:n_molecule-1;
[S_k_p___{1+nmolecule}] = get_template_0(verbose,n_k_p_r,k_p_r_,k_p_r_max,weight_k_p_r_,l_max_,a_k_Y_quad__(:,1+nmolecule),viewing_k_eq_d,template_k_eq_d);
end;%for nmolecule=0:n_molecule-1;
n_S = n_viewing_all; %<-- number of templates stored in S_k_p__. ;
%%%%%%%%;
% Use templates S_k_p__ to estimate typical signal strength (per unit area in k_p). ;
%%%%%%%%;
template_area_element_ = weight_2d_k_all_*(2*pi)^2; %<-- now sum(template_area_element_)==pi*k_p_r_max^2;
S_k_p_l1_avg_ = zeros(n_molecule,1);
for nmolecule=0:n_molecule-1;
S_k_p_l1_avg_(1+nmolecule) = mean(transpose(template_area_element_)*abs(S_k_p___{1+nmolecule}.*CTF_avg_k_p_))/(pi*k_p_r_max^2);
end;%for nmolecule=0:n_molecule-1;
S_k_p_l1_avg = sum(molecule_density_.*S_k_p_l1_avg_);
if (syn_snr<=0); syn_sigma = 0; end;
if (syn_snr> 0); syn_sigma = S_k_p_l1_avg/syn_snr; end;
%%%%%%%%;
% First generate synthetic images from templates S_k_p__. ;
%%%%%%%%;
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
tmp_M_k_p_ = S_k_p___{1+syn_M_molecule_type_(1+nM)}(:,1+syn_S_permutation_(1+nS)).*CTF_avg_k_p_; %<-- extract random template. ;
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
tmp_UX_M_k_p_ = tmp_UX_M_k_p_ + UX_(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_M_k_p__(:,1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
syn_UX_M_k_p___(:,1+nM,1+nUX_rank) = tmp_UX_M_k_p_;
syn_UX_M_k_q___(:,1+nM,1+nUX_rank) = interp_p_to_q(1,n_w_max,n_w_max,tmp_UX_M_k_p_);
end;%for nUX_rank=0:syn_n_UX_rank-1;
nS=nS+1; if (nS>=n_S); nS=0; end;
end;%for nM=0:syn_n_M-1;
clear S_k_p___;
for nmolecule=0:n_molecule-1;
tmp_d = length(find(syn_M_molecule_type_==nmolecule));
disp(sprintf(' %% molecule_type: %d <-- %d/%d = %0.2f',nmolecule,tmp_d,syn_n_M,tmp_d/syn_n_M));
end;%for nmolecule=0:n_molecule-1;
%%%%%%%%;
% use true euler-angles to solve for 'true' model (across all shells). ;
%%%%%%%%;
syn_a_k_Y_0lsq_ = cg_lsq_3(syn_n_order,n_k_p_r,l_max_,n_w_,syn_n_M,syn_M_k_p__,1,CTF_avg_k_p_,true_viewing_polar_a_,true_viewing_azimu_b_,true_viewing_gamma_z_);
%%%%%%%%;
% Compare current model (across all shells) to a_k_Y_quad_mavg_ (across all shells). ;
%%%%%%%%;
a_k_Y_quad_mavg_ = zeros(n_lm_sum,1);
for nmolecule=0:n_molecule-1;
a_k_Y_quad_mavg_ = a_k_Y_quad_mavg_ + molecule_density_(1+nmolecule)*a_k_Y_quad__(:,1+nmolecule);
end;%for nmolecule=0:n_molecule-1;
[syn_a_k_Y_0lsq_X_best,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_k_p_r_,0,l_max_,a_k_Y_quad_mavg_,syn_a_k_Y_0lsq_);
disp(sprintf(' %% [across all shells] a_k_Y_quad_mavg_ vs syn_a_k_Y_0lsq_: correlation %+0.6f',syn_a_k_Y_0lsq_X_best));
%%%%%%%%;
save(fname_mat ...
     ,'n_lm_','n_lm_max','n_lm_sum','n_lm_csum_','l_max_max','m_max_','n_m_max'...
     ,'n_w_','weight_2d_k_p_r_','weight_2d_k_all_','n_viewing_all','viewing_azimu_b_all_','viewing_polar_a_all_','n_viewing_polar_a','viewing_polar_a_','n_viewing_azimu_b_','template_k_c_0__','template_k_c_1__','template_k_c_2__'...
     ,'n_w_max','n_w_sum','n_w_csum'...
     ,'n_S','S_k_p_l1_avg','syn_sigma' ...
     ,'true_viewing_polar_a_' ...
     ,'true_viewing_azimu_b_' ...
     ,'true_viewing_gamma_z_' ...
     ,'molecule_density_' ...
     ,'syn_M_molecule_type_' ...
     ,'syn_M_k_p__' ...
     ,'syn_UX_M_k_p___' ...
     ,'syn_UX_M_k_q___' ...
     ,'syn_S_permutation_' ...
     ,'syn_a_k_Y_0lsq_' ...
     ,'a_k_Y_quad_mavg_' ...
     ,'syn_a_k_Y_0lsq_X_best' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
fname_fig = sprintf('%s_mat/dir_jpg/tpmhham_synthetic_5_UX_%s_M_k_p___n%.3ds%.4dr%.3d',dir_trunk,syn_infix,syn_n_M,floor(1000*syn_snr),syn_n_UX_rank);
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
fname_fig = sprintf('%s_mat/dir_jpg/tpmhham_synthetic_5_UX_%s_M_k_p___A_n%.3ds%.4dr%.3d',dir_trunk,syn_infix,syn_n_M,floor(1000*syn_snr),syn_n_UX_rank);
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

%%%%%%%%;
% Now set up alternating minimization for successive combinations of principled-image-rings. ;
%%%%%%%%;
for nUX_rank=0:syn_n_UX_rank-1;
fname_pre = sprintf('%s_mat/tpmhham_synthetic_5_UX_%s_Y_n%.3ds%.4dr%.3d_allshell_rng%.3dnUX%.3d',dir_trunk,syn_infix,syn_n_M,floor(1000*syn_snr),syn_n_UX_rank,syn_rseed,nUX_rank);
fname_mat = sprintf('%s.mat',fname_pre);
fname_tmp = sprintf('%s.tmp',fname_pre);
if (~exist(fname_mat,'file') & ~exist(fname_tmp,'file'));
disp(sprintf(' %% %s not found, creating',fname_pre));
save(fname_tmp,'fname_mat');
%%%%%%%%;
tmp_verbose=0;
tmp_n_k_p_r = 1+nUX_rank;
tmp_k_p_r_ = ones(1+nUX_rank,1); tmp_k_p_r_max = 1;
tmp_n_w_ = n_w_max*ones(1+nUX_rank,1); tmp_n_w_max = max(tmp_n_w_); tmp_n_w_sum = sum(tmp_n_w_); 
tmp_weight_k_p_r_ = ones(1+nUX_rank,1); tmp_weight_2d_k_p_r_ = ones(1+nUX_rank,1);
tmp_l_max_ = l_max_max*ones(1+nUX_rank,1);
tmp_n_lm_ = (1+tmp_l_max_).^2; tmp_n_lm_sum = sum(tmp_n_lm_);
tmp_a_k_Y_quad_mavg_ = reshape(a_UX_Y_quad_mavg__(:,1+(0:nUX_rank)),[tmp_n_lm_sum,1]); %<-- use as ground truth for this set of principled-image-rings. ;
tmp_M_k_p__ = reshape(permute(syn_UX_M_k_p___(:,:,1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,syn_n_M]);
tmp_M_k_q__ = reshape(permute(syn_UX_M_k_q___(:,:,1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,syn_n_M]);
%%%%%%%%;
tmp_rseed=syn_rseed;tmp_n_iteration=syn_n_iteration;tmp_n_iteration_register=1;tmp_viewing_k_eq_d = 1/(2*pi);tmp_n_order = syn_n_order;
[syn_X_best_allshell_,syn_a_k_Y_0lsq__,euler_polar_a__,euler_azimu_b__,euler_gamma_z__] = am_2(tmp_rseed,tmp_n_iteration,tmp_n_iteration_register,tmp_viewing_k_eq_d,tmp_n_order,tmp_n_k_p_r,tmp_k_p_r_,tmp_k_p_r_max,tmp_weight_k_p_r_,tmp_weight_2d_k_p_r_,tmp_n_w_,syn_n_M,tmp_M_k_p__,tmp_M_k_q__,1,ones(tmp_n_w_sum,1),tmp_l_max_,tmp_a_k_Y_quad_mavg_,[],[],[],1,[]);
save(fname_mat ...
     ,'syn_X_best_allshell_' ...
     ,'syn_a_k_Y_0lsq__' ...
     ,'euler_polar_a__','euler_azimu_b__','euler_gamma_z__'...
     );
%%%%%%%%;
delete(fname_tmp);
end;%if (~exist(fname_mat,'file') & ~exist(fname_tmp,'file'));
if ( exist(fname_mat,'file') | exist(fname_tmp,'file'));
disp(sprintf(' %% %s found, not creating',fname_pre));
end;%if ( exist(fname_mat,'file') | exist(fname_tmp,'file'));
end;%for nUX_rank=0:syn_n_UX_rank-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;





