function tpmham_experimental_4b(dat_infix,dat_rseed,dat_n_M,UX_M_k_p___,UX_M_k_q___,dat_n_UX_rank,dat_n_iteration,dat_n_order,dir_trunk,n_k_p_r,n_w_,l_max_,a_UX_Y_quad__);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% test_principled_marching_helper_alternating_minimization_experimental_4b. ;
% uses all shells. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

rng(0);
assert(dat_n_M<=size(UX_M_k_p___,2));
n_UX_rank = size(UX_M_k_p___,3);
l_max_max = max(l_max_);
n_w_max = max(n_w_); assert(n_w_max==size(UX_M_k_p___,1));
n_w_sum = sum(n_w_);
n_w_csum = cumsum([0;n_w_]);

%%%%%%%%;
fname_fig = sprintf('%s_mat/dir_jpg/tpmham_experimental_4_UX_%s_M_k_p___n%.3d',dir_trunk,dat_infix,dat_n_M);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
colormap(colormap_beach());
nUX_rank=1-1; clim_ = 1.5*std(real(UX_M_k_p___(:,1+(0:dat_n_M-1),1+nUX_rank)),1,'all')*[-1,+1];
nUX_rank= 1-1; subplot(2,3,1); imagesc(real(squeeze(UX_M_k_p___(:,1+(0:dat_n_M-1),1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank= 2-1; subplot(2,3,2); imagesc(real(squeeze(UX_M_k_p___(:,1+(0:dat_n_M-1),1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank= 3-1; subplot(2,3,3); imagesc(real(squeeze(UX_M_k_p___(:,1+(0:dat_n_M-1),1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank= 4-1; subplot(2,3,4); imagesc(real(squeeze(UX_M_k_p___(:,1+(0:dat_n_M-1),1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank= 5-1; subplot(2,3,5); imagesc(real(squeeze(UX_M_k_p___(:,1+(0:dat_n_M-1),1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=dat_n_UX_rank-1; subplot(2,3,6); imagesc(real(squeeze(UX_M_k_p___(:,1+(0:dat_n_M-1),1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_mat/dir_jpg/tpmham_experimental_4_UX_%s_M_k_p___A_n%.3d',dir_trunk,dat_infix,dat_n_M);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1);
for nUX_rank=0:min(15,dat_n_UX_rank)-1;
subplot(3,5,1+nUX_rank); hold on;
for nM=0:32:dat_n_M-1;
nc = max(0,min(n_c-1,floor(n_c*nM/dat_n_M)));
plot(2*pi*(0:n_w_max-1)/n_w_max,real(UX_M_k_p___(:,1+nM,1+nUX_rank)),'.','Color',c_(1+nc,:)); 
end;%for nM=0:dat_n_M-1;
xlim([0,2*pi]);
ylim(9*[-1,+1]); 
xlabel('gamma_z','Interpreter','none');
ylabel('real(UX_M_k_p)','Interpreter','none');
title(sprintf('rank %d',nUX_rank));
end; %for nUX_rank=0:min(15,dat_n_UX_rank)-1;
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

rng(dat_rseed);

%{
%%%%%%%%;
% Now set up alternating minimization for each single principled-image-ring. ;
%%%%%%%%;
for nUX_rank=0:dat_n_UX_rank-1;
fname_pre = sprintf('%s_mat/tpmham_experimental_4_UX_%s_Y_n%.3d_singlepc_rng%.3dnUX%.3d',dir_trunk,dat_infix,dat_n_M,dat_rseed,nUX_rank);
fname_mat = sprintf('%s.mat',fname_pre);
fname_tmp = sprintf('%s.tmp',fname_pre);
if (~exist(fname_mat,'file') & ~exist(fname_tmp,'file'));
disp(sprintf(' %% %s not found, creating',fname_pre));
save(fname_tmp,'fname_mat');
%%%%%%%%;
tmp_verbose=0;
tmp_n_k_p_r = 1;
tmp_l_max = l_max_max;
tmp_a_k_Y_quad_ = a_UX_Y_quad__(:,1+nUX_rank); %<-- use as ground truth for this principled-image-ring. ;
tmp_k_p_r = 1; tmp_k_p_r_max = 1;
tmp_n_w = n_w_max; tmp_n_w_max = tmp_n_w; tmp_n_w_sum = tmp_n_w; 
tmp_weight_k_p_r = 1; tmp_weight_2d_k_p_r = 1;
tmp_M_k_p__ = squeeze(UX_M_k_p___(:,1+(0:dat_n_M-1),1+nUX_rank));
tmp_M_k_q__ = squeeze(UX_M_k_q___(:,1+(0:dat_n_M-1),1+nUX_rank));
%%%%%%%%;
tmp_rseed=dat_rseed;tmp_n_iteration=dat_n_iteration;tmp_n_iteration_register=1;tmp_viewing_k_eq_d = 1/(2*pi);tmp_n_order = dat_n_order;
[dat_X_best_singlepc_,~,euler_polar_a__,euler_azimu_b__,euler_gamma_z__] = am_2(tmp_rseed,tmp_n_iteration,tmp_n_iteration_register,tmp_viewing_k_eq_d,tmp_n_order,tmp_n_k_p_r,tmp_k_p_r,tmp_k_p_r_max,tmp_weight_k_p_r,tmp_weight_2d_k_p_r,tmp_n_w,dat_n_M,tmp_M_k_p__,tmp_M_k_q__,1,ones(tmp_n_w_sum,1),tmp_l_max,tmp_a_k_Y_quad_,[],[],[],1,[]);
save(fname_mat ...
     ,'dat_X_best_singlepc_' ...
     ,'euler_polar_a__','euler_azimu_b__','euler_gamma_z__' ...
     );
%%%%%%%%;
delete(fname_tmp);
end;%if (~exist(fname_mat,'file') & ~exist(fname_tmp,'file'));
if ( exist(fname_mat,'file') | exist(fname_tmp,'file'));
disp(sprintf(' %% %s found, not creating',fname_pre));
end;%if ( exist(fname_mat,'file') | exist(fname_tmp,'file'));
end;%for nUX_rank=0:dat_n_UX_rank-1;
 %}

%%%%%%%%;
% Now set up alternating minimization for 'MS-phase' of successive combinations of principled-image-rings. ;
%%%%%%%%;
for nUX_rank=0:dat_n_UX_rank-1;
fname_pre = sprintf('%s_mat/tpmham_experimental_4_UX_%s_Y_n%.3d_allshell_rng%.3dnUX%.3d',dir_trunk,dat_infix,dat_n_M,dat_rseed,nUX_rank);
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
tmp_a_k_Y_quad_ = reshape(a_UX_Y_quad__(:,1+(0:nUX_rank)),[tmp_n_lm_sum,1]); %<-- use as ground truth for this set of principled-image-rings. ;
tmp_M_k_p__ = reshape(permute(UX_M_k_p___(:,1+(0:dat_n_M-1),1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,dat_n_M]);
tmp_M_k_q__ = reshape(permute(UX_M_k_q___(:,1+(0:dat_n_M-1),1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,dat_n_M]);
flag_MS_vs_SM = 1;
%%%%%%%%;
tmp_rseed=dat_rseed;tmp_n_iteration=dat_n_iteration;tmp_n_iteration_register=1;tmp_viewing_k_eq_d = 1/(2*pi);tmp_n_order = dat_n_order;
[dat_X_best_allshell_,~,euler_polar_a__,euler_azimu_b__,euler_gamma_z__] = am_2(tmp_rseed,tmp_n_iteration,tmp_n_iteration_register,tmp_viewing_k_eq_d,tmp_n_order,tmp_n_k_p_r,tmp_k_p_r_,tmp_k_p_r_max,tmp_weight_k_p_r_,tmp_weight_2d_k_p_r_,tmp_n_w_,dat_n_M,tmp_M_k_p__,tmp_M_k_q__,1,ones(tmp_n_w_sum,1),tmp_l_max_,tmp_a_k_Y_quad_,[],[],[],flag_MS_vs_SM,[]);
save(fname_mat ...
     ,'dat_X_best_allshell_' ...
     ,'euler_polar_a__','euler_azimu_b__','euler_gamma_z__' ...
     );
%%%%%%%%;
delete(fname_tmp);
end;%if (~exist(fname_mat,'file') & ~exist(fname_tmp,'file'));
if ( exist(fname_mat,'file') | exist(fname_tmp,'file'));
disp(sprintf(' %% %s found, not creating',fname_pre));
end;%if ( exist(fname_mat,'file') | exist(fname_tmp,'file'));
end;%for nUX_rank=0:dat_n_UX_rank-1;

%%%%%%%%;
% Now set up alternating minimization for 'SM-phase' of successive combinations of principled-image-rings. ;
%%%%%%%%;
for nUX_rank=0:dat_n_UX_rank-1;
fname_pre = sprintf('%s_mat/tpmham_experimental_4_UX_%s_Y_n%.3d_allshell_rng%.3dnUX%.3d',dir_trunk,dat_infix,dat_n_M,dat_rseed,nUX_rank);
fname_mat = sprintf('%s.mat',fname_pre);
fname_tmp = sprintf('%s.tmp',fname_pre);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, skipping',fname_pre));
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, loading',fname_pre));
MS_tmp_ = load(fname_mat);
%%%%%%%%;

SM_fname_pre = sprintf('%s_mat/tpmham_experimental_4_UX_%s_Y_n%.3d_allshell_SM_rng%.3dnUX%.3d',dir_trunk,dat_infix,dat_n_M,dat_rseed,nUX_rank);
SM_fname_mat = sprintf('%s.mat',SM_fname_pre);
SM_fname_tmp = sprintf('%s.tmp',SM_fname_pre);
if (~exist(SM_fname_mat,'file') & ~exist(SM_fname_tmp,'file'));
disp(sprintf(' %% %s not found, creating',SM_fname_pre));
save(SM_fname_tmp,'SM_fname_mat');
%%%%%%%%;
tmp_verbose=0;
tmp_n_k_p_r = 1+nUX_rank;
tmp_k_p_r_ = ones(1+nUX_rank,1); tmp_k_p_r_max = 1;
tmp_n_w_ = n_w_max*ones(1+nUX_rank,1); tmp_n_w_max = max(tmp_n_w_); tmp_n_w_sum = sum(tmp_n_w_); 
tmp_weight_k_p_r_ = ones(1+nUX_rank,1); tmp_weight_2d_k_p_r_ = ones(1+nUX_rank,1);
tmp_l_max_ = l_max_max*ones(1+nUX_rank,1);
tmp_n_lm_ = (1+tmp_l_max_).^2; tmp_n_lm_sum = sum(tmp_n_lm_);
tmp_a_k_Y_quad_ = reshape(a_UX_Y_quad__(:,1+(0:nUX_rank)),[tmp_n_lm_sum,1]); %<-- use as ground truth for this set of principled-image-rings. ;
tmp_M_k_p__ = reshape(permute(UX_M_k_p___(:,1+(0:dat_n_M-1),1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,dat_n_M]);
tmp_M_k_q__ = reshape(permute(UX_M_k_q___(:,1+(0:dat_n_M-1),1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,dat_n_M]);
tmp_euler_polar_a_MS_ = MS_tmp_.euler_polar_a__(:,end);
tmp_euler_azimu_b_MS_ = MS_tmp_.euler_azimu_b__(:,end);
tmp_euler_gamma_z_MS_ = MS_tmp_.euler_gamma_z__(:,end);
flag_MS_vs_SM = 0; f_rand = 0.05;
%%%%%%%%;
tmp_rseed=dat_rseed;tmp_n_iteration=dat_n_iteration;tmp_n_iteration_register=1;tmp_viewing_k_eq_d = 1/(2*pi);tmp_n_order = dat_n_order;
[dat_X_best_allshell_SM_,~,euler_polar_a_SM__,euler_azimu_b_SM__,euler_gamma_z_SM__] = am_2(tmp_rseed,tmp_n_iteration,tmp_n_iteration_register,tmp_viewing_k_eq_d,tmp_n_order,tmp_n_k_p_r,tmp_k_p_r_,tmp_k_p_r_max,tmp_weight_k_p_r_,tmp_weight_2d_k_p_r_,tmp_n_w_,dat_n_M,tmp_M_k_p__,tmp_M_k_q__,1,ones(tmp_n_w_sum,1),tmp_l_max_,tmp_a_k_Y_quad_,tmp_euler_polar_a_MS_,tmp_euler_azimu_b_MS_,tmp_euler_gamma_z_MS_,flag_MS_vs_SM,f_rand);
save(SM_fname_mat ...
     ,'dat_X_best_allshell_SM_' ...
     ,'euler_polar_a_SM__','euler_azimu_b_SM__','euler_gamma_z_SM__' ...
     );
%%%%%%%%;
delete(SM_fname_tmp);
end;%if (~exist(SM_fname_mat,'file') & ~exist(SM_fname_tmp,'file'));
if ( exist(SM_fname_mat,'file') | exist(SM_fname_tmp,'file'));
disp(sprintf(' %% %s found, not creating',SM_fname_pre));
end;%if ( exist(SM_fname_mat,'file') | exist(SM_fname_tmp,'file'));
%%%%%%%%;
end;%if ( exist(fname_mat,'file'));
end;%for nUX_rank=0:dat_n_UX_rank-1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;





