load('ampmh_X_wSM___5_debug.mat');
if (verbose); disp(sprintf(' %% [entering ampmh_X_wSM___5_experiment_3]')); end;

template_tree = get_template_tree_0(viewing_k_eq_d,5);

%%%%%%%%;
% Pull out CTF-modes for the principal-templates. ;
%%%%%%%%;
[n_S] = sample_shell_5(pm_k_p_r_max,viewing_k_eq_d,'L') ; %<-- obtain viewing angles on outer shell. ;
tmp_t = tic();
tmp_verbose=0;
UCTF_UX_S_k_p_wSc___ = zeros(pm_n_w_sum,n_S,n_CTF_rank);
for nCTF_rank=0:n_CTF_rank-1;
[ ...
 UCTF_UX_S_k_p_wSc___(:,:,1+nCTF_rank) ...
,~ ...
,~ ...
,~ ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
] = ...
get_template_1( ...
 tmp_verbose ...
,pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_k_p_r_max ...
,pm_weight_k_p_r_ ...
,pm_l_max_ ...
,a_UCTF_UX_Y_0lsq_ync__(:,1+nCTF_rank) ...
,viewing_k_eq_d ...
,-1 ...
,pm_n_w_ ...
);
assert(n_S==n_viewing_all);
end;%for nCTF_rank=0:n_CTF_rank-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% get_template_1: %0.3fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% group images by micrograph (i.e., sort by CTF_index_). ;
%%%%%%%%;
u_CTF_index_ = unique(CTF_index_(1:n_M)); n_u_CTF_index = numel(u_CTF_index_);
index_M_CTF_index__ = cell(n_u_CTF_index,1);
n_u_CTF_index_ = zeros(n_u_CTF_index,1);
for nu_CTF_index=0:n_u_CTF_index-1;
u_CTF_index = u_CTF_index_(1+nu_CTF_index);
index_M_CTF_index__{1+nu_CTF_index} = efind(CTF_index_(1:n_M)==u_CTF_index);
n_u_CTF_index_(1+nu_CTF_index) = numel(index_M_CTF_index__{1+nu_CTF_index});
end;%for nu_CTF_index=0:n_u_CTF_index-1;
if (verbose); disp(sprintf(' %% n_u_CTF_index %d',n_u_CTF_index)); end;
%%%%%%%%;
if (verbose);
for nu_CTF_index=0:n_u_CTF_index-1;
disp(sprintf(' %% nu_CTF_index %.3d/%.3d (%.8d) <-- n %.3d',nu_CTF_index,n_u_CTF_index,u_CTF_index_(1+nu_CTF_index),n_u_CTF_index_(1+nu_CTF_index)));
end;%for nu_CTF_index=0:n_u_CTF_index-1;
disp(sprintf(' %% n_M %.4d sum(n_u_CTF_index_) = %.4d',n_M,sum(n_u_CTF_index_)));
end;%if (verbose);

%%%%%%%%;
% Pick just one micrograph for testing. ;
% Calculate the templates associated with that particular CTF-function, ;
% and then calculate innerproducts between those templates and all the images. ;
%%%%%%%%;
nu_CTF_index = 0;%for nu_CTF_index=0:n_u_CTF_index-1;
tmp_index_M_ = index_M_CTF_index__{1+nu_CTF_index};
tmp_n_M = n_u_CTF_index_(1+nu_CTF_index);
if (verbose); disp(sprintf(' %% nu_CTF_index %d/%d --> tmp_n_M %d [%d,..,%d] ',nu_CTF_index,n_u_CTF_index,tmp_n_M,tmp_index_M_(0+1),tmp_index_M_(tmp_n_M-1+1))); end;
%%%%%%%%;
% combine CTF-modes to form templates for that particular CTF-function. ;
%%%%%%%%;
tmp_t = tic();
VSCTF_avg_ = mean(VSCTF_Mc__(1+tmp_index_M_,:),1);
VSCTF_std_ = std(VSCTF_Mc__(1+tmp_index_M_,:),1,1);
assert(max(VSCTF_std_./max(1e-12,abs(VSCTF_avg_)))<1e-6); %<-- consider lowering threshold to 1e-3. ;
CTF_UX_S_k_p_wnS__ = zeros(pm_n_w_sum,n_S);
for nCTF_rank=0:n_CTF_rank-1;
CTF_UX_S_k_p_wnS__ = CTF_UX_S_k_p_wnS__ + UCTF_UX_S_k_p_wSc___(:,:,1+nCTF_rank) * VSCTF_avg_(1+nCTF_rank);
end;%for nCTF_rank=0:n_CTF_rank-1;
%%%%%%%%;
CTF_UX_S_l2_ = zeros(n_S,1);
for nS=0:n_S-1;
CTF_UX_S_l2_(1+nS) = ...
innerproduct_p_quad( ...
 pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_weight_2d_k_p_r_/(2*pi) ...
,pm_n_w_ ...
,pm_n_w_sum ...
,CTF_UX_S_k_p_wnS__(:,1+nS) ...
,CTF_UX_S_k_p_wnS__(:,1+nS) ...
);
end;%for nS=0:n_S-1;
%%%%%%%%;
CTF_UX_S_k_q_wnS__ = zeros(pm_n_w_sum,n_S);
for nS=0:n_S-1;
CTF_UX_S_k_q_wnS__(:,1+nS) = ...
interp_p_to_q( ...
 pm_n_k_p_r ...
,pm_n_w_ ...
,pm_n_w_sum ...
,CTF_UX_S_k_p_wnS__(:,1+nS) ...
); 
end;%for nS=0:n_S-1; 
%%%%%%%%;

parameter = struct('type','parameter');
parameter.svd_eps_use = 0.01;
parameter.n_svd_l_use = 0;
parameter.n_delta_v_use = 48;
parameter.pm_n_UX_rank_use = pm_n_UX_rank;
parameter.n_w_max_use = 80;
parameter.flag_optimize_over_gamma_z = 1;

%%%%%%%%;
% Now test out template alignment for that particular tmp_n_M. ;
%%%%%%%%;
tmp_t = tic();
[ ...
 X_SM_full__ ...
,delta_x_SM_full__ ...
,delta_y_SM_full__ ...
,gamma_z_SM_full__ ...
,I_value_SM_full__ ...
] = ...
ampmh_X_local_SM__7( ...
 FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,n_S_per_Sbatch ...
,CTF_UX_S_k_p_wnS__ ...
,CTF_UX_S_l2_ ...
,tmp_n_M ...
,n_M_per_Mbatch ...
,svd_VUXM_lwnM____(:,:,:,1+tmp_index_M_) ...
,UX_M_l2_dM__(:,1+tmp_index_M_) ...
,parameter ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM_full__: %0.3fs',tmp_t)); end;
t_full = tmp_t;
f_full = n_S*tmp_n_M;

parameter.n_neighborhood_retain = 4;
%%%%%%%%;
% Compare to local search. ;
%%%%%%%%;
%%%%%%%%;
% Now test out template alignment for that particular tmp_n_M. ;
%%%%%%%%;
%profile on;
tmp_t = tic();
[ ...
 X_SM_tree__ ...
,delta_x_SM_tree__ ...
,delta_y_SM_tree__ ...
,gamma_z_SM_tree__ ...
,I_value_SM_tree__ ...
,t_level_SM_tree__ ...
] = ...
ampmh_X_local_SM__7( ...
 FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,n_S_per_Sbatch ...
,CTF_UX_S_k_p_wnS__ ...
,CTF_UX_S_l2_ ...
,tmp_n_M ...
,n_M_per_Mbatch ...
,svd_VUXM_lwnM____(:,:,:,1+tmp_index_M_) ...
,UX_M_l2_dM__(:,1+tmp_index_M_) ...
,parameter ...
,template_tree ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM_tree__: %0.3fs',tmp_t)); end;
%profile viewer;
%profile off;

flag_plot=0;
if flag_plot;
nM=0;
X_lim_ = [min(X_SM_full__(:,1+nM)),max(X_SM_full__(:,1+nM))];
p_row = 2; p_col = 3; np = 0;
figure(1);clf;figbig;figbeach();
%%%%%%%%;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_polar_a_azimu_b_0(template_tree.viewing_polar_a_all__{n_level},template_tree.viewing_azimu_b_all__{n_level},t_level_SM_tree__(:,1+nM),[-1,n_level],colormap_beach(),1);
xlim([0,2*pi]); ylim([0,1*pi]); axisnotick;
title(sprintf('nM %d t_level',nM),'Interpreter','none');
%%%%%%%%;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_polar_a_azimu_b_0(template_tree.viewing_polar_a_all__{n_level},template_tree.viewing_azimu_b_all__{n_level},X_SM_tree__(:,1+nM),X_lim_,colormap_beach(),1);
xlim([0,2*pi]); ylim([0,1*pi]); axisnotick;
title(sprintf('nM %d X_tree',nM),'Interpreter','none');
%%%%%%%%;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_polar_a_azimu_b_0(template_tree.viewing_polar_a_all__{n_level},template_tree.viewing_azimu_b_all__{n_level},X_SM_full__(:,1+nM),X_lim_,colormap_beach(),1);
xlim([0,2*pi]); ylim([0,1*pi]); axisnotick;
title(sprintf('nM %d X_full',nM),'Interpreter','none');
%%%%%%%%;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_polar_a_azimu_b_0(template_tree.viewing_polar_a_all__{n_level},template_tree.viewing_azimu_b_all__{n_level},gamma_z_SM_tree__(:,1+nM),[0,2*pi],colormap_beach(),1);
xlim([0,2*pi]); ylim([0,1*pi]); axisnotick;
title(sprintf('nM %d gamma_z_tree',nM),'Interpreter','none');
%%%%%%%%;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_polar_a_azimu_b_0(template_tree.viewing_polar_a_all__{n_level},template_tree.viewing_azimu_b_all__{n_level},gamma_z_SM_full__(:,1+nM),[0,2*pi],colormap_beach(),1);
xlim([0,2*pi]); ylim([0,1*pi]); axisnotick;
title(sprintf('nM %d gamma_z_full',nM),'Interpreter','none');
%%%%%%%%;
end;%if flag_plot;

%%%%%%%%;
% Compare to local search. ;
%%%%%%%%;
n_neighbor_retain_ = 2.^[0:5]; n_n_neighbor_retain = numel(n_neighbor_retain_);
X_p_tree__ = zeros(tmp_n_M,n_n_neighbor_retain);
X_p_full__ = zeros(tmp_n_M,n_n_neighbor_retain);
t_tree_ = zeros(n_n_neighbor_retain,1);
f_tree_ = zeros(n_n_neighbor_retain,1);
for nn_neighbor_retain=0:n_n_neighbor_retain-1;
n_neighbor_retain = n_neighbor_retain_(1+nn_neighbor_retain);
parameter.n_neighborhood_retain = n_neighbor_retain;
tmp_t = tic();
[ ...
 X_SM_tree__ ...
,delta_x_SM_tree__ ...
,delta_y_SM_tree__ ...
,gamma_z_SM_tree__ ...
,I_value_SM_tree__ ...
,t_level_SM_tree__ ...
] = ...
ampmh_X_local_SM__7( ...
 FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,n_S_per_Sbatch ...
,CTF_UX_S_k_p_wnS__ ...
,CTF_UX_S_l2_ ...
,tmp_n_M ...
,n_M_per_Mbatch ...
,svd_VUXM_lwnM____(:,:,:,1+tmp_index_M_) ...
,UX_M_l2_dM__(:,1+tmp_index_M_) ...
,parameter ...
,template_tree ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM_tree__: %0.3fs',tmp_t)); end;
t_tree_(1+nn_neighbor_retain) = tmp_t;
f_tree_(1+nn_neighbor_retain) = numel(efind(t_level_SM_tree__(:)>=0));
for nM=0:tmp_n_M-1;
[~,tmp_ij] = max(X_SM_tree__(:,1+nM));
tmp_p = numel(efind(X_SM_full__(:,1+nM)<=X_SM_full__(tmp_ij,1+nM)))/n_S;
X_p_tree__(1+nM,1+nn_neighbor_retain) = tmp_p;
tmp_index_ = efind(t_level_SM_tree__(:,1+nM)>=0);
[~,tmp_ij] = max(X_SM_tree__(1+tmp_index_,1+nM));
assert(fnorm(X_SM_tree__(1+tmp_index_(tmp_ij),1+nM)-X_SM_full__(1+tmp_index_(tmp_ij),1+nM))<1e-3);
tmp_p = numel(efind(X_SM_full__(:,1+nM)<=X_SM_full__(1+tmp_index_(tmp_ij),1+nM)))/n_S;
X_p_full__(1+nM,1+nn_neighbor_retain) = tmp_p;
end;%for nM=0:tmp_n_M-1;
end;%for nn_neighbor_retain=0:n_n_neighbor_retain-1;
%%%%%%%%;
figure(1);figbig;clf;
p_row=1;p_col=3; np=0;
%%%%%%%%;
subplot(p_row,p_col,1+np);np=np+1;
plot(0,f_full/t_full,'ko',1:n_n_neighbor_retain,f_tree_./t_tree_,'ro-');
xlabel('strategy');ylabel('comparisons/s');
%%%%%%%%;
subplot(p_row,p_col,1+np);np=np+1;
plot(t_tree_,mean(X_p_tree__,1),'ro-',t_full+[0,0],[0,1],'k-');
xlabel('time'); ylabel('p_tree','Interpreter','none');
xlim([0,2*t_full]);ylim([0,1]);
%%%%%%%%;
subplot(p_row,p_col,1+np);np=np+1;
plot(t_tree_,mean(X_p_full__,1),'ro-',t_full+[0,0],[0,1],'k-');
xlabel('time'); ylabel('p_full','Interpreter','none');
xlim([0,2*t_full]);ylim([0,1]);
%%%%%%%%;


%%%%%%%%;
%end;%for nu_CTF_index=0:n_u_CTF_index-1;
%%%%%%%%;



disp('returning'); return;
