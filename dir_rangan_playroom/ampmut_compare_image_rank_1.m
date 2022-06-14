function ...
[ ...
 parameter ...
] = ...
ampmut_compare_image_rank_1( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,n_CTF ...
,n_CTF_rank ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,l_max_ ...
,a_k_Y_true_ ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
,X_TM_ ...
,X_fname_mat ...
,X_fname_compare_image_rank_pre ...
);

verbose=1;
if (verbose); disp(sprintf(' %% [entering ampmut_compare_image_rank_1]')); end;

if isempty(parameter); parameter = struct('type','parameter'); end;
if (~isfield(parameter,'flag_replot')); parameter.flag_replot = 0; end;
if (~isfield(parameter,'compare_image_rank_flag_nrm')); parameter.compare_image_rank_flag_nrm = 0; end;
flag_replot = parameter.flag_replot;
compare_image_rank_flag_nrm = parameter.compare_image_rank_flag_nrm;

if (~exist(X_fname_mat,'file'));
disp(sprintf(' %% %s not found, skipping',X_fname_mat));
end;%if (~exist(X_fname_mat,'file'));
if ( exist(X_fname_mat,'file'));
disp(sprintf(' %% %s found, not skipping',X_fname_mat));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
fname_fig = X_fname_compare_image_rank_pre;
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
l_max_max = max(l_max_);
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);

%%%%%%%%;
% construct CTF of same size as images. ;
%%%%%%%%
if (size(CTF_k_p_wkC__,1)==n_k_p_r);
CTF_k_p_r_kC__ = CTF_k_p_wkC__;
CTF_k_p_wkC__ = zeros(n_w_sum,n_CTF);
for nCTF=0:n_CTF-1;
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
tmp_index_ = n_w_csum_(1+nk_p_r) + [0:n_w-1];
CTF_k_p_wkC__(1+tmp_index_,1+nCTF) = CTF_k_p_r_kC__(1+nk_p_r,1+nCTF);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nCTF=0:n_CTF-1;
end;%if (size(CTF_k_p_wkC__,1)==n_k_p_r);

if isempty(X_TM_);
[ ...
 parameter ...
,X_TM_ ...
] = ...
ampmut_X_TM_0( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,n_CTF ...
,n_CTF_rank ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,l_max_ ...
,a_k_Y_true_ ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
);
end;%if isempty(X_TM_);

compare_image_rank_flag_nrm = 0;
tmp_t = tic();
M_k_p_l2_ = zeros(n_M,1);
M_k_p_nrm__ = M_k_p_wkM__;
for nM=0:n_M-1;
M_k_p_ = M_k_p_wkM__(:,1+nM);
tmp_l2 = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,M_k_p_,M_k_p_);
if compare_image_rank_flag_nrm; M_k_p_nrm__(:,1+nM) = M_k_p_/sqrt(tmp_l2); end;
M_k_p_l2_(1+nM) = tmp_l2;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% M_k_p_nrm__ time %0.2fs',tmp_t)); end;

%%%%%%%%;
% Load one of the reconstructed molecules. ;
%%%%%%%%;
tmp_t = tic();
tmp_ = load(X_fname_mat);
image_X_value_Mi__ = tmp_.image_X_value_Mi__;
n_iteration_use = size(image_X_value_Mi__,2)-1;
R_k_p_l2_Mi__ = zeros(n_M,n_iteration_use);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for niteration_use=0:n_iteration_use-1;
tmp_euler_polar_a_true_ = tmp_.euler_polar_a_Mi__(:,1+niteration_use+1);
tmp_euler_azimu_b_true_ = tmp_.euler_azimu_b_Mi__(:,1+niteration_use+1);
tmp_euler_gamma_z_true_ = tmp_.euler_gamma_z_Mi__(:,1+niteration_use+1);
tmp_image_delta_x_true_ = tmp_.image_delta_x_acc_Mi__(:,1+niteration_use+1) + tmp_.image_delta_x_upd_Mi__(:,1+niteration_use+1);
tmp_image_delta_y_true_ = tmp_.image_delta_y_acc_Mi__(:,1+niteration_use+1) + tmp_.image_delta_y_upd_Mi__(:,1+niteration_use+1);
tmp_image_X_value_true_ = tmp_.image_X_value_Mi__(:,1+niteration_use+0);
parameter.flag_loading_svd_vs_iterate = 1;
parameter.flag_loading_skip_loading = 1;
[ ...
 parameter ...
,~ ...
,~ ...
,tmp_T_k_p_wnM__ ...
,~ ...
,tmp_S_k_p_wnM__ ...
,tmp_R_k_p_wnM__ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
] = ...
get_loading_qbp_1( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_max*ones(n_k_p_r,1) ...
,n_w_ ...
,n_M ...
,M_k_p_nrm__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,tmp_euler_polar_a_true_ ...
,tmp_euler_azimu_b_true_ ...
,tmp_euler_gamma_z_true_ ...
,tmp_image_delta_x_true_ ...
,tmp_image_delta_y_true_ ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% get_loading_qbp_1 time %0.2fs',tmp_t)); end;
%%%%%%%%;
% Measure residual norm. ;
%%%%%%%%;
tmp_t = tic();
tmp_T_k_p_l2_ = zeros(n_M,1);
tmp_S_k_p_l2_ = zeros(n_M,1);
tmp_R_k_p_l2_ = zeros(n_M,1);
tmp_X_ = zeros(n_M,1);
for nM=0:n_M-1;
tmp_T_k_p_ = tmp_T_k_p_wnM__(:,1+nM);
tmp_S_k_p_ = tmp_S_k_p_wnM__(:,1+nM);
tmp_R_k_p_ = tmp_R_k_p_wnM__(:,1+nM);
tmp_l2 = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,tmp_T_k_p_,tmp_T_k_p_);
tmp_T_k_p_l2_(1+nM) = tmp_l2;
tmp_l2 = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,tmp_S_k_p_,tmp_S_k_p_);
tmp_S_k_p_l2_(1+nM) = tmp_l2;
tmp_l2 = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,tmp_R_k_p_,tmp_R_k_p_);
tmp_R_k_p_l2_(1+nM) = tmp_l2;
tmp_TT = tmp_T_k_p_l2_(1+nM);
tmp_SS = tmp_S_k_p_l2_(1+nM);
tmp_TS = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,tmp_T_k_p_,tmp_S_k_p_);
tmp_X_(1+nM) = tmp_TS/sqrt(tmp_TT*tmp_SS);
end;%for nM=0:n_M-1;
%%%%%%%%;
% Naively, we expect that (oriented) image T and template S should align as: ;
% X = <T,S> / ( |T| * |S| ) ;
% whereas: ;
% R = T - S ;
% |R|^2 = |T|^2 + |S|^2 - 2<T,S> ;
% so: ;
% |T|^2 + |S|^2 - 2*X|T||S| = |R|^2. ;
%%%%%%%%;
if (verbose); disp(sprintf(' %% niteration_use %d/%d: check X against R: %0.16f',niteration_use,n_iteration_use,fnorm(tmp_T_k_p_l2_ + tmp_S_k_p_l2_ - 2*real(tmp_X_).*sqrt(tmp_T_k_p_l2_).*sqrt(tmp_S_k_p_l2_) - tmp_R_k_p_l2_)/fnorm(tmp_R_k_p_wnM__))); end;
%%%%%%%%;
R_k_p_l2_Mi__(:,1+niteration_use) = tmp_R_k_p_l2_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for niteration_use=0:n_iteration_use-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% R_k_p_l2_Mi__ time %0.2fs',tmp_t)); end;

%%%%%%%%;
tmp_t = tic();
tmp_ij_min_threshold_vs_X_TM_avg_Mi__ = zeros(n_M,n_iteration_use);
for niteration_use=0:n_iteration_use-1;
tmp_image_X_value_ = image_X_value_Mi__(:,1+niteration_use);
tmp_R_k_p_l2_ = R_k_p_l2_Mi__(:,1+niteration_use);
[~,tmp_ij_X_value_sort_] = sort(tmp_image_X_value_,'ascend'); [~,tmp_ij_X_value_sort_] = sort(tmp_ij_X_value_sort_,'ascend');
[~,tmp_ij_R_k_p_l2_sort_] = sort(tmp_R_k_p_l2_,'descend'); [~,tmp_ij_R_k_p_l2_sort_] = sort(tmp_ij_R_k_p_l2_sort_,'ascend');
tmp_ij_min_sort_ = min(tmp_ij_X_value_sort_,tmp_ij_R_k_p_l2_sort_);
tmp_ij_max_sort_ = max(tmp_ij_X_value_sort_,tmp_ij_R_k_p_l2_sort_);
%%%%;
[~,tmp_ij_X_TM_sort_] = sort(X_TM_,'ascend'); [~,tmp_ij_X_TM_sort_] = sort(tmp_ij_X_TM_sort_,'ascend');
%%%%;
tmp_ij_X_value_vs_ij_X_TM__ = sparse(tmp_ij_X_value_sort_,tmp_ij_X_TM_sort_,1,n_M,n_M);
tmp_ij_X_value_threshold_vs_X_TM_avg_ = sum(bsxfun(@times,cumsum(tmp_ij_X_value_vs_ij_X_TM__,1,'reverse'),reshape(1:n_M,[1,n_M])),2)./sum(cumsum(tmp_ij_X_value_vs_ij_X_TM__,1,'reverse'),2);
tmp_ij_R_k_p_l2_vs_ij_X_TM__ = sparse(tmp_ij_R_k_p_l2_sort_,tmp_ij_X_TM_sort_,1,n_M,n_M);
tmp_ij_R_k_p_l2_threshold_vs_X_TM_avg_ = sum(bsxfun(@times,cumsum(tmp_ij_R_k_p_l2_vs_ij_X_TM__,1,'reverse'),reshape(1:n_M,[1,n_M])),2)./sum(cumsum(tmp_ij_R_k_p_l2_vs_ij_X_TM__,1,'reverse'),2);
tmp_ij_min_vs_ij_X_TM__ = sparse(tmp_ij_min_sort_,tmp_ij_X_TM_sort_,1,n_M,n_M);
tmp_ij_min_threshold_vs_X_TM_avg_ = sum(bsxfun(@times,cumsum(tmp_ij_min_vs_ij_X_TM__,1,'reverse'),reshape(1:n_M,[1,n_M])),2)./sum(cumsum(tmp_ij_min_vs_ij_X_TM__,1,'reverse'),2);
tmp_ij_min_threshold_vs_X_TM_avg_Mi__(:,1+niteration_use) = tmp_ij_min_threshold_vs_X_TM_avg_;
end;%for niteration_use=0:n_iteration_use-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% tmp_ij_min_threshold_vs_X_TM_avg_Mi__ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_image_X_value_ = image_X_value_Mi__(:,n_iteration_use);
tmp_R_k_p_l2_ = R_k_p_l2_Mi__(:,n_iteration_use);
[~,tmp_ij_X_value_sort_] = sort(tmp_image_X_value_,'ascend'); [~,tmp_ij_X_value_sort_] = sort(tmp_ij_X_value_sort_,'ascend');
[~,tmp_ij_R_k_p_l2_sort_] = sort(tmp_R_k_p_l2_,'descend'); [~,tmp_ij_R_k_p_l2_sort_] = sort(tmp_ij_R_k_p_l2_sort_,'ascend');
tmp_ij_min_sort_ = min(tmp_ij_X_value_sort_,tmp_ij_R_k_p_l2_sort_);
tmp_ij_max_sort_ = max(tmp_ij_X_value_sort_,tmp_ij_R_k_p_l2_sort_);
%%%%;
[~,tmp_ij_X_TM_sort_] = sort(X_TM_,'ascend'); [~,tmp_ij_X_TM_sort_] = sort(tmp_ij_X_TM_sort_,'ascend');
%%%%;
tmp_ij_X_value_vs_ij_X_TM__ = sparse(tmp_ij_X_value_sort_,tmp_ij_X_TM_sort_,1,n_M,n_M);
tmp_ij_X_value_threshold_vs_X_TM_avg_ = sum(bsxfun(@times,cumsum(tmp_ij_X_value_vs_ij_X_TM__,1,'reverse'),reshape(1:n_M,[1,n_M])),2)./sum(cumsum(tmp_ij_X_value_vs_ij_X_TM__,1,'reverse'),2);
tmp_ij_R_k_p_l2_vs_ij_X_TM__ = sparse(tmp_ij_R_k_p_l2_sort_,tmp_ij_X_TM_sort_,1,n_M,n_M);
tmp_ij_R_k_p_l2_threshold_vs_X_TM_avg_ = sum(bsxfun(@times,cumsum(tmp_ij_R_k_p_l2_vs_ij_X_TM__,1,'reverse'),reshape(1:n_M,[1,n_M])),2)./sum(cumsum(tmp_ij_R_k_p_l2_vs_ij_X_TM__,1,'reverse'),2);
tmp_ij_min_vs_ij_X_TM__ = sparse(tmp_ij_min_sort_,tmp_ij_X_TM_sort_,1,n_M,n_M);
tmp_ij_min_threshold_vs_X_TM_avg_ = sum(bsxfun(@times,cumsum(tmp_ij_min_vs_ij_X_TM__,1,'reverse'),reshape(1:n_M,[1,n_M])),2)./sum(cumsum(tmp_ij_min_vs_ij_X_TM__,1,'reverse'),2);
%%%%%%%%;

%%%%%%%%;
% save to mat-file. ;
%%%%%%%%;
tmp_fname_mat = sprintf('%s.mat',X_fname_compare_image_rank_pre);
save(tmp_fname_mat ...
     ,'n_M' ...
     ,'n_iteration_use' ...
     ,'image_X_value_Mi__' ...
     ,'R_k_p_l2_Mi__' ...
     ,'tmp_image_X_value_' ...
     ,'tmp_R_k_p_l2_' ...
     ,'X_TM_' ...
     ,'tmp_ij_min_threshold_vs_X_TM_avg_Mi__' ...
     ,'tmp_ij_min_threshold_vs_X_TM_avg_' ...
     );
%%%%%%%%;
figure(1);clf;figbig;figbeach();ns=0;
markersize_use = 6;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
c_80s__ = colormap_80s(); n_c_80s = size(c_80s__,1);
X_value_lim_ = prctile(tmp_image_X_value_,[ 1,99]);
X_value_lim_ = mean(X_value_lim_) + 1.25*0.5*diff(X_value_lim_)*[-1,+1];
R_k_p_l2_lim_ = prctile(tmp_R_k_p_l2_,[ 1,99]);
R_k_p_l2_lim_ = mean(R_k_p_l2_lim_) + 1.25*0.5*diff(R_k_p_l2_lim_)*[-1,+1];
X_TM_lim_ = prctile(X_TM_,[ 1,99]);
X_TM_lim_ = mean(X_TM_lim_) + 1.25*0.5*diff(X_TM_lim_)*[-1,+1];
%%%%;
subplot(2,3,1+ns);ns=ns+1;
hold on;
for nM=0:n_M-1;
nc = max(0,min(n_c_80s-1,floor(n_c_80s*tmp_ij_X_TM_sort_(1+nM)/n_M)));
plot(tmp_image_X_value_(1+nM),tmp_R_k_p_l2_(1+nM),'o','MarkerSize',markersize_use,'MarkerFaceColor',c_80s__(1+nc,:),'MarkerEdgeColor',0.65*[1,1,1]);
end;%for nM=0:n_M-1;
hold off;
xlim(X_value_lim_);
ylim(R_k_p_l2_lim_);
xlabel('image_X_value_','Interpreter','none');
ylabel('tmp_R_k_p_l2_','Interpreter','none');
title('colored by true-rank','Interpreter','none');
%%%%;
subplot(2,3,1+ns);ns=ns+2;
hold on;
for nM=0:n_M-1;
nc = max(0,min(n_c_beach-1,floor(n_c_beach*tmp_ij_min_sort_(1+nM)/n_M)));
plot(tmp_image_X_value_(1+nM),tmp_R_k_p_l2_(1+nM),'o','MarkerSize',markersize_use,'MarkerFaceColor',c_beach__(1+nc,:),'MarkerEdgeColor',0.65*[1,1,1]);
end;%for nM=0:n_M-1;
hold off;
xlim(X_value_lim_);
ylim(R_k_p_l2_lim_);
xlabel('image_X_value_','Interpreter','none');
ylabel('tmp_R_k_p_l2_','Interpreter','none');
title('colored by estimated-minimum-rank','Interpreter','none');
%%%%;
subplot(2,3,1+ns);ns=ns+1;
hold on;
for nM=0:n_M-1;
nc = max(0,min(n_c_beach-1,floor(n_c_beach*tmp_ij_min_sort_(1+nM)/n_M)));
plot(tmp_image_X_value_(1+nM),X_TM_(1+nM),'o','MarkerSize',markersize_use,'MarkerFaceColor',c_beach__(1+nc,:),'MarkerEdgeColor',0.65*[1,1,1]);
end;%for nM=0:n_M-1;
hold off;
xlim(X_value_lim_);
ylim(X_TM_lim_);
xlabel('image_X_value_','Interpreter','none');
ylabel('X_TM_','Interpreter','none');
title('colored by estimated-minimum-rank','Interpreter','none');
%%%%;
subplot(2,3,1+ns);ns=ns+1;
hold on;
for nM=0:n_M-1;
nc = max(0,min(n_c_beach-1,floor(n_c_beach*tmp_ij_min_sort_(1+nM)/n_M)));
plot(tmp_R_k_p_l2_(1+nM),X_TM_(1+nM),'o','MarkerSize',markersize_use,'MarkerFaceColor',c_beach__(1+nc,:),'MarkerEdgeColor',0.65*[1,1,1]);
end;%for nM=0:n_M-1;
hold off;
xlim(R_k_p_l2_lim_);
ylim(X_TM_lim_);
xlabel('tmp_R_k_p_l2_','Interpreter','none');
ylabel('X_TM_','Interpreter','none');
title('colored by estimated-minimum-rank','Interpreter','none');
%%%%;
subplot(2,3,3);ns=ns+1;
hold on;
plot(1:n_M,tmp_ij_X_value_threshold_vs_X_TM_avg_,'r.-');
plot(1:n_M,tmp_ij_R_k_p_l2_threshold_vs_X_TM_avg_,'g.-');
plot(1:n_M,tmp_ij_min_threshold_vs_X_TM_avg_,'k.-');
plot(1:n_M,n_M/2*ones(1,n_M),'k-');
hold off;
xlim([0,n_M]);
ylim([n_M/2,n_M]);
grid on;
xlabel('threshold for estimated-minimum-rank');
ylabel('average true-rank exceeding threshold');
legend({'X value','R k p l2','estimated-minimum'},'Location','NorthWest');
title('average-true-rank for final-iteration','Interpreter','none');
%%%%;
subplot(2,3,6);ns=ns+1;
figbeach();
imagesc(transpose(tmp_ij_min_threshold_vs_X_TM_avg_Mi__),[n_M/2,n_M]);
colorbar;
ylabel('iteration number (ampm)');
xlabel('threshold for estimated-minimum-rank');
title('average true-rank exceeding threshold (across iterations)','Interpreter','none');
%%%%%%%%;
sgtitle(X_fname_mat,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ( exist(X_fname_mat,'file'));

if (verbose); disp(sprintf(' %% [finished ampmut_compare_image_rank_1]')); end;




