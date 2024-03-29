
%%%%%%%%%%%%%%%%;
% Now calculate the bayesian likelihood associated with true, pole and pol3, as applied to the 'expanded' set of images. ;
%%%%%%%%%%%%%%%%;
S_p123_k_q_wkS__ = S_k_q__;
%%%%%%%%;
% First estimate distances between templates and images. ;
%%%%%%%%;
ncluster = 0; %<-- pick one of the CTF clusters. ;
UX_kn__ = UX_knc___(:,:,1+ncluster);
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
T_p123_orig_k_p_wkS__ = reshape(bsxfun(@times,reshape(S_p123_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_p123_orig_k_q_wkS__ = reshape(bsxfun(@times,reshape(S_p123_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_p123_garb_k_p_wkS__ = reshape(bsxfun(@times,reshape(p123_pole_from_R_S_reco_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_p123_garb_k_q_wkS__ = reshape(bsxfun(@times,reshape(p123_pole_from_R_S_reco_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_p123_pole_k_p_wkS__ = reshape(bsxfun(@times,reshape(S_p123_pole_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_p123_pole_k_q_wkS__ = reshape(bsxfun(@times,reshape(S_p123_pole_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_p231_pole_k_p_wkS__ = reshape(bsxfun(@times,reshape(S_p231_pole_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_p231_pole_k_q_wkS__ = reshape(bsxfun(@times,reshape(S_p231_pole_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_p312_pole_k_p_wkS__ = reshape(bsxfun(@times,reshape(S_p312_pole_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_p312_pole_k_q_wkS__ = reshape(bsxfun(@times,reshape(S_p312_pole_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
tmp_t = tic();
svd_VUXT_p123_orig_lwnS____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_S,T_p123_orig_k_q_wkS__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXT_p123_orig_lwnS____: %0.3fs',tmp_t)); end;
tmp_t = tic();
svd_VUXT_p123_garb_lwnS____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_S,T_p123_garb_k_q_wkS__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXT_p123_garb_lwnS____: %0.3fs',tmp_t)); end;
tmp_t = tic();
svd_VUXT_p123_pole_lwnS____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_S,T_p123_pole_k_q_wkS__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXT_p123_pole_lwnS____: %0.3fs',tmp_t)); end;
tmp_t = tic();
svd_VUXT_p231_pole_lwnS____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_S,T_p231_pole_k_q_wkS__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXT_p231_pole_lwnS____: %0.3fs',tmp_t)); end;
tmp_t = tic();
svd_VUXT_p312_pole_lwnS____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_S,T_p312_pole_k_q_wkS__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXT_p312_pole_lwnS____: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[UX_T_p123_orig_k_q_wnS___,UX_T_p123_orig_k_p_wnS___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,n_UX_rank,n_S,svd_VUXT_p123_orig_lwnS____,zeros(n_S,1),zeros(n_S,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_p123_orig_k_q_wnS___: %0.6fs',tmp_t)); end;
UX_T_p123_orig_k_p_wnS__ = reshape(UX_T_p123_orig_k_p_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
UX_T_p123_orig_k_q_wnS__ = reshape(UX_T_p123_orig_k_q_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
%%%%%%%%;
tmp_t = tic();
[UX_T_p123_garb_k_q_wnS___,UX_T_p123_garb_k_p_wnS___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,n_UX_rank,n_S,svd_VUXT_p123_garb_lwnS____,zeros(n_S,1),zeros(n_S,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_p123_garb_k_q_wnS___: %0.6fs',tmp_t)); end;
UX_T_p123_garb_k_p_wnS__ = reshape(UX_T_p123_garb_k_p_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
UX_T_p123_garb_k_q_wnS__ = reshape(UX_T_p123_garb_k_q_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
%%%%%%%%;
tmp_t = tic();
[UX_T_p123_pole_k_q_wnS___,UX_T_p123_pole_k_p_wnS___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,n_UX_rank,n_S,svd_VUXT_p123_pole_lwnS____,zeros(n_S,1),zeros(n_S,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_p123_pole_k_q_wnS___: %0.6fs',tmp_t)); end;
UX_T_p123_pole_k_p_wnS__ = reshape(UX_T_p123_pole_k_p_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
UX_T_p123_pole_k_q_wnS__ = reshape(UX_T_p123_pole_k_q_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
%%%%%%%%;
tmp_t = tic();
[UX_T_p231_pole_k_q_wnS___,UX_T_p231_pole_k_p_wnS___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,n_UX_rank,n_S,svd_VUXT_p231_pole_lwnS____,zeros(n_S,1),zeros(n_S,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_p231_pole_k_q_wnS___: %0.6fs',tmp_t)); end;
UX_T_p231_pole_k_p_wnS__ = reshape(UX_T_p231_pole_k_p_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
UX_T_p231_pole_k_q_wnS__ = reshape(UX_T_p231_pole_k_q_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
%%%%%%%%;
tmp_t = tic();
[UX_T_p312_pole_k_q_wnS___,UX_T_p312_pole_k_p_wnS___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,n_UX_rank,n_S,svd_VUXT_p312_pole_lwnS____,zeros(n_S,1),zeros(n_S,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_p312_pole_k_q_wnS___: %0.6fs',tmp_t)); end;
UX_T_p312_pole_k_p_wnS__ = reshape(UX_T_p312_pole_k_p_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
UX_T_p312_pole_k_q_wnS__ = reshape(UX_T_p312_pole_k_q_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
%%%%%%%%;
% Now calculate norms of the templates. ;
%%%%%%%%;
tmp_p123_orig_TT_S_ = zeros(n_S,1);
tmp_p123_garb_TT_S_ = zeros(n_S,1);
tmp_p123_pole_TT_S_ = zeros(n_S,1);
tmp_p231_pole_TT_S_ = zeros(n_S,1);
tmp_p312_pole_TT_S_ = zeros(n_S,1);
for nS=0:n_S-1;
tmp_p123_orig_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_p123_orig_k_p_wkS__(:,1+nS),T_p123_orig_k_p_wkS__(:,1+nS))/(2*pi);
tmp_p123_orig_TT_S_(1+nS) = tmp_p123_orig_TT;
tmp_p123_garb_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_p123_garb_k_p_wkS__(:,1+nS),T_p123_garb_k_p_wkS__(:,1+nS))/(2*pi);
tmp_p123_garb_TT_S_(1+nS) = tmp_p123_garb_TT;
tmp_p123_pole_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_p123_pole_k_p_wkS__(:,1+nS),T_p123_pole_k_p_wkS__(:,1+nS))/(2*pi);
tmp_p123_pole_TT_S_(1+nS) = tmp_p123_pole_TT;
tmp_p231_pole_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_p231_pole_k_p_wkS__(:,1+nS),T_p231_pole_k_p_wkS__(:,1+nS))/(2*pi);
tmp_p231_pole_TT_S_(1+nS) = tmp_p231_pole_TT;
tmp_p312_pole_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_p312_pole_k_p_wkS__(:,1+nS),T_p312_pole_k_p_wkS__(:,1+nS))/(2*pi);
tmp_p312_pole_TT_S_(1+nS) = tmp_p312_pole_TT;
end;%for nS=0:n_S-1;
UX_T_p123_orig_l2_S_ = tmp_p123_orig_TT_S_;
UX_T_p123_garb_l2_S_ = tmp_p123_garb_TT_S_;
UX_T_p123_pole_l2_S_ = tmp_p123_pole_TT_S_;
UX_T_p231_pole_l2_S_ = tmp_p231_pole_TT_S_;
UX_T_p312_pole_l2_S_ = tmp_p312_pole_TT_S_;
%%%%%%%%;
% Now generate svd_VUXM2_lwnM2____. ;
%%%%%%%%;
tmp_t = tic();
svd_VUXM2_lwnM2____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M2,M2_k_q__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXM2_lwnM2____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the images. ;
%%%%%%%%;
tmp_t = tic();
UX_M2_l2_dM2__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M2,n_UX_rank,svd_VUXM2_lwnM2____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M2_l2_dM2__: %0.3fs',tmp_t)); end;
if (verbose); disp(sprintf(' %% average l2-norm of images: %0.16f',mean(UX_M2_l2_dM2__(:))/(pi*k_p_r_max^2))); end;
tmp_M2M2_M2_ = zeros(n_M2,1);
for nM2=0:n_M2-1;
tmp_M2M2 = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,M2_k_p__(:,1+nM2),M2_k_p__(:,1+nM2))/(2*pi);
tmp_M2M2_M2_(1+nM2) = tmp_M2M2;
end;%for nM2=0:n_M2-1;
tmp_index = efind((FTK.delta_x_.^2 + FTK.delta_y_.^2)<1e-6);
UX_M2_l2_M2_ = transpose(UX_M2_l2_dM2__(1+tmp_index,:));
if (verbose); disp(sprintf(' %% tmp_M2M2_M2_ vs UX_M2_l2_M2_: %0.16f',fnorm(tmp_M2M2_M2_ - UX_M2_l2_M2_)/fnorm(tmp_M2M2_M2_))); end;
flag_plot=0;
if flag_plot;
tmp_index = efind((FTK.delta_x_.^2 + FTK.delta_y_.^2)<1e-6);
subplot(1,2,1); hold on; 
plot(0:n_M2-1,UX_M2_l2_M2_/(pi*k_p_r_max^2),'rx'); xlabel('nM2'); ylabel('l2');
plot(0:n_M2-1,tmp_M2M2_M2_/(pi*k_p_r_max^2),'bo'); xlabel('nM2'); ylabel('l2');
hold off;
subplot(1,2,2); plot(UX_M2_l2_M2_/(pi*k_p_r_max^2),tmp_M2M2_M2_/(pi*k_p_r_max^2),'g.'); xlabel('l2_A'); ylabel('l2_B');
end;%if flag_plot;
%%%%%%%%;
tmp_t = tic();
[UX_M2_k_q_wnM2___,UX_M2_k_p_wnM2___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,n_UX_rank,n_M2,svd_VUXM2_lwnM2____,zeros(n_M2,1),zeros(n_M2,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M2_k_q_wnM2___: %0.6fs',tmp_t)); end;
UX_M2_k_p_wnM2__ = reshape(UX_M2_k_p_wnM2___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_M2]);
UX_M2_k_q_wnM2__ = reshape(UX_M2_k_q_wnM2___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_M2]);
%%%%%%%%;
% Calculate ampmh_X_wSM2___8. ;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_optimize_over_gamma_z = 1;
parameter.flag_compute_I_value = 0;
parameter.tolerance_master = tolerance_master;
parameter.pm_n_UX_rank_use = n_UX_rank;
tmp_t = tic();
[ ...
 parameter ...
,X_p123_orig_TM2__ ...
,delta_x_p123_orig_TM2__ ...
,delta_y_p123_orig_TM2__ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK ...
,n_w_max ...
,n_UX_rank ...
,n_S ...
,UX_T_p123_orig_k_q_wnS__ ...
,UX_T_p123_orig_l2_S_ ...
,n_M2 ...
,svd_VUXM2_lwnM2____ ...
,UX_M2_l2_dM2__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_p123_orig_TM2__: %0.3fs',tmp_t)); end;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_optimize_over_gamma_z = 1;
parameter.flag_compute_I_value = 0;
parameter.tolerance_master = tolerance_master;
parameter.pm_n_UX_rank_use = n_UX_rank;
tmp_t = tic();
[ ...
 parameter ...
,X_p123_garb_TM2__ ...
,delta_x_p123_garb_TM2__ ...
,delta_y_p123_garb_TM2__ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK ...
,n_w_max ...
,n_UX_rank ...
,n_S ...
,UX_T_p123_garb_k_q_wnS__ ...
,UX_T_p123_garb_l2_S_ ...
,n_M2 ...
,svd_VUXM2_lwnM2____ ...
,UX_M2_l2_dM2__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_p123_garb_TM2__: %0.3fs',tmp_t)); end;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_optimize_over_gamma_z = 1;
parameter.flag_compute_I_value = 0;
parameter.tolerance_master = tolerance_master;
parameter.pm_n_UX_rank_use = n_UX_rank;
tmp_t = tic();
[ ...
 parameter ...
,X_p123_pole_TM2__ ...
,delta_x_p123_pole_TM2__ ...
,delta_y_p123_pole_TM2__ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK ...
,n_w_max ...
,n_UX_rank ...
,n_S ...
,UX_T_p123_pole_k_q_wnS__ ...
,UX_T_p123_pole_l2_S_ ...
,n_M2 ...
,svd_VUXM2_lwnM2____ ...
,UX_M2_l2_dM2__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_p123_pole_TM2__: %0.3fs',tmp_t)); end;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_optimize_over_gamma_z = 1;
parameter.flag_compute_I_value = 0;
parameter.tolerance_master = tolerance_master;
parameter.pm_n_UX_rank_use = n_UX_rank;
tmp_t = tic();
[ ...
 parameter ...
,X_p231_pole_TM2__ ...
,delta_x_p231_pole_TM2__ ...
,delta_y_p231_pole_TM2__ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK ...
,n_w_max ...
,n_UX_rank ...
,n_S ...
,UX_T_p231_pole_k_q_wnS__ ...
,UX_T_p231_pole_l2_S_ ...
,n_M2 ...
,svd_VUXM2_lwnM2____ ...
,UX_M2_l2_dM2__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_p231_pole_TM__: %0.3fs',tmp_t)); end;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_optimize_over_gamma_z = 1;
parameter.flag_compute_I_value = 0;
parameter.tolerance_master = tolerance_master;
parameter.pm_n_UX_rank_use = n_UX_rank;
tmp_t = tic();
[ ...
 parameter ...
,X_p312_pole_TM2__ ...
,delta_x_p312_pole_TM2__ ...
,delta_y_p312_pole_TM2__ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK ...
,n_w_max ...
,n_UX_rank ...
,n_S ...
,UX_T_p312_pole_k_q_wnS__ ...
,UX_T_p312_pole_l2_S_ ...
,n_M2 ...
,svd_VUXM2_lwnM2____ ...
,UX_M2_l2_dM2__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_p312_pole_TM__: %0.3fs',tmp_t)); end;
%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% calculate log-unlikelihood ratio for expanded set of n_M2 images. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%;
sigma_bayesian_ = transpose([0,2.^[-12:0.125:0]]); n_sigma_bayesian = numel(sigma_bayesian_);
%%%%%%%%;
[ ...
 ~ ...
,ssnll2_uni_tru_s_ ...
,ssnll2_uni_tru_M2s__ ...
,inten2_uni_tru_vs__ ...
] = ...
ssnll_from_X_1( ...
 [] ...
,1 ... 
,[n_S] ...
,n_M2 ...
,{X_p123_orig_TM2__} ...
,{UX_T_p123_orig_l2_S_} ...
,UX_M2_l2_M2_ ...
,{viewing_weight_all_} ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
);
%%%%;
[ ...
 ~ ...
,ssnll2_uni_tr0_s_ ...
,ssnll2_uni_tr0_M2s__ ...
,inten2_uni_tr0_vs__ ...
] = ...
ssnll_from_X_1( ...
 [] ...
,2 ... 
,[n_S,n_S] ...
,n_M2 ...
,{X_p123_orig_TM2__,zeros(size(X_p123_orig_TM2__))} ...
,{UX_T_p123_orig_l2_S_,ones(size(UX_T_p123_orig_l2_S_))} ...
,UX_M2_l2_M2_ ...
,{viewing_weight_all_,viewing_weight_all_} ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
,repmat(log(inten2_uni_tru_vs__),[2,1]) ...
);
%%%%;
[ ...
 ~ ...
,ssnll2_uni_tr1_s_ ...
,ssnll2_uni_tr1_M2s__ ...
,inten2_uni_tr1_vs__ ...
] = ...
ssnll_from_X_1( ...
 [] ...
,2 ... 
,[n_S,n_S] ...
,n_M2 ...
,{X_p123_orig_TM2__,X_p123_garb_TM2__} ...
,{UX_T_p123_orig_l2_S_,UX_T_p123_garb_l2_S_} ...
,UX_M2_l2_M2_ ...
,{viewing_weight_all_,viewing_weight_all_} ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
,repmat(log(inten2_uni_tru_vs__),[2,1]) ...
);
%%%%;
[ ...
 ~ ...
,ssnll2_uni_ppp_s_ ...
,ssnll2_uni_ppp_M2s__ ...
,inten2_uni_ppp_vs__ ...
] = ...
ssnll_from_X_1( ...
 [] ...
,3 ... 
,[n_S,n_S,n_S] ...
,n_M2 ...
,{X_p123_pole_TM2__,X_p231_pole_TM2__,X_p312_pole_TM2__} ...
,{UX_T_p123_pole_l2_S_,UX_T_p231_pole_l2_S_,UX_T_p312_pole_l2_S_} ...
,UX_M2_l2_M2_ ...
,{viewing_weight_all_,viewing_weight_all_,viewing_weight_all_} ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
,repmat(log(inten2_uni_tru_vs__),[3,1]) ...
);
%%%%;
ssnll2r_uni_tru_vs_ppp_M2s__ = ssnll2_uni_tru_M2s__ - ssnll2_uni_ppp_M2s__; %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
ssnll2r_uni_tru_vs_ppp_s_ = sum(ssnll2r_uni_tru_vs_ppp_M2s__,1); %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
ssnll2r_uni_tr0_vs_ppp_M2s__ = ssnll2_uni_tr0_M2s__ - ssnll2_uni_ppp_M2s__; %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
ssnll2r_uni_tr0_vs_ppp_s_ = sum(ssnll2r_uni_tr0_vs_ppp_M2s__,1); %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
ssnll2r_uni_tr1_vs_ppp_M2s__ = ssnll2_uni_tr1_M2s__ - ssnll2_uni_ppp_M2s__; %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
ssnll2r_uni_tr1_vs_ppp_s_ = sum(ssnll2r_uni_tr1_vs_ppp_M2s__,1); %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
%%%%%%%%;

%%%%%%%%;
[ ...
 ~ ...
,ssnll2_emp_tru_s_ ...
,ssnll2_emp_tru_M2s__ ...
,inten2_emp_tru_vs__ ...
] = ...
ssnll_from_X_1( ...
 [] ...
,1 ... 
,[n_S] ...
,n_M2 ...
,{X_p123_orig_TM2__} ...
,{UX_T_p123_orig_l2_S_} ...
,UX_M2_l2_M2_ ...
,{viewing_weight_p123_} ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
);
%%%%;
[ ...
 ~ ...
,ssnll2_emp_tr0_s_ ...
,ssnll2_emp_tr0_M2s__ ...
,inten2_emp_tr0_vs__ ...
] = ...
ssnll_from_X_1( ...
 [] ...
,2 ... 
,[n_S,n_S] ...
,n_M2 ...
,{X_p123_orig_TM2__,zeros(size(X_p123_orig_TM2__))} ...
,{UX_T_p123_orig_l2_S_,ones(size(UX_T_p123_orig_l2_S_))} ...
,UX_M2_l2_M2_ ...
,{viewing_weight_p123_,viewing_weight_p123_} ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
,repmat(log(inten2_emp_tru_vs__),[2,1]) ...
);
%%%%;
[ ...
 ~ ...
,ssnll2_emp_tr1_s_ ...
,ssnll2_emp_tr1_M2s__ ...
,inten2_emp_tr1_vs__ ...
] = ...
ssnll_from_X_1( ...
 [] ...
,2 ... 
,[n_S,n_S] ...
,n_M2 ...
,{X_p123_orig_TM2__,X_p123_garb_TM2__} ...
,{UX_T_p123_orig_l2_S_,UX_T_p123_garb_l2_S_} ...
,UX_M2_l2_M2_ ...
,{viewing_weight_p123_,viewing_weight_p123_} ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
,repmat(log(inten2_emp_tru_vs__),[2,1]) ...
);
%%%%;
[ ...
 ~ ...
,ssnll2_emp_ppp_s_ ...
,ssnll2_emp_ppp_M2s__ ...
,inten2_emp_ppp_vs__ ...
] = ...
ssnll_from_X_1( ...
 [] ...
,3 ... 
,[n_S,n_S,n_S] ...
,n_M2 ...
,{X_p123_pole_TM2__,X_p231_pole_TM2__,X_p312_pole_TM2__} ...
,{UX_T_p123_pole_l2_S_,UX_T_p231_pole_l2_S_,UX_T_p312_pole_l2_S_} ...
,UX_M2_l2_M2_ ...
,{viewing_weight_p123_,viewing_weight_p123_,viewing_weight_p123_} ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
,repmat(log(inten2_emp_tru_vs__),[3,1]) ...
);
%%%%;
ssnll2r_emp_tru_vs_ppp_M2s__ = ssnll2_emp_tru_M2s__ - ssnll2_emp_ppp_M2s__; %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
ssnll2r_emp_tru_vs_ppp_s_ = sum(ssnll2r_emp_tru_vs_ppp_M2s__,1); %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
ssnll2r_emp_tr0_vs_ppp_M2s__ = ssnll2_emp_tr0_M2s__ - ssnll2_emp_ppp_M2s__; %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
ssnll2r_emp_tr0_vs_ppp_s_ = sum(ssnll2r_emp_tr0_vs_ppp_M2s__,1); %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
ssnll2r_emp_tr1_vs_ppp_M2s__ = ssnll2_emp_tr1_M2s__ - ssnll2_emp_ppp_M2s__; %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
ssnll2r_emp_tr1_vs_ppp_s_ = sum(ssnll2r_emp_tr1_vs_ppp_M2s__,1); %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
%%%%%%%%;

%%%%%%%%;
save(fname_ssnll_pole_mat ...
,'viewing_weight_all_' ...
,'viewing_weight_p123_' ...
,'n_M' ...
,'n_R' ...
,'n_M2' ...
,'sigma_bayesian_' ...
,'n_sigma_bayesian' ...
,'ssnll2_uni_tru_s_' ...
,'ssnll2_uni_tru_M2s__' ...
,'inten2_uni_tru_vs__' ...
,'ssnll2_uni_tr0_s_' ...
,'ssnll2_uni_tr0_M2s__' ...
,'inten2_uni_tr0_vs__' ...
,'ssnll2_uni_tr1_s_' ...
,'ssnll2_uni_tr1_M2s__' ...
,'inten2_uni_tr1_vs__' ...
,'ssnll2_uni_ppp_s_' ...
,'ssnll2_uni_ppp_M2s__' ...
,'inten2_uni_ppp_vs__' ...
,'ssnll2r_uni_tru_vs_ppp_M2s__' ...
,'ssnll2r_uni_tru_vs_ppp_s_' ...
,'ssnll2r_uni_tr0_vs_ppp_M2s__' ...
,'ssnll2r_uni_tr0_vs_ppp_s_' ...
,'ssnll2r_uni_tr1_vs_ppp_M2s__' ...
,'ssnll2r_uni_tr1_vs_ppp_s_' ...
,'ssnll2_emp_tru_s_' ...
,'ssnll2_emp_tru_M2s__' ...
,'inten2_emp_tru_vs__' ...
,'ssnll2_emp_tr0_s_' ...
,'ssnll2_emp_tr0_M2s__' ...
,'inten2_emp_tr0_vs__' ...
,'ssnll2_emp_tr1_s_' ...
,'ssnll2_emp_tr1_M2s__' ...
,'inten2_emp_tr1_vs__' ...
,'ssnll2_emp_ppp_s_' ...
,'ssnll2_emp_ppp_M2s__' ...
,'inten2_emp_ppp_vs__' ...
,'ssnll2r_emp_tru_vs_ppp_M2s__' ...
,'ssnll2r_emp_tru_vs_ppp_s_' ...
,'ssnll2r_emp_tr0_vs_ppp_M2s__' ...
,'ssnll2r_emp_tr0_vs_ppp_s_' ...
,'ssnll2r_emp_tr1_vs_ppp_M2s__' ...
,'ssnll2r_emp_tr1_vs_ppp_s_' ...
);
close_fname_tmp(fname_ssnll_pole_pre);
%%%%%%%%%%%%%%%%;
