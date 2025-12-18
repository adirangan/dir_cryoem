flag_verbose=1; flag_disp=1;nf=0;
flag_speed_vs_error = 0;

disp(sprintf(' %% testing tfpmh_UX_T_M_l2_dM__1: flag_speed_vs_error %d',flag_speed_vs_error));
n_1 = 1; n_2 = 2;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi); str_L = 'L';
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 0*flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_L ...
);
%%%%;
n_w_int = 1;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
n_w_max = n_w_int*2*(l_max_max+1); n_w_0in_ = n_w_max*ones(n_k_p_r,1);
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_p_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 0*flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_0in_ ...
,weight_3d_k_p_r_ ...
);
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%;
if flag_speed_vs_error==0; n_M = 17; end;
if flag_speed_vs_error==1; n_M = 1024; end;
n_source = 3;
delta_2sM___ = permute(reshape(linspace(-0.05,+0.05,n_2*n_source*n_M),[n_M,n_2,n_source]),1+[1,2,0]);
M_k_p_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
for nsource=0:n_source-1;
delta_0 = delta_2sM___(1+0,1+nsource,1+nM); delta_1 = delta_2sM___(1+1,1+nsource,1+nM);
M_k_p_wkM__(:,1+nM) = M_k_p_wkM__(:,1+nM) + 1*exp(+i*2*pi*(k_c_0_wk_*delta_0 + k_c_1_wk_*delta_1));
end;%for nsource=0:n_source-1;
end;%for nM=0:n_M-1;
%%%%;
n_CTF = 1;
CTF_k_p_r_k_ = 0.5 + 1.0*cos(2*pi*k_p_r_);
%%%%;
[ ...
 X_kk__ ...
,X_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_1( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
);
if flag_speed_vs_error==0; pm_n_UX_rank = min(2,n_k_p_r-1); end; %<-- should be accurate even when radial-compression is lossy. ;
if flag_speed_vs_error==1; pm_n_UX_rank = 17; end;
[tmp_UX_kn__,tmp_SX_k__,tmp_VX_kn__] = svds(X_kk__,pm_n_UX_rank);
pm_UX_kn__(:,1:pm_n_UX_rank) = tmp_UX_kn__(:,1:pm_n_UX_rank);
pm_X_weight_r_ = sqrt(weight_2d_k_p_r_);
%%%%;
delta_r_max = 0.5/max(1e-12,k_p_r_max); 
if flag_speed_vs_error==0; svd_eps = 1e-9; n_delta_v_requested =   9; end;
if flag_speed_vs_error==1; svd_eps = 1e-5; n_delta_v_requested = 512; end;

%FTK = tfh_FTK_2(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
parameter = struct('type','parameter');
parameter.flag_verbose = max(0,flag_verbose-1);
parameter.delta_r_max = delta_r_max;
parameter.svd_eps = svd_eps;
parameter.n_delta_v_requested = n_delta_v_requested;
[parameter,FTK] = tfh_FTK_3(parameter,n_k_p_r,k_p_r_,k_p_r_max);
n_delta_v = FTK.n_delta_v;
n_svd_l = FTK.n_svd_l;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% n_k_p_r %d n_w_max %d pm_n_UX_rank %d n_delta_v %d n_svd_l %d n_M %d',n_k_p_r,n_w_max,pm_n_UX_rank,n_delta_v,n_svd_l,n_M)); end;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
plot(FTK.delta_x_,FTK.delta_y_,'ro','MarkerFaceColor','g');
axis(1.25*FTK.delta_r_max*[-1,+1,-1,+1]); axis image; axisnotick;
drawnow();
end;%if flag_disp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_speed_vs_error==1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_t = tic();
M_k_q_wkM__ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wkM__);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% M_k_q_wkM__: %0.6fs',tmp_t)); end;
tmp_t = tic();
svd_V_UX_M_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M,M_k_q_wkM__,pm_n_UX_rank,pm_UX_kn__,pm_X_weight_r_);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% svd_V_UX_M_lwnM____: time %0.6fs',tmp_t)); end;
tmp_t = tic();
UX_T_M_l2_dM__ = tfpmh_UX_T_M_l2_dM__1(FTK,n_w_,n_M,pm_n_UX_rank,svd_V_UX_M_lwnM____);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% tfpmh_UX_T_M_l2_dM__1: time %0.6fs',tmp_t)); end;
tmp_t = tic();
UX_T_M_l2_dM__ = tfpmh_UX_T_M_l2_dM__3(FTK,n_w_,n_M,pm_n_UX_rank,svd_V_UX_M_lwnM____);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% tfpmh_UX_T_M_l2_dM__3: time %0.6fs',tmp_t)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_speed_vs_error==1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_speed_vs_error==0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_t = tic();
M_k_q_wkM__ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wkM__);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% M_k_q_wkM__: %0.6fs',tmp_t)); end;
tmp_t = tic();
svd_V_UX_M_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M,M_k_q_wkM__,pm_n_UX_rank,pm_UX_kn__,pm_X_weight_r_);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% svd_V_UX_M_lwnM____: time %0.6fs',tmp_t)); end;
tmp_t = tic();
UX_T_M_l2_dM__ = tfpmh_UX_T_M_l2_dM__1(FTK,n_w_,n_M,pm_n_UX_rank,svd_V_UX_M_lwnM____);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% tfpmh_UX_T_M_l2_dM__1: time %0.6fs',tmp_t)); end;
tmp_t = tic();
UX_T_M_l2_dM__ = tfpmh_UX_T_M_l2_dM__3(FTK,n_w_,n_M,pm_n_UX_rank,svd_V_UX_M_lwnM____);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% tfpmh_UX_T_M_l2_dM__3: time %0.6fs',tmp_t)); end;
%%%%%%%%;
% check. ;
%%%%%%%%;
gamma_z_ = linspace(0,2*pi,1+n_w_max); gamma_z_ = transpose(gamma_z_(1:n_w_max));
pm_n_k_p_r = pm_n_UX_rank; pm_n_w_max = n_w_max;
pm_n_w_ = pm_n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_sum = pm_n_k_p_r*pm_n_w_max;
pm_wUX_kn__ = diag(pm_X_weight_r_)*pm_UX_kn__;
CTF_k_p_wk_ = reshape(repmat(reshape(CTF_k_p_r_k_,[1,n_k_p_r]),[n_w_max,1]),[n_w_sum,1]);
tmp_UX_T_M_l2_errrel = 0.0;
for nM=0:n_M-1;
for ndelta_v=0:n_delta_v-1;
UX_T_M_l2_tfpm = UX_T_M_l2_dM__(1+ndelta_v,1+nM);
delta_x = FTK.delta_x_(1+ndelta_v); delta_y = FTK.delta_y_(1+ndelta_v);
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
T_M_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wk_,+delta_x,+delta_y);
UX_T_M_k_p_wn_ = reshape(reshape(T_M_k_p_wk_,[n_w_max,n_k_p_r])*pm_wUX_kn__,[pm_n_w_sum,1]);
UX_T_M_l2_quad = sum(conj(UX_T_M_k_p_wn_).*UX_T_M_k_p_wn_) / max(1,pm_n_w_max) ;
tmp_UX_T_M_l2_errrel = tmp_UX_T_M_l2_errrel + fnorm_disp(max(0,flag_verbose-1),'UX_T_M_l2_quad',UX_T_M_l2_quad,'UX_T_M_l2_tfpm',UX_T_M_l2_tfpm,' %%<-- should be zero');
end;%for ndelta_v=0:n_delta_v-1;
end;%for nM=0:n_M-1;
if (flag_verbose>0); disp(sprintf(' %% tmp_UX_T_M_l2_errrel: %0.16f',tmp_UX_T_M_l2_errrel)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_speed_vs_error==0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
