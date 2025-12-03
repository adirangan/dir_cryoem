flag_verbose=1; flag_disp=0;nf=0;
disp(sprintf(' %% testing tpmh_VUXM_lwnM___3'));
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
n_M = 4; n_2 = 2;
delta_2M__ = reshape(linspace(-0.05,+0.05,n_2*n_M),[n_2,n_M]);
M_k_p_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
delta_0 = delta_2M__(1+0,1+nM); delta_1 = delta_2M__(1+1,1+nM);
M_k_p_wkM__(:,1+nM) = exp(+i*2*pi*(k_c_0_wk_*delta_0 + k_c_1_wk_*delta_1));
end;%for nM=0:n_M-1;
%%%%;
CTF_phi = pi*2/3;
CTF_k_p_wk_ = 2*k_p_r_wk_.*cos(k_p_w_wk_ - CTF_phi);
%%%%%%%%;

%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions.; 
pm_n_UX_rank = n_UX_rank; %<-- just to check dimension. ;
UX0_kn__ = zeros(n_k_p_r,n_UX_rank);
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
[tmp_UX0_kn__,tmp_SX0_k__,tmp_VX0_kn__] = svds(X_kk__,n_UX_rank); tmp_SX0_k_ = diag(tmp_SX0_k__);
UX0_kn__(:,1:n_UX_rank) = tmp_UX0_kn__(:,1:n_UX_rank);
SX0_n_ = tmp_SX0_k_(1:n_UX_rank);
%%%%%%%%;

%%%%%%%%;
delta_r_max = 0.5/max(1e-12,k_p_r_max); svd_eps = 1e-6; n_delta_v_requested = 128;
FTK = tfh_FTK_2(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
CTF_M_k_p_wkM__ = bsxfun(@times,CTF_k_p_wk_,M_k_p_wkM__);
CTF_M_k_q_wkM__ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,CTF_M_k_p_wkM__);
svd_VUXCTFM_ori_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M,CTF_M_k_q_wkM__,pm_n_UX_rank,UX0_kn__,X_weight_r_);
%%%%%%%%;

%%%%;
% compare to ../dir_rangan_python/dir_ascii/. ;
%%%%;
dir_ascii = '../dir_rangan_python/dir_ascii';
n_d = 2;
for nd=0:n_d-1;
na=0;
if nd==na; tmp_str_ascii = 'UX_kn__'; tmp_str_local = 'UX0_kn__'; tmp_local = UX0_kn__; tab_p=0; end; na=na+1;
if nd==na; tmp_str_ascii = 'svd_VUXCTFM_lwnM____'; tmp_str_local = 'svd_VUXCTFM_ori_lwnM____'; tmp_local = svd_VUXCTFM_ori_lwnM____; tab_p=1; end; na=na+1;
fname_ascii = sprintf('%s/%s.ascii',dir_ascii,tmp_str_ascii);
if ~exist(fname_ascii,'file'); disp(sprintf(' %% Warning, %s not found',fname_ascii)); end;
if  exist(fname_ascii,'file');
if (flag_verbose>1); disp(sprintf(' %% %s found, loading',fname_ascii)); end;
fid = fopen(fname_ascii,'r');
if tab_p==0; tmp_ascii = reshape(cell2mat(textscan(fid,'%f')),size(tmp_local)); end;
if tab_p==1; tmp_ascii = reshape(cell2mat(textscan(fid,'(%f)')),size(tmp_local)); end;
fclose(fid);
fnorm_disp(flag_verbose,tmp_str_local,tmp_local,tmp_str_ascii,tmp_ascii,' %% <-- can be large');
na=0;
if nd==na; UX1_kn__ = reshape(tmp_ascii,size(UX0_kn__)); end; na=na+1;
if nd==na; svd_VUXCTFM_xpu_lwnM____ = reshape(tmp_ascii,size(svd_VUXCTFM_ori_lwnM____)); end; na=na+1;
end;%if  exist(fname_ascii,'file');
end;%for nd=0:n_d-1;
%%%%%%%%;

%%%%%%%%;
% Now repeat with alternative UX1_kn__. ;
%%%%%%%%;
svd_VUXCTFM_upd_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M,CTF_M_k_q_wkM__,pm_n_UX_rank,UX1_kn__,X_weight_r_);
fnorm_disp(flag_verbose,'svd_VUXCTFM_upd_lwnM____',svd_VUXCTFM_upd_lwnM____,'svd_VUXCTFM_xpu_lwnM____',svd_VUXCTFM_xpu_lwnM____,' %<-- should be <1e-6');
%%%%%%%%;

tmp_str_ascii = 'svd_VUXCTFM_cpu_lwnM____'; tmp_str_local = 'svd_VUXCTFM_upd_lwnM____'; tmp_local = svd_VUXCTFM_upd_lwnM____; tab_p=1;
fname_ascii = sprintf('%s/%s.ascii',dir_ascii,tmp_str_ascii);
if ~exist(fname_ascii,'file'); disp(sprintf(' %% Warning, %s not found',fname_ascii)); end;
if  exist(fname_ascii,'file');
if (flag_verbose>1); disp(sprintf(' %% %s found, loading',fname_ascii)); end;
fid = fopen(fname_ascii,'r');
if tab_p==0; tmp_ascii = reshape(cell2mat(textscan(fid,'%f')),size(tmp_local)); end;
if tab_p==1; tmp_ascii = reshape(cell2mat(textscan(fid,'(%f)')),size(tmp_local)); end;
fclose(fid);
fnorm_disp(flag_verbose,tmp_str_local,tmp_local,tmp_str_ascii,tmp_ascii,' %% <-- should be <1e-6');
svd_VUXCTFM_cpu_lwnM____ = reshape(tmp_ascii,size(svd_VUXCTFM_upd_lwnM____));
end;%if  exist(fname_ascii,'file');

tmp_str_ascii = 'svd_VUXCTFM_gpu_lwnM____'; tmp_str_local = 'svd_VUXCTFM_upd_lwnM____'; tmp_local = svd_VUXCTFM_upd_lwnM____; tab_p=1;
fname_ascii = sprintf('%s/%s.ascii',dir_ascii,tmp_str_ascii);
if ~exist(fname_ascii,'file'); disp(sprintf(' %% Warning, %s not found',fname_ascii)); end;
if  exist(fname_ascii,'file');
if (flag_verbose>1); disp(sprintf(' %% %s found, loading',fname_ascii)); end;
fid = fopen(fname_ascii,'r');
if tab_p==0; tmp_ascii = reshape(cell2mat(textscan(fid,'%f')),size(tmp_local)); end;
if tab_p==1; tmp_ascii = reshape(cell2mat(textscan(fid,'(%f)')),size(tmp_local)); end;
fclose(fid);
fnorm_disp(flag_verbose,tmp_str_local,tmp_local,tmp_str_ascii,tmp_ascii,' %% <-- should be <1e-6');
svd_VUXCTFM_gpu_lwnM____ = reshape(tmp_ascii,size(svd_VUXCTFM_upd_lwnM____));
end;%if  exist(fname_ascii,'file');

%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;fig81s;
p_row = 1; p_col = 2; np=0;
llim_ = [-7,0];
%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(log10(SX0_n_),log10(sum(abs(UX0_kn__-UX1_kn__).^2,1)),'ko-');
xlabel('log10(SX0_n_)','Interpreter','none');
ylabel('log10(sum(abs(UX0_kn__-UX1_kn__).^2,1))','Interpreter','none');
title('UX_kn__ 0 vs 1','Interpreter','none');
%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
tmp_n_ = squeeze(sum(abs(svd_VUXCTFM_xpu_lwnM____ - svd_VUXCTFM_upd_lwnM____).^2,1+[0,1,3]));
plot(log10(SX0_n_),log10(tmp_n_),'ko-');
xlabel('log10(SX0_n_)','Interpreter','none');
ylabel('log10(sum(abs(svd_VUXCTFM_xpu_lwnM____ - svd_VUXCTFM_upd_lwnM____).^2,1+[0,2,3]))','Interpreter','none');
title('svd_VUXCTFM_lwnM____ 0 vs 1','Interpreter','none');
%%;
end;%if flag_disp;
%%%%%%%%;
