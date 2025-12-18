function ...
[ ...
 parameter ...
,FTK ...
] = ...
tfh_FTK_3( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
);
%%%%%%%%;
% prepends a zero-translation to the list produced by get_delta_3.m ;
%%%%%%%%;
str_thisfunction = 'tfh_FTK_3';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'delta_r_max'); parameter.delta_r_max = 0; end;
delta_r_max = parameter.delta_r_max;
if ~isfield(parameter,'svd_eps'); parameter.svd_eps = 1e-12; end;
svd_eps = parameter.svd_eps;
if ~isfield(parameter,'n_delta_v_requested'); parameter.n_delta_v_requested = 0; end;
n_delta_v_requested = parameter.n_delta_v_requested;
if ~isfield(parameter,'l_max'); parameter.l_max = []; end;
l_max = parameter.l_max;
if ~isfield(parameter,'n_a_degree'); parameter.n_a_degree = []; end;
n_a_degree = parameter.n_a_degree;
if ~isfield(parameter,'n_b_degree'); parameter.n_b_degree = []; end;
n_b_degree = parameter.n_b_degree;

if isempty(delta_r_max); delta_r_max = 0; end;
if isempty(svd_eps); svd_eps = 1e-12; end;

pm_N_pixel = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2)) ; 
FTK = gen_Jsvd_FTK_8(k_p_r_max,pm_N_pixel,svd_eps,l_max,n_a_degree,n_b_degree);
%%%%;
if isempty(n_delta_v_requested); n_delta_v_requested = 2*FTK.n_svd_l; end;
if (n_delta_v_requested<=0); n_delta_v_requested = 2*FTK.n_svd_l; end;
FTK.delta_r_max = delta_r_max;
%%%%;
if delta_r_max==0; 
FTK.n_delta_p_r = 0;
FTK.delta_p_r_ = zeros(FTK.n_delta_p_r,1);
FTK.weight_2d_p_r_ = zeros(FTK.n_delta_p_r,1);
FTK.n_delta_p_w = 0;
FTK.delta_p_w_ = zeros(FTK.n_delta_p_w,1);
FTK.delta_p_0_wr__ = zeros(FTK.n_delta_p_w,FTK.n_delta_p_r);
FTK.delta_p_1_wr__ = zeros(FTK.n_delta_p_w,FTK.n_delta_p_r);
FTK.weight_2d_p_wr__ = zeros(FTK.n_delta_p_w,FTK.n_delta_p_r);
end;%if delta_r_max==0; 
%%%%;
if delta_r_max> 0;
[ ...
 parameter ...
,FTK.n_delta_p_r ...
,FTK.delta_p_r_ ...
,FTK.weight_2d_p_r_ ...
,FTK.n_delta_p_w ...
,FTK.delta_p_w_ ...
,FTK.delta_p_0_wr__ ...
,FTK.delta_p_1_wr__ ...
,FTK.weight_2d_p_wr__ ...
] = ....
get_delta_3( ....
 parameter ....
,delta_r_max ...
,n_delta_v_requested ...
);
end;%if delta_r_max> 0;
%%%%;
FTK.n_delta_v = 1 + FTK.n_delta_p_w*FTK.n_delta_p_r;
FTK.delta_x_ = [0;FTK.delta_p_0_wr__(:)];
FTK.delta_y_ = [0;FTK.delta_p_1_wr__(:)];
%%%%;
FTK.svd_d_max = delta_r_max;
FTK.svd_chebval_U_d_ = reshape(get_svd_chebval_U_d_0(FTK.svd_d_max,FTK.n_svd_d,FTK.svd_d_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_U_d_chebcoef_,FTK.n_delta_v,FTK.delta_x_,FTK.delta_y_),[FTK.n_svd_l*FTK.n_delta_v,1]);
FTK.svd_r_max = 2*pi*k_p_r_max;
FTK.svd_chebval_V_r_ = reshape(get_svd_chebval_V_r_0(FTK.svd_r_max,FTK.n_svd_r,FTK.svd_r_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_V_r_chebcoef_,n_k_p_r,2*pi*k_p_r_),[FTK.n_svd_l*n_k_p_r,1]);
FTK.svd_expiw__ = reshape(exp(-i*(pi/2 - reshape(atan2(FTK.delta_y_,FTK.delta_x_),[FTK.n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l])),[FTK.n_delta_v,FTK.n_svd_l]);
FTK.svd_U_d_expiw_s__ = (permute(reshape(FTK.svd_chebval_U_d_,[FTK.n_svd_l,FTK.n_delta_v]),1+[1,0]).*FTK.svd_expiw__)*diag(FTK.svd_s_);
%%%%;
FTK.svd_chebval_U_d_l0_ = reshape(get_svd_chebval_U_d_0(FTK.svd_d_max,FTK.n_svd_d,FTK.svd_d_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_U_d_chebcoef_,1,0,0),[FTK.n_svd_l*1,1]);
%FTK.svd_chebval_U_d_rl__ = permute(reshape(get_svd_chebval_U_d_0(FTK.svd_d_max,FTK.n_svd_d,FTK.svd_d_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_U_d_chebcoef_,FTK.n_delta_p_r,FTK.delta_p_r_,zeros(FTK.n_delta_p_r,1)),[FTK.n_svd_l,FTK.n_delta_p_r]),1+[1,0]);
FTK.svd_chebval_U_d_lr__ = reshape(get_svd_chebval_U_d_0(FTK.svd_d_max,FTK.n_svd_d,FTK.svd_d_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_U_d_chebcoef_,FTK.n_delta_p_r,FTK.delta_p_r_,zeros(FTK.n_delta_p_r,1)),[FTK.n_svd_l,FTK.n_delta_p_r]);
%%%%;
FTK.svd_expil_l_ = reshape(exp(-i*(pi/2)*reshape(FTK.svd_l_,[FTK.n_svd_l,1])),[FTK.n_svd_l,1]);
FTK.svd_l_from_nl__ = sparse(1+periodize(FTK.svd_l_,0,FTK.n_delta_p_w),1:FTK.n_svd_l,1,FTK.n_delta_p_w,FTK.n_svd_l);
FTK.svd_expiwl_wl__ = reshape(exp(+i*reshape(FTK.delta_p_w_,[FTK.n_delta_p_w,1])*reshape(FTK.svd_l_,[1,FTK.n_svd_l])),[FTK.n_delta_p_w,FTK.n_svd_l]);
%%%%%%%%;
