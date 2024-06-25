function ...
FTK = ...
ampmh_FTK_1( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,delta_r_max ...
,svd_eps ...
,n_delta_v_requested ...
,delta_x_0in_ ...
,delta_y_0in_ ...
);

na=0;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); delta_r_max=[]; end; na=na+1;
if (nargin<1+na); svd_eps=[]; end; na=na+1;
if (nargin<1+na); n_delta_v_requested=[]; end; na=na+1;
if (nargin<1+na); delta_x_0in_=[]; end; na=na+1;
if (nargin<1+na); delta_y_0in_=[]; end; na=na+1;

if isempty(delta_r_max); delta_r_max = 0; end;
if isempty(svd_eps); svd_eps = 1e-12; end;

pm_N_pixel = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2)) ; 
FTK = gen_Jsvd_FTK_7(k_p_r_max,pm_N_pixel,svd_eps);
%%%%;
if isempty(n_delta_v_requested); n_delta_v_requested = 2*FTK.n_svd_l; end;
if (n_delta_v_requested<=0); n_delta_v_requested = 2*FTK.n_svd_l; end;
%%%%;
if delta_r_max==0; [FTK.n_delta_v,FTK.delta_x_,FTK.delta_y_] = get_delta_2(delta_r_max,1); end;
%%%%;
if delta_r_max> 0;
if ~isempty(delta_x_0in_) & ~isempty(delta_y_0in_);
assert(numel(delta_x_0in_)==n_delta_v_requested);
assert(numel(delta_y_0in_)==n_delta_v_requested);
FTK.n_delta_v = n_delta_v_requested;
FTK.delta_x_ = delta_x_0in_(:);
FTK.delta_y_ = delta_y_0in_(:);
end;%if ~isempty(delta_x_0in_) & ~isempty(delta_y_0in_);
if  isempty(delta_x_0in_) |  isempty(delta_y_0in_);
[FTK.n_delta_v,FTK.delta_x_,FTK.delta_y_] = get_delta_2(delta_r_max,n_delta_v_requested);
end;%if  isempty(delta_x_0in_) |  isempty(delta_y_0in_);
end;%if delta_r_max> 0;
%%%%;
FTK.svd_d_max = delta_r_max;
FTK.svd_polyval_U_d_ = get_svd_polyval_U_d_0(FTK.svd_d_max,FTK.n_svd_d,FTK.svd_d_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_U_d_,FTK.n_delta_v,FTK.delta_x_,FTK.delta_y_);
FTK.svd_r_max = 2*pi*k_p_r_max;
FTK.svd_polyval_V_r_ = get_svd_polyval_V_r_0(FTK.svd_r_max,FTK.n_svd_r,FTK.svd_r_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_V_r_,n_k_p_r,2*pi*k_p_r_);
FTK.svd_expiw__ = exp(-i*(pi/2 - reshape(atan2(FTK.delta_y_,FTK.delta_x_),[FTK.n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l]));
FTK.svd_U_d_expiw_s__ = (permute(reshape(FTK.svd_polyval_U_d_,[FTK.n_svd_l,FTK.n_delta_v]),[2,1]).*FTK.svd_expiw__)*diag(FTK.svd_s_);
