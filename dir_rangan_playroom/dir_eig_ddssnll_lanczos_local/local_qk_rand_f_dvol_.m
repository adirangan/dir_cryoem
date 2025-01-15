function [f_dvol_qk_,f_dvol_qk__] = local_qk_rand_f_dvol_(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,rseed);
na=0;
if nargin<1+na; n_q=[]; end; na=na+1;
if nargin<1+na; n_k_p_r=[]; end; na=na+1;
if nargin<1+na; weight_3d_riesz_k_p_qk_=[]; end; na=na+1;
if nargin<1+na; rseed=[]; end; na=na+1;

if isempty(rseed); rseed=0; end;
rng(rseed);

n_qk = n_q*n_k_p_r;
scaling_volumetric = (4*pi)^2 * sqrt(pi/2);
f_dvol_qk_ = randn(n_qk,1); %<-- does not yet include imaginary component. Fix later! ;
tmp_ff_dvol = local_qk_f_dvol_bar_dot_g_dvol_(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,f_dvol_qk_,f_dvol_qk_);
f_dvol_qk_ = f_dvol_qk_/max(1e-12,sqrt(tmp_ff_dvol));
f_dvol_qk__ = local_qk__from_qk_(n_q,n_k_p_r,f_dvol_qk_);
