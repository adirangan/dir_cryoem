function ...
[ ...
 T_M_k_p_wkdS___ ...
] = ....
transf_p_to_p_uniform_over_n_k_p_r_0( ...
 n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,n_S ...
,S_k_p_wkS__ ...
,n_delta_v ...
,delta_x_dS__ ...
,delta_y_dS__ ...
);
%%%%%%%%;
% Assumes that n_w_ = n_w_max*ones(n_k_p_r,1);
% This 'per-template' translation is appropriate ;
% for the construction of an array of the form M_k_p_wkdM___ ;
% (as might be required for the FTK). ;
%%%%%%%%;

str_thisfunction = 'transf_p_to_p_uniform_over_n_k_p_r_0';

na=0;
if nargin<1+na; n_k_p_r=[]; end; na=na+1;
if nargin<1+na; k_p_r_=[]; end; na=na+1;
if nargin<1+na; n_w_=[]; end; na=na+1;
if nargin<1+na; n_S=[]; end; na=na+1;
if nargin<1+na; S_k_p_wkS__=[]; end; na=na+1;
if nargin<1+na; n_delta_v=[]; end; na=na+1;
if nargin<1+na; delta_x_dS__=[]; end; na=na+1;
if nargin<1+na; delta_y_dS__=[]; end; na=na+1;

n_w_ = n_w_(:); n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
if numel(unique(n_w_))> 1; disp(sprintf(' %% Warning, n_w_ not uniform in %s',str_thisfunction)); end;

%%%%%%%%%%%%%%%%;
gamma_z_ = linspace(0,2*pi,n_w_max+1); gamma_z_ = transpose(gamma_z_(1:n_w_max));
k_c_0_wk_ = reshape(+cos(gamma_z_)*transpose(k_p_r_),[n_w_sum,1]);
k_c_1_wk_ = reshape(+sin(gamma_z_)*transpose(k_p_r_),[n_w_sum,1]);
L_c_wkdS___ = bsxfun(@times,k_c_0_wk_,reshape(delta_x_dS__,[1,n_delta_v,n_S])) + bsxfun(@times,k_c_1_wk_,reshape(delta_y_dS__,[1,n_delta_v,n_S])) ;
C_c_wkdS___ = exp(-i*2*pi*L_c_wkdS___);
T_M_k_p_wkdS___ = reshape(bsxfun(@times,C_c_wkdS___,reshape(S_k_p_wkS__,[n_w_sum,1,n_S])),[n_w_sum,n_delta_v,n_S]);
%%%%%%%%%%%%%%%%;
