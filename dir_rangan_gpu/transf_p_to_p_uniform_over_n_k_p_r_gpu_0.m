function ...
[ ...
 T_M_k_p_gpu_wkdS___ ...
] = ....
transf_p_to_p_uniform_over_n_k_p_r_gpu_0( ...
 n_k_p_r ...
,k_p_gpu_r_ ...
,n_w_ ...
,n_S ...
,S_k_p_gpu_wkS__ ...
,n_delta_v ...
,delta_x_gpu_dS__ ...
,delta_y_gpu_dS__ ...
);
%%%%%%%%;
% Assumes that n_w_ = n_w_max*ones(n_k_p_r,1);
% This 'per-template' translation is appropriate ;
% for the construction of an array of the form M_k_p_wkdM___ ;
% (as might be required for the FTK). ;
%%%%%%%%;

str_thisfunction = 'transf_p_to_p_uniform_over_n_k_p_r_gpu_0';

na=0;
if nargin<1+na; n_k_p_r=[]; end; na=na+1;
if nargin<1+na; k_p_gpu_r_=[]; end; na=na+1;
if nargin<1+na; n_w_=[]; end; na=na+1;
if nargin<1+na; n_S=[]; end; na=na+1;
if nargin<1+na; S_k_p_gpu_wkS__=[]; end; na=na+1;
if nargin<1+na; n_delta_v=[]; end; na=na+1;
if nargin<1+na; delta_x_gpu_dS__=[]; end; na=na+1;
if nargin<1+na; delta_y_gpu_dS__=[]; end; na=na+1;

f_zero = gpuArray( single(0.0));
if ~strcmp(class(k_p_gpu_r_),'gpuArray'); tmp_t = tic(); k_p_gpu_r_ = gpuArray( (k_p_gpu_r_)); tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% k_p_gpu_r_: time %0.6fs',tmp_t)); end; end;
if ~strcmp(class(S_k_p_gpu_wkS__),'gpuArray'); tmp_t = tic(); S_k_p_gpu_wkS__ = gpuArray( (S_k_p_gpu_wkS__)); tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% S_k_p_gpu_wkS__: time %0.6fs',tmp_t)); end; end;
if ~strcmp(class(delta_x_gpu_dS__),'gpuArray'); tmp_t = tic(); delta_x_gpu_dS__ = gpuArray( (delta_x_gpu_dS__)); tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% delta_x_gpu_dS__: time %0.6fs',tmp_t)); end; end;
if ~strcmp(class(delta_y_gpu_dS__),'gpuArray'); tmp_t = tic(); delta_y_gpu_dS__ = gpuArray( (delta_y_gpu_dS__)); tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% delta_y_gpu_dS__: time %0.6fs',tmp_t)); end; end;

n_w_ = n_w_(:); n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
if numel(unique(n_w_))> 1; disp(sprintf(' %% Warning, n_w_ not uniform in %s',str_thisfunction)); end;

%%%%%%%%%%%%%%%%;
gamma_gpu_z_ = gpuArray( linspace(0,2*pi,n_w_max+1) ); gamma_gpu_z_ = gpuArray( transpose(gamma_gpu_z_(1:n_w_max)) );
k_c_0_gpu_wk_ = reshape(+cos(gamma_gpu_z_)*transpose(k_p_gpu_r_),[n_w_sum,1]);
k_c_1_gpu_wk_ = reshape(+sin(gamma_gpu_z_)*transpose(k_p_gpu_r_),[n_w_sum,1]);
L_c_gpu_wkdS___ = bsxfun(@times,k_c_0_gpu_wk_,reshape(delta_x_gpu_dS__,[1,n_delta_v,n_S])) + bsxfun(@times,k_c_1_gpu_wk_,reshape(delta_y_gpu_dS__,[1,n_delta_v,n_S])) ;
C_c_gpu_wkdS___ = exp(-i*2*pi*L_c_gpu_wkdS___);
T_M_k_p_gpu_wkdS___ = reshape(bsxfun(@times,C_c_gpu_wkdS___,reshape(S_k_p_gpu_wkS__,[n_w_sum,1,n_S])),[n_w_sum,n_delta_v,n_S]);
%%%%%%%%%%%%%%%%;
