function ...
S_x_p_tx_ = ...
interp_k_p_to_x_p_xxnufft( ...
 n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,n_t_ ...
,weight_2d_x_p_tx_ ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_2d_k_p_wk_ ...
,S_k_p_wk_ ...
);
%%%%%%%%;
str_thisfunction = 'interp_k_p_to_x_p_xxnufft';
flag_verbose=0;
if (flag_verbose>0); sprintf(' %% [entering %s]',str_thisfunction); end;

n_t_sum = sum(n_t_);
%%%%%%%%;
if numel(unique(n_t_)==1);
n_t_max = n_t_(1+0);
theta_ = 2*pi*transpose([0:n_t_max-1])/max(1,n_t_max);
x_p_t_tx_ = reshape(repmat(theta_,[1,n_x_p_r]),[n_t_sum,1]);
x_p_r_tx_ = reshape(transpose(repmat(x_p_r_,[1,n_t_max])),[n_t_sum,1]);
x_c_0_tx_ = x_p_r_tx_ .* cos(x_p_t_tx_);
x_c_1_tx_ = x_p_r_tx_ .* sin(x_p_t_tx_);
end;%if numel(unique(n_t_)==1);
%%%%%%%%
if numel(unique(n_t_)> 1);
n_t_sum = sum(n_t_(1:n_x_p_r));
x_c_0_tx_ = zeros(n_t_sum,1);
x_c_1_tx_ = zeros(n_t_sum,1);
nt_sum=0;
for nx_p_r=0:n_x_p_r-1;
x_p_r = x_p_r_(1+nx_p_r);
n_t = n_t_(1+nx_p_r);
for nt=0:n_t-1;
theta = (2.0d0*pi*nt)/max(1,n_t);
x_c_0_tx_(1+nt_sum) = x_p_r*cos(theta);
x_c_1_tx_(1+nt_sum) = x_p_r*sin(theta);
nt_sum = nt_sum+1;
end;%for nt=0:n_t-1;
end;%for nx_p_r=0:n_x_p_r-1;
if (nt_sum~=n_t_sum); disp(sprintf(' %% Warning, nt_sum %d vs %d',nt_sum,n_t_sum)); end;
end;%if numel(unique(n_t_)> 1);
%%%%%%%%;

n_w_sum = sum(n_w_);
%%%%%%%%;
if numel(unique(n_w_)==1);
n_w_max = n_w_(1+0);
omega_ = 2*pi*transpose([0:n_w_max-1])/max(1,n_w_max);
k_p_w_wk_ = reshape(repmat(omega_,[1,n_k_p_r]),[n_w_sum,1]);
k_p_r_wk_ = reshape(transpose(repmat(k_p_r_,[1,n_w_max])),[n_w_sum,1]);
k_c_0_wk_ = k_p_r_wk_ .* cos(k_p_w_wk_);
k_c_1_wk_ = k_p_r_wk_ .* sin(k_p_w_wk_);
end;%if numel(unique(n_w_)==1);
%%%%%%%%
if numel(unique(n_w_)> 1);
n_w_sum = sum(n_w_(1:n_k_p_r));
k_c_0_wk_ = zeros(n_w_sum,1);
k_c_1_wk_ = zeros(n_w_sum,1);
nw_sum=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
for nw=0:n_w-1;
omega = (2.0d0*pi*nw)/max(1,n_w);
k_c_0_wk_(1+nw_sum) = k_p_r*cos(omega);
k_c_1_wk_(1+nw_sum) = k_p_r*sin(omega);
nw_sum = nw_sum+1;
end;%for nw=0:n_w-1;
end;%for nk_p_r=0:n_k_p_r-1;
if (nw_sum~=n_w_sum); disp(sprintf(' %% Warning, nw_sum %d vs %d',nw_sum,n_w_sum)); end;
end;%if numel(unique(n_w_)> 1);
%%%%%%%%;

S_x_p_tx_ = zeros(n_t_sum,1);
eta = 1.0; %eta = pi/max(1e-12,k_p_r_max);
S_x_p_tx_ = xxnufft2d3(n_w_sum,2*pi*k_c_0_wk_*eta,2*pi*k_c_1_wk_*eta,S_k_p_wk_.*weight_2d_k_p_wk_,+1,1e-12,n_t_sum,x_c_0_tx_/eta,x_c_1_tx_/eta);

if (flag_verbose>0); sprintf(' %% [finished %s]',str_thisfunction); end;
