function M_p_ = transf_p_to_p(n_r,grid_p_,n_w_,n_w_sum,S_p_,delta_x,delta_y,flag_each);
%%%%%%%%;
% Note that if delta_x and delta_y are given as arrays, ;
% then we expect n_S==n_delta_v. ;
% Thus, this 'per-template' translation is not appropriate ;
% for the construction of an array of the form M_wkdM___ ;
% (as might be required for the FTK). ;
%%%%%%%%;

str_thisfunction = 'transf_p_to_p';

%%%%%%%%;
if nargin<1;
%%%%%%%%;
flag_verbose=1;
n_w_int = 2;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi); str_L = 'L';
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_L ...
);
%%%%;
if (flag_verbose>0); disp(sprintf(' %% testing not adaptive')); end;
template_k_eq_d = -1;
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
,template_k_eq_d ...
,n_w_0in_ ...
,weight_3d_k_p_r_ ...
);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
if (flag_verbose>0); disp(sprintf(' %% n_w_ [ %d .. %d ]',n_w_(1),n_w_(end))); end;
n_S = 5;
S_k_p_wkS__ = randn(n_w_sum,n_S) + i*randn(n_w_sum,n_S);
delta_x = +0.1; delta_y = -0.2;
R_k_p_wkS__ = reshape(transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,S_k_p_wkS__,delta_x,delta_y),[n_w_sum,n_S]);
T_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
T_k_p_wkS__(:,1+nS) = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS),delta_x,delta_y,1);
end;%for nS=0:n_S-1;
fnorm_disp(flag_verbose,'T_k_p_wkS__',T_k_p_wkS__,'R_k_p_wkS__',R_k_p_wkS__,' %%<-- should be zero');
n_delta_v = n_S; delta_x_ = +0.1*randn(n_delta_v,1); delta_y_ = +0.1*randn(n_delta_v,1);
R_k_p_wkS__ = reshape(transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,S_k_p_wkS__,delta_x_,delta_y_),[n_w_sum,n_S]);
T_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
T_k_p_wkS__(:,1+nS) = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS),delta_x_(1+nS),delta_y_(1+nS),1);
end;%for nS=0:n_S-1;
fnorm_disp(flag_verbose,'T_k_p_wkS__',T_k_p_wkS__,'R_k_p_wkS__',R_k_p_wkS__,' %%<-- should be zero');
%%%%;
if (flag_verbose>0); disp(sprintf(' %% testing yes adaptive')); end;
template_k_eq_d = -1;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
l_max_ = round(2*pi*k_p_r_);
n_w_0in_ = n_w_int*2*(l_max_+1);
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
,template_k_eq_d ...
,n_w_0in_ ...
,weight_3d_k_p_r_ ...
);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
if (flag_verbose>0); disp(sprintf(' %% n_w_ [ %d .. %d ]',n_w_(1),n_w_(end))); end;
n_S = 5;
S_k_p_wkS__ = randn(n_w_sum,n_S) + i*randn(n_w_sum,n_S);
delta_x = +0.1; delta_y = -0.2;
R_k_p_wkS__ = reshape(transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,S_k_p_wkS__,delta_x,delta_y),[n_w_sum,n_S]);
T_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
T_k_p_wkS__(:,1+nS) = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS),delta_x,delta_y,1);
end;%for nS=0:n_S-1;
fnorm_disp(flag_verbose,'T_k_p_wkS__',T_k_p_wkS__,'R_k_p_wkS__',R_k_p_wkS__,' %%<-- should be zero');
n_delta_v = n_S; delta_x_ = +0.1*randn(n_delta_v,1); delta_y_ = +0.1*randn(n_delta_v,1);
R_k_p_wkS__ = reshape(transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,S_k_p_wkS__,delta_x_,delta_y_),[n_w_sum,n_S]);
T_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
T_k_p_wkS__(:,1+nS) = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS),delta_x_(1+nS),delta_y_(1+nS),1);
end;%for nS=0:n_S-1;
fnorm_disp(flag_verbose,'T_k_p_wkS__',T_k_p_wkS__,'R_k_p_wkS__',R_k_p_wkS__,' %%<-- should be zero');
%%%%;
disp('returning'); return;
%%%%%%%%;
end;%if nargin<1;
%%%%%%%%;

na=0;
if nargin<1+na; n_r=[]; end; na=na+1;
if nargin<1+na; grid_p_=[]; end; na=na+1;
if nargin<1+na; n_w_=[]; end; na=na+1;
if nargin<1+na; n_w_sum=[]; end; na=na+1;
if nargin<1+na; S_p_=[]; end; na=na+1;
if nargin<1+na; delta_x=[]; end; na=na+1;
if nargin<1+na; delta_y=[]; end; na=na+1;
if nargin<1+na; flag_each=[]; end; na=na+1;

if isempty(flag_each); flag_each = 0; end;

n_S = idivide(cast(numel(S_p_),'int64'),cast(n_w_sum,'int64'));
if numel(S_p_)~=n_w_sum*n_S;
disp(sprintf(' %% Warning, n_w_sum %d n_S %d in %s',n_w_sum,n_S,str_thisfunction));
end;%if numel(S_p_)~=n_w_sum*n_S;

assert(numel(delta_y)==numel(delta_x));
n_delta_v = numel(delta_x);
delta_x_ = reshape(delta_x,[n_delta_v,1]);
delta_y_ = reshape(delta_y,[n_delta_v,1]);
if (n_delta_v> 1) & (n_S~=n_delta_v);
disp(sprintf(' %% Warning, n_S %d n_delta_v %d in %s',n_S,n_delta_v,str_thisfunction));
end;%if (n_delta_v> 1) & (n_S~=n_delta_v);

%%%%%%%%%%%%%%%%;
if  flag_each | (numel(unique(n_w_))> 1);
M_p_ = zeros(n_w_sum*n_S,1);
%%%%%%%%;
if n_S==1;
%%%%;
% single transformation on adaptive grid. ;
%%%%;
M_p_ = zeros(size(S_p_));
nw_sum=0;
for nr=0:n_r-1;
R_c = grid_p_(1+nr);
n_w = n_w_(1+nr);
for nw=0:n_w-1;
W_c = 0.0 + nw*(2*pi)/max(1,n_w);
X_c = R_c*cos(W_c);
Y_c = R_c*sin(W_c);
L_c = (X_c * delta_x) + (Y_c * delta_y);
C_c = +cos(2*pi*L_c) - i*sin(2*pi*L_c);
M_p_(1+nw_sum) = C_c*S_p_(1+nw_sum);
nw_sum = nw_sum + 1;
end;%for nw=0:n_w-1;
end;%for nr=0:n_r-1;
if (nw_sum~=n_w_sum); disp(sprintf(' %% Warning: nw_sum %d vs %d',nw_sum,n_w_sum)); end;
%%%%;
end;%if n_S==1;
%%%%%%%%;
if n_S> 1;
for nS=0:n_S-1;
tmp_index_ = cast(nS*n_w_sum,'double') + [0:n_w_sum-1];
if n_delta_v==1;
M_p_(1+tmp_index_) = transf_p_to_p(n_r,grid_p_,n_w_,n_w_sum,S_p_(1+tmp_index_),delta_x,delta_y);
end;%if n_delta_v==1;
if n_delta_v> 1;
M_p_(1+tmp_index_) = transf_p_to_p(n_r,grid_p_,n_w_,n_w_sum,S_p_(1+tmp_index_),delta_x_(1+nS),delta_y_(1+nS));
end;%if n_delta_v> 1;
end;%for nS=0:n_S-1;
M_p_ = reshape(M_p_,size(S_p_));
end;%if n_S> 1;
%%%%%%%%;
end;%if  flag_each | (numel(unique(n_w_))> 1);
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%;
if ~flag_each & (numel(unique(n_w_))==1);
n_w = n_w_(1+0); n_w_max = n_w;
gamma_z_ = linspace(0,2*pi,n_w_max+1); gamma_z_ = transpose(gamma_z_(1:n_w_max));
k_c_0_wk_ = reshape(+cos(gamma_z_)*transpose(grid_p_),[n_w_sum,1]);
k_c_1_wk_ = reshape(+sin(gamma_z_)*transpose(grid_p_),[n_w_sum,1]);
L_c_wkv__ = k_c_0_wk_*reshape(delta_x_,[1,n_delta_v]) + k_c_1_wk_*reshape(delta_y_,[1,n_delta_v]) ;
C_c_wkv__ = exp(-i*2*pi*L_c_wkv__);
M_p_ = reshape(C_c_wkv__.*reshape(S_p_,[n_w_sum,n_S]),[n_w_sum*n_S,1]);
M_p_ = reshape(M_p_,size(S_p_));
end;%if ~flag_each & (numel(unique(n_w_))==1);
%%%%%%%%%%%%%%%%;

%%%%%%%%;
% older code for reference. ;
%%%%%%%%;
%{
% Assumes that M_p_ is the same size and dimensions as S_p_;
% Multiplication performed in place;
M_p_ = zeros(size(S_p_));
nw_sum=0;
for nr=0:n_r-1;
R_c = grid_p_(1+nr);
for nw=0:n_w_(1+nr)-1;
W_c = 0.0 + nw*(2*pi)/(n_w_(1+nr));
X_c = R_c*cos(W_c);
Y_c = R_c*sin(W_c);
L_c = (X_c * delta_x) + (Y_c * delta_y);
C_c = +cos(2*pi*L_c) - i*sin(2*pi*L_c);
M_p_(1+nw_sum) = C_c*S_p_(1+nw_sum);
nw_sum = nw_sum + 1;
end;%for;
end;%for;
if (nw_sum~=n_w_sum); disp(sprintf(' %% Warning: nw_sum %d vs %d',nw_sum,n_w_sum)); end;
%}
