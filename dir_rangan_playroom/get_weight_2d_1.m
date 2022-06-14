function ...
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
] = ...
get_weight_2d_1( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
);

na=0;
if (nargin<1+na); verbose=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); template_k_eq_d=[]; end; na=na+1;
if (nargin<1+na); n_w_0in_=[]; end; na=na+1;

if isempty(template_k_eq_d); template_k_eq_d=-1; end;

n_w_ = zeros(n_k_p_r,1);
if (template_k_eq_d>0);
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_equator = 3+round(2*pi*k_p_r/template_k_eq_d);
n_polar_a = 3+round(n_equator/2);
n_w_(1+nk_p_r) = 2*n_polar_a;
end;%for nk_p_r=0:n_k_p_r-1;
end;%if (template_k_eq_d>0);
if (template_k_eq_d<=0);
n_w_ = n_w_0in_;
assert(numel(n_w_)==n_k_p_r); assert(min(n_w_)>0);
end;%if (template_k_eq_d<=0);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
if (verbose); disp(sprintf(' %% n_w_max %d n_w_sum %d',n_w_max,n_w_sum)); end;
%%%%%%%%;
if (verbose); disp(sprintf(' %% Set up integration weights for the templates.')); end;
%%%%%%%%;
tmp_P_ = zeros(n_k_p_r,n_k_p_r); %<-- polynomials of order 0:n_k_p_r-1 evaluated on k_p_r_/k_p_r_max. ;
tmp_I_ = zeros(n_k_p_r,1); %<-- integrals of those polynomials on the 2d-disc of radius 1. ;
for nk_p_r=0:n_k_p_r-1;
tmp_x = @(x) x.^nk_p_r;
tmp_P_(1+nk_p_r,:) = tmp_x(k_p_r_/k_p_r_max);
tmp_I_(1+nk_p_r) = 2*pi*1/(nk_p_r+2);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_W_ = pinv(tmp_P_,1e-6)*tmp_I_;
if (verbose>1); disp(sprintf(' %% weight error: %0.16f',fnorm(tmp_P_*tmp_W_ - tmp_I_)/fnorm(tmp_I_))); end;
weight_2d_k_p_r_ = tmp_W_*k_p_r_max^2;
weight_2d_k_all_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
weight_2d_k_all_(1+tmp_ij_) = weight_2d_k_p_r_(1+nk_p_r) / max(1,n_w_(1+nk_p_r)) / (2*pi)^2;
end;%for nk_p_r=0:n_k_p_r-1;
