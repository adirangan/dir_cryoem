import numpy as np
from matlab_scalar_round import matlab_scalar_round

'''
function ...
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
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
,weight_3d_k_p_r_ ...
);

str_thisfunction = 'get_weight_2d_2';

na=0;
if (nargin<1+na); flag_verbose=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); template_k_eq_d=[]; end; na=na+1;
if (nargin<1+na); n_w_0in_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;

if isempty(flag_verbose); flag_verbose = 0; end;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

flag_calc = 0; 
if isempty(n_k_p_r) | isempty(k_p_r_) | isempty(weight_3d_k_p_r_); flag_calc=1; end;
if flag_calc;
if (flag_verbose>0); disp(sprintf(' %% Warning, precomputation required')); end;
end;%if flag_calc;
if isempty(template_k_eq_d); template_k_eq_d=-1; end;
if isempty(n_w_0in_);
if (flag_verbose>0); disp(sprintf(' %% calculating n_w_0in_')); end;
l_max_upb = round(2*pi*k_p_r_max); %<-- typically sufficient for 2-3 digits of precision. ;
n_w_max = 2*(l_max_upb+1);
n_w_0in_ = n_w_max*ones(n_k_p_r,1);
end;%if isempty(n_w_0in_);
if  isempty(weight_3d_k_p_r_); flag_pinv=1; end;
if ~isempty(weight_3d_k_p_r_); flag_pinv=0; end;

n_w_ = zeros(n_k_p_r,1);
if (template_k_eq_d>0);
if (flag_verbose>0); disp(sprintf(' %% template_k_eq_d> 0')); end;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_equator = 3+round(2*pi*k_p_r/template_k_eq_d);
n_polar_a = 3+round(n_equator/2);
n_w_(1+nk_p_r) = 2*n_polar_a;
end;%for nk_p_r=0:n_k_p_r-1;
end;%if (template_k_eq_d>0);
if (template_k_eq_d<=0);
if (flag_verbose>0); disp(sprintf(' %% template_k_eq_d<=0')); end;
n_w_ = n_w_0in_;
assert(numel(n_w_)==n_k_p_r); assert(min(n_w_)>0);
end;%if (template_k_eq_d<=0);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
if (flag_verbose>0); disp(sprintf(' %% n_w_max %d n_w_sum %d',n_w_max,n_w_sum)); end;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% Set up integration weights for the templates.')); end;
%%%%%%%%;
if flag_pinv==1;
if (flag_verbose>0); disp(sprintf(' %% Warning, pinv not implemented')); end;
end;%if flag_pinv==1;
if flag_pinv==0;
if (flag_verbose>0); disp(sprintf(' %% rescaling weight_3d_k_p_r_')); end;
weight_2d_k_p_r_ = 2*pi*reshape(weight_3d_k_p_r_,[n_k_p_r,1])./max(1e-12,k_p_r_);
if flag_verbose>0;
tmp_P_ = zeros(n_k_p_r,n_k_p_r); %<-- polynomials of order 0:n_k_p_r-1 evaluated on k_p_r_/k_p_r_max. ;
tmp_I_ = zeros(n_k_p_r,1); %<-- integrals of those polynomials on the 2d-disc of radius 1. ;
for nk_p_r=0:n_k_p_r-1;
tmp_x = @(x) x.^nk_p_r;
tmp_P_(1+nk_p_r,:) = tmp_x(k_p_r_/k_p_r_max);
tmp_I_(1+nk_p_r) = 2*pi*1/(nk_p_r+2);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_PW_ = tmp_P_*weight_2d_k_p_r_/k_p_r_max^2;
disp(sprintf(' %% tmp_I_ vs tmp_PW_: %0.16f',fnorm(tmp_I_-tmp_PW_)/max(1e-12,fnorm(tmp_I_))));
end;%if flag_verbose>0;
end;%if flag_pinv==0;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% constructing weight_2d_k_p_wk_')); end;
weight_2d_k_p_wk_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
weight_2d_k_p_wk_(1+tmp_index_) = weight_2d_k_p_r_(1+nk_p_r) / max(1,n_w_(1+nk_p_r)) / (2*pi)^2;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% Set up integration nodes for the templates.')); end;
%%%%%%%%;
k_c_0_wk_ = zeros(n_w_sum,1);
k_c_1_wk_ = zeros(n_w_sum,1);
k_p_r_wk_ = zeros(n_w_sum,1);
k_p_w_wk_ = zeros(n_w_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
k_p_r = k_p_r_(1+nk_p_r);
for nw=0:n_w-1;
gamma_z = (2*pi)*(1.0*nw)/max(1,n_w); cc = cos(gamma_z); sc = sin(gamma_z);
k_c_0_wk_(1+na) = k_p_r*cc;
k_c_1_wk_(1+na) = k_p_r*sc;
k_p_r_wk_(1+na) = k_p_r;
k_p_w_wk_(1+na) = gamma_z;
na=na+1;
end;%for nw=0:n_w-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
'''
def get_weight_2d_2(
        flag_verbose,
        n_k_p_r,
        k_p_r_,
        k_p_r_max,
        template_k_eq_d,
        n_w_0in_,
        weight_3d_k_p_r_,
):

    if flag_verbose is None:
        flag_verbose = 0

    if flag_verbose > 0:
        print(f" %% [entering get_weight_2d_2]")

    flag_calc = 0
    if n_k_p_r is None or k_p_r_ is None or weight_3d_k_p_r_ is None:
        flag_calc = 1

    if flag_calc:
        print(f" %% Warning, precomputation required")

    if template_k_eq_d is None:
        template_k_eq_d = -1
    if n_w_0in_ is None:
        if flag_verbose > 0:
            print(f" %% calculating n_w_0in_")
        l_max_upb = matlab_scalar_round(2 * np.pi * k_p_r_max)
        n_w_max = 2 * (l_max_upb + 1)
        n_w_0in_ = n_w_max * np.ones(n_k_p_r, dtype=int)

    n_w_ = np.zeros(n_k_p_r, dtype=int)
    if template_k_eq_d > 0:
        if flag_verbose > 0:
            print(f" %% template_k_eq_d > 0")
        for nk_p_r in range(n_k_p_r):
            k_p_r = k_p_r_[nk_p_r]
            n_equator = 3 + matlab_scalar_round(2 * np.pi * k_p_r / template_k_eq_d)
            n_polar_a = 3 + matlab_scalar_round(n_equator / 2)
            n_w_[nk_p_r] = 2 * n_polar_a
    else:
        if flag_verbose > 0:
            print(f" %% template_k_eq_d <= 0")
        n_w_ = n_w_0in_
        assert len(n_w_) == n_k_p_r and np.min(n_w_) > 0

    n_w_ = np.array(n_w_, dtype=int)
    n_w_max = np.max(n_w_)
    n_w_sum = np.sum(n_w_)
    n_w_csum_ = np.cumsum(np.concatenate(([0], n_w_)))

    if flag_verbose > 0:
        print(f" %% n_w_max {n_w_max} n_w_sum {n_w_sum}")

    weight_2d_k_p_r_ = 2 * np.pi * np.reshape(weight_3d_k_p_r_, n_k_p_r) / np.maximum(1e-12,k_p_r_)

    weight_2d_k_p_wk_ = np.zeros(n_w_sum)
    for nk_p_r in range(n_k_p_r):
        tmp_index_ = n_w_csum_[nk_p_r] + np.arange(n_w_[nk_p_r])
        weight_2d_k_p_wk_[tmp_index_] = weight_2d_k_p_r_[nk_p_r] / np.maximum(1, n_w_[nk_p_r]) / (2 * np.pi) ** 2

    k_c_0_wk_ = np.zeros(n_w_sum)
    k_c_1_wk_ = np.zeros(n_w_sum)
    k_p_r_wk_ = np.zeros(n_w_sum)
    k_p_w_wk_ = np.zeros(n_w_sum)

    na = 0
    for nk_p_r in range(n_k_p_r):
        n_w = n_w_[nk_p_r]
        k_p_r = k_p_r_[nk_p_r]
        for nw in range(n_w):
            gamma_z = (2 * np.pi) * nw / np.maximum(1, n_w)
            cc = np.cos(gamma_z)
            sc = np.sin(gamma_z)
            k_c_0_wk_[na] = k_p_r * cc
            k_c_1_wk_[na] = k_p_r * sc
            k_p_r_wk_[na] = k_p_r
            k_p_w_wk_[na] = gamma_z
            na += 1

    if flag_verbose > 0:
        print(f" %% [finished get_weight_2d_2]")

    return (
        n_w_,
        weight_2d_k_p_r_,
        weight_2d_k_p_wk_,
        k_p_r_wk_,
        k_p_w_wk_,
        k_c_0_wk_,
        k_c_1_wk_,
    )
