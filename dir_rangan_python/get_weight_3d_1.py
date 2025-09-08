import numpy as np
from scipy.special import roots_jacobi

'''
function ...
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_T_vs_L ...
);

na=0;
if (nargin<1+na); flag_verbose=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); k_eq_d=[]; end; na=na+1;
if (nargin<1+na); str_T_vs_L=[]; end; na=na+1;

if isempty(flag_verbose); flag_verbose=0; end;
if isempty(k_p_r_max); k_p_r_max = 1.0; end;
if isempty(k_eq_d); k_eq_d = 1.0/(2*pi); end;
if isempty(str_T_vs_L); str_T_vs_L = 'L'; end;

n_k_p_r = 1+ceil(k_p_r_max/k_eq_d); [a_jx_,a_jw_] = jacpts(n_k_p_r,0,2);
k_p_r_ = (a_jx_+1.0)*k_p_r_max/2; weight_3d_k_p_r_ = a_jw_*(k_p_r_max/2)^3;
'''
def get_weight_3d_1(flag_verbose=None, k_p_r_max=None, k_eq_d=None, str_T_vs_L=None):
    if flag_verbose is None:
        flag_verbose = 0
    if k_p_r_max is None:
        k_p_r_max = 1.0
    if k_eq_d is None:
        k_eq_d = 1.0 / (2 * np.pi)
    if str_T_vs_L is None:
        str_T_vs_L = 'L'

    n_k_p_r = 1 + int(np.ceil(k_p_r_max / k_eq_d))
    a_jx_, a_jw_ = roots_jacobi(n_k_p_r, 0, 2)
    k_p_r_ = np.reshape((a_jx_ + 1.0) * k_p_r_max / 2,-1)
    weight_3d_k_p_r_ = np.reshape(a_jw_ * (k_p_r_max / 2) ** 3,-1)

    return (
        n_k_p_r,
        k_p_r_, 
        weight_3d_k_p_r_,
    )
