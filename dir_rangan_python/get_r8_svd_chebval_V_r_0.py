import numpy as np ; pi = np.pi ; i = 1j ; import torch ; import timeit ;
from matlab_index_2d_0 import matlab_index_2d_0 ;
from matlab_index_3d_0 import matlab_index_3d_0 ;
from matlab_index_4d_0 import matlab_index_4d_0 ;
from matlab_scalar_round import matlab_scalar_round ;
from r8_chebval_0 import r8_chebval_0 ;
mtr = lambda a : tuple(reversed(a)) ; #<-- matlab-arranged size (i.e., tuple(reversed(...))). ;
msr = lambda str : str[::-1] ; #<-- for einsum (i.e., string reversed (...)). ;
mts = lambda a : tuple(len(a) - x - 1 for x in a) ; #<-- for permute (i.e., tuple subtract (...)). ;
n_byte_per_float32 = 4;
n_byte_per_complex64 = 8;

def get_r8_svd_chebval_V_r_0(
        r8_svd_r_max=None,
        n_svd_r=None,
        r8_svd_r_=None,
        n_svd_l=None,
        i4_svd_l_=None,
        r8_svd_V_r_chebcoef_=None,
        n_r=None,
        r8_grid_p_=None,
):
    str_thisfunction = 'get_r8_svd_chebval_V_r_0' ;
    flag_verbose=0; flag_warning=1;
    if (flag_verbose>0): print(f' %% [entering {str_thisfunction}]');
    r8_svd_r_m = r8_svd_r_max / 2.0;
    r8_svd_r_c = r8_svd_r_m;
    #%%%%%%%%;
    #% vect version. ;
    #%%%%%%%%;
    tolerance_margin = 1e-6;
    if torch.sum(r8_grid_p_>r8_svd_r_max+tolerance_margin).item()>0 and flag_warning:
        print(f' %% Warning, r8_grid_p_ > r8_svd_r_max {r8_svd_r_max:.2f}');
    #end;%if sum(r8_grid_p_>r8_svd_r_max+tolerance_margin & flag_warning);
    tmp_r8_svd_r_ = (r8_grid_p_ - r8_svd_r_m)/np.maximum(1e-12,r8_svd_r_c);
    r8_svd_chebval_V_r_lr__ = torch.zeros(mtr((n_svd_l,n_r))).to(dtype=torch.float64);
    for nl in range(n_svd_l):
        tmp_index_rhs_ = 0+nl*n_svd_r+torch.arange(n_svd_r).to(dtype=torch.int32);
        r8_svd_chebval_V_r_lr__[:,nl] = r8_chebval_0(n_svd_r,r8_svd_V_r_chebcoef_[tmp_index_rhs_],n_r,tmp_r8_svd_r_);
    #end;%for nl=0:n_svd_l-1;
    r8_svd_chebval_V_r_ = r8_svd_chebval_V_r_lr__.ravel();
    assert(r8_svd_chebval_V_r_.numel()==n_svd_l*n_r);

    if (flag_verbose>0): print(f' %% [finished {str_thisfunction}]');
    return(r8_svd_chebval_V_r_);

