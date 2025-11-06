import numpy as np ; pi = np.pi ; i = 1j ; import torch ; import timeit ;
from matlab_index_2d_0 import matlab_index_2d_0 ;
from matlab_index_3d_0 import matlab_index_3d_0 ;
from matlab_index_4d_0 import matlab_index_4d_0 ;
from matlab_scalar_round import matlab_scalar_round ;
mtr = lambda a : tuple(reversed(a)) ; #<-- matlab-arranged size (i.e., tuple(reversed(...))). ;
msr = lambda str : str[::-1] ; #<-- for einsum (i.e., string reversed (...)). ;
mts = lambda a : tuple(len(a) - x - 1 for x in a) ; #<-- for permute (i.e., tuple subtract (...)). ;
n_byte_per_float32 = 4;
n_byte_per_complex64 = 8;

def r8_chebval_0(
        n_t=None,
        r8_t_chebcoef_=None,
        n_x=None,
        r8_x_=None,
):
    r8_1p = torch.tensor(+1.0).to(dtype=torch.float64);
    r8_1n = torch.tensor(-1.0).to(dtype=torch.float64);
    r8_v_ = torch.zeros(n_x).to(dtype=torch.float64);
    r8_nt_ = torch.arange(n_t).to(dtype=torch.float64); #<-- integer valued floats. ;
    r8_x__,r8_nt__ = torch.meshgrid(r8_x_.ravel(),r8_nt_,indexing='ij'); #<-- reversed to match matlab. ;
    r8_v_ = torch.reshape(torch.sum(torch.reshape(r8_t_chebcoef_,mtr((n_t,1))) * torch.cos(r8_nt__ * torch.acos(torch.maximum(r8_1n,torch.minimum(r8_1p,r8_x__)))),1-0),mtr(r8_x_.size())) ;
    return(r8_v_);

