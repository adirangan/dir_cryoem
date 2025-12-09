import numpy as np ; pi = np.pi ; import torch ; import timeit ;
from matlab_index_2d_0 import matlab_index_2d_0 ;
from matlab_index_3d_0 import matlab_index_3d_0 ;
from matlab_index_4d_0 import matlab_index_4d_0 ;
from matlab_scalar_round import matlab_scalar_round ;
cumsum_0 = lambda a : torch.cumsum(torch.concatenate((torch.tensor([0]),a)) , 0).to(torch.int32) ;
mtr = lambda a : tuple(reversed(a)) ; #<-- matlab-arranged size (i.e., tuple(reversed(...))). ;

def local_yk__from_yk_(
        n_k_p_r,
        l_max_,
        tmp_yk_,
):
    n_y_ = (l_max_+1)**2;
    n_y_max = int(torch.max(n_y_).item());
    n_y_sum = int(torch.sum(n_y_).item());
    n_y_csum_ = cumsum_0(n_y_);
    l_max_max = int(torch.max(l_max_).item());
    tmp_yk__ = torch.zeros(mtr((n_y_max,n_k_p_r))).to(dtype=torch.complex64);
    for nk_p_r in range(n_k_p_r):
        n_y = int(n_y_[nk_p_r].item());
        tmp_i8_index_rhs_ = int(n_y_csum_[nk_p_r].item()) + torch.arange(n_y);
        tmp_i8_index_lhs_ = matlab_index_2d_0(n_y_max,torch.arange(n_y),n_k_p_r,nk_p_r);
        tmp_yk__.ravel()[tmp_i8_index_lhs_] = tmp_yk_[tmp_i8_index_rhs_].to(dtype=torch.complex64);
    #end;%for nk_p_r=0:n_k_p_r-1;
    return(tmp_yk__) ;
