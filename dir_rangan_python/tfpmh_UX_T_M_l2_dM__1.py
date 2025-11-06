import numpy as np ; pi = np.pi ; i = 1j ; import torch ; import timeit ;
from matlab_index_2d_0 import matlab_index_2d_0 ;
from matlab_index_3d_0 import matlab_index_3d_0 ;
from matlab_index_4d_0 import matlab_index_4d_0 ;
cumsum_0 = lambda a : torch.cumsum(torch.concatenate((torch.tensor([0]),a)) , 0).to(torch.int32) ;
fnorm = lambda a : torch.linalg.norm(a).item() ;
mtr = lambda a : tuple(reversed(a)) ; #<-- matlab-arranged size (i.e., tuple(reversed(...))). ;
msr = lambda str : str[::-1] ; #<-- for einsum (i.e., string reversed (...)). ;
mts = lambda a : tuple(len(a) - x - 1 for x in a) ; #<-- for permute (i.e., tuple subtract (...)). ;
tic = lambda : timeit.default_timer() ;
toc = lambda a : tic() - a ;
mmmm = lambda A , B : torch.einsum( msr('ab') + ',' + msr('bc') + '->' + msr('ac') , A , B ) ; #<-- matlab matrix matrix multiplication. ;
mmvm = lambda A , B : torch.einsum( msr('ab') + ',' +  msr('b') + '->' +  msr('a') , A , B ) ; #<-- matlab matrix vector multiplication. ;
mvmm = lambda A , B : torch.einsum(  msr('b') + ',' + msr('bc') + '->' +  msr('c') , A , B ) ; #<-- matlab vector matrix multiplication. ;
mvvm = lambda A , B : torch.einsum(  msr('b') + ',' +  msr('b') + '->' +   msr('') , A , B ) ; #<-- matlab vector vector multiplication. ;
efind = lambda a : torch.where(a)[0] ;
n_1 = int(1); n_2 = int(2); n_3 = int(3);

def tfpmh_UX_T_M_l2_dM__1(
        FTK=None,
        n_w_=None,
        n_M=None,
        pm_n_UX_rank=None,
        svd_VUXM_lwnM____=None,
):
    flag_verbose=0;
    str_thisfunction = 'tfpmh_UX_T_M_l2_dM__1';
    if flag_verbose > 0: print(f' %% [entering {str_thisfunction}]');
    UX_T_M_l2_dM__ = torch.zeros(mtr((FTK['n_delta_v'],n_M))).to(dtype=torch.float32);
    n_w_max = int(torch.max(n_w_).item());
    UX_T_M_k_q_dwnM____ = torch.reshape( mmmm( torch.reshape(FTK['c16_svd_U_d_expiw_s__'].to(dtype=torch.complex64),mtr((FTK['n_delta_v'],FTK['n_svd_l']))) , torch.reshape(svd_VUXM_lwnM____.to(dtype=torch.complex64),mtr((FTK['n_svd_l'],n_w_max*pm_n_UX_rank*n_M))) ) , mtr((FTK['n_delta_v'],n_w_max,pm_n_UX_rank,n_M)));
    UX_T_M_l2_dM__ = torch.reshape(torch.real(torch.sum(torch.reshape(torch.abs(torch.permute(UX_T_M_k_q_dwnM____,mtr(mts((1,2,0,3)))))**2,mtr((n_w_max*pm_n_UX_rank,FTK['n_delta_v'],n_M))),2-0) * n_w_max),mtr((FTK['n_delta_v'],n_M)));
    if flag_verbose > 0: print(f' %% [finished {str_thisfunction}]');
    return(UX_T_M_l2_dM__);
