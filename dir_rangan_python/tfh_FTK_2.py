import numpy as np ; pi = np.pi ; i = 1j ; import torch ; import timeit ;
from matlab_index_2d_0 import matlab_index_2d_0 ;
from matlab_index_3d_0 import matlab_index_3d_0 ;
from matlab_index_4d_0 import matlab_index_4d_0 ;
from matlab_scalar_round import matlab_scalar_round ;
from gen_Jsvd_FTK_8 import gen_Jsvd_FTK_8 ;
from get_r8_delta_2 import get_r8_delta_2 ;
from get_r8_svd_chebval_U_d_0 import get_r8_svd_chebval_U_d_0 ;
from get_r8_svd_chebval_V_r_0 import get_r8_svd_chebval_V_r_0 ;
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
n_1 = int(1); n_2 = int(2); n_3 = int(3);
n_byte_per_float32 = 4;
n_byte_per_complex64 = 8;

def ampmh_FTK_2(
        n_k_p_r=None,
        r8_k_p_r_=None,
        r8_k_p_r_max=None,
        r8_delta_r_max=None,
        r8_svd_eps=None,
        n_delta_v_requested=None,
        r8_delta_x_0in_=None,
        r8_delta_y_0in_=None,
        l_max=None,
        n_a_degree=None,
        n_b_degree=None,
):
    flag_verbose=0;
    str_thisfunction = 'ampmh_FTK_2' ;

    #%%%%%%%%;
    if (flag_verbose): print(f' %% [entering {str_thisfunction}]');
    #%%%%%%%%;

    if r8_delta_r_max is None: r8_delta_r_max = 0.0;
    if r8_svd_eps is None: r8_svd_eps = 1e-4;
    if l_max is None: l_max = int(24);
    if n_a_degree is None: n_a_degree = int(64);
    if n_b_degree is None: n_b_degree = int(65);

    pm_N_pixel = r8_delta_r_max * (2*pi*r8_k_p_r_max) / (pi*np.sqrt(2)) ; 
    FTK = gen_Jsvd_FTK_8(r8_k_p_r_max,pm_N_pixel,r8_svd_eps,l_max,n_a_degree,n_b_degree);
    #%%%%;
    if n_delta_v_requested is None: n_delta_v_requested = 2*FTK['n_svd_l'];
    if (n_delta_v_requested<=0): n_delta_v_requested = 2*FTK['n_svd_l'];
    #%%%%;
    if r8_delta_r_max==0: FTK['n_delta_v'],FTK['r8_delta_x_'],FTK['r8_delta_y_'] = get_r8_delta_2(r8_delta_r_max,1);
    #%%%%;
    if r8_delta_r_max> 0:
        if (r8_delta_x_0in_ is not None) and (r8_delta_y_0in_ is not None):
            assert(r8_delta_x_0in_.numel()==n_delta_v_requested);
            assert(r8_delta_y_0in_.numel()==n_delta_v_requested);
            FTK['n_delta_v'] = n_delta_v_requested;
            FTK['r8_delta_x_'] = r8_delta_x_0in_.ravel();
            FTK['r8_delta_y_'] = r8_delta_y_0in_.ravel();
        #end;%if ~isempty(r8_delta_x_0in_) & ~isempty(r8_delta_y_0in_);
        if (r8_delta_x_0in_ is     None)  or (r8_delta_y_0in_ is     None):
            FTK['n_delta_v'],FTK['r8_delta_x_'],FTK['r8_delta_y_'] = get_r8_delta_2(r8_delta_r_max,n_delta_v_requested);
        #end;%if  isempty(r8_delta_x_0in_) |  isempty(r8_delta_y_0in_);
    #end;%if r8_delta_r_max> 0;
    #%%%%;

    FTK['r8_svd_d_max'] = r8_delta_r_max;
    FTK['r8_svd_chebval_U_d_'] = get_r8_svd_chebval_U_d_0(FTK['r8_svd_d_max'],FTK['n_svd_d'],FTK['r8_svd_d_'],FTK['n_svd_l'],FTK['i4_svd_l_'],FTK['r8_svd_U_d_chebcoef_'],FTK['n_delta_v'],FTK['r8_delta_x_'],FTK['r8_delta_y_']).ravel();
    assert(FTK['r8_svd_chebval_U_d_'].numel()==FTK['n_svd_l']*FTK['n_delta_v']);
    FTK['r8_svd_r_max'] = 2*pi*r8_k_p_r_max;
    FTK['r8_svd_chebval_V_r_'] = get_r8_svd_chebval_V_r_0(FTK['r8_svd_r_max'],FTK['n_svd_r'],FTK['r8_svd_r_'],FTK['n_svd_l'],FTK['i4_svd_l_'],FTK['r8_svd_V_r_chebcoef_'],n_k_p_r,2*pi*r8_k_p_r_).ravel();
    assert(FTK['r8_svd_chebval_V_r_'].numel()==FTK['n_svd_l']*n_k_p_r);
    FTK['c16_svd_expiw__'] = torch.reshape(torch.exp(-i*(pi/2 - torch.reshape(torch.atan2(FTK['r8_delta_y_'],FTK['r8_delta_x_']),mtr((FTK['n_delta_v'],1))))*torch.reshape(FTK['i4_svd_l_'],mtr((1,FTK['n_svd_l'])))),mtr((FTK['n_delta_v'],FTK['n_svd_l'])));
    FTK['c16_svd_U_d_expiw_s__'] = torch.permute(torch.reshape(FTK['r8_svd_chebval_U_d_'],mtr((FTK['n_svd_l'],FTK['n_delta_v']))), mtr(mts((1,0)))) * mmmm( FTK['c16_svd_expiw__'] , torch.diagflat(FTK['r8_svd_s_']).to(dtype=torch.complex128) );

    #%%%%%%%%;
    if (flag_verbose): print(f' %% [finished {str_thisfunction}]');
    #%%%%%%%%;
    return(FTK);

