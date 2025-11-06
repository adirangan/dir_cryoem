import numpy as np ; pi = np.pi ; i = 1j ; import torch ; import timeit ;
from matlab_index_2d_0 import matlab_index_2d_0 ;
from matlab_index_3d_0 import matlab_index_3d_0 ;
from matlab_index_4d_0 import matlab_index_4d_0 ;
from matlab_scalar_round import matlab_scalar_round ;
numel_unique = lambda a : np.unique(a.numpy().ravel()).size ;
mtr = lambda a : tuple(reversed(a)) ; #<-- matlab-arranged size (i.e., tuple(reversed(...))). ;
msr = lambda str : str[::-1] ; #<-- for einsum (i.e., string reversed (...)). ;
mts = lambda a : tuple(len(a) - x - 1 for x in a) ; #<-- for permute (i.e., tuple subtract (...)). ;
n_byte_per_float32 = 4;
n_byte_per_complex64 = 8;

def interp_p_to_q(
        n_r=None,
        n_w_=None,
        n_A=None,
        S_p_=None,
):
    flag_verbose = 0;
    str_thisfunction = 'interp_p_to_q' ;

    #%%%%%%%%;
    if (flag_verbose): print(f' %% [entering {str_thisfunction}]');
    #%%%%%%%%;

    n_S = int( int(S_p_.numel()) / int(n_A) );
    if S_p_.numel()!=n_A*n_S: print(f' %% Warning, n_A {n_A} n_S {n_S} in {str_thisfunction}');

    #%%%%%%%%%%%%%%%%;
    if numel_unique(n_w_) > 1:
        S_q_ = torch.zeros(n_A*n_S).to(dtype=torch.complex64);
        #%%%%%%%%;
        if n_S==1:
            #%%%%;
            #% single transformation on adaptive grid. ;
            #%%%%;
            n_w_max = int(n_w_[n_r-1].item());
            ic=0;
            for nr in range(n_r):
                n_w = int(n_w_[nr]).item();
                if (n_w>0):
                    tmp_index_ = int(ic) + torch.arange(n_w).to(dtype=torch.int32);
                    S_q_[tmp_index_] = torch.tensor(np.fft.fft(S_p_[tmp_index_].numpy())/np.maximum(1,np.sqrt(n_w))).to(dtype=torch.complex64);
                    ic = ic + n_w;
                #end;%if (n_w>0);
            #end;%for nr=0:n_r-1;
            assert(ic==n_A);
            #%%%%;
        #end;%if n_S==1;
        #%%%%%%%%;
        if n_S> 1:
            for nS in range(n_S):
                tmp_index_ = int(nS*n_A) + torch.arange(n_A).to(dtype=torch.int32);
                S_q_[tmp_index_] = interp_p_to_q(n_r,n_w_,n_A,S_p_[tmp_index_]);
            #end;%for nS=0:n_S-1;
        #end;%if n_S> 1;
        #%%%%%%%%;
    #end;%if (numel(unique(n_w_))> 1);
    #%%%%%%%%%%%%%%%%;

    #%%%%%%%%%%%%%%%%;
    if numel_unique(n_w_)==1:
        n_w = int(n_w_[0].item());
        S_q_ = torch.tensor(np.fft.fft(torch.reshape(S_p_,mtr((n_w,n_r,n_S))).numpy(),axis=2-0)/np.maximum(1,np.sqrt(n_w))).to(dtype=torch.complex64).ravel();
        assert(S_q_.numel()==S_p_.numel());
    #end;%if (numel(unique(n_w_))==1);
    #%%%%%%%%%%%%%%%%;

    #%%%%%%%%;
    if (flag_verbose): print(f' %% [finished {str_thisfunction}]');
    #%%%%%%%%;
    return(S_q_);


