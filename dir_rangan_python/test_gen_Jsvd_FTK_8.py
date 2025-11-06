import numpy as np ; pi = np.pi ; i = 1j ; import torch ; import timeit ;
from matlab_index_2d_0 import matlab_index_2d_0 ;
from matlab_index_3d_0 import matlab_index_3d_0 ;
from matlab_index_4d_0 import matlab_index_4d_0 ;
from matlab_scalar_round import matlab_scalar_round ;
from fnorm_disp import fnorm_disp ;
from get_r8_weight_3d_1 import get_r8_weight_3d_1 ;
from gen_Jsvd_FTK_8 import gen_Jsvd_FTK_8 ;
from get_r8_delta_2 import get_r8_delta_2 ;
from get_r8_svd_chebval_U_d_0 import get_r8_svd_chebval_U_d_0 ;
from get_r8_svd_chebval_V_r_0 import get_r8_svd_chebval_V_r_0 ;
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

flag_verbose=1;
print(f' %% testing gen_Jsvd_FTK_8');
r8_k_p_r_max = 48.0/(2*pi); r8_k_eq_d = 1.0/(2*pi); str_L = 'L';
(
    n_k_p_r,
    r8_k_p_r_,
    r8_weight_3d_k_p_r_,
) = get_r8_weight_3d_1(
    flag_verbose,
    r8_k_p_r_max,
    r8_k_eq_d,
    str_L,
);
pm_N_pixel = 1.0; r8_delta_r_max = pm_N_pixel / (2*pi*r8_k_p_r_max) * (pi*np.sqrt(2)) ;
r8_svd_eps = 1e-4; l_max = 24; n_a_degree = 64; n_b_degree = 65;
FTK = gen_Jsvd_FTK_8(r8_k_p_r_max,pm_N_pixel,r8_svd_eps,l_max,n_a_degree,n_b_degree);
n_delta_v_requested = 2*FTK['n_svd_l'];
FTK['n_delta_v'],FTK['r8_delta_x_'],FTK['r8_delta_y_'] = get_r8_delta_2(r8_delta_r_max,n_delta_v_requested);
FTK['r8_svd_d_max'] = r8_delta_r_max;
FTK['r8_svd_chebval_U_d_'] = get_r8_svd_chebval_U_d_0(FTK['r8_svd_d_max'],FTK['n_svd_d'],FTK['r8_svd_d_'],FTK['n_svd_l'],FTK['i4_svd_l_'],FTK['r8_svd_U_d_chebcoef_'],FTK['n_delta_v'],FTK['r8_delta_x_'],FTK['r8_delta_y_']).ravel();
assert(FTK['r8_svd_chebval_U_d_'].numel()==FTK['n_svd_l']*FTK['n_delta_v']);
FTK['r8_svd_r_max'] = 2*pi*r8_k_p_r_max;
FTK['r8_svd_chebval_V_r_'] = get_r8_svd_chebval_V_r_0(FTK['r8_svd_r_max'],FTK['n_svd_r'],FTK['r8_svd_r_'],FTK['n_svd_l'],FTK['i4_svd_l_'],FTK['r8_svd_V_r_chebcoef_'],n_k_p_r,2*pi*r8_k_p_r_).ravel();
assert(FTK['r8_svd_chebval_V_r_'].numel()==FTK['n_svd_l']*n_k_p_r);
FTK['c16_svd_expiw__'] = torch.reshape(torch.exp(-i*(pi/2 - torch.reshape(torch.atan2(FTK['r8_delta_y_'],FTK['r8_delta_x_']),mtr((FTK['n_delta_v'],1))))*torch.reshape(FTK['i4_svd_l_'],mtr((1,FTK['n_svd_l'])))),mtr((FTK['n_delta_v'],FTK['n_svd_l'])));
FTK['c16_svd_U_d_expiw_s__'] = torch.permute(torch.reshape(FTK['r8_svd_chebval_U_d_'],mtr((FTK['n_svd_l'],FTK['n_delta_v']))), mtr(mts((1,0)))) * mmmm( FTK['c16_svd_expiw__'] , torch.diagflat(FTK['r8_svd_s_']).to(dtype=torch.complex128) );

dir_base = '/data/rangan' ;
dir_ascii = dir_base + '/dir_cryoem/dir_rangan_python/dir_ascii' ;

#fname_ascii = dir_ascii + '/FTK_n_r_degree.ascii' 
#print(f' %% writing fname_ascii: {fname_ascii}');
#np.savetxt(fname_ascii,FTK['n_r_degree'].numpy().ravel()); 
#fname_ascii = dir_ascii + '/FTK_n_d_degree.ascii'
#print(f' %% writing fname_ascii: {fname_ascii}');
#np.savetxt(fname_ascii,FTK['n_d_degree'].numpy().ravel());
#fname_ascii = dir_ascii + '/FTK_n_svd_r.ascii'
#print(f' %% writing fname_ascii: {fname_ascii}');
#np.savetxt(fname_ascii,FTK['n_svd_r'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_svd_r_.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['r8_svd_r_'].numpy().ravel());
#fname_ascii = dir_ascii + '/FTK_svd_r_m.ascii'
#print(f' %% writing fname_ascii: {fname_ascii}');
#np.savetxt(fname_ascii,FTK['r8_svd_r_m'].numpy().ravel());
#fname_ascii = dir_ascii + '/FTK_svd_r_c.ascii'
#print(f' %% writing fname_ascii: {fname_ascii}');
#np.savetxt(fname_ascii,FTK['r8_svd_r_c'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_svd_r_w_.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['r8_svd_r_w_'].numpy().ravel());
#fname_ascii = dir_ascii + '/FTK_n_svd_d.ascii'
#print(f' %% writing fname_ascii: {fname_ascii}');
#np.savetxt(fname_ascii,FTK['n_svd_d'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_svd_d_.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['r8_svd_d_'].numpy().ravel());
#fname_ascii = dir_ascii + '/FTK_svd_d_m.ascii'
#print(f' %% writing fname_ascii: {fname_ascii}');
#np.savetxt(fname_ascii,FTK['r8_svd_d_m'].numpy().ravel());
#fname_ascii = dir_ascii + '/FTK_svd_d_c.ascii'
#print(f' %% writing fname_ascii: {fname_ascii}');
#np.savetxt(fname_ascii,FTK['r8_svd_d_c'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_svd_d_w_.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['r8_svd_d_w_'].numpy().ravel());
#fname_ascii = dir_ascii + '/FTK_n_svd_l.ascii'
#print(f' %% writing fname_ascii: {fname_ascii}');
#np.savetxt(fname_ascii,FTK['n_svd_l'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_svd_l_.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['i4_svd_l_'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_svd_U_d_jacocoef_.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['r8_svd_U_d_jacocoef_'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_svd_s_.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['r8_svd_s_'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_svd_V_r_jacocoef_.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['r8_svd_V_r_jacocoef_'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_svd_U_d_chebcoef_.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['r8_svd_U_d_chebcoef_'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_svd_V_r_chebcoef_.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['r8_svd_V_r_chebcoef_'].numpy().ravel());
#fname_ascii = dir_ascii + '/FTK_n_delta_v.ascii'
#print(f' %% writing fname_ascii: {fname_ascii}');
#np.savetxt(fname_ascii,FTK['n_delta_v'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_delta_x_.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['r8_delta_x_'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_delta_y_.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['r8_delta_y_'].numpy().ravel());
#fname_ascii = dir_ascii + '/FTK_svd_d_max.ascii'
#print(f' %% writing fname_ascii: {fname_ascii}');
#np.savetxt(fname_ascii,FTK['r8_svd_d_max'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_svd_chebval_U_d_.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['r8_svd_chebval_U_d_'].numpy().ravel());
#fname_ascii = dir_ascii + '/FTK_svd_r_max.ascii'
#print(f' %% writing fname_ascii: {fname_ascii}');
#np.savetxt(fname_ascii,FTK['r8_svd_r_max'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_svd_chebval_V_r_.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['r8_svd_chebval_V_r_'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_svd_expiw__.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['c16_svd_expiw__'].numpy().ravel());
fname_ascii = dir_ascii + '/FTK_svd_U_d_expiw_s__.ascii'
print(f' %% writing fname_ascii: {fname_ascii}');
np.savetxt(fname_ascii,FTK['c16_svd_U_d_expiw_s__'].numpy().ravel());



