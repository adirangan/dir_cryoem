# M-x reb-change-syntax --> string ;
# M-x query-replace-regexp ;
# (\([^(),]*\),\([^(),]*\),\([^(),]*\)) ;
# [\3,\2,\1] ;

import numpy as np ; pi = np.pi ; import torch ;
from ylgndr_2 import ylgndr_2 ;
from fnorm_disp import fnorm_disp ;

flag_verbose=2; flag_disp=0; nf=0;
print(' %% testing ylgndr_2')
l_max = 4; n_x = 25;
x_ = torch.linspace(-0.9,+0.9,n_x).reshape((n_x,));
x_.numel()
(
    _,
    tmp_y_ylgndr_jlm___,
    _,
    _,
    _,
    _,
    _,
    _,
    _,
) = ylgndr_2(
    [],
    l_max,
    x_,
);
dir_base = '/data/rangan' ;
dir_ascii = dir_base + '/dir_cryoem/dir_rangan_python/dir_ascii' ;
fname_ascii = dir_ascii + '/tmp_y_ylgndr_jlm___.ascii' ;
print(f' %% writing fname_ascii: {fname_ascii}')
np.savetxt(fname_ascii,tmp_y_ylgndr_jlm___.numpy().ravel());

l_max = 4; n_x = 25; x_ = torch.linspace(-0.9,+0.9,n_x).reshape((n_x,));
x_ = torch.cos(torch.acos(x_));
(
    _,
    d0y_mid_jlm___,
    d1y_mid_jlm___,
    d2y_mid_jlm___,
    sqrt_2lp1_,
    sqrt_2mp1_,
    sqrt_rat0_m_,
    sqrt_rat3_lm__,
    sqrt_rat4_lm__,
) = ylgndr_2(
    [],
    l_max,
    x_,
);
fname_ascii = dir_ascii + '/d0y_mid_jlm___.ascii' ;
print(f' %% writing fname_ascii: {fname_ascii}')
np.savetxt(fname_ascii,d0y_mid_jlm___.numpy().ravel());
fname_ascii = dir_ascii + '/d1y_mid_jlm___.ascii' ;
print(f' %% writing fname_ascii: {fname_ascii}')
np.savetxt(fname_ascii,d1y_mid_jlm___.numpy().ravel());
fname_ascii = dir_ascii + '/d2y_mid_jlm___.ascii' ;
print(f' %% writing fname_ascii: {fname_ascii}')
np.savetxt(fname_ascii,d2y_mid_jlm___.numpy().ravel());
fname_ascii = dir_ascii + '/sqrt_2lp1_.ascii' ;
print(f' %% writing fname_ascii: {fname_ascii}')
np.savetxt(fname_ascii,sqrt_2lp1_.numpy().ravel());
fname_ascii = dir_ascii + '/sqrt_2mp1_.ascii' ;
print(f' %% writing fname_ascii: {fname_ascii}')
np.savetxt(fname_ascii,sqrt_2mp1_.numpy().ravel());
fname_ascii = dir_ascii + '/sqrt_rat0_m_.ascii' ;
print(f' %% writing fname_ascii: {fname_ascii}')
np.savetxt(fname_ascii,sqrt_rat0_m_.numpy().ravel());
fname_ascii = dir_ascii + '/sqrt_rat3_lm__.ascii' ;
print(f' %% writing fname_ascii: {fname_ascii}')
np.savetxt(fname_ascii,sqrt_rat3_lm__.numpy().ravel());
fname_ascii = dir_ascii + '/sqrt_rat4_lm__.ascii' ;
print(f' %% writing fname_ascii: {fname_ascii}')
np.savetxt(fname_ascii,sqrt_rat4_lm__.numpy().ravel());

