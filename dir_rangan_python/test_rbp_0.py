import numpy as np
from matlab_scalar_round import matlab_scalar_round
from fnorm_disp import fnorm_disp
from cg_rhs_2 import cg_rhs_2
from sample_sphere_7 import sample_sphere_7
from sample_shell_6 import sample_shell_6
from get_weight_2d_2 import get_weight_2d_2
from rotate_p_to_p_fftw import rotate_p_to_p_fftw
from scipy.sparse import csr_matrix

def test_rbp_0(
        k_eq_d_double = None,
        k_int = None,
        viewing_k_eq_d_double = None,
):
    str_thisfunction = 'test_rbp_0';
    flag_verbose=1;
    if (flag_verbose> 0): print(f" %% [entering {str_thisfunction}]");
    np.random.seed(0);
    if k_eq_d_double is None: k_eq_d_double = 1.0;
    if k_int is None: k_int = 48;
    if viewing_k_eq_d_double is None: viewing_k_eq_d_double = 1.0;

    if flag_verbose > 0:
        print(" %% %% %% %% %% %% %% %% ");
        print(f" %% setting k_eq_d_double = {k_eq_d_double:.6f}");
        print(f" %% setting k_int = {k_int} ");
        print(f" %% setting viewing_k_eq_d_double = {viewing_k_eq_d_double:.6f}");

    k_p_r_max = k_int/(2*np.pi); k_eq_d = k_eq_d_double/(2*np.pi); str_T_vs_L = 'T';
    flag_uniform_over_n_k_p_r = 1; flag_uniform_over_polar_a = 0;
    (
        n_qk,
        n_qk_csum_,
        k_p_r_qk_,
        k_p_azimu_b_qk_,
        k_p_polar_a_qk_,
        weight_3d_k_p_qk_,
        weight_shell_k_,
        n_k_p_r,
        k_p_r_,
        weight_3d_k_p_r_,
        k_c_0_qk_,
        k_c_1_qk_,
        k_c_2_qk_,
        J_node_,
        J_weight_,
        J_chebfun_,
        J_polyval_,
        n_polar_a_k_,
        polar_a_ka__,
        n_azimu_b_ka__,
    ) = sample_sphere_7(
        flag_verbose,
        k_p_r_max,
        k_eq_d,
        str_T_vs_L,
        flag_uniform_over_n_k_p_r,
        flag_uniform_over_polar_a,
        );

    if flag_verbose > 0:
        print(f" %% setting k_p_r_max = {k_p_r_max:.6f}");
        print(f" %% setting k_eq_d = {k_eq_d:.6f}");
        print(f" %% n_qk = {n_qk}");
        
    index_outer_shell_ = range(n_qk_csum_[n_k_p_r-1],n_qk_csum_[n_k_p_r-1+1]);
    n_outer_shell = len(index_outer_shell_);
    k_outer_shell_r_q_ = k_p_r_qk_[index_outer_shell_];
    assert(np.std(k_outer_shell_r_q_,ddof=1)< 1e-6);
    k_outer_shell_azimu_b_q_ = k_p_azimu_b_qk_[index_outer_shell_];
    k_outer_shell_polar_a_q_ = k_p_polar_a_qk_[index_outer_shell_];

    if flag_verbose> 2:
        np.set_printoptions(threshold=10);
        print('index_outer_shell_');print(index_outer_shell_);
        print('n_outer_shell');print(n_outer_shell);
        print('k_outer_shell_r_q_');print(k_outer_shell_r_q_);
        print('k_outer_shell_azimu_b_q_');print(k_outer_shell_azimu_b_q_);
        print('k_outer_shell_polar_a_q_');print(k_outer_shell_polar_a_q_);

    n_w_int = 1;
    l_max_upb = matlab_scalar_round(2*np.pi*k_p_r_max);
    n_w_max = n_w_int*2*(l_max_upb+1);
    template_k_eq_d = -1;
    n_w_0in_ = n_w_max*np.ones(n_k_p_r);
    (
        n_w_,
        weight_2d_k_p_r_,
        weight_2d_k_p_wk_,
        k_p_r_wk_,
        k_p_w_wk_,
        k_c_0_wk_,
        k_c_1_wk_,
    ) = get_weight_2d_2(
        flag_verbose,
        n_k_p_r,
        k_p_r_,
        k_p_r_max,
        template_k_eq_d,
        n_w_0in_,
        weight_3d_k_p_r_,
    );
    n_w_max = np.max(n_w_);
    n_w_sum = np.sum(n_w_);
    n_w_csum_ = np.cumsum([0] + n_w_);

    if flag_verbose> 2:
        np.set_printoptions(threshold=10);
        print('n_w_');print(n_w_);
        print('weight_2d_k_p_r_');print(weight_2d_k_p_r_);
        print('weight_2d_k_p_wk_');print(weight_2d_k_p_wk_);
        print('k_p_r_wk_');print(k_p_r_wk_);
        print('k_p_w_wk_');print(k_p_w_wk_);
        print('k_c_0_wk_');print(k_c_0_wk_);
        print('k_c_1_wk_');print(k_c_1_wk_);

    delta_a_c_3s__ = np.array([
        [ +1.5 , -0.5 , +0.3 ] ,
        [ -0.5 , -1.5 , +2.0 ] ,
        ]) / 2.0 / max(1e-12,k_p_r_max);
    n_source = delta_a_c_3s__.shape[0];

    if flag_verbose> 2:
        np.set_printoptions(threshold=10);
        print('delta_a_c_3s__');print(delta_a_c_3s__);
        print('n_source');print(n_source);

    a_k_p_form_ = np.zeros(n_qk,dtype=complex);
    for nsource in range(n_source):
        delta_a_c_ = delta_a_c_3s__[nsource,:]
        a_k_p_form_ += np.exp(1j*2*np.pi*(
            k_c_0_qk_*delta_a_c_[0] +
            k_c_1_qk_*delta_a_c_[1] +
            k_c_2_qk_*delta_a_c_[2]
        ));

    if flag_verbose> 2:
        np.set_printoptions(threshold=10);
        print('a_k_p_form_');print(a_k_p_form_);

    viewing_k_eq_d = viewing_k_eq_d_double*1.0/max(1e-12,k_p_r_max);
    flag_uniform_over_polar_a = 0;
    (
        n_viewing_S,
        viewing_azimu_b_S_,
        viewing_polar_a_S_,
        viewing_weight_S_,
        _,
        _,
        _,
        _,
        _,
        _,
    )= sample_shell_6(
        1.0,
        viewing_k_eq_d,
        'L',
        flag_uniform_over_polar_a,
    );
    n_S = n_viewing_S;
    if flag_verbose > 0:
        print(f"%% viewing_k_eq_d_double: {viewing_k_eq_d_double:.2f}, n_S: {n_S}");

    def R2(gamma_z):
        return np.array([
            [ +np.cos(gamma_z) , -np.sin(gamma_z) ],
            [ +np.sin(gamma_z) , +np.cos(gamma_z) ],
        ]);
    def Rz(azimu_b):
        return np.array([
            [+np.cos(azimu_b) , -np.sin(azimu_b) , 0 ],
            [+np.sin(azimu_b) , +np.cos(azimu_b) , 0 ],
            [               0 ,                0 , 1 ],
        ]);
    def Ry(polar_a):
        return np.array([
            [ +np.cos(polar_a) , 0 , +np.sin(polar_a) ],
            [                0 , 1 ,                0 ],
            [ -np.sin(polar_a) , 0 , +np.cos(polar_a) ],
        ]);

    S_k_p_wkS__ = np.zeros((n_S,n_w_sum),dtype=complex);
    for nS in range(n_S):
        tmp_azimu_b = viewing_azimu_b_S_[nS];
        tmp_polar_a = viewing_polar_a_S_[nS];
        tmp_gamma_z = 0.0;
        tmp_R__ = Rz(-tmp_gamma_z) @ Ry(-tmp_polar_a) @ Rz(-tmp_azimu_b);
        S_k_p_wk_ = np.zeros(n_w_sum,dtype=complex);
        for nsource in range(n_source):
            tmp_delta_ = tmp_R__ @ delta_a_c_3s__[nsource,:];
            S_k_p_wk_ += np.exp(1j*2*np.pi*(
                k_c_0_wk_*tmp_delta_[0] + k_c_1_wk_*tmp_delta_[1]
            ));
        S_k_p_wkS__[nS,:] = S_k_p_wk_.flatten();

    if flag_verbose> 2:
        np.set_printoptions(threshold=10);
        print('S_k_p_wkS__');print(S_k_p_wkS__);

    n_M = n_S*2
    index_nS_from_nM_ = np.maximum(0,np.minimum(n_S-1,np.floor(np.mod(np.arange(n_M),n_S)))).astype(int);
    euler_polar_a_M_ = viewing_polar_a_S_[index_nS_from_nM_];
    euler_azimu_b_M_ = viewing_azimu_b_S_[index_nS_from_nM_];
    euler_gamma_z_M_ = 2*np.pi*np.arange(n_M) / max(1,n_M) + np.pi/2.0;

    if flag_verbose> 2:
        np.set_printoptions(threshold=10);
        print('n_M');print(n_M);
        print('index_nS_from_nM_');print(index_nS_from_nM_);
        print('euler_polar_a_M_');print(euler_polar_a_M_);
        print('euler_azimu_b_M_');print(euler_azimu_b_M_);
        print('euler_gamma_z_M_');print(euler_gamma_z_M_);

    n_CTF = 8  # <-- some number of CTF-functions. ;
    CTF_phi_C_ = np.zeros(n_CTF);
    CTF_k_p_wkC__ = np.zeros((n_CTF,n_w_sum));
    for nCTF in range(n_CTF):
        CTF_phi = 2*np.pi*nCTF / max(1,n_CTF) + np.pi/2.0;
        CTF_phi_C_[nCTF] = CTF_phi;
        CTF_k_p_wk_ = 2*k_p_r_wk_ * np.cos(k_p_w_wk_-CTF_phi);
        CTF_k_p_wkC__[nCTF,:] = CTF_k_p_wk_;
    index_nCTF_from_nM_ = np.mod(np.arange(n_M),n_CTF).astype(int);

    if flag_verbose> 2:
        np.set_printoptions(threshold=10);
        print('n_CTF');print(n_CTF);
        print('index_nCTF_from_nM_');print(index_nCTF_from_nM_);
        print('CTF_phi_C_');print(CTF_phi_C_);
        print('CTF_k_p_wkC__');print(CTF_k_p_wkC__);

    M_k_p_wkM__ = np.zeros((n_M,n_w_sum),dtype=complex);
    for nM in range(n_M):
        nCTF = index_nCTF_from_nM_[nM];
        CTF_k_p_wk_ = CTF_k_p_wkC__[nCTF,:];
        nS = index_nS_from_nM_[nM];
        gamma_z = euler_gamma_z_M_[nM];
        S_k_p_wk_ = S_k_p_wkS__[nS,:];
        M_k_p_wk_ = CTF_k_p_wk_ * rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,+gamma_z);
        M_k_p_wkM__[nM,:] = M_k_p_wk_;

    if flag_verbose> 2:
        np.set_printoptions(threshold=10);
        print('M_k_p_wkM__');print(M_k_p_wkM__);

    CTF_k_p_wkM__ = CTF_k_p_wkC__[index_nCTF_from_nM_,:]
    M0C2_k_p_wkM__ = np.abs(CTF_k_p_wkM__)**2;
    M1C1_k_p_wkM__ = M_k_p_wkM__ * CTF_k_p_wkM__;
    
    M0C2_k_p_wMk__ = np.reshape(
        np.transpose(
            np.reshape(M0C2_k_p_wkM__,(n_M,n_k_p_r,n_w_max)),(1,0,2)
        ),
        (n_k_p_r,n_w_max*n_M)
    );
    
    M1C1_k_p_wMk__ = np.reshape(
        np.transpose(
            np.reshape(M1C1_k_p_wkM__,(n_M,n_k_p_r,n_w_max)),(1,0,2)
        ),
        (n_k_p_r,n_w_max*n_M,)
    );

    if flag_verbose> 2:
        np.set_printoptions(threshold=10);
        print('CTF_k_p_wkM__');print(CTF_k_p_wkM__);
        tmp_index_ = (1e6*np.arange(10)).astype(int);
        print('M0C2_k_p_wkM__[tmp_index_]');print(M0C2_k_p_wkM__.flatten()[tmp_index_]);
        print('M0C2_k_p_wMk__[tmp_index_]');print(M0C2_k_p_wMk__.flatten()[tmp_index_]);
        print('M1C1_k_p_wkM__[tmp_index_]');print(M1C1_k_p_wkM__.flatten()[tmp_index_]);
        print('M1C1_k_p_wMk__[tmp_index_]');print(M1C1_k_p_wMk__.flatten()[tmp_index_]);

    n_3 = 3;
    (
        k_p_polar_a_wM__,
        k_p_azimu_b_wM__,
        k_c_0_wM__,
        k_c_1_wM__,
        k_c_2_wM__,
        k_p_r01_wM__,
    )= cg_rhs_2(
        n_M,
        n_w_max,
        euler_polar_a_M_,
        euler_azimu_b_M_,
        euler_gamma_z_M_,
    );

    if flag_verbose> 2:
        np.set_printoptions(threshold=10);
        print('k_p_polar_a_wM__');print(k_p_polar_a_wM__);
        print('k_p_azimu_b_wM__');print(k_p_azimu_b_wM__);
        print('k_c_0_wM__');print(k_c_0_wM__);
        print('k_c_1_wM__');print(k_c_1_wM__);
        print('k_c_2_wM__');print(k_c_2_wM__);
        print('k_p_r01_wM__');print(k_p_r01_wM__);

    polar_a_a_ = polar_a_ka__[0];
    n_polar_a = len(polar_a_a_);
    assert np.std(np.diff(polar_a_a_),ddof=1) < 1e-12;
    dpolar_a = np.mean(np.diff(polar_a_a_)) #<-- will be negative. ;
    n_azimu_b_a_ = n_azimu_b_ka__[0];
    dazimu_b_a_ = (2*np.pi) / np.maximum(1,n_azimu_b_a_);
    n_azimu_b_csum_ = np.cumsum(np.concatenate(([0],n_azimu_b_a_)));
    assert n_outer_shell == n_azimu_b_csum_[n_polar_a];

    npolar_a_from_data_wM_ = np.reshape(
        np.clip(
            np.round((np.pi - k_p_polar_a_wM__) / max(1e-12,-dpolar_a)), #<-- the banker-round should not impact the results much. ;
            0,
            n_polar_a - 1,
        ),
        (n_w_max*n_M),
    );
        
    nazimu_b_from_data_wMa__ = np.zeros((n_polar_a,n_w_max*n_M));
    for npolar_a in range(n_polar_a):
        n_azimu_b = n_azimu_b_a_[npolar_a];
        dazimu_b = dazimu_b_a_[npolar_a];
        nazimu_b_from_data_wM_ = np.reshape(
            np.mod(
                np.round(k_p_azimu_b_wM__ / max(1e-12,dazimu_b)), #<-- the banker-round should not impact the results much. ;
                n_azimu_b,
            ),
            (n_w_max*n_M),
        );
        nazimu_b_from_data_wMa__[npolar_a,:] = nazimu_b_from_data_wM_.flatten();

    index_babin_shell_from_data_wM_ = np.zeros((n_w_max*n_M),dtype=int);
    for npolar_a in range(n_polar_a):
        tmp_index_ = np.where(npolar_a_from_data_wM_.flatten() == npolar_a)[0];
        index_babin_shell_from_data_wM_[tmp_index_] = (
            n_azimu_b_csum_[npolar_a]
            + nazimu_b_from_data_wMa__[npolar_a,tmp_index_]
        );

    babin_shell_from_data_qwM__ = csr_matrix(
        (
            np.ones(n_w_max*n_M),
            (np.arange(n_w_max*n_M),index_babin_shell_from_data_wM_),
        ),
        shape=(n_w_max*n_M,n_outer_shell),
    );

    if flag_verbose> 2:
        np.set_printoptions(threshold=10);
        print('polar_a_a_');print(polar_a_a_);
        print('n_polar_a');print(n_polar_a);
        print('dpolar_a');print(dpolar_a);
        print('n_azimu_b_a_');print(n_azimu_b_a_);
        print('dazimu_b_a_');print(dazimu_b_a_);
        print('n_azimu_b_csum_');print(n_azimu_b_csum_);
        print('npolar_a_from_data_wM_');print(npolar_a_from_data_wM_);
        print('nazimu_b_from_data_wMa__');print(nazimu_b_from_data_wMa__);
        print('index_babin_shell_from_data_wM_');print(index_babin_shell_from_data_wM_);
        print('babin_shell_from_data_qwM__');print(babin_shell_from_data_qwM__);

    a_k_p_babin_ = np.reshape((babin_shell_from_data_qwM__.T.dot(M1C1_k_p_wMk__.T).T)/np.maximum(1e-12,babin_shell_from_data_qwM__.T.dot(M0C2_k_p_wMk__.T).T),n_qk); #<-- all these non-conjugate transposes so that this matches the matlab code in line 307 of test_CryoLike_rbp_20250730.m. ;

    _=fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_babin_',a_k_p_babin_,' %<-- total error can be 1e-1 or so');
        
    if (flag_verbose> 0): print(f" %% [finished {str_thisfunction}]");

if __name__ == '__main__':
    print('running test_rbp_0');
    test_rbp_0();
    print('returning');



