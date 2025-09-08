import numpy as np
from fnorm_disp import fnorm_disp
from periodize import periodize

def cg_rhs_2(
        n_M = None,
        n_w = None,
        viewing_polar_a_M_ = None,
        viewing_azimu_b_M_ = None,
        viewing_gamma_z_M_ = None,
        flag_cleanup = None,
):
    str_thisfunction = 'cg_rhs_2';
    flag_verbose=0;
    ################################################################;
    if (flag_verbose>0): print(f" %% [entering {str_thisfunction}]");
    ################################################################;

    if viewing_gamma_z_M_ is None or len(viewing_gamma_z_M_) == 0: viewing_gamma_z_M_ = np.zeros(n_M);
    if flag_cleanup is None: flag_cleanup = 1;
    flag_1 = 1;
    flag_d = 0
    flag_dd = 0
    tolerance_cg_rhs = 1e-16;
    tolerance_pole = np.sqrt(tolerance_cg_rhs);
    tolerance_cg_rhs_upb = 1.0/np.sqrt(1e-16);

    if flag_1:
        sa_wM__ = np.tile(np.sin(viewing_polar_a_M_).reshape(n_M,1), (1,n_w))
        ca_wM__ = np.tile(np.cos(viewing_polar_a_M_).reshape(n_M,1), (1,n_w))
        sb_wM__ = np.tile(np.sin(viewing_azimu_b_M_).reshape(n_M,1), (1,n_w))
        cb_wM__ = np.tile(np.cos(viewing_azimu_b_M_).reshape(n_M,1), (1,n_w))
        viewing_gamma_z_M_ = viewing_gamma_z_M_.reshape(n_M,1)
        inplane_gamma_z_w_ = (2*np.pi*np.arange(n_w) / max(1,n_w)).reshape(1,n_w)
        combine_gamma_z_wM__ = np.tile(inplane_gamma_z_w_, (n_M,1)) - np.tile(viewing_gamma_z_M_, (1,n_w))
        sc_wM__ = np.sin(combine_gamma_z_wM__)
        cc_wM__ = np.cos(combine_gamma_z_wM__)
        k_c_0_wM__ = cb_wM__*ca_wM__*cc_wM__ - sb_wM__*sc_wM__
        k_c_1_wM__ = sb_wM__*ca_wM__*cc_wM__ + cb_wM__*sc_wM__
        k_c_2_wM__ = -sa_wM__*cc_wM__
        k_p_r01_bkp_wM__ = np.sqrt(k_c_0_wM__**2 + k_c_1_wM__**2)
        k_p_r01_wM__ = np.sqrt((ca_wM__*cc_wM__)**2 + sc_wM__**2)
        _ = fnorm_disp(flag_verbose,'k_p_r01_wM__',k_p_r01_wM__,'k_p_r01_bkp_wM__',k_p_r01_bkp_wM__,' %<-- should be <1e-12');
        k_p_polar_a_wM__ = np.arctan2(k_p_r01_wM__, k_c_2_wM__)
        k_p_azimu_b_wM__ = np.arctan2(k_c_1_wM__, k_c_0_wM__)
        for nM in range(n_M):
            if np.abs(periodize(viewing_polar_a_M_[nM] - np.pi / 2, -np.pi / 2, np.pi / 2)) < tolerance_pole:
                tmp_index_0_pos_ = np.where(k_c_0_wM__[nM,:] > tolerance_pole)[0]
                tmp_index_0_neg_ = np.where(k_c_0_wM__[nM,:] < -tolerance_pole)[0]
                flag_0_use = np.linalg.norm(k_c_0_wM__[nM,:]) >= 0.25
                tmp_index_1_pos_ = np.where(k_c_1_wM__[nM,:] > tolerance_pole)[0]
                tmp_index_1_neg_ = np.where(k_c_1_wM__[nM,:] < -tolerance_pole)[0]
                flag_1_use = np.linalg.norm(k_c_1_wM__[nM,:]) >= 0.25
                tmp_index_2_pos_ = np.where(np.abs(k_c_2_wM__[nM,:] - 1) < tolerance_pole)[0]
                tmp_index_2_neg_ = np.where(np.abs(k_c_2_wM__[nM,:] + 1) < tolerance_pole)[0]
                assert flag_0_use or flag_1_use
                if flag_0_use:
                    tmp_azimu_b_pos = np.mean(k_p_azimu_b_wM__[nM,tmp_index_0_pos_])
                    tmp_azimu_b_neg = np.mean(k_p_azimu_b_wM__[nM,tmp_index_0_neg_])
                else:
                    tmp_azimu_b_pos = np.mean(k_p_azimu_b_wM__[nM,tmp_index_1_pos_])
                    tmp_azimu_b_neg = np.mean(k_p_azimu_b_wM__[nM,tmp_index_1_neg_])
                    tmp_azimu_b_mid = periodize(0.5*(tmp_azimu_b_pos + tmp_azimu_b_neg), 0, 2*np.pi)
                    k_p_azimu_b_wM__[nM,tmp_index_2_pos_] = periodize(tmp_azimu_b_mid + 0*np.pi, 0, 2*np.pi)
                    k_p_azimu_b_wM__[nM,tmp_index_2_neg_] = periodize(tmp_azimu_b_mid + 1*np.pi, 0, 2*np.pi)

    if flag_1:
        if not np.isfinite(k_p_polar_a_wM__).all():
            print(f"Warning, k_p_polar_a_wM__ not finite in {str_thisfunction}")
        if not np.isfinite(k_p_azimu_b_wM__).all():
            print(f"Warning, k_p_azimu_b_wM__ not finite in {str_thisfunction}")
        if not np.isfinite(k_c_0_wM__).all():
            print(f"Warning, k_c_0_wM__ not finite in {str_thisfunction}")
        if not np.isfinite(k_c_1_wM__).all():
            print(f"Warning, k_c_1_wM__ not finite in {str_thisfunction}")
        if not np.isfinite(k_c_2_wM__).all():
            print(f"Warning, k_c_2_wM__ not finite in {str_thisfunction}")
        if not np.isfinite(k_p_r01_wM__).all():
            print(f"Warning, k_p_r01_wM__ not finite in {str_thisfunction}")

    if flag_1:
        return(
            k_p_polar_a_wM__,
            k_p_azimu_b_wM__,
            k_c_0_wM__,
            k_c_1_wM__,
            k_c_2_wM__,
            k_p_r01_wM__,
        );
    
    ################################################################;
    if (flag_verbose>0): print(f" %% [finished {str_thisfunction}]");
    ################################################################;
