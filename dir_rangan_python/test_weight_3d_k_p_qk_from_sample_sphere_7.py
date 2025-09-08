import numpy as np ; import warnings ;

warnings.filterwarnings("ignore", message=".*__array_wrap__.*", category=DeprecationWarning)

from sample_sphere_7 import sample_sphere_7

def h3d(kd: float = 0.0):
    '''
    This function provides a formula for the bandlimited integral of a plane-wave in 3d.
    It can also be used to calculate the bandlimited-inner-product of two plane-waves in 3d.

    More specifically, this function calculates the following integral:
    \[ \int_{k=0}^{k=K} \int_{\epolara=0}^{\epolara=\pi} \int_{\eazimub=0}^{\eazimub=2\pi} \fT(\vec{\kappa}) d\eazimub \sin(\epolara)d\epolara k^{2}\dk \]
    normalized so that the domain has size $4\pi/3$.
    In this integrand $\fT(\vec{\kappa})$ is a plane-wave:
    \[ \fT(\vec{\kappa}) = \exp(+i2\pi \vec{\kappa}\cdot\vec{\delta}_{T}), \]
    defined by the displacement-vector $\vec{\delta}_{T} \in \Real^{3}$.

    The output is:
    \[ (4\pi) \frac{ \sin(\kd) - \kd\cos(\kd) }{ \kd^{3} }, \]
    where the intermediate quantity $\kd$ is:
    \[ \kd = 2pi K \delta_{T}, \]
    where $\delta_{T}$ is the magnitude of the displacement $\vec{\delta}_{T}$.

    This function can also be used to calculate the following bandlimited inner-product:
    \[ \int_{k=0}^{k=K} \int_{\epolara=0}^{\epolara=\pi} \int_{\eazimub=0}^{\eazimub=2\pi} \fS(\vec{\kappa})^{\dagger} \fM(\vec{\kappa}) d\eazimub \sin(\epolara)d\epolara k^{2}\dk \]
    where $\fS(\vec{\kappa})$ and $\fM(\vec{\kappa})$ are plane-waves defined via:
    \[ \fS(\vec{\kappa}) = \exp(+i2\pi \vec{\kappa}\cdot\vec{\delta}_{S}), \]
    \[ \fM(\vec{\kappa}) = \exp(+i2\pi \vec{\kappa}\cdot\vec{\delta}_{M}), \]
    each defined by their own displacement-vector.
    Here we simply multiply the two plane-waves together to get $\fT$ defined via:
    \[ \vec{\delta}_{T} & = & \vec{\delta}_{S} - \vec{\delta}_{M}. \]
    '''   
    output = (4*np.pi) * (np.sin(kd) - kd*np.cos(kd))/max(1e-12,kd**3)
    return output

def test_weight_3d_k_p_qk_from_sample_sphere_7(
        k_int: int = 48, #<-- highest frequency (2*pi*k_p_r_max), watch out for aliasing! 
        k_eq_d_double: np.float64 = 1.0, #<-- prefactor for k_eq_d, determines density of sampling in frequency-space. 
        flag_uniform_over_n_k_p_r: int= 0, #<-- SamplingStrategy = SamplingStrategy.ADAPTIVE,
        flag_uniform_over_polar_a: int= 0, #<-- SamplingStrategy = SamplingStrategy.ADAPTIVE,
        n_frequency: int = 1+3*8, #<-- 8 plane-waves tested per unit-interval in k-space.
        frequency_max: np.float64 = 3.0, #<-- test delta_T in [0,3].
):
    '''
    This function tests the spherical-quadrature-grid in 3-dimensions.
    '''   
    device = 'cpu'
    ####
    flag_verbose = 2
    if flag_verbose: print('flag_verbose: ',flag_verbose)
    ####
    k_p_r_max = k_int/(2.0*np.pi); k_eq_d = k_eq_d_double/(2*np.pi); str_T_vs_L = 'L';
    ####
    # define quadrature grid. 
    ####
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
    ) = sample_sphere_7(
        0*flag_verbose,
        k_p_r_max,
        k_eq_d,
        str_T_vs_L,
        flag_uniform_over_n_k_p_r,
        flag_uniform_over_polar_a,
    )[:13]
    if flag_verbose>1: print('2*np.pi*k_p_r_max',2*np.pi*k_p_r_max,'k_p_r_max',k_p_r_max)
    if flag_verbose>1: print('n_k_p_r ',n_k_p_r)
    if flag_verbose>1: print('k_p_r_.shape',k_p_r_.shape)
    if flag_verbose>2: print('weight_3d_k_p_r_.shape',weight_3d_k_p_r_.shape)
    if flag_verbose>1: print('n_qk ',n_qk)
    ####
    # Here we just check the basic integral for various plane-waves of varying frequency: 
    ####
    frequency_ = np.linspace(0, frequency_max, n_frequency).reshape(-1, 1)
    I_T_errrel_max = 0.0
    for nfrequency in range(n_frequency):
        frequency = frequency_[nfrequency]
        delta_T = frequency
        tmp_polar_a = 1 * np.pi * np.random.rand()
        tmp_azimu_b = 2 * np.pi * np.random.rand()
        delta_T_0 = delta_T * np.cos(tmp_azimu_b) * np.sin(tmp_polar_a)
        delta_T_1 = delta_T * np.sin(tmp_azimu_b) * np.sin(tmp_polar_a)
        delta_T_2 = delta_T * np.cos(tmp_polar_a)
        T_k_p_qk_ = np.exp( 2*np.pi*1j* (k_c_0_qk_*delta_T_0 + k_c_1_qk_*delta_T_1 + k_c_2_qk_*delta_T_2) )
        tmp_T_kd = 2*np.pi*k_p_r_max*delta_T
        I_T_quad = np.sum(T_k_p_qk_ * np.array(weight_3d_k_p_qk_))
        if np.abs(tmp_T_kd)< 1e-6: I_T_form = (4*np.pi/3) * k_p_r_max**3
        if np.abs(tmp_T_kd)>=1e-6: I_T_form = h3d(tmp_T_kd) * k_p_r_max**3
        I_T_errrel = np.abs(I_T_form - I_T_quad)/max(1e-12,np.abs(I_T_form))
        if flag_verbose>1: print(' %%',' nfrequency ',nfrequency,'/',n_frequency,' frequency ',frequency)
        if flag_verbose>1: print(' %%',' I_T_form: ',I_T_form,' I_T_quad: ',I_T_quad,' I_T_errrel: ',I_T_errrel)
        I_T_errrel_max = max(I_T_errrel_max,I_T_errrel)
    if flag_verbose>1: print(' %%',' I_T_errrel_max: ',I_T_errrel_max)
    assert(np.abs(I_T_errrel_max) < 1e-2)

if __name__ == '__main__':
    print('running test_weight_3d_k_p_qk_from_sample_sphere_7')
    test_weight_3d_k_p_qk_from_sample_sphere_7()
    print('returning')
