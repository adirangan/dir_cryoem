import numpy as np ; import warnings ;

warnings.filterwarnings("ignore", message=".*__array_wrap__.*", category=DeprecationWarning)

from sample_shell_6 import sample_shell_6

def hshell(kd: float = 0.0):
    '''
    This function provides a formula for the integral of a plane-wave on a spherical shell in 3d.
    It can also be used to calculate the inner-product of two plane-waves on that same shell.

    More specifically, this function calculates the following integral:
    \[ \int_{\epolara=0}^{\epolara=\pi} \int_{\eazimub=0}^{\eazimub=2\pi} \fT(\vec{\kappa}) d\eazimub \sin(\epolara)d\epolara \]
    for a fixed value of $|\vec{\kappa}|$, all normalized so that the domain has area $4\pi$.
    In this integrand $\fT(\vec{\kappa})$ is a plane-wave:
    \[ \fT(\vec{\kappa}) = \exp(+i2\pi \vec{\kappa}\cdot\vec{\delta}_{T}), \]
    defined by the displacement-vector $\vec{\delta}_{T} \in \Real^{3}$.

    The output is:
    \[ (4\pi) \frac{ \sin(\kd) }{ \kd^{1} }, \]
    where the intermediate quantity $\kd$ is:
    \[ \kd = 2pi K \delta_{T}, \]
    where $\delta_{T}$ is the magnitude of the displacement $\vec{\delta}_{T}$.

    This function can also be used to calculate the following shell-limited inner-product (on the same shell specified by $|\vec{\kappa}|$ as above):
    \[ \int_{\epolara=0}^{\epolara=\pi} \int_{\eazimub=0}^{\eazimub=2\pi} \fS(\vec{\kappa})^{\dagger} \fM(\vec{\kappa}) d\eazimub \sin(\epolara)d\epolara \]
    where $\fS(\vec{\kappa})$ and $\fM(\vec{\kappa})$ are plane-waves defined via:
    \[ \fS(\vec{\kappa}) = \exp(+i2\pi \vec{\kappa}\cdot\vec{\delta}_{S}), \]
    \[ \fM(\vec{\kappa}) = \exp(+i2\pi \vec{\kappa}\cdot\vec{\delta}_{M}), \]
    each defined by their own displacement-vector.
    Here we simply multiply the two plane-waves together to get $\fT$ defined via:
    \[ \vec{\delta}_{T} & = & \vec{\delta}_{S} - \vec{\delta}_{M}. \]
    '''   
    output = (4*np.pi) * (np.sin(kd))/max(1e-12,kd)
    return output

def test_weight_shell_q_from_sample_shell_6(
        k_int: int = 48, #<-- highest frequency (2*pi*k_p_r_max), watch out for aliasing! 
        k_eq_d_double: np.float64 = 1.0, #<-- prefactor for k_eq_d, determines density of sampling in frequency-space. 
        flag_uniform_over_polar_a: int = 0, #<-- SamplingStrategy = SamplingStrategy.ADAPTIVE,
        n_frequency: int = 1+3*8, #<-- 8 plane-waves tested per unit-interval in k-space.
        frequency_max: np.float64 = 3.0, #<-- test delta_T in [0,3].
):
    '''
    This function tests the spherical-quadrature-grid on a single shell in 3-dimensions.
    '''   
    device = 'cpu'
    ####
    flag_verbose = 2
    if flag_verbose: print('flag_verbose: ',flag_verbose)
    ####
    k_p_r_max = k_int/(2.0*np.pi); k_eq_d = k_eq_d_double/(2*np.pi)
    ####
    # define quadrature grid. 
    ####
    (
        n_q,
        azimu_b_q_,
        polar_a_q_,
        weight_shell_q_,
        k_c_0_q_,
        k_c_1_q_,
        k_c_2_q_,
    ) = sample_shell_6(
        k_p_r_max = k_p_r_max,
        k_eq_d = k_eq_d,
        str_T_vs_L = 'L',
        flag_uniform_over_polar_a = flag_uniform_over_polar_a,
    )[:7]
    if flag_verbose>1: print('2*np.pi*k_p_r_max',2*np.pi*k_p_r_max,'k_p_r_max',k_p_r_max)
    if flag_verbose>1: print('n_q ',n_q)
    if flag_verbose>2: print('polar_a_q_.shape',polar_a_q_.shape)
    if flag_verbose>2: print('azimu_b_q_.shape',azimu_b_q_.shape)
    if flag_verbose>2: print('weight_shell_q_.shape',weight_shell_q_.shape)
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
        T_k_p_q_ = np.exp(2 * np.pi * 1j * (k_c_0_q_ * delta_T_0 + k_c_1_q_ * delta_T_1 + k_c_2_q_ * delta_T_2))
        tmp_T_kd = 2*np.pi*k_p_r_max*delta_T
        I_T_quad = np.sum(T_k_p_q_ * np.array(weight_shell_q_))
        if np.abs(tmp_T_kd)< 1e-6: I_T_form = (4*np.pi) * k_p_r_max**2
        if np.abs(tmp_T_kd)>=1e-6: I_T_form = hshell(tmp_T_kd) * k_p_r_max**2
        I_T_errrel = np.abs(I_T_form - I_T_quad)/max(1e-12,np.abs(I_T_form))
        if flag_verbose>1: print(' %%',' nfrequency ',nfrequency,'/',n_frequency,' frequency ',frequency)
        if flag_verbose>1: print(' %%',' I_T_form: ',I_T_form,' I_T_quad: ',I_T_quad,' I_T_errrel: ',I_T_errrel)
        I_T_errrel_max = max(I_T_errrel_max,I_T_errrel)
    if flag_verbose>1: print(' %%',' I_T_errrel_max: ',I_T_errrel_max)
    assert(np.abs(I_T_errrel_max) < 1e-2)

if __name__ == '__main__':
    print('running test_weight_shell_q_from_sample_shell_6')
    test_weight_shell_q_from_sample_shell_6()
    print('returning')
