import numpy as np ; pi = np.pi; i = 1j ; import torch ;
from scipy.special import jv ;
from hshell import hshell ;

from EMPM.grids import SphereShell
from EMPM.util import (
    CrossCorrelationReturnType,
    FloatArrayType,
    QuadratureType,
    SamplingStrategy,
    Precision,
    to_torch,
)

def test_weight_shell_q(
        k_int: int = 48, #<-- highest frequency (2*pi*k_p_r_max), watch out for aliasing! 
        k_eq_d_double: np.float64 = 1.0, #<-- prefactor for k_eq_d, determines density of sampling in frequency-space. 
        sampling_strategy_azimuthal_sampling: SamplingStrategy = SamplingStrategy.ADAPTIVE,
        n_frequency: int = 1+3*8, #<-- 8 plane-waves tested per unit-interval in k-space.
        frequency_max: np.float64 = 3.0, #<-- test delta_T in [0,3].
):
    '''
    This function tests the spherical-quadrature-grid on a single shell in 3-dimensions.
    '''   
    device = 'cpu' ;
    #### ;
    flag_verbose = 2 ;
    if flag_verbose: print('flag_verbose: ',flag_verbose) ;
    precision = Precision.DOUBLE ;
    (torch_float_type, torch_complex_type, _) = precision.get_dtypes(Precision.DOUBLE) ;
    #### ;
    k_p_r_max = k_int/(2.0*pi); k_eq_d = k_eq_d_double/(2*pi) ;
    #### ;
    # define quadrature grid.  ;
    #### ;
    shell_grid = SphereShell(
        radius = k_p_r_max,
        dist_eq = k_eq_d,
        azimuthal_sampling = sampling_strategy_azimuthal_sampling,
    ) ;
    n_q = shell_grid.n_points ;
    k_p_polar_a_q_ = torch.tensor(shell_grid.polar_points).to(dtype=torch.float32) ;
    k_p_azimu_b_q_ = torch.tensor(shell_grid.azimu_points).to(dtype=torch.float32) ;
    weight_shell_q_ = torch.tensor(shell_grid.weight_points).to(dtype=torch.float32) ;
    if flag_verbose>1: print('2*pi*k_p_r_max',2*pi*k_p_r_max,'k_p_r_max',k_p_r_max) ;
    if flag_verbose>1: print('n_q ',n_q) ;
    if flag_verbose>2: print('k_p_polar_a_q_.shape',k_p_polar_a_q_.shape) ;
    if flag_verbose>2: print('k_p_azimu_b_q_.shape',k_p_azimu_b_q_.shape) ;
    if flag_verbose>2: print('weight_shell_q_.shape',weight_shell_q_.shape) ;
    k_p_0_q_ = k_p_r_max * torch.cos(k_p_azimu_b_q_) * torch.sin(k_p_polar_a_q_) ;
    k_p_1_q_ = k_p_r_max * torch.sin(k_p_azimu_b_q_) * torch.sin(k_p_polar_a_q_) ;
    k_p_2_q_ = k_p_r_max * torch.cos(k_p_polar_a_q_) ;
    #### ;
    # Here we just check the basic integral for various plane-waves of varying frequency:  ;
    #### ;
    frequency_ = torch.linspace(0, frequency_max, n_frequency).to(dtype=torch.float32) ;
    I_T_errrel_max = 0.0 ;
    for nfrequency in range(n_frequency):
        frequency = frequency_[nfrequency].item() ;
        delta_T = frequency ;
        tmp_polar_a = 1 * pi * np.random.rand() ;
        tmp_azimu_b = 2 * pi * np.random.rand() ;
        delta_T_0 = delta_T * np.cos(tmp_azimu_b) * np.sin(tmp_polar_a) ;
        delta_T_1 = delta_T * np.sin(tmp_azimu_b) * np.sin(tmp_polar_a) ;
        delta_T_2 = delta_T * np.cos(tmp_polar_a) ;
        T_k_p_q_ = torch.exp(2*pi*i* (k_p_0_q_*delta_T_0 + k_p_1_q_*delta_T_1 + k_p_2_q_*delta_T_2) ) ;
        tmp_T_kd = 2*pi*k_p_r_max*delta_T ;
        # For now this test explicitly references the vector weight_shell_q_. ;
        # In the future we should modify this test to call the actual function in EMPM ;
        # that integrates an arbitrary collection of function-values against the quadrature-weights.  ;
        I_T_quad = torch.sum(T_k_p_q_ * weight_shell_q_).item() ; 
        if np.abs(tmp_T_kd)< 1e-6: I_T_form = (4*pi) * k_p_r_max**2 ;
        if np.abs(tmp_T_kd)>=1e-6: I_T_form = (hshell(tmp_T_kd) * k_p_r_max**2).item() ;
        I_T_errrel = np.abs(I_T_form - I_T_quad) / np.maximum(1e-12,np.abs(I_T_form)) ;
        if flag_verbose>1: print(' %%',' nfrequency ',nfrequency,'/',n_frequency,' frequency ',frequency) ;
        if flag_verbose>1: print(' %%',' I_T_form: ',I_T_form,' I_T_quad: ',I_T_quad,' I_T_errrel: ',I_T_errrel) ;
        I_T_errrel_max = np.maximum(I_T_errrel_max,I_T_errrel) ;
    if flag_verbose>1: print(' %%',' I_T_errrel_max: ',I_T_errrel_max) ;
    assert(np.abs(I_T_errrel_max) < 1e-2) ;

if __name__ == '__main__':
    print('running test_weight_shell_q') ;
    test_weight_shell_q() ;
    print('returning') ;
