import numpy as np ; pi = np.pi; i = 1j ; import torch ; from scipy.special import jv ; import warnings ;
from h2d import h2d ; from h3d import h3d ;
warnings.filterwarnings("ignore", message=".*__array_wrap__.*", category=DeprecationWarning)

from EMPM.grids import PolarGrid, SphereGrid
from EMPM.util import (
    CrossCorrelationReturnType,
    FloatArrayType,
    QuadratureType,
    SamplingStrategy,
    Precision,
    to_torch,
)

def test_weight_3d_k_p_qk(
        k_int: int = 48, #<-- highest frequency (2*pi*k_p_r_max), watch out for aliasing! 
        k_eq_d_double: np.float64 = 1.0, #<-- prefactor for k_eq_d, determines density of sampling in frequency-space. 
        sampling_strategy_dist_type: SamplingStrategy = SamplingStrategy.ADAPTIVE,
        sampling_strategy_azimuthal_sampling: SamplingStrategy = SamplingStrategy.ADAPTIVE,
        n_frequency: int = 1+3*8, #<-- 8 plane-waves tested per unit-interval in k-space.
        frequency_max: np.float64 = 3.0, #<-- test delta_T in [0,3].
):
    '''
    This function tests the spherical-quadrature-grid in 3-dimensions.
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
    sphere_grid = SphereGrid(
        radius_max = k_p_r_max,
        dist_eq = k_eq_d,
        dist_type = sampling_strategy_dist_type,
        azimuthal_sampling = sampling_strategy_azimuthal_sampling,
    ) ;
    k_p_r_ = sphere_grid.radius_shells ;
    k_p_r_k_ = k_p_r_ ;
    n_k_p_r = sphere_grid.n_shells ;
    weight_3d_k_p_r_ = sphere_grid.weights_radius_shells ;
    n_qk = sphere_grid.n_points ;
    k_p_r_qk_ = sphere_grid.radius_points ;
    k_p_polar_a_qk_ = sphere_grid.polar_points ;
    k_p_azimu_b_qk_ = sphere_grid.azimu_points ;
    weight_3d_k_p_qk_ = sphere_grid.weight_points ;
    if flag_verbose>1: print('2*pi*k_p_r_max',2*pi*k_p_r_max,'k_p_r_max',k_p_r_max) ;
    if flag_verbose>1: print('n_k_p_r ',n_k_p_r) ;
    if flag_verbose>1: print('k_p_r_.shape',k_p_r_.shape) ;
    if flag_verbose>2: print('weight_3d_k_p_r_.shape',weight_3d_k_p_r_.shape) ;
    if flag_verbose>1: print('n_qk ',n_qk) ;
    if flag_verbose>2: print('k_p_r_qk_.shape',k_p_r_qk_.shape) ;
    if flag_verbose>2: print('k_p_polar_a_qk_.shape',k_p_polar_a_qk_.shape) ;
    if flag_verbose>2: print('k_p_azimu_b_qk_.shape',k_p_azimu_b_qk_.shape) ;
    if flag_verbose>2: print('weight_3d_k_p_qk_.shape',weight_3d_k_p_qk_.shape) ;
    k_p_0_qk_ = k_p_r_qk_ * np.cos(k_p_azimu_b_qk_) * np.sin(k_p_polar_a_qk_) ;
    k_p_1_qk_ = k_p_r_qk_ * np.sin(k_p_azimu_b_qk_) * np.sin(k_p_polar_a_qk_) ;
    k_p_2_qk_ = k_p_r_qk_ * np.cos(k_p_polar_a_qk_) ;
    #### ;
    # Here we just check the basic integral for various plane-waves of varying frequency:  ;
    #### ;
    frequency_ = np.linspace(0, frequency_max, n_frequency).reshape(-1, 1) ;
    I_T_errrel_max = 0.0 ;
    for nfrequency in range(n_frequency):
        frequency = frequency_[nfrequency].item() ;
        delta_T = frequency ;
        tmp_polar_a = 1 * pi * np.random.rand() ;
        tmp_azimu_b = 2 * pi * np.random.rand() ;
        delta_T_0 = delta_T * np.cos(tmp_azimu_b) * np.sin(tmp_polar_a) ;
        delta_T_1 = delta_T * np.sin(tmp_azimu_b) * np.sin(tmp_polar_a) ;
        delta_T_2 = delta_T * np.cos(tmp_polar_a) ;
        T_k_p_qk_ = torch.exp(torch.tensor(2*pi*i* (k_p_0_qk_*delta_T_0 + k_p_1_qk_*delta_T_1 + k_p_2_qk_*delta_T_2) )) ;
        tmp_T_kd = 2*pi*k_p_r_max*delta_T ;
        # For now this test explicitly references the vector weight_3d_k_p_qk_. ;
        # In the future we should modify this test to call the actual function in EMPM ;
        # that integrates an arbitrary collection of function-values against the quadrature-weights.  ;
        I_T_quad = torch.sum(T_k_p_qk_ * torch.tensor(weight_3d_k_p_qk_)).item() ; 
        if np.abs(tmp_T_kd)< 1e-6: I_T_form = 4*pi/3 * k_p_r_max**3 ;
        if np.abs(tmp_T_kd)>=1e-6: I_T_form = h3d(tmp_T_kd).item() * k_p_r_max**3 ;
        I_T_errrel = np.abs(I_T_form - I_T_quad) / np.maximum(1e-12,np.abs(I_T_form)) ;
        if flag_verbose>1: print(' %%',' nfrequency ',nfrequency,'/',n_frequency,' frequency ',frequency) ;
        if flag_verbose>1: print(' %%',' I_T_form: ',I_T_form,' I_T_quad: ',I_T_quad,' I_T_errrel: ',I_T_errrel) ;
        I_T_errrel_max = max(I_T_errrel_max,I_T_errrel) ;
    #end;%for;
    if flag_verbose>1: print(' %%',' I_T_errrel_max: ',I_T_errrel_max) ;
    assert(np.abs(I_T_errrel_max) < 1e-2) ;

if __name__ == '__main__':
    print('running test_weight_3d_k_p_qk') ;
    test_weight_3d_k_p_qk() ;
    print('returning') ;
