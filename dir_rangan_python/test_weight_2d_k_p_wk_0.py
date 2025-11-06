import numpy as np ; pi = np.pi ; i = 1j ; import torch ;
from scipy.special import jv ; import warnings ;
from h2d import h2d ;

from EMPM.grids import PolarGrid
from EMPM.util import (
    CrossCorrelationReturnType,
    FloatArrayType,
    QuadratureType,
    Precision,
    to_torch,
)

def test_weight_2d_k_p_wk(
        k_int: int = 48, #<-- highest frequency (2*pi*k_p_r_max), watch out for aliasing! 
        k_eq_d_double: np.float64 = 0.5, #<-- prefactor for k_eq_d, determines density of sampling in frequency-space. 
        n_w_int: int = 2, #<-- prefactor for n_w_max, determines the number of distinct angles (i.e., n_gamma_z) used in frequency-space 2d-polar-grid. 
        quadrature_requested: QuadratureType = QuadratureType.GAUSS_JACOBI_BETA_1, #<-- quadrature style.
        n_frequency: int = 1+3*8, #<-- 8 plane-waves tested per unit-interval in k-space.
        frequency_max: np.float64 = 3.0, #<-- test delta_T in [0,3].
):
    '''
    This function tests the polar-quadrature-grid in 2-dimensions.
    Here we assume that the grid is `non-adaptive', in the sense that the each image-ring has the same number of points (i.e., the angular-discretization does not depend on $k$).
    '''   
    device = 'cpu' ;
    #### ;
    flag_verbose = 2 ;
    if flag_verbose: print('flag_verbose: ',flag_verbose) ;
    precision = Precision.DOUBLE ;
    (torch_float_type, torch_complex_type, _) = precision.get_dtypes(Precision.DOUBLE) ;
    #### ;
    k_p_r_max = k_int/(2.0*pi); k_eq_d = k_eq_d_double/(2*pi) ;
    l_max_upb = np.round(2*pi*k_p_r_max) ;
    n_w_max = n_w_int*2*(l_max_upb+1) ;
    #### ;
    # define quadrature grid.  ;
    #### ;
    polar_grid = PolarGrid(
        radius_max = k_p_r_max,
        dist_radii = k_eq_d,
        n_inplanes = n_w_max,
        uniform = True,
        quadrature = quadrature_requested,
    ) ;
    k_p_r_ = polar_grid.radius_shells ;
    k_p_r_k_ = k_p_r_ ;
    n_k_p_r = polar_grid.n_shells ;
    gamma_z_ = polar_grid.theta_shell ;
    n_w_max = polar_grid.n_inplanes ;
    k_p_r_wk_ = polar_grid.radius_points ;
    k_p_w_wk_ = polar_grid.theta_points ;
    weight_2d_k_p_wk_ = polar_grid.weight_points ;
    if flag_verbose>1: print('2*pi*k_p_r_max',2*pi*k_p_r_max,'k_p_r_max',k_p_r_max) ;
    if flag_verbose>1: print('k_p_r_.shape',k_p_r_.shape) ;
    if flag_verbose>2: print('k_p_r_ ',k_p_r_) ;
    if flag_verbose>2: print('np.max(k_p_r_)',np.max(k_p_r_)) ;
    if flag_verbose>2: print('gamma_z_.shape',gamma_z_.shape) ;
    if flag_verbose>2: print('gamma_z_ ',gamma_z_) ;
    if flag_verbose>1: print('n_w_max ',n_w_max) ;
    if flag_verbose>2: print('polar_grid.radius_points.shape',polar_grid.radius_points.shape) ;
    if flag_verbose>2: print('polar_grid.theta_points.shape',polar_grid.theta_points.shape) ;
    k_p_0_wk_ = k_p_r_wk_ * np.cos(k_p_w_wk_) ;
    k_p_1_wk_ = k_p_r_wk_ * np.sin(k_p_w_wk_) ;
    #### ;
    # Here we just check the basic integral for various plane-waves of varying frequency:  ;
    #### ;
    frequency_ = np.linspace(0, frequency_max, n_frequency).reshape(-1, 1) ;
    I_T_errrel_max = 0.0 ;
    for nfrequency in range(n_frequency):
        frequency = frequency_[nfrequency] ;
        delta_T = frequency ;
        omega = 2 * pi * np.random.rand() ;
        delta_T_0 = delta_T * np.cos(omega) ;
        delta_T_1 = delta_T * np.sin(omega) ;
        T_k_p_wk_ = torch.exp(torch.tensor(2*pi*i* (k_p_0_wk_*delta_T_0 + k_p_1_wk_*delta_T_1) )) ;
        tmp_T_kd = 2*pi*k_p_r_max*delta_T ;
        # For now this test explicitly references the vector weight_2d_k_p_wk_. ;
        # In the future we should modify this test to call the actual function in EMPM ;
        # that integrates an arbitrary collection of function-values against the quadrature-weights.  ;
        I_T_quad = torch.sum(T_k_p_wk_ * torch.tensor(weight_2d_k_p_wk_)).item() * (2*pi)**2 ; 
        I_T_form = h2d(tmp_T_kd).item() / (2*pi)**2 * (pi * k_p_r_max**2) ;
        I_T_errrel = np.abs(I_T_form - I_T_quad) / np.maximum(1e-12,np.abs(I_T_form)) ;
        if flag_verbose>1: print(' %%',' nfrequency ',nfrequency,'/',n_frequency,' frequency ',frequency) ;
        if flag_verbose>1: print(' %%',' I_T_form: ',I_T_form,' I_T_quad: ',I_T_quad,' I_T_errrel: ',I_T_errrel) ;
        I_T_errrel_max = np.maximum(I_T_errrel_max,I_T_errrel) ;
    #end;%for;
    if flag_verbose>1: print(' %%',' I_T_errrel_max: ',I_T_errrel_max) ;
    assert(np.abs(I_T_errrel_max) < 1e-2) ;

if __name__ == '__main__':
    print('running test_weight_2d_k_p_wk') ;
    test_weight_2d_k_p_wk() ;
    print('returning') ;
