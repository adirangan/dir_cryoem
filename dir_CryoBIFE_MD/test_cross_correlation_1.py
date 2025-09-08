import numpy as np
import torch
from scipy.special import jv

from EMPM.grids import PolarGrid, FourierImages
from EMPM.stacks import Templates, Images
from EMPM.microscopy import CTF
from EMPM.metadata import ViewingAngles
from EMPM.cross_correlation_likelihood import CrossCorrelationLikelihood, conform_ctf
from EMPM.util import (
    CrossCorrelationReturnType,
    FloatArrayType,
    Precision,
    to_torch,
)

def h2d(kd: float = 0.0):
    output = 4*np.pi**2 * (jv(0,kd) + jv(2,kd));
    return output;

def test_cross_correlation_PP():
    device = 'cpu';
    ####;
    flag_verbose = 1;
    if flag_verbose: print('flag_verbose: ',flag_verbose);
    n_viewings = 1;
    box_size = 2.0;
    n_pixels = 128;
    pixel_size = box_size / n_pixels;
    ####;
    precision = Precision.DOUBLE;
    (torch_float_type, torch_complex_type, _) = precision.get_dtypes(Precision.DOUBLE);
    radius_max = n_pixels / (2.0 * np.pi) * np.pi / 2.0;
    if flag_verbose>1: print('radius_max: ',radius_max);
    dist_radii = 0.5 / (2.0 * np.pi) * np.pi / 2.0;
    if flag_verbose>1: print('dist_radii: ',dist_radii);
    n_inplanes = n_pixels * 4;
    if flag_verbose>1: print('n_pixels: ',n_pixels);
    if flag_verbose>1: print('n_inplanes: ',n_inplanes);
    max_displacement = 0.08;
    n_displacements_x = 3;
    n_displacements_y = 5;
    ####;
    # define quadrature grid. ;
    ####;
    polar_grid = PolarGrid(
        radius_max = radius_max,
        dist_radii = dist_radii,
        n_inplanes = n_inplanes,
        uniform = True
    );
    k_p_r_max = radius_max;
    k_p_r_ = polar_grid.radius_shells;
    k_p_r_k_ = k_p_r_;
    n_k_p_r = polar_grid.n_shells;
    gamma_z_ = polar_grid.theta_shell;
    n_w_max = polar_grid.n_inplanes;
    k_p_r_wk_ = polar_grid.radius_points;
    k_p_w_wk_ = polar_grid.theta_points;
    weight_2d_wk_ = polar_grid.weight_points;
    assert(n_w_max==n_inplanes);
    if flag_verbose>1: print('k_p_r_.shape',k_p_r_.shape);
    if flag_verbose>1: print('k_p_r_ ',k_p_r_);
    if flag_verbose>1: print('gamma_z_.shape',gamma_z_.shape);
    if flag_verbose>1: print('gamma_z_ ',gamma_z_);
    if flag_verbose>1: print('n_w_max ',n_w_max);
    if flag_verbose>1: print('polar_grid.radius_points.shape',polar_grid.radius_points.shape);
    if flag_verbose>1: print('polar_grid.theta_points.shape',polar_grid.theta_points.shape);
    ####;
    # Here we just check the basic integral. ;
    ####;
    k_p_0_wk_ = k_p_r_wk_ * np.cos(k_p_w_wk_);
    k_p_1_wk_ = k_p_r_wk_ * np.sin(k_p_w_wk_);
    delta_S_0 = +0.05; delta_S_1 = -0.13; delta_S_fnorm = np.sqrt(delta_S_0**2 + delta_S_1**2);
    S_k_p_wk_ = torch.exp(torch.tensor(2*np.pi*1j* (k_p_0_wk_*delta_S_0 + k_p_1_wk_*delta_S_1) ));
    delta_S_S_fnorm = 0.0;
    tmp_S_S_kd = 2*np.pi*k_p_r_max*delta_S_S_fnorm;
    I_S_S_form = h2d(tmp_S_S_kd) / (2*np.pi)**2 * (np.pi * k_p_r_max**2);
    I_S_S_quad = torch.sum(np.conj(S_k_p_wk_) * S_k_p_wk_ * torch.tensor(weight_2d_wk_)) * (2*np.pi)**2 ; 
    I_S_S_quad = I_S_S_quad.numpy();
    I_S_S_errrel = np.abs(I_S_S_form - I_S_S_quad)/max(1e-12,np.abs(I_S_S_form));
    if flag_verbose>1: print(' I_S_S_form: ',I_S_S_form,' I_S_S_quad: ',I_S_S_quad,' I_S_S_errrel: ',I_S_S_errrel);
    delta_M_0 = -0.17; delta_M_1 = -0.03; delta_M_fnorm = np.sqrt(delta_M_0**2 + delta_M_1**2);
    M_k_p_wk_ = torch.exp(torch.tensor(2*np.pi*1j* (k_p_0_wk_*delta_M_0 + k_p_1_wk_*delta_M_1) ));
    delta_M_M_fnorm = 0.0;
    tmp_M_M_kd = 2*np.pi*k_p_r_max*delta_M_M_fnorm;
    I_M_M_form = h2d(tmp_M_M_kd) / (2*np.pi)**2 * (np.pi * k_p_r_max**2);
    I_M_M_quad = torch.sum(np.conj(M_k_p_wk_) * M_k_p_wk_ * torch.tensor(weight_2d_wk_)) * (2*np.pi)**2 ; 
    I_M_M_quad = I_M_M_quad.numpy();
    I_M_M_errrel = np.abs(I_M_M_form - I_M_M_quad)/max(1e-12,np.abs(I_M_M_form));
    if flag_verbose>1: print(' I_M_M_form: ',I_M_M_form,' I_M_M_quad: ',I_M_M_quad,' I_M_M_errrel: ',I_M_M_errrel);
    delta_S_M_fnorm = np.sqrt( (delta_S_0-delta_M_0)**2 + (delta_S_1-delta_M_1)**2 );
    tmp_S_M_kd = 2*np.pi*k_p_r_max*delta_S_M_fnorm;
    I_S_M_form = h2d(tmp_S_M_kd) / (2*np.pi)**2 * (np.pi * k_p_r_max**2);
    I_S_M_quad = torch.sum(np.conj(S_k_p_wk_) * M_k_p_wk_ * torch.tensor(weight_2d_wk_)) * (2*np.pi)**2 ; 
    I_S_M_quad = I_S_M_quad.numpy();
    I_S_M_errrel = np.abs(I_S_M_form - I_S_M_quad)/max(1e-12,np.abs(I_S_M_form));
    if flag_verbose>1: print(' I_S_M_form: ',I_S_M_form,' I_S_M_quad: ',I_S_M_quad,' I_S_M_errrel: ',I_S_M_errrel,' #<-- should be <1e-6');
    ####;
    # Now we try and recalculate this integral using the image and template structures. ;
    ####;
    n_viewing_S = 1;
    S_k_p_wkS___ = S_k_p_wk_.reshape(n_viewing_S,n_k_p_r,n_w_max);
    S_image = Images( fourier_data = FourierImages( images_fourier = S_k_p_wkS___ , polar_grid = polar_grid ) ) ;
    M_k_p_wkS___ = M_k_p_wk_.reshape(n_viewing_S,n_k_p_r,n_w_max);
    M_image = Images( fourier_data = FourierImages( images_fourier = M_k_p_wkS___ , polar_grid = polar_grid ) ) ;
    S_template = Templates( 
        fourier_data = FourierImages( images_fourier = S_k_p_wkS___ , polar_grid = polar_grid ) ,
        viewing_angles = ViewingAngles([0],[0],[0]) ,
    ) ;
    cc = CrossCorrelationLikelihood(
        templates = S_template,
        max_displacement = max_displacement,
        n_displacements_x = n_displacements_x,
        n_displacements_y = n_displacements_y,
        precision = precision,
        device = device,
        verbose = False
    ) ;
    if flag_verbose>2: print(cc.n_displacements);
    if flag_verbose>2: print(cc.x_displacements_expt_scale);
    if flag_verbose>2: print(cc.y_displacements_expt_scale);
    ctf_one_wk_ = np.ones(k_p_r_wk_.shape);
    ctf_one_wk1___ = ctf_one_wk_.reshape(1, n_k_p_r, n_w_max).astype(np.float64)
    ctf_one = CTF(
        polar_grid = polar_grid,
        box_size = box_size, # is this really in Angstroms ? ;
        anisotropy = True,
        ctf_descriptor = ctf_one_wk1___
    )
    torch_ctf_one_wk1___ = conform_ctf(to_torch(ctf_one.ctf,precision,device),ctf_one.anisotropy);
    res = cc._compute_cross_correlation_likelihood(
        device=torch.device(device),
        images_fourier = M_image.images_fourier,
        ctf=torch_ctf_one_wk1___,
        n_pixels_phys = n_pixels*n_pixels,
        n_templates_per_batch=n_viewing_S,
        n_images_per_batch=n_viewing_S,
        return_type=CrossCorrelationReturnType.FULL_TENSOR,
        return_integrated_likelihood=False,
    ) ;
    I_x_errrel_sum = 0.0; ntotal=0;
    for nw in range(n_w_max):
        gamma_z = (2*np.pi)*nw/max(1,n_w_max);
        cz = np.cos(+gamma_z);
        sz = np.sin(+gamma_z);
        for ny in range(n_displacements_y): #<-- may not be interchangeable with nx. ;
            for nx in range(n_displacements_x): #<-- may not be interchangeable with ny. ;
                nxy = nx + ny*n_displacements_x;
                delta_D_0 = cc.x_displacements_expt_scale[0+nxy];
                delta_D_1 = cc.y_displacements_expt_scale[0+nxy];
                if flag_verbose>2: print(' gamma_z: ',gamma_z,' delta_D_0: ',delta_D_0,' delta_D_1: ',delta_D_1);
                # Here we... ;
                # (a) Rotate delta_M_ by R(+gamma_z), implying we rotate M_k_p_ (as a function) by R(+gamma_z). ;
                # (b) Add delta_D_ to delta_RM_, implying we multiply RM_k_p_ (as a function) by exp(+2*pi*i*k_p_*delta_D_). ;
                delta_RM_0 = + cz * delta_M_0 - sz * delta_M_1 ; #<-- rotate delta_M by R(+gamma_z). ;
                delta_RM_1 = + sz * delta_M_0 + cz * delta_M_1 ; #<-- rotate delta_M by R(+gamma_z). ;
                delta_DRM_0 = delta_RM_0 + delta_D_0; delta_DRM_1 = delta_RM_1 + delta_D_1;
                delta_DRM_DRM_fnorm = 0.0;
                tmp_DRM_DRM_kd = 2*np.pi*k_p_r_max*delta_DRM_DRM_fnorm;
                I_DRM_DRM_form = h2d(tmp_DRM_DRM_kd) / (2*np.pi)**2 * (np.pi * k_p_r_max**2);
                delta_S_DRM_fnorm = np.sqrt( (delta_S_0-delta_DRM_0)**2 + (delta_S_1-delta_DRM_1)**2 );
                tmp_S_DRM_kd = 2*np.pi*k_p_r_max*delta_S_DRM_fnorm;
                I_S_DRM_form = h2d(tmp_S_DRM_kd) / (2*np.pi)**2 * (np.pi * k_p_r_max**2);
                I_x_form = I_S_DRM_form/max(1e-12,np.sqrt(I_S_S_form*I_DRM_DRM_form));
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_form: ',I_x_form);
                I_x_quad = res.cross_correlation_SMdw[:,:,0+nxy,0+nw].numpy();
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_quad: ',I_x_quad);
                I_x_errrel = np.abs(I_x_form - I_x_quad)/max(1e-12,np.abs(I_x_form));
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_errrel: ',I_x_errrel);
                I_x_errrel_sum += I_x_errrel; ntotal += 1;
                # Of course this should be equivalent to... ;
                # (a) subtract delta_D_ from delta_S_, implying we multiply S_k_p_ by exp(-2*pi*i*k_p_*delta_D_). ;
                # (b) Rotate delta_DS_ by R(-gamma_z), implying we rotate DS_k_p_ (as a function) by R(-gamma_z). ;
                delta_DS_0 = delta_S_0 - delta_D_0; delta_DS_1 = delta_S_1 - delta_D_1;
                delta_RDS_0 = + cz * delta_DS_0 + sz * delta_DS_1 ; #<-- rotate delta_DS by R(-gamma_z). ;
                delta_RDS_1 = - sz * delta_DS_0 + cz * delta_DS_1 ; #<-- rotate delta_DS by R(-gamma_z). ;
                delta_RDS_RDS_fnorm = 0.0;
                tmp_RDS_RDS_kd = 2*np.pi*k_p_r_max*delta_RDS_RDS_fnorm;
                I_RDS_RDS_form = h2d(tmp_RDS_RDS_kd) / (2*np.pi)**2 * (np.pi * k_p_r_max**2);
                delta_M_M_fnorm = 0.0;
                tmp_M_M_kd = 2*np.pi*k_p_r_max*delta_M_M_fnorm;
                I_M_M_form = h2d(tmp_M_M_kd) / (2*np.pi)**2 * (np.pi * k_p_r_max**2);
                delta_RDS_M_fnorm = np.sqrt( (delta_RDS_0-delta_M_0)**2 + (delta_RDS_1-delta_M_1)**2 );
                tmp_RDS_M_kd = 2*np.pi*k_p_r_max*delta_RDS_M_fnorm;
                I_RDS_M_form = h2d(tmp_RDS_M_kd) / (2*np.pi)**2 * (np.pi * k_p_r_max**2);
                I_x_form = I_RDS_M_form/max(1e-12,np.sqrt(I_RDS_RDS_form*I_M_M_form));
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_form: ',I_x_form);
                I_x_quad = res.cross_correlation_SMdw[:,:,0+nxy,0+nw].numpy();
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_quad: ',I_x_quad);
                I_x_errrel = np.abs(I_x_form - I_x_quad)/max(1e-12,np.abs(I_x_form));
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_errrel: ',I_x_errrel);
                I_x_errrel_sum += I_x_errrel; ntotal += 1;
    n_total = ntotal;
    I_x_errrel_avg = I_x_errrel_sum / max(1,n_total)
    if flag_verbose>0: print(' I_x_errrel_avg: ',I_x_errrel_avg,' #<-- should be <1e-6');

def planewave__bessel0_planewave_integral_0(
    k_p_r_max: float = 48.0 / (2.0 * np.pi),
    delta_S: float | FloatArrayType = 0.1,
    omega_S: float | FloatArrayType = 0.0,
    delta_M: float | FloatArrayType = 0.1,
    omega_M: float | FloatArrayType = 0.0,
    alpha : float = 1.0
    ):
    delta_S_0 = delta_S * np.cos(omega_S); delta_S_1 = delta_S * np.sin(omega_S);
    delta_M_0 = delta_M * np.cos(omega_M); delta_M_1 = delta_M * np.sin(omega_M);
    delta_T_0 = delta_M_0 - delta_S_0; delta_T_1 = delta_M_1 - delta_S_1;
    delta_T = np.sqrt(delta_T_0**2 + delta_T_1**2);
    omega_T = np.arctan2(delta_T_1,delta_T_0);
    a = alpha * k_p_r_max
    b = 2 * np.pi * k_p_r_max * delta_T
    c = np.maximum(a, b)
    d = np.minimum(a, b)
    I = 2 * np.pi * k_p_r_max ** 2 * (d * jv(-1, d) * jv(0, c) - c * jv(-1, c) * jv(0, d)) / np.maximum(1e-12, c ** 2 - d ** 2)
    return I

def planewave__bessel0_planewave_integral_1(
    k_p_r_max: float = 48.0 / (2.0 * np.pi),
    delta_S_0: float | FloatArrayType = 0.0,
    delta_S_1: float | FloatArrayType = 0.0,
    delta_M_0: float | FloatArrayType = 0.0,
    delta_M_1: float | FloatArrayType = 0.0,
    alpha : float = 1.0
    ):
    delta_T_0 = delta_M_0 - delta_S_0; delta_T_1 = delta_M_1 - delta_S_1;
    delta_T = np.sqrt(delta_T_0**2 + delta_T_1**2);
    omega_T = np.arctan2(delta_T_1,delta_T_0);
    a = alpha * k_p_r_max
    b = 2 * np.pi * k_p_r_max * delta_T
    c = np.maximum(a, b)
    d = np.minimum(a, b)
    I = 2 * np.pi * k_p_r_max ** 2 * (d * jv(-1, d) * jv(0, c) - c * jv(-1, c) * jv(0, d)) / np.maximum(1e-12, c ** 2 - d ** 2)
    return I

def test_cross_correlation_PJ0P():
    device = 'cpu';
    ####;
    flag_verbose = 1;
    if flag_verbose: print('flag_verbose: ',flag_verbose);
    n_viewings = 1;
    box_size = 2.0;
    n_pixels = 128;
    pixel_size = box_size / n_pixels;
    ####;
    precision = Precision.DOUBLE;
    (torch_float_type, torch_complex_type, _) = precision.get_dtypes(Precision.DOUBLE);
    radius_max = n_pixels / (2.0 * np.pi) * np.pi / 2.0;
    if flag_verbose>1: print('radius_max: ',radius_max);
    dist_radii = 0.5 / (2.0 * np.pi) * np.pi / 2.0;
    if flag_verbose>1: print('dist_radii: ',dist_radii);
    n_inplanes = n_pixels * 4;
    if flag_verbose>1: print('n_pixels: ',n_pixels);
    if flag_verbose>1: print('n_inplanes: ',n_inplanes);
    max_displacement = 0.08;
    n_displacements_x = 3;
    n_displacements_y = 5;
    ####;
    # define quadrature grid. ;
    ####;
    polar_grid = PolarGrid(
        radius_max = radius_max,
        dist_radii = dist_radii,
        n_inplanes = n_inplanes,
        uniform = True
    );
    k_p_r_max = radius_max;
    k_p_r_ = polar_grid.radius_shells;
    k_p_r_k_ = k_p_r_;
    n_k_p_r = polar_grid.n_shells;
    gamma_z_ = polar_grid.theta_shell;
    n_w_max = polar_grid.n_inplanes;
    k_p_r_wk_ = polar_grid.radius_points;
    k_p_w_wk_ = polar_grid.theta_points;
    weight_2d_wk_ = polar_grid.weight_points;
    assert(n_w_max==n_inplanes);
    if flag_verbose>1: print('k_p_r_.shape',k_p_r_.shape);
    if flag_verbose>1: print('k_p_r_ ',k_p_r_);
    if flag_verbose>1: print('gamma_z_.shape',gamma_z_.shape);
    if flag_verbose>1: print('gamma_z_ ',gamma_z_);
    if flag_verbose>1: print('n_w_max ',n_w_max);
    if flag_verbose>1: print('polar_grid.radius_points.shape',polar_grid.radius_points.shape);
    if flag_verbose>1: print('polar_grid.theta_points.shape',polar_grid.theta_points.shape);
    ####;
    # Here we just check the basic integral. ;
    ####;
    k_p_0_wk_ = k_p_r_wk_ * np.cos(k_p_w_wk_);
    k_p_1_wk_ = k_p_r_wk_ * np.sin(k_p_w_wk_);
    delta_S_0 = +0.05; delta_S_1 = -0.13; delta_S_fnorm = np.sqrt(delta_S_0**2 + delta_S_1**2);
    S_k_p_wk_ = torch.exp(torch.tensor(2*np.pi*1j* (k_p_0_wk_*delta_S_0 + k_p_1_wk_*delta_S_1) ));
    delta_S_S_fnorm = 0.0;
    tmp_S_S_kd = 2*np.pi*k_p_r_max*delta_S_S_fnorm;
    I_S_S_form = h2d(tmp_S_S_kd) / (2*np.pi)**2 * (np.pi * k_p_r_max**2);
    I_S_S_quad = torch.sum(np.conj(S_k_p_wk_) * S_k_p_wk_ * torch.tensor(weight_2d_wk_)) * (2*np.pi)**2 ; 
    I_S_S_quad = I_S_S_quad.numpy();
    I_S_S_errrel = np.abs(I_S_S_form - I_S_S_quad)/max(1e-12,np.abs(I_S_S_form));
    if flag_verbose>1: print(' I_S_S_form: ',I_S_S_form,' I_S_S_quad: ',I_S_S_quad,' I_S_S_errrel: ',I_S_S_errrel);
    delta_M_0 = -0.17; delta_M_1 = -0.03; delta_M_fnorm = np.sqrt(delta_M_0**2 + delta_M_1**2);
    M_k_p_wk_ = torch.exp(torch.tensor(2*np.pi*1j* (k_p_0_wk_*delta_M_0 + k_p_1_wk_*delta_M_1) ));
    delta_M_M_fnorm = 0.0;
    tmp_M_M_kd = 2*np.pi*k_p_r_max*delta_M_M_fnorm;
    I_M_M_form = h2d(tmp_M_M_kd) / (2*np.pi)**2 * (np.pi * k_p_r_max**2);
    I_M_M_quad = torch.sum(np.conj(M_k_p_wk_) * M_k_p_wk_ * torch.tensor(weight_2d_wk_)) * (2*np.pi)**2 ; 
    I_M_M_quad = I_M_M_quad.numpy();
    I_M_M_errrel = np.abs(I_M_M_form - I_M_M_quad)/max(1e-12,np.abs(I_M_M_form));
    if flag_verbose>1: print(' I_M_M_form: ',I_M_M_form,' I_M_M_quad: ',I_M_M_quad,' I_M_M_errrel: ',I_M_M_errrel);
    delta_S_M_fnorm = np.sqrt( (delta_S_0-delta_M_0)**2 + (delta_S_1-delta_M_1)**2 );
    tmp_S_M_kd = 2*np.pi*k_p_r_max*delta_S_M_fnorm;
    I_S_M_form = h2d(tmp_S_M_kd) / (2*np.pi)**2 * (np.pi * k_p_r_max**2);
    I_S_M_quad = torch.sum(np.conj(S_k_p_wk_) * M_k_p_wk_ * torch.tensor(weight_2d_wk_)) * (2*np.pi)**2 ; 
    I_S_M_quad = I_S_M_quad.numpy();
    I_S_M_errrel = np.abs(I_S_M_form - I_S_M_quad)/max(1e-12,np.abs(I_S_M_form));
    if flag_verbose>1: print(' I_S_M_form: ',I_S_M_form,' I_S_M_quad: ',I_S_M_quad,' I_S_M_errrel: ',I_S_M_errrel);
    ctf_alpha = 0.6;
    ctf_J0_wk_ = jv(0,ctf_alpha*k_p_r_wk_);
    CTF_S_k_p_wk_ = S_k_p_wk_ * ctf_J0_wk_ ;
    I_CTF_S_CTF_S_quad = torch.sum(np.conj(CTF_S_k_p_wk_) * CTF_S_k_p_wk_ * torch.tensor(weight_2d_wk_)) * (2*np.pi)**2 ; 
    I_CTF_S_CTF_S_quad = I_CTF_S_CTF_S_quad.numpy();
    I_CTF_S_M_form = planewave__bessel0_planewave_integral_1( k_p_r_max , delta_S_0,delta_S_1,delta_M_0,delta_M_1,ctf_alpha );
    I_CTF_S_M_quad = torch.sum(np.conj(CTF_S_k_p_wk_) * M_k_p_wk_ * torch.tensor(weight_2d_wk_)) * (2*np.pi)**2 ; 
    I_CTF_S_M_errrel = np.abs(I_CTF_S_M_form - I_CTF_S_M_quad)/max(1e-12,np.abs(I_CTF_S_M_form));
    if flag_verbose>1: print(' I_CTF_S_M_form: ',I_CTF_S_M_form,' I_CTF_S_M_quad: ',I_CTF_S_M_quad,' I_CTF_S_M_errrel: ',I_CTF_S_M_errrel,' #<-- should be <1e-6');
    ####;
    # Now we try and recalculate this integral using the image and template structures. ;
    ####;
    n_viewing_S = 1;
    S_k_p_wkS___ = S_k_p_wk_.reshape(n_viewing_S,n_k_p_r,n_w_max);
    S_image = Images( fourier_data = FourierImages( images_fourier = S_k_p_wkS___ , polar_grid = polar_grid ) ) ;
    M_k_p_wkS___ = M_k_p_wk_.reshape(n_viewing_S,n_k_p_r,n_w_max);
    M_image = Images( fourier_data = FourierImages( images_fourier = M_k_p_wkS___ , polar_grid = polar_grid ) ) ;
    S_template = Templates( 
        fourier_data = FourierImages( images_fourier = S_k_p_wkS___ , polar_grid = polar_grid ) ,
        viewing_angles = ViewingAngles([0],[0],[0]) ,
    ) ;
    cc = CrossCorrelationLikelihood(
        templates = S_template,
        max_displacement = max_displacement,
        n_displacements_x = n_displacements_x,
        n_displacements_y = n_displacements_y,
        precision = precision,
        device = device,
        verbose = False
    ) ;
    if flag_verbose>2: print(cc.n_displacements);
    if flag_verbose>2: print(cc.x_displacements_expt_scale);
    if flag_verbose>2: print(cc.y_displacements_expt_scale);
    ctf_J0_wk1___ = ctf_J0_wk_.reshape(1, n_k_p_r, n_w_max).astype(np.float64)
    ctf_J0 = CTF(
        polar_grid = polar_grid,
        box_size = box_size, # is this really in Angstroms ? ;
        anisotropy = True,
        ctf_descriptor = ctf_J0_wk1___
    )
    torch_ctf_J0_wk1___ = conform_ctf(to_torch(ctf_J0.ctf,precision,device),ctf_J0.anisotropy);
    res = cc._compute_cross_correlation_likelihood(
        device=torch.device(device),
        images_fourier = M_image.images_fourier,
        ctf=torch_ctf_J0_wk1___,
        n_pixels_phys = n_pixels*n_pixels,
        n_templates_per_batch=n_viewing_S,
        n_images_per_batch=n_viewing_S,
        return_type=CrossCorrelationReturnType.FULL_TENSOR,
        return_integrated_likelihood=False,
    ) ;
    I_x_errrel_sum = 0.0; ntotal = 0;
    #for nw in range(n_w_max);
    nw_min = min(120,n_w_max); nw_max = min(124,n_w_max);
    for nw in range(nw_min,nw_max): #<-- restricted range. ;
        gamma_z = (2*np.pi)*nw/max(1,n_w_max);
        cz = np.cos(+gamma_z);
        sz = np.sin(+gamma_z);
        for ny in range(n_displacements_y): #<-- may not be interchangeable with nx. ;
            for nx in range(n_displacements_x): #<-- may not be interchangeable with ny. ;
                nxy = nx + ny*n_displacements_x;
                delta_D_0 = cc.x_displacements_expt_scale[0+nxy];
                delta_D_1 = cc.y_displacements_expt_scale[0+nxy];
                if flag_verbose>2: print(' gamma_z: ',gamma_z,' delta_D_0: ',delta_D_0,' delta_D_1: ',delta_D_1);
                # Here we... ;
                # (a) Rotate delta_M_ by R(+gamma_z), implying we rotate M_k_p_ (as a function) by R(+gamma_z). ;
                # (b) Add delta_D_ to delta_RM_, implying we multiply RM_k_p_ (as a function) by exp(+2*pi*i*k_p_*delta_D_). ;
                delta_RM_0 = + cz * delta_M_0 - sz * delta_M_1 ; #<-- rotate delta_M by R(+gamma_z). ;
                delta_RM_1 = + sz * delta_M_0 + cz * delta_M_1 ; #<-- rotate delta_M by R(+gamma_z). ;
                delta_DRM_0 = delta_RM_0 + delta_D_0; delta_DRM_1 = delta_RM_1 + delta_D_1;
                delta_DRM_DRM_fnorm = 0.0;
                tmp_DRM_DRM_kd = 2*np.pi*k_p_r_max*delta_DRM_DRM_fnorm;
                I_DRM_DRM_form = h2d(tmp_DRM_DRM_kd) / (2*np.pi)**2 * (np.pi * k_p_r_max**2);
                I_CTF_S_DRM_form = planewave__bessel0_planewave_integral_1( k_p_r_max , delta_S_0,delta_S_1,delta_DRM_0,delta_DRM_1,ctf_alpha );
                I_x_form = I_CTF_S_DRM_form/max(1e-12,np.sqrt(I_CTF_S_CTF_S_quad*I_DRM_DRM_form));
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_form: ',I_x_form);
                I_x_quad = res.cross_correlation_SMdw[:,:,0+nxy,0+nw].numpy();
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_quad: ',I_x_quad);
                I_x_errrel = np.abs(I_x_form - I_x_quad)/max(1e-12,np.abs(I_x_form));
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_errrel: ',I_x_errrel);
                I_x_errrel_sum += I_x_errrel; ntotal += 1;
                # Of course this should be equivalent to... ;
                # (a) subtract delta_D_ from delta_S_, implying we multiply S_k_p_ by exp(-2*pi*i*k_p_*delta_D_). ;
                # (b) Rotate delta_DS_ by R(-gamma_z), implying we rotate DS_k_p_ (as a function) by R(-gamma_z). ;
                delta_DS_0 = delta_S_0 - delta_D_0; delta_DS_1 = delta_S_1 - delta_D_1;
                delta_RDS_0 = + cz * delta_DS_0 + sz * delta_DS_1 ; #<-- rotate delta_DS by R(-gamma_z). ;
                delta_RDS_1 = - sz * delta_DS_0 + cz * delta_DS_1 ; #<-- rotate delta_DS by R(-gamma_z). ;
                RDS_k_p_wk_ = torch.exp(torch.tensor(2*np.pi*1j* (k_p_0_wk_*delta_RDS_0.numpy() + k_p_1_wk_*delta_RDS_1.numpy()) )); #<-- this takes a long time. ;
                CTF_RDS_k_p_wk_ = RDS_k_p_wk_ * ctf_J0_wk_ ; #<-- this takes a long time. ;
                I_CTF_RDS_CTF_RDS_quad = torch.sum(np.conj(CTF_RDS_k_p_wk_) * CTF_RDS_k_p_wk_ * torch.tensor(weight_2d_wk_)) * (2*np.pi)**2 ; #<-- this takes a long time. ;
                I_CTF_RDS_CTF_RDS_quad = I_CTF_RDS_CTF_RDS_quad.numpy();
                delta_M_M_fnorm = 0.0;
                tmp_M_M_kd = 2*np.pi*k_p_r_max*delta_M_M_fnorm;
                I_M_M_form = h2d(tmp_M_M_kd) / (2*np.pi)**2 * (np.pi * k_p_r_max**2);
                I_CTF_RDS_M_form = planewave__bessel0_planewave_integral_1( k_p_r_max , delta_RDS_0,delta_RDS_1,delta_M_0,delta_M_1,ctf_alpha );
                I_x_form = I_CTF_RDS_M_form/max(1e-12,np.sqrt(I_CTF_RDS_CTF_RDS_quad*I_M_M_form));
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_form: ',I_x_form);
                I_x_quad = res.cross_correlation_SMdw[:,:,0+nxy,0+nw].numpy();
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_quad: ',I_x_quad);
                I_x_errrel = np.abs(I_x_form - I_x_quad)/max(1e-12,np.abs(I_x_form));
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_errrel: ',I_x_errrel);
                I_x_errrel_sum += I_x_errrel; ntotal += 1;
    n_total = ntotal;
    I_x_errrel_avg = I_x_errrel_sum / max(1,n_total)
    if flag_verbose>0: print(' I_x_errrel_avg: ',I_x_errrel_avg,' #<-- should be <1e-6');

def planewave_planar_planar_planewave(
    k_p_r_max: float = 48.0 / (2.0 * np.pi),
    phi_S: float | FloatArrayType = 0.0,
    delta_S: float | FloatArrayType = 0.0,
    omega_S: float | FloatArrayType = 0.0,
    phi_M: float | FloatArrayType = 0.0,
    delta_M: float | FloatArrayType = 0.0,
    omega_M: float | FloatArrayType = 0.0,
    ):
    # calculates the integral: ;
    # \int_{k=0}^{K} \int_{psi=0}^{2*pi} conj[CTF * S] * M dpsi * kdk ;
    # where: ;
    # CTF = 2*k*cos(psi - phi_S) ;
    # S   = exp( +i*2*pi*k*delta_S * cos(psi - omega_S) ) ;
    # M   = 2*k*cos(psi - phi_M) * exp( +i*2*pi*k*delta_M * cos(psi - omega_M) ) ;
    delta_M_0 = delta_M * np.cos(omega_M); delta_M_1 = delta_M * np.sin(omega_M);
    delta_S_0 = delta_S * np.cos(omega_S); delta_S_1 = delta_S * np.sin(omega_S);
    delta_T_0 = delta_M_0 - delta_S_0; delta_T_1 = delta_M_1 - delta_S_1;
    delta_T = np.sqrt(delta_T_0**2 + delta_T_1**2);
    omega_T = np.arctan2(delta_T_1,delta_T_0);
    phi_pos = phi_S + phi_M;
    phi_neg = phi_S - phi_M;
    t = 2*np.pi*delta_T;
    tK = t*k_p_r_max;
    b0 = jv(0,tK);
    b2 = jv(2,tK);
    b4 = jv(4,tK);
    mode_2 = -4*np.pi*np.cos(phi_pos - 2*omega_T) * k_p_r_max**4 * ( b2 + b4 )/6 ;
    mode_0 = +4*np.pi*np.cos(phi_neg) * k_p_r_max**4 * ( b0/4 + b2/6 - b4/12 ) ;
    output = mode_0 + mode_2 ;
    return output

def test_cross_correlation_PxxP():
    device = 'cpu';
    ####;
    flag_verbose = 1;
    if flag_verbose: print('flag_verbose: ',flag_verbose);
    n_viewings = 1;
    box_size = 2.0;
    n_pixels = 128;
    pixel_size = box_size / n_pixels;
    ####;
    precision = Precision.DOUBLE;
    (torch_float_type, torch_complex_type, _) = precision.get_dtypes(Precision.DOUBLE);
    radius_max = n_pixels / (2.0 * np.pi) * np.pi / 2.0;
    if flag_verbose>1: print('radius_max: ',radius_max);
    dist_radii = 0.5 / (2.0 * np.pi) * np.pi / 2.0;
    if flag_verbose>1: print('dist_radii: ',dist_radii);
    n_inplanes = n_pixels * 4;
    if flag_verbose>1: print('n_pixels: ',n_pixels);
    if flag_verbose>1: print('n_inplanes: ',n_inplanes);
    max_displacement = 0.08;
    n_displacements_x = 3;
    n_displacements_y = 5;
    ####;
    # define quadrature grid. ;
    ####;
    polar_grid = PolarGrid(
        radius_max = radius_max,
        dist_radii = dist_radii,
        n_inplanes = n_inplanes,
        uniform = True
    );
    k_p_r_max = radius_max;
    k_p_r_ = polar_grid.radius_shells;
    k_p_r_k_ = k_p_r_;
    n_k_p_r = polar_grid.n_shells;
    gamma_z_ = polar_grid.theta_shell;
    n_w_max = polar_grid.n_inplanes;
    k_p_r_wk_ = polar_grid.radius_points;
    k_p_w_wk_ = polar_grid.theta_points;
    weight_2d_wk_ = polar_grid.weight_points;
    assert(n_w_max==n_inplanes);
    if flag_verbose>1: print('k_p_r_.shape',k_p_r_.shape);
    if flag_verbose>1: print('k_p_r_ ',k_p_r_);
    if flag_verbose>1: print('gamma_z_.shape',gamma_z_.shape);
    if flag_verbose>1: print('gamma_z_ ',gamma_z_);
    if flag_verbose>1: print('n_w_max ',n_w_max);
    if flag_verbose>1: print('polar_grid.radius_points.shape',polar_grid.radius_points.shape);
    if flag_verbose>1: print('polar_grid.theta_points.shape',polar_grid.theta_points.shape);
    ####;
    # Here we just check the basic integral. ;
    ####;
    k_p_0_wk_ = k_p_r_wk_ * np.cos(k_p_w_wk_);
    k_p_1_wk_ = k_p_r_wk_ * np.sin(k_p_w_wk_);
    delta_S_0 = +0.05; delta_S_1 = -0.13; delta_S = np.sqrt(delta_S_0**2 + delta_S_1**2); omega_S = np.arctan2(delta_S_1,delta_S_0);
    S_k_p_wk_ = torch.exp(torch.tensor(2*np.pi*1j* (k_p_0_wk_*delta_S_0 + k_p_1_wk_*delta_S_1) ));
    phi_S = -np.pi/5.0;
    CTF_k_p_wk_ = 2.0*k_p_r_wk_*np.cos(k_p_w_wk_ - phi_S);
    CTF_S_k_p_wk_ = torch.tensor(CTF_k_p_wk_) * S_k_p_wk_;
    delta_M_0 = -0.17; delta_M_1 = -0.03; delta_M = np.sqrt(delta_M_0**2 + delta_M_1**2); omega_M = np.arctan2(delta_M_1,delta_M_0);
    M_pwv_k_p_wk_ = torch.exp(torch.tensor(2*np.pi*1j* (k_p_0_wk_*delta_M_0 + k_p_1_wk_*delta_M_1) )) ;
    phi_M = +np.pi/3.0;
    M_pla_k_p_wk_ = 2.0*k_p_r_wk_*np.cos(k_p_w_wk_ - phi_M);
    M_k_p_wk_ = torch.tensor(M_pla_k_p_wk_) * M_pwv_k_p_wk_;
    I_CTF_S_M_form = planewave_planar_planar_planewave(k_p_r_max,phi_S,delta_S,omega_S,phi_M,delta_M,omega_M);
    I_CTF_S_M_quad = torch.sum(np.conj(CTF_S_k_p_wk_) * M_k_p_wk_ * torch.tensor(weight_2d_wk_)) * (2*np.pi)**2 ; 
    I_CTF_S_M_quad = I_CTF_S_M_quad.numpy();
    I_CTF_S_M_errrel = np.abs(I_CTF_S_M_form - I_CTF_S_M_quad)/max(1e-12,np.abs(I_CTF_S_M_form));
    if flag_verbose>1: print(' I_CTF_S_M_form: ',I_CTF_S_M_form,' I_CTF_S_M_quad: ',I_CTF_S_M_quad,' I_CTF_S_M_errrel: ',I_CTF_S_M_errrel,' #<-- should be <1e-6');
    ####;
    # Now we try and recalculate this integral using the image and template structures. ;
    ####;
    n_viewing_S = 1;
    S_k_p_wkS___ = S_k_p_wk_.reshape(n_viewing_S,n_k_p_r,n_w_max);
    S_image = Images( fourier_data = FourierImages( images_fourier = S_k_p_wkS___ , polar_grid = polar_grid ) ) ;
    M_k_p_wkS___ = M_k_p_wk_.reshape(n_viewing_S,n_k_p_r,n_w_max);
    M_image = Images( fourier_data = FourierImages( images_fourier = M_k_p_wkS___ , polar_grid = polar_grid ) ) ;
    S_template = Templates( 
        fourier_data = FourierImages( images_fourier = S_k_p_wkS___ , polar_grid = polar_grid ) ,
        viewing_angles = ViewingAngles([0],[0],[0]) ,
    ) ;
    cc = CrossCorrelationLikelihood(
        templates = S_template,
        max_displacement = max_displacement,
        n_displacements_x = n_displacements_x,
        n_displacements_y = n_displacements_y,
        precision = precision,
        device = device,
        verbose = False
    ) ;
    CTF_k_p_wk1___ = CTF_k_p_wk_.reshape(1, n_k_p_r, n_w_max).astype(np.float64)
    CTF_k_p = CTF(
        polar_grid = polar_grid,
        box_size = box_size, # is this really in Angstroms ? ;
        anisotropy = True,
        ctf_descriptor = CTF_k_p_wk1___
    )
    torch_CTF_k_p_wk1___ = conform_ctf(to_torch(CTF_k_p.ctf,precision,device),CTF_k_p.anisotropy);
    res = cc._compute_cross_correlation_likelihood(
        device=torch.device(device),
        images_fourier = M_image.images_fourier,
        ctf=torch_CTF_k_p_wk1___,
        n_pixels_phys = n_pixels*n_pixels,
        n_templates_per_batch=n_viewing_S,
        n_images_per_batch=n_viewing_S,
        return_type=CrossCorrelationReturnType.FULL_TENSOR,
        return_integrated_likelihood=False,
    ) ;
    I_x_errrel_sum = 0.0; ntotal = 0;
    for nw in range(n_w_max):
        gamma_z = (2*np.pi)*nw/max(1,n_w_max);
        cz = np.cos(+gamma_z);
        sz = np.sin(+gamma_z);
        for ny in range(n_displacements_y): #<-- may not be interchangeable with nx. ;
            for nx in range(n_displacements_x): #<-- may not be interchangeable with ny. ;
                nxy = nx + ny*n_displacements_x;
                delta_D_0 = cc.x_displacements_expt_scale[0+nxy];
                delta_D_1 = cc.y_displacements_expt_scale[0+nxy];
                if flag_verbose>2: print(' gamma_z: ',gamma_z,' delta_D_0: ',delta_D_0,' delta_D_1: ',delta_D_1);
                # Here we... ;
                # (a) Rotate delta_M_ by R(+gamma_z), implying we rotate M_k_p_ (as a function) by R(+gamma_z). ;
                # (b) Add delta_D_ to delta_RM_, implying we multiply RM_k_p_ (as a function) by exp(+2*pi*i*k_p_*delta_D_). ;
                delta_RM_0 = + cz * delta_M_0 - sz * delta_M_1 ; #<-- rotate delta_M by R(+gamma_z). ;
                delta_RM_1 = + sz * delta_M_0 + cz * delta_M_1 ; #<-- rotate delta_M by R(+gamma_z). ;
                delta_DRM_0 = delta_RM_0 + delta_D_0; delta_DRM_1 = delta_RM_1 + delta_D_1;
                phi_RM = phi_M + gamma_z; #<-- rotate M_k_p_ by R(+gamma_z). ;
                delta_DRM = np.sqrt(delta_DRM_0**2 + delta_DRM_1**2); omega_DRM = np.arctan2(delta_DRM_1,delta_DRM_0);
                I_DRM_DRM_form = planewave_planar_planar_planewave(k_p_r_max,phi_RM,delta_DRM,omega_DRM,phi_RM,delta_DRM,omega_DRM);
                phi_RS = phi_S + gamma_z; #<-- also rotate CTF_k_p_ by R(+gamma_z). ;
                I_CTF_S_CTF_S_form = planewave_planar_planar_planewave(k_p_r_max,phi_RS,delta_S,omega_S,phi_RS,delta_S,omega_S);
                I_CTF_S_DRM_form = planewave_planar_planar_planewave(k_p_r_max,phi_RS,delta_S,omega_S,phi_RM,delta_DRM,omega_DRM);
                I_x_form = I_CTF_S_DRM_form/max(1e-12,np.sqrt(I_CTF_S_CTF_S_form*I_DRM_DRM_form));
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_form: ',I_x_form);
                I_x_quad = res.cross_correlation_SMdw[:,:,0+nxy,0+nw].numpy();
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_quad: ',I_x_quad);
                I_x_errrel = np.abs(I_x_form - I_x_quad)/max(1e-12,np.abs(I_x_form));
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_errrel: ',I_x_errrel);
                I_x_errrel_sum += I_x_errrel; ntotal += 1;
                # Of course this should be equivalent to... ;
                # (a) subtract delta_D_ from delta_S_, implying we multiply S_k_p_ by exp(-2*pi*i*k_p_*delta_D_). ;
                # (b) Rotate delta_DS_ by R(-gamma_z), implying we rotate DS_k_p_ (as a function) by R(-gamma_z). ;
                delta_DS_0 = delta_S_0 - delta_D_0; delta_DS_1 = delta_S_1 - delta_D_1;
                delta_RDS_0 = + cz * delta_DS_0 + sz * delta_DS_1 ; #<-- rotate delta_DS by R(-gamma_z). ;
                delta_RDS_1 = - sz * delta_DS_0 + cz * delta_DS_1 ; #<-- rotate delta_DS by R(-gamma_z). ;
                #phi_S = phi_S; #<-- do not rotate CTF_k_p_. ;
                delta_RDS = np.sqrt(delta_RDS_0**2 + delta_RDS_1**2); omega_RDS = np.arctan2(delta_RDS_1,delta_RDS_0);
                I_CTF_RDS_CTF_RDS_form = planewave_planar_planar_planewave(k_p_r_max,phi_S,delta_RDS,omega_RDS,phi_S,delta_RDS,omega_RDS);
                I_M_M_form = planewave_planar_planar_planewave(k_p_r_max,phi_M,delta_M,omega_M,phi_M,delta_M,omega_M);
                I_CTF_RDS_M_form = planewave_planar_planar_planewave(k_p_r_max,phi_S,delta_RDS,omega_RDS,phi_M,delta_M,omega_M);
                I_x_form = I_CTF_RDS_M_form/max(1e-12,np.sqrt(I_CTF_RDS_CTF_RDS_form*I_M_M_form));
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_form: ',I_x_form);
                I_x_quad = res.cross_correlation_SMdw[:,:,0+nxy,0+nw].numpy();
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_quad: ',I_x_quad);
                I_x_errrel = np.abs(I_x_form - I_x_quad)/max(1e-12,np.abs(I_x_form));
                if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_errrel: ',I_x_errrel);
                I_x_errrel_sum += I_x_errrel; ntotal += 1;
    n_total = ntotal;
    I_x_errrel_avg = I_x_errrel_sum / max(1,n_total)
    if flag_verbose>0: print(' I_x_errrel_avg: ',I_x_errrel_avg,' #<-- should be <1e-6');

if __name__ == '__main__':
    print('running test_cross_correlation_PP');
    test_cross_correlation_PP();
    print('running test_cross_correlation_PJ0P');
    test_cross_correlation_PJ0P();
    print('running test_cross_correlation_PxxP');
    test_cross_correlation_PxxP();
    print('returning');
