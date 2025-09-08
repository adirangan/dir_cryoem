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

def test_cross_correlation_PxxP_from_a_k_p():
    device = 'cpu';
    ####;
    flag_verbose = 1;
    str_thisfunction = 'test_cross_correlation_PxxP_from_a_k_p';
    if flag_verbose>0: print('[entering ',str_thisfunction,']');
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
    n_w_sum = n_w_max*n_k_p_r;
    k_p_r_wk_ = polar_grid.radius_points;
    k_p_w_wk_ = polar_grid.theta_points;
    weight_2d_wk_ = polar_grid.weight_points;
    k_p_0_wk_ = k_p_r_wk_ * np.cos(k_p_w_wk_);
    k_p_1_wk_ = k_p_r_wk_ * np.sin(k_p_w_wk_);
    assert(n_w_max==n_inplanes);
    if flag_verbose>1: print('k_p_r_.shape',k_p_r_.shape);
    if flag_verbose>1: print('k_p_r_ ',k_p_r_);
    if flag_verbose>1: print('gamma_z_.shape',gamma_z_.shape);
    if flag_verbose>1: print('gamma_z_ ',gamma_z_);
    if flag_verbose>1: print('n_w_max ',n_w_max);
    if flag_verbose>1: print('n_w_sum ',n_w_sum);
    if flag_verbose>1: print('polar_grid.radius_points.shape',polar_grid.radius_points.shape);
    if flag_verbose>1: print('polar_grid.theta_points.shape',polar_grid.theta_points.shape);
    ####;
    # Here we define a_k_p, ;
    # from which we 'slice' a collection of templates. ;
    ####;
    delta_a_0 = +0.13; delta_a_1 = -0.11; delta_a_2 = +0.07; #<-- wave-vector associated with a_k_p. ;
    delta_a_ = np.array([delta_a_0,delta_a_1,delta_a_2],np.float64);
    def a_k_p(k_c_3wkS___):
        tmp_wkS__ = delta_a_0 * k_c_3wkS___.to(device)[:,:,0+0] + delta_a_1 * k_c_3wkS___.to(device)[:,:,0+1] + delta_a_2 * k_c_3wkS___.to(device)[:,:,0+2] ;
        S_k_p_wkS__ = torch.exp( 2*np.pi*1j* tmp_wkS__ ) ;
        return S_k_p_wkS__.to(torch_complex_type);
    viewing_polar_a_S_ = 1*np.pi*torch.tensor([ +0.28 , +0.09 , +0.72 ],dtype=torch.float64);
    viewing_azimu_b_S_ = 2*np.pi*torch.tensor([ +0.10 , +0.32 , +0.85 ],dtype=torch.float64);
    viewing_gamma_z_S_ = 2*np.pi*torch.tensor([ +0.71 , +0.14 , +0.48 ],dtype=torch.float64); #<-- Here we can test with arbitrary viewing_gamma_z_S_. ;
    #viewing_gamma_z_S_ = 2*np.pi*torch.tensor([ +0.00 , +0.00 , +0.00 ],dtype=torch.float64); #<-- in pm_template_2 and sph_template_3 we assume viewing_gamma_z_S_ is 0. ;
    n_viewing_S = len(viewing_polar_a_S_);
    viewing_angles = ViewingAngles(viewing_azimu_b_S_, viewing_polar_a_S_, viewing_gamma_z_S_);
    assert(viewing_angles.n_angles==n_viewing_S);
    S_k_p_wkS__ = Templates.generate_from_function(a_k_p, viewing_angles, polar_grid, precision=precision, use_cuda=False)
    if flag_verbose>1: print('S_k_p_wkS__.images_fourier.shape',S_k_p_wkS__.images_fourier.shape);
    ####;
    # Now we reconstitute each template. ;
    # The general formula used here is as follows. ;
    # let sa and ca be sin(polar_a) and cos(polar_a), respectively. ;
    # let sb and cb be sin(azimu_b) and cos(azimu_b), respectively. ;
    # let sc and cc be sin(gamma_z) and cos(gamma_z), respectively. ;
    # And rotation by azimu_b about the +z-axis is represented as: ;
    # Rz(azimu_b) = ;
    # [ +cb -sb 0 ] ;
    # [ +sb +cb 0 ] ;
    # [  0   0  1 ] ;
    # And rotation by polar_a about the +y-axis is represented as: ;
    # Ry(polar_a) = ;
    # [ +ca 0 +sa ] ;
    # [  0  1  0  ] ;
    # [ -sa 0 +ca ] ;
    # And rotation by gamma_z about the +z-axis is represented as: ;
    # Rz(gamma_z) = ;
    # [ +cc -sc 0 ] ;
    # [ +sc +cc 0 ] ;
    # [  0   0  1 ] ;
    # Which, collectively, implies that under the transform: ;
    # Rz(azimu_b) * Ry(polar_a) * Rz(gamma_z), ;
    # Which is the same as: ;
    # [ +cb -sb 0 ] [ +ca*cc -ca*sc +sa ]   [ +cb*ca*cc - sb*sc , -cb*ca*sc -sb*cc , +cb*sa ];
    # [ +sb +cb 0 ] [ +sc    +cc    0   ] = [ +sb*ca*cc + cb*sc , -sb*ca*sc +cb*cc , +sb*sa ];
    # [  0   0  1 ] [ -sa*cc +sa*sc +ca ]   [ -sa*cc            , +sa*sc           , +ca    ];
    # the point [1;0;0] is mapped to: ;
    # [ template_k_c_0 ; template_k_c_1 ; template_k_c_2 ] = [ +cb*ca*cc - sb*sc ; +sb*ca*cc + cb*sc ; -sa*cc ];
    ####;
    # Given the convention above, we associated the template-value ;
    # S_k_p_wkS__(1+na,1+nS) with the indices: ;
    # nS = template-index determining viewing angles. ;
    # na = multi-index nw + n_w_csum_(1+nk_p_r)
    # nw =angle psi = 2*pi*nw/max(1,n_w_(1+nk_p_r)), ;
    # nk_p_r = radius k = k_p_r_(1+nk_p_r). ;
    # with the location Rz(+psi)*[1;0;0] under the above map in 3-dimensional k_p_ space. ;
    ####;
    # Similarly, we can consider reconstructing the template S_k_p_wkS__(:,1+nS) ;
    # by first rotating the volume a_k_p by the transformation: ;
    # inv( Rz(azimu_b) * Ry(polar_a) * Rz(gamma_z) ) = Rz(-gamma_z) * Ry(-polar_a) * Rz(-azimu_b) ;
    # and then taking the equatorial-slice: (k,psi) = [k*cos(psi);k*sin(psi);0] in 3-dimensions. ;
    ####;
    # Thus, if: ;
    # a_k_p = exp(+2*pi*i * transpose(k_p_) * delta_a_ ) ;
    # we can see that: ;
    # R__ * a_k_p = exp(+2*pi*i * transpose((inv(R) * k_p_)) * delta_a_ ) ;
    # R__ * a_k_p = exp(+2*pi*i * transpose(k_p_) * R__ * delta_a_ ) ;
    # and thus we can reconstitute the template: ;
    # S_k_p_wkS__(:,1+nS) ;
    # by applying the transformation: ;
    # Rz(-gamma_z) * Ry(-polar_a) * Rz(-azimu_b) ;
    # to the wave-vector delta_a_, ;
    # and then restricting the result to the equatorial-plane (i.e., ignoring the final coordinate). ;
    ####;
    def Rz(azimu_b):
        R__ = np.zeros((3,3),np.float64);
        R__[0+0,0+0] = +np.cos(azimu_b) ; R__[0+0,0+1] = -np.sin(azimu_b) ; R__[0+0,0+2] = 0.0 ;
        R__[0+1,0+0] = +np.sin(azimu_b) ; R__[0+1,0+1] = +np.cos(azimu_b) ; R__[0+1,0+2] = 0.0 ;
        R__[0+2,0+0] = 0.0              ; R__[0+2,0+1] =  0.0             ; R__[0+2,0+2] = 1.0 ;
        return R__;
    def Ry(polar_a):
        R__ = np.zeros((3,3),np.float64);
        R__[0+0,0+0] = +np.cos(polar_a) ; R__[0+0,0+1] = 0.0 ; R__[0+0,0+2] = +np.sin(polar_a) ;
        R__[0+1,0+0] = 0.0              ; R__[0+1,0+1] = 1.0 ; R__[0+1,0+2] = 0.0              ;
        R__[0+2,0+0] = -np.sin(polar_a) ; R__[0+2,0+1] = 0.0 ; R__[0+2,0+2] = +np.cos(polar_a) ;
        return R__;
    delta_S_2S__ = np.zeros((n_viewing_S,2),np.float64);
    for nS in range(n_viewing_S):
        tmp_polar_a = viewing_polar_a_S_[0+nS]; tmp_polar_a = tmp_polar_a.item();
        tmp_azimu_b = viewing_azimu_b_S_[0+nS]; tmp_azimu_b = tmp_azimu_b.item();
        tmp_gamma_z = viewing_gamma_z_S_[0+nS]; tmp_gamma_z = tmp_gamma_z.item();
        if flag_verbose>1: print(' nS: ',nS,' tmp_polar_a: ',tmp_polar_a,' tmp_azimu_b: ',tmp_azimu_b,' tmp_gamma_z: ',tmp_gamma_z);
        Rz_z = Rz(-tmp_gamma_z); Ry_a = Ry(-tmp_polar_a); Rz_b = Rz(-tmp_azimu_b);
        tmp_R__ = Rz_z.dot(Ry_a.dot(Rz_b)) ;
        tmp_delta_a_ = tmp_R__.dot(delta_a_); tmp_delta_a_0 = tmp_delta_a_[0+0]; tmp_delta_a_1 = tmp_delta_a_[0+1];
        delta_S_2S__[0+nS,0+0] = tmp_delta_a_0; delta_S_2S__[0+nS,0+1] = tmp_delta_a_1;
        T_k_p_wk_ = torch.exp(torch.tensor(2*np.pi*1j* (k_p_0_wk_*tmp_delta_a_0 + k_p_1_wk_*tmp_delta_a_1) ));
        S_k_p_wk_ = torch.reshape(S_k_p_wkS__.images_fourier[0+nS,:,:],(1,n_w_sum));
        I_S_S_quad = torch.sum(np.conj(S_k_p_wk_) * S_k_p_wk_ * torch.tensor(weight_2d_wk_)) * (2*np.pi)**2 ; 
        I_T_T_quad = torch.sum(np.conj(T_k_p_wk_) * T_k_p_wk_ * torch.tensor(weight_2d_wk_)) * (2*np.pi)**2 ; 
        I_S_vs_T_quad = torch.sum(np.conj(S_k_p_wk_ - T_k_p_wk_) * (S_k_p_wk_ - T_k_p_wk_) * torch.tensor(weight_2d_wk_)) * (2*np.pi)**2 ; 
        I_S_vs_T_errrel_quad = np.abs(I_S_vs_T_quad)/max(1e-12,np.abs(I_S_vs_T_quad));
        if flag_verbose>0: print(' nS: ',nS,' I_S_S_quad: ',I_S_S_quad,' #<-- should be large');
        if flag_verbose>0: print(' nS: ',nS,' I_S_vs_T_errrel_quad: ',I_S_vs_T_errrel_quad,' #<-- should be small (e.g., <1e-6)');
    ####;
    # Now we define the CTF-function that we will use later. ;
    ####;
    phi_S = -np.pi/5.0;
    CTF_k_p_wk_ = 2.0*k_p_r_wk_*np.cos(k_p_w_wk_ - phi_S);
    ####;
    # Now we define the image we will use later. ;
    ####;
    delta_M_0 = -0.17; delta_M_1 = -0.03; delta_M = np.sqrt(delta_M_0**2 + delta_M_1**2); omega_M = np.arctan2(delta_M_1,delta_M_0);
    M_pwv_k_p_wk_ = torch.exp(torch.tensor(2*np.pi*1j* (k_p_0_wk_*delta_M_0 + k_p_1_wk_*delta_M_1) )) ;
    phi_M = +np.pi/3.0;
    M_pla_k_p_wk_ = 2.0*k_p_r_wk_*np.cos(k_p_w_wk_ - phi_M);
    M_k_p_wk_ = torch.tensor(M_pla_k_p_wk_) * M_pwv_k_p_wk_;
    ####;
    # Now we check the basic integral for each template. ;
    ####;
    for nS in range(n_viewing_S):
        delta_S_0 = delta_S_2S__[0+nS,0+0]; delta_S_1 = delta_S_2S__[0+nS,0+1];
        delta_S = np.sqrt(delta_S_0**2 + delta_S_1**2); omega_S = np.arctan2(delta_S_1,delta_S_0);
        S_k_p_wk_ = torch.exp(torch.tensor(2*np.pi*1j* (k_p_0_wk_*delta_S_0 + k_p_1_wk_*delta_S_1) ));
        CTF_S_k_p_wk_ = torch.tensor(CTF_k_p_wk_) * S_k_p_wk_;
        I_CTF_S_M_form = planewave_planar_planar_planewave(k_p_r_max,phi_S,delta_S,omega_S,phi_M,delta_M,omega_M);
        I_CTF_S_M_quad = torch.sum(np.conj(CTF_S_k_p_wk_) * M_k_p_wk_ * torch.tensor(weight_2d_wk_)) * (2*np.pi)**2 ; 
        I_CTF_S_M_quad = I_CTF_S_M_quad.numpy();
        I_CTF_S_M_errrel = np.abs(I_CTF_S_M_form - I_CTF_S_M_quad)/max(1e-12,np.abs(I_CTF_S_M_form));
        if flag_verbose>0: print(' nS: ',nS,' I_CTF_S_M_form: ',I_CTF_S_M_form,' I_CTF_S_M_quad: ',I_CTF_S_M_quad,' I_CTF_S_M_errrel: ',I_CTF_S_M_errrel,' #<-- should be <1e-6');
    ####;
    # Now we try and recalculate this integral using the image and template structures. ;
    ####;
    n_M = 1;
    M_k_p_wkM___ = M_k_p_wk_.reshape(n_M,n_k_p_r,n_w_max);
    M_image = Images( fourier_data = FourierImages( images_fourier = M_k_p_wkM___ , polar_grid = polar_grid ) ) ;
    cc = CrossCorrelationLikelihood(
        templates = S_k_p_wkS__,
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
    if flag_verbose>-1: print(' res.cross_correlation_SMdw.shape: ',res.cross_correlation_SMdw.shape); #<-- Why is this of size [n_w_max,n_xy,n_S,n_M] ? ;
    I_x_errrel_sum = 0.0; ntotal = 0;
    nM = 0;
    for nS in range(n_viewing_S):
        delta_S_0 = delta_S_2S__[0+nS,0+0]; delta_S_1 = delta_S_2S__[0+nS,0+1];
        delta_S = np.sqrt(delta_S_0**2 + delta_S_1**2); omega_S = np.arctan2(delta_S_1,delta_S_0);
        for nw in range(n_w_max):
            gamma_z = (2*np.pi)*nw/max(1,n_w_max);
            cz = np.cos(+gamma_z);
            sz = np.sin(+gamma_z);
            for ny in range(n_displacements_y): #<-- not interchangeable with nx. ;
                for nx in range(n_displacements_x): #<-- not interchangeable with ny. ;
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
                    I_x_quad = res.cross_correlation_SMdw[0+nM,0+nS,0+nxy,0+nw].numpy();
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
                    I_x_quad = res.cross_correlation_SMdw[0+nM,0+nS,0+nxy,0+nw].numpy();
                    if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_quad: ',I_x_quad);
                    I_x_errrel = np.abs(I_x_form - I_x_quad)/max(1e-12,np.abs(I_x_form));
                    if flag_verbose>2: print(' nw: ',nw,' nx: ',nx,' ny: ',ny,' I_x_errrel: ',I_x_errrel);
                    I_x_errrel_sum += I_x_errrel; ntotal += 1;
    n_total = ntotal;
    I_x_errrel_avg = I_x_errrel_sum / max(1,n_total)
    if flag_verbose>0: print(' I_x_errrel_avg: ',I_x_errrel_avg,' #<-- should be <1e-6');
    if flag_verbose>0: print('[finished ',str_thisfunction,']');
    return;

if __name__ == '__main__':
    print('running test_cross_correlation_PxxP_from_a_k_p');
    test_cross_correlation_PxxP_from_a_k_p();
    print('returning');
