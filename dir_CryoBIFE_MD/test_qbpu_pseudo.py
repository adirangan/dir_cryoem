# This is pseudo code for the cartesian quadrature back-propagation (via uniform sampling). ;
# I imagine something like this code starting from line 73 of test_likelihood.py. ;

# After loading the 1uao.pdb atomic positions, let's generate a collection of ctf-functions ;

n_images = tp.n_templates
defocus = np.random.uniform(300.0, 900.0, n_images) #<-- I want to make sure that there are multiple ctfs included. ;
#defocus = np.ones(n_images, dtype=np.float64) * 300.0
Astigmatism = np.random.uniform(0, 20, n_images) #<-- this isn't too important, as we can safely ignore astigmatism at low-resolution. ;
defocusU = defocus + Astigmatism / 2
defocusV = defocus - Astigmatism / 2
defocusAng = np.zeros(n_images, dtype=np.float64)#np.random.uniform(-90, 90, n_images)
phaseShift = np.zeros(n_images, dtype=np.float64)
ctf_desc = LensDescriptor(
    defocusU = defocusU, # in Angstrom
    defocusV = defocusV,  # in Angstrom
    defocusAng = defocusAng, # in degrees, defocus angle
    sphericalAberration = 2.7,  # in mm, spherical aberration
    voltage = 300,  # in kV
    amplitudeContrast = 0.1,    # amplitude contrast
    phaseShift = phaseShift,   # in degrees
)
ctf = CTF(
    ctf_descriptor=ctf_desc,
    polar_grid = polar_grid,
    box_size = box_size, # in Angstrom
    anisotropy = True#True#
)

# Now let's generate the images from the templates. ;
# We can set the noise to 0 if you like (e.g., snr = None or snr = float('inf')). ;
# Note that we need to update add_noise_phys to allow for zero noise (or infinite snr). ; 

image = Images.from_templates(templates = tp)
image.apply_ctf(ctf)
image.displace_images_fourier(
    x_displacements = true_displacement_x,
    y_displacements = true_displacement_y,
    precision = precision,
)
image.rotate_images_fourier(true_rotation)
image.transform_to_spatial(grid=(n_pixels, pixel_size), use_cuda = use_cuda)

assert image.images_phys is not None
image_true = deepcopy(image.images_phys.real)
_, sigma_noise = image.add_noise_phys(snr = snr) #<-- remember to update add_noise_phys to allow for zero noise. ;
image.transform_to_fourier(polar_grid, use_cuda = use_cuda)
mean, std = image.center_physical_image_signal()
image.transform_to_fourier(polar_grid, use_cuda = use_cuda)
image.transform_to_spatial(grid=(n_pixels, pixel_size), use_cuda = use_cuda)
image_true -= mean
image_true /= std
assert image.images_phys is not None
sigma_noise = torch.sqrt(torch.mean((image.images_phys - image_true) ** 2, dim = (1, 2)))

# Now let's calculate the cross-correlations. ;
# For now we'll only look at the optimal pose, ;
# as this is required for the 'maximum-likelihood' reconstruction. ;
# However, later on we'll also need to perform a 'maximum-entropy' reconsctruction. ;
# For the maximum-entropy reconstruction we'll need the cross-correlation array X_SM__ ;
# (i.e., X_dwSM____ maximized over d and w). ;

cclik = CrossCorrelationLikelihood(
    templates = tp,
    max_displacement = max_displacement,
    n_displacements_x = n_displacements_x,
    n_displacements_y = n_displacements_y,
    precision = precision,
    device = device,
    verbose = True
)

_imgs = image.images_fourier
_ctf = to_torch(ctf.ctf, precision, "cpu")
assert _imgs is not None

res = cclik._compute_cross_correlation_likelihood(
    device=torch.device("cuda"),
    images_fourier = _imgs,
    ctf = _ctf,
    n_pixels_phys = image.phys_grid.n_pixels[0] * image.phys_grid.n_pixels[1],
    n_images_per_batch=64,#128,#256,
    n_templates_per_batch=16,
    return_type=CrossCorrelationReturnType.OPTIMAL_POSE, #<-- remember to test a more general version later. ;
    return_integrated_likelihood=True
)
assert len(res) == 2
(optimal_pose, log_likelihood_fourier_integrated) = res
assert isinstance(optimal_pose, OptimalPoseReturn)

# Now let's set the 'estimated' pose parameters to be the optimal ones. ;
n_M = n_images
image_delta_x_est_M_ = optimal_pose.optimal_displacement_x_S
image_delta_y_est_M_ = optimal_pose.optimal_displacement_y_S
image_gamma_z_est_M_ = optimal_pose.optimal_inplane_rotation_S

# And now let's transform the image-stack so that we 'undo' the effect of the optimal displacement and rotation. ;
# This is very similar to the manipulations within: ;
# calc_distance_optimal_templates_vs_physical_image ;
# within likelihood.py. ;
# Importantly, the function calc_distance_optimal_templates_vs_physical_image constructs the templates: ;
# R(-gamma_z) T(+delta_x,+delta_y) S, ;
# which are then compared against the images. ;
# Thus, for this function, we will calculate: ;
# T(-delta_x,-delta_y) R(+gamma_z) M. ;
# Note that this particular implementation enforces an 'on-grid' (i.e., discretized) in-plane rotation gamma_z. ;
# This should be fixed later (although it is probably fine for low-resolution reconstruction). ;
# Note also that, if we have an anisotropic CTF, we'll need to rotate each CTF-function as well.

def M_trn_k_p_wkM___from_M_ori_k_p_wkM___(
    M_ori_k_p_ : Images, #<-- these is the original image-stack. ;
    image_delta_x_est_M_ : torch.Tensor | None = None, #<-- estimated delta_0. ;
    image_delta_y_est_M_ : torch.Tensor | None = None, #<-- estimated delta_1. ;
    image_gamma_z_est_M_ : torch.Tensor | None = None, #<-- estimated gamma_z. ;
    CTF_ori_k_p_ : CTF | None = None,
    precision : Precision = Precision.SINGLE,
    use_cuda : bool = True
):
    device = 'cuda' if use_cuda and torch.cuda.is_available() else 'cpu'
    assert M_ori_k_p_.images_fourier is not None
    n_M = M_ori_k_p_.n_images
    image_gamma_z_dsc_M_ = None
    if image_gamma_z_est_M_ is not None:
        dgamma_z = 2 * np.pi / M_ori_k_p_.polar_grid.n_inplanes
        image_gamma_z_dsc_M_ = - torch.round(image_gamma_z_est_M_ / dgamma_z).to(torch.int64)
    M_trn_k_p_wkM___ = to_torch(M_ori_k_p_.images_fourier, precision, device) #<-- This is supposed to extract the images_fourier array and store it as M_trn_k_p_wkM___. ;
    if image_gamma_z_est_M_ is not None:
        assert image_gamma_z_dsc_M_ is not None
        for i in range(M_trn_k_p_wkM___.shape[0]):
            inplane_rotation_discrete = int(image_gamma_z_dsc_M_[i].item())
            # Note that we reverse the direction of the rotation, as these are being applied to the images and not the templates. ;
            M_trn_k_p_wkM___[i] = torch.roll(M_trn_k_p_wkM___[i], shifts = -inplane_rotation_discrete, dims = 1)
    translation_kernel___ = None
    if image_delta_x_est_M_ is not None and image_delta_y_est_M_ is not None:
        image_delta_x_est_M_ *= 2.0 / M_ori_k_p_.box_size[0]
        image_delta_y_est_M_ *= 2.0 / M_ori_k_p_.box_size[1]
        assert image_delta_x_est_M_ is not None
        assert image_delta_y_est_M_ is not None
        # Note that we reverse the direction of the translation, as these are being applied to the images and not the templates. ;
        translation_kernel___ = translation_kernel_fourier(M_ori_k_p_.polar_grid,-image_delta_x_est_M_,-image_delta_y_est_M_, precision, device) 
    if translation_kernel___ is not None:
        M_trn_k_p_wkM___ = M_trn_k_p_wkM___ * translation_kernel___
    if CTF_ori_k_p_ is not None:
        if type(M_ori_k_p_.images_fourier) is np.ndarray:
            CTF_ori_k_p_wkM___ = np.array(CTF_ori_k_p_.ctf, dtype = M_ori_k_p_.images_fourier.dtype)
        if type(M_ori_k_p_.images_fourier) is torch.Tensor:
            CTF_ori_k_p_wkM___ = torch.tensor(CTF_ori_k_p_.ctf, dtype = M_ori_k_p_.images_fourier.dtype, device = M_ori_k_p_.images_fourier.device)
        CTF_trn_k_p_wkM___ = to_torch(CTF_ori_k_p_wkM___, precision, device) #<-- this is supposed to copy the CTF_ori_k_p_wkM___ array. ;
        if image_gamma_z_est_M_ is not None:
            assert image_gamma_z_dsc_M_ is not None
            for i in range(CTF_trn_k_p_wkM___.shape[0]):
                inplane_rotation_discrete = int(image_gamma_z_dsc_M_[i].item())
                # Note that we rotate the CTF functions the same way that we rotated the images. ;
                CTF_trn_k_p_wkM___[i] = torch.roll(CTF_trn_k_p_wkM___[i], shifts = -inplane_rotation_discrete, dims = 1)
    return M_trn_k_p_wkM___ , CTF_trn_k_p_wkM___

M_trn_k_p_wkM___ , CTF_trn_k_p_wkM___ = M_trn_k_p_wkM___from_M_ori_k_p_wkM___(
    M_ori_k_p_ = Images,
    image_delta_x_est_M_= image_delta_x_est_M_,
    image_delta_y_est_M_= image_delta_y_est_M_,
    image_gamma_z_est_M_= image_gamma_z_est_M_,
    CTF_ori_k_p_ = ctf,
    precision = precision,
    use_cuda = True
)

# Now we'll perform the volumetric reconstruction (i.e., cartesian quadrature back-propagation onto a uniform grid). ;

k_p_r_max = radius_max #<-- I'm not sure exactly where to draw this from, but it should depend on the Images object. ;
n_x_u = 64 #<-- we can pick this to be any resolution. ;
n_k_u = n_x_u #<-- can be smaller if there is not much data. ;
half_diameter_k_c = k_p_r_max
diameter_k_c = 2.0 * half_diameter_k_c
# Here we recreate the matlab ndgrid for k-space. ;
k_u_0_ = np.linspace(-k_p_r_max, +k_p_r_max, n_k_u).reshape(-1, 1)
d_k_0 = np.mean(np.diff(k_u_0_))
k_u_1_ = np.linspace(-k_p_r_max, +k_p_r_max, n_k_u).reshape(-1, 1)
d_k_1 = np.mean(np.diff(k_u_1_))
k_u_2_ = np.linspace(-k_p_r_max, +k_p_r_max, n_k_u).reshape(-1, 1)
d_k_2 = np.mean(np.diff(k_u_2_))
k_u_0___, k_u_1___, k_u_2___ = np.meshgrid(k_u_0_, k_u_1_, k_u_2_, indexing='ij')
n_kkk_u = n_k_u ** 3
weight_kkk_u_ = d_k_0 * d_k_1 * d_k_2
# Here we need to extract the cartesian coordinates for the image-data in k-space (at essentially arbitrary points). ;
index_nS_from_nM_ = optimal_pose.optimal_template_S #<-- I want to extract the template-index associated with each image-index. ;
k_c_0_wkS___ , k_c_1_wkS___ , k_c_2_wkS___ = do_something_like(rearrange(tp.templates_fourier.fourier_slices_)) #<-- I'm not sure how to access the actual cartesian-positions 
# For memory reasons we can think about delaying the index access until later. ;
k_c_0_wkM___ = k_c_0_wkS___[:, :, index_nS_from_nM_] #<-- this should be the k_c_0 coordinate for the various image-data points. ;
k_c_1_wkM___ = k_c_1_wkS___[:, :, index_nS_from_nM_] #<-- this should be the k_c_1 coordinate for the various image-data points. ;
k_c_2_wkM___ = k_c_2_wkS___[:, :, index_nS_from_nM_] #<-- this should be the k_c_2 coordinate for the various image-data points. ;
# And now I want to define the (uniform-grid) indices associated with each cartesian point for the image-data. ;
# These indices will link the original data (at essentially arbitrary points) to the uniform-grid defined by k_u_0_, k_u_1_ and k_u_2_. ;
index_k_u_0_wkM___ = np.clip(np.round((k_c_0_wkM___ + half_diameter_k_c) / max(1e-12, d_k_0)), 0, n_k_u - 1).astype(int)
index_k_u_1_wkM___ = np.clip(np.round((k_c_1_wkM___ + half_diameter_k_c) / max(1e-12, d_k_1)), 0, n_k_u - 1).astype(int)
index_k_u_2_wkM___ = np.clip(np.round((k_c_2_wkM___ + half_diameter_k_c) / max(1e-12, d_k_2)), 0, n_k_u - 1).astype(int)
index_k_u_wkM___ = index_k_u_0_wkM___ + (index_k_u_1_wkM___ + index_k_u_2_wkM___ * n_k_u) * n_k_u
# And now we use the sparse matrix structure as a way to sum up the moments. ;
M0C2_k_u_ = csr_matrix((np.abs(CTF_trn_k_p_wkM___) ** 2, (index_k_u_wkM___.ravel(), np.zeros(n_kkk_u))), shape=(n_kkk_u, 1))
# Note that, in the following code, we divide and then multiply by M_upb. ;
# This is in an attempt to ensure that very large values of M do not overwhelm the thresholding of the CTF in the denominator. ;
# There is probably a better way to do this. ;
# My current thinking is as follows: ;
# We want the resulting a_k_u_qpbu_ to be accurate to, say 1e-6. ;
# Thus, for any frequency-voxel associated with only a single small CTF-value (say, less than 1e-6), ;
# we want to prevent the division M1C1/M0C2 from being greater than M_upb*1e6. ;
# To do this we can simply threshold the denominator at 1e-6. ;
# Now any CTF values of 1e-6 and smaller will produce M1C1/M0C2 around M_upb, which should be acceptable. ;
M_upb = np.max(np.abs(M_trn_k_p_wkM__)) #<-- upper bound for abs(M). ;
M1C1_k_u_ = csr_matrix(((CTF_trn_k_p_wkM___ * (M_trn_k_p_wkM___ / max(1e-12, M_upb)), (index_k_u_wkM___.ravel(), np.zeros(n_kkk_u))), shape=(n_kkk_u, 1))
# Now we reconstruct the volume a_k_u_qbpu_ in frequency-space. ;
# Note the thresholding of the denominator. ;
a_k_u_qbpu_ = (M1C1_k_u_ / max(1e-6, M0C2_k_u_)
# Now we define a uniform-spatial-grid. ;
diameter_x_c = 2.0
half_diameter_x_c = 0.5*diameter_x_c
x_p_r_max = half_diameter_x_c
x_u_0_ = np.linspace(-x_p_r_max, +x_p_r_max, n_x_u).reshape(-1, 1)
d_x_0 = np.mean(np.diff(x_u_0_))
x_u_1_ = np.linspace(-x_p_r_max, +x_p_r_max, n_x_u).reshape(-1, 1)
d_x_1 = np.mean(np.diff(x_u_1_))
x_u_2_ = np.linspace(-x_p_r_max, +x_p_r_max, n_x_u).reshape(-1, 1)
d_x_2 = np.mean(np.diff(x_u_2_))
x_u_0___, x_u_1___, x_u_2___ = np.meshgrid(x_u_0_, x_u_1_, x_u_2_, indexing='ij')
# Now we use the fft to reconstruct the volume a_x_u_qbpu_ in physical space from the frequency-space a_k_u_qbpu_. ;
# This particular implementation uses the slower (but more intuitive) nufft3d3. ;
# This can of course be replaced by an appropriately padded fft. ;
eta = np.pi / k_p_r_max
tolerance_nufft = 1e-6
a_x_u_qbpu_ = finufft.nufft3d3(
    2.0*np.pi*k_c_0_wkM___*eta,
    2.0*np.pi*k_c_1_wkM___*eta,
    2.0*np.pi*k_c_2_wkM___*eta,
    a_k_u_qbpu_.*weight_kkk_u_,
    +1,
    tolerance_nufft,
    x_u_0___(:)/eta,
    x_u_1___(:)/eta,
    x_u_2___(:)/eta
)
# Now we finally rescale by M_upb. ;
a_x_u_qbpu_ = M_upb * a_x_u_qbpu_

# Update 20241203: Chatting with Wai-Shing: ;
# We discussed a regime where the number of images n_M is much larger than the number of templates. ;
# In this regime memory is a consideration. ;
# For this large-image limit we can consider the following version of the quadrature back-propagation. ;
#         
# Step 0: The first thing to consider is that images_fourier (i.e., M_ori_k_p_wkM___) already takes up a lot of space. ;
#         More specifically, M_ori_k_p_wkM___ takes up n_w_max*n_k_p_r*n_M*8 bytes. ;
#         If n_w_max=128, n_k_p_r=64 and n_M = 1024^2 (for a million images at low resolution), ;
#         then M_ori_k_p_wkM___ already takes up 8*8GB=64GB of space. ;
#         By contrast, n_S might only be 1024 or so, ;
#         implying that S_ori_k_p_wkS___ (i.e., templates_fourier) might only take up 64MB of space. ;
#         As a result, we don't actually want to create and store M_trn_k_p_wkM___ and CTF_trn_k_p_wkM___. ;
#         Instead, we want to use M_ori_k_p_wkM___ and CTF_ori_k_p_wkM___ (i.e., whatever is already stored in images_fourier and the ctf). ;
#         
# Step 1: Now let's note that rotations and translations commute as follows: ;
#         R(gamma_ori)*T(delta_ori_) = T(delta_rot_)*R(gamma_ori), ;
#         where delta_rot_ = R(+gamma_ori)*delta_ori_ corresponds to the rotated delta-vector. ;
#         Equivalently, we can write: ;
#         R(gamma_ori)*T(delta_inv_) = T(delta_ori_)*R(gamma_ori), ;
#         where delta_inv_ = R(-gamma_ori)*delta_ori_ corresponds to the inversely rotated delta-vector. ;
#         With this in mind, we can revisit the original definition of M_trn_k_p_wkM___. ;
#         For the moment, let's denote the nM-image via M_ori_k_p_wk__ := M_ori_k_p_wkM___(:,:,nM). ;
#         Now M_trn_k_p_wk__ is defined as follows: ;
#         M_trn_k_p_wk__ = T(-delta_est_) R(+gamma_est) M_ori_k_p_wk__ ;
#         and we can rewrite it as: ;
#         M_trn_k_p_wk__ = R(+gamma_est) T(-delta_inv_) M_ori_k_p_wk__, ;
#         where -delta_inv_ = R(-gamma_est)*-delta_est_. ;
#         This observation suggests the following strategy: ;
#         Each time we determine the estimated translations and rotations for each image, ;
#         we can replace (i.e., overwrite) each M_ori_k_p_wk__ with the appropriately centered version of that image: ;
#         M_cen_k_p_wk__ := T(-delta_inv_)*M_ori_k_p_wk__. ;
#         In this way we only have to track the estimated rotations gamma_est for each image (denoted by image_gamma_z_est_M_ above). ;
#         Note also that, when centering (i.e., translating) the images in this way, ;
#         the associated CTF-functions do not change (even if they happen to be anisotropic). ;
#         
# Step 2: Okay, so now let's assume that we've maintained a centered stack of images M_cen_k_p_wkM__ ;
#         (overwriting images_fourier each time we obtain an estimate of image_delta_x_est_M_, image_delta_y_est_M_ and image_gamma_z_est_M_). ;
#         Now, for the volumetric-reconstruction, we previously calculated each transformed-image: ;
#         M_trn_k_p_wk__ = R(+gamma_est) M_cen_k_p_wk__ , ;
#         CTF_trn_k_p_wk__ = R(+gamma_est) CTF_ori_k_p_wk__ , ;
#         where, for each image-index nM, gamma_est := index_gamma_z_est_M_(nM). ;
#         We then associated that transformed image- and ctf-data with a particular set of locations k_c_?_wk__ in frequency-space. ;
#         Thus, for example, if the image M_trn_k_p_wk__ is associated with the template S_ori_k_p_wk__ with viewing-angles (polar_a,azimu_b), ;
#         then we would associate M_trn_k_p_wk__ and CTF_trn_k_p_wk__ with the 'fourier-circles' generated with viewing-angles (polar_a,azimu_b), ;
#         which are referred to as k_c_?_wkS___(:,:,index_nS_from_nM_(nM)) above. ;
#         A straightforward observation at this point is that, ;
#         instead of associating k_c_?_wkS___(:,:,index_nS_from_nM_(nM)) with the transformed-image nM, ;
#         we can instead associate the transformed-version of the fourier-data with the original-image. ;
#         That is to say, we can associate the original-image M_cen_k_p_wk__ with the 'fourier-circle' generated using a non-zero in-plane-rotation. ;
#         You already have the code to do this (see the '_fourier_circles' function in 'template.py'). ;
#         However, when you construct the 'ViewingAngles' object (see 'viewing_angles.py') you typically set 'gammas=None'. ;
#         Instead, you could set 'gammas= pick_the_correct_sign_to_use * image_gamma_z_est_M_', ;
#         and generate the appropriately arranged fourier-data directly. ;
#         
# Step 3: Now, as suggested above, the generation of all this in-plane-rotated fourier-data will still take up quite a bit of space. ;
#         For example, the 'xyz_template_points' array will typically contain n_w_max*n_k_p_r*n_M*4*3 bytes, which is 96GB of space. ;
#         This is also too expensive when dealing with many images. ;
#         However, there is a straightforward work-around, which is never to actually store the xyz_template_points, ;
#         but instead to only store their 'bins' (i.e., corresponding indices within the uniform fourier-grid denoted by k_u_?___ above). ;
#         These bins are referred to above as index_k_u_?_wkM___. ;
#         To perform this efficiently, let's imagine that we fix (i.e., hardcode) the resolution n_k_u of the uniform fourier-grid to be <=256. ;
#         For presentation, let's just assume that n_k_u is 256. ;
#         This means that, for each particular value of k_c_0_wkS___ (i.e., denoted by x_template_points) ;
#         there is an associated unsigned 8-bit integer: ;
#         uint8_k_u_0_wkS___ = np.clip(np.round((k_c_0_wkS___ + half_diameter_k_c) / max(1e-12, d_k_0)), 0, n_k_u - 1).astype(uint8)
#         with similar definitions for uint8_k_u_1_wkS___ and uint8_k_u_2_wkS___. ;
#         Now, as long as gamma_est is an 'on-grid' rotation (i.e., a multiple of (2*pi)/number_of_thetas) then an in-plane rotation of the ;
#         fourier-data corresponds to a circular-shift by some number of steps (in either the positive- or negative-direction). ;
#         Thus, we can generate the bin-indices for the appropriately-oriented fourier-data for the nM-image just by: ;
#         3a. determining the template-index nS corresponding to the nM-image. ;
#             uint8_ori_k_u_0_wk__ = uint8_k_u_0_wkS___(:,:,nS) ;
#             uint8_ori_k_u_1_wk__ = uint8_k_u_1_wkS___(:,:,nS) ;
#             uint8_ori_k_u_2_wk__ = uint8_k_u_2_wkS___(:,:,nS) ;
#         3b. determining the circular-shift corresponding to gamma_est for the nM-image, ;
#             ngamma_z_est_M_ = torch.round(image_gamma_z_est_M_ / dgamma_z).to(torch.int64)
#             ngamma_z_est = int(ngamma_z_est_M_[i].item())
#         3c. determining the bin-indices for the nM-image, ;
#             uint8_trn_k_u_0_wk__ = torch.roll(uint8_ori_k_u_0_wk__k_p_wk__, shifts = pick_the_right_sign * ngamma_z_est, dims = 1)
#             uint8_trn_k_u_1_wk__ = torch.roll(uint8_ori_k_u_1_wk__k_p_wk__, shifts = pick_the_right_sign * ngamma_z_est, dims = 1)
#             uint8_trn_k_u_2_wk__ = torch.roll(uint8_ori_k_u_2_wk__k_p_wk__, shifts = pick_the_right_sign * ngamma_z_est, dims = 1)
#         For example, the 'xyz_template_points' array will typically contain n_w_max*n_k_p_r*n_M*4*3 bytes, which is 96GB of space. ;
#         And now the collection of uint8_trn_k_u_0_wkM___, uint8_trn_k_u_1_wkM___ and uint8_trn_k_u_2_wkM___ ;
#         should take only n_w_max*n_k_p_r*n_M*1*3 = 24GB of space, which is roughly 3/8 of the storage for M_cen_k_p_wkM___. ;
#         Moreover, we don't actually need to store all of the uint8_trn_k_u_?_wkM___ simultaneously (as these can be generated per image-batch). ;
#         
# Step 4: Given the bin-indices uint8_trn_k_u_?_wkM___ above, we could once again combine them via: ;
#         uint32_trn_k_u_wkM___ = uint8_trn_k_u_0_wkM___ + (uint8_trn_k_u_1_wkM___ + uint8_trn_k_u_2_wkM___ * n_k_u) * n_k_u, ;
#         but this now produces a new array of size n_w_max*n_k_p_r*n_M*4 = 32GB of space. ;
#         Instead, there is probably a better way to directly bucket/bin the data using the uint8_trn_k_u_?_wkM___ indices. ;
#         If that's too hard, remember that ultimately all we need are the two binned-arrays on the uniform-grid: M0C2_k_u_ and M1C1_k_u_. ;
#         Thus, sub-arrays of the form M0C2_sub_k_u_ and M1C1_sub_k_u_ can be created for each image-batch, ;
#         and then summed up (across batches) to form the final M0C2_k_u_ and M1C1_k_u_. ;
#         Such a batching process can avoid creating uint32_trn_k_u_wkM___ across the entire image-pool. ;




