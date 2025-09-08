from numpy import pi
import torch
from torch.testing import assert_close
import pytest

from EMPM.likelihoods.cross_correlation_likelihood import CrossCorrelationLikelihood, conform_ctf
from EMPM.util import (
    CrossCorrelationReturnType,
    Precision,
    to_torch,
)

from cross_correlation_fixtures import (
    make_polar_grid,
    make_viewing_angles,
    make_planewave_templates,
    viewing_angles_to_cartesian_displacements,
    get_planar_ctf,
    # p_xx_p_multiwave_iterative,
    p_xx_p_multiwave_vectorized,
)


class parameters():
    device: str
    n_pixels: int
    precision: Precision
    wavevector: torch.Tensor
    max_displacement: float  # this is per axis; real max will be this * root-2
    template_ctf_angle: float
    image_ctf_angle: float

    def __init__(
        self,
        device: str,
        n_pixels: int,
        precision: Precision,
        wavevector: torch.Tensor,
        template_ctf_angle: float,
        image_ctf_angle: float,
        max_displacement: float,
    ):
        self.device = device
        self.n_pixels = n_pixels
        self.precision = precision
        self.wavevector = wavevector
        self.template_ctf_angle = template_ctf_angle
        self.image_ctf_angle = image_ctf_angle
        self.max_displacement = max_displacement


    def duplicate(self, *,
        device: str | None = None,
        n_pixels: int | None = None,
        precision: Precision | None = None,
        wavevector: torch.Tensor | None = None,
        template_ctf_angle: float | None = None,
        image_ctf_angle: float | None = None,
        max_displacement: float | None = None,
    ):
        (float_type, _, _) = self.precision.get_dtypes(default=Precision.DOUBLE)
        if precision is not None:
            (float_type, _, _) = precision.get_dtypes(default=Precision.DOUBLE)
        wv = wavevector if wavevector is not None else self.wavevector

        return parameters(
            device=device if device is not None else self.device,
            n_pixels=n_pixels if n_pixels is not None else self.n_pixels,
            precision=precision if precision is not None else self.precision,
            wavevector=wv.to(float_type),
            template_ctf_angle=template_ctf_angle if template_ctf_angle is not None else self.template_ctf_angle,
            image_ctf_angle=image_ctf_angle if image_ctf_angle is not None else self.image_ctf_angle,
            max_displacement=max_displacement if max_displacement is not None else self.max_displacement,
        )


    @staticmethod
    def default():
        (default_float, _, _) = Precision.DOUBLE.get_dtypes(default=Precision.DOUBLE)
        return parameters(
            device='cpu',
            n_pixels=128,
            precision=Precision.DOUBLE,
            wavevector=torch.tensor([[-0.17, -0.03, 0.07]], dtype=default_float),
            template_ctf_angle= -pi / 5.,
            image_ctf_angle=pi / 3.,
            max_displacement=0.08,
        )


# TODO: Parameterize based on: max_displacement,
    # this basically simulates the effect of a displaced atom
    # So maybe also try -.8 to +.8
    # try 3 to 5 possible values in that range? Further away = higher values = lower accuracy

def make_cases() -> list[parameters]:
    cases = [parameters.default()]
    with_single_precision = [x.duplicate(precision=Precision.SINGLE) for x in cases]
    cases.extend(with_single_precision)
    with_cuda = [x.duplicate(device='cuda') for x in cases]
    cases.extend(with_cuda)
    low_pixel = [x.duplicate(n_pixels=64) for x in cases]
    high_pixel = [x.duplicate(n_pixels=256) for x in cases]
    cases.extend(low_pixel)
    cases.extend(high_pixel)
    multiwave_sources = [torch.tensor([[-0.17, -0.03, 0.07], [0.00, 0.00, 0.00]]),
                         torch.tensor([[-0.15, -0.07, 0.07], [0.12, 0.08, 0.02]]),
                         torch.tensor([[0.1, 0.03, -0.02], [0.12, 0.03, -0.02], [-0.02, 0.14, -0.08]])]
    multiwave = [x.duplicate(wavevector=y) for y in multiwave_sources for x in cases]
    cases.extend(multiwave)
    return cases


param_matrix = make_cases()
@pytest.mark.parametrize("params", param_matrix)
def test_cross_correlation_PxxP_from_a_k_p(params: parameters):
    if (params.device == 'cuda' and not torch.cuda.is_available()):
        pytest.skip("Test cannot run because CUDA is not present.")

    (torch_float_type, _, _) = params.precision.get_dtypes(default=Precision.DOUBLE)
    with torch.device(params.device):
        torch.set_default_dtype(torch_float_type)
        ## More stuff to probably parametrize over?
        box_size = 2.0
        n_displacements_x = 3
        n_displacements_y = 5
        wavevector_planewave = params.wavevector
        template_ctf_angle = params.template_ctf_angle  # phi_S
        image_ctf_angle = params.image_ctf_angle        # phi_M
        displacement_planewave_image = torch.tensor([0.2, -0.18])
        ####

        polar_grid = make_polar_grid(params.n_pixels)
        viewing_angles = make_viewing_angles(params.device, torch_float_type)

        ### Set up numeric (EMPM-native) objects
        # CTFs are defined by the angles (and numerically also depend on the grid)
        planar_ctf_template = get_planar_ctf(polar_grid, template_ctf_angle, box_size, params.precision, params.device)
        planar_ctf_image = get_planar_ctf(polar_grid, image_ctf_angle, box_size, params.precision, params.device)
        templates = make_planewave_templates(wavevector_planewave, viewing_angles, polar_grid, params.precision)
        ## Images are a copy of the templates, with an applied displacement vector
        images = templates.to_images()
        images.displace_images_fourier(
            x_displacements = displacement_planewave_image[0].item(),
            y_displacements = displacement_planewave_image[1].item()
        )
        images.apply_ctf(planar_ctf_image)

        cc = CrossCorrelationLikelihood(
            templates = templates,
            max_displacement = params.max_displacement,
            n_displacements_x = n_displacements_x,
            n_displacements_y = n_displacements_y,
            precision = params.precision,
            device = params.device,
            verbose = False
        )
        ctf_tensor = conform_ctf(to_torch(planar_ctf_template.ctf, params.precision, params.device), planar_ctf_template.anisotropy)
        res = cc._compute_cross_correlation_likelihood(
            device=torch.device(params.device),
            images_fourier = images.images_fourier,
            ctf=ctf_tensor,
            n_pixels_phys = params.n_pixels*params.n_pixels,
            n_templates_per_batch=viewing_angles.n_angles,
            n_images_per_batch=viewing_angles.n_angles,
            return_type=CrossCorrelationReturnType.FULL_TENSOR,
            return_integrated_likelihood=False,
        )

        # Analytic objects
        wavevector_templates = viewing_angles_to_cartesian_displacements(viewing_angles, wavevector_planewave).to(params.device)
        wavevector_images = wavevector_templates.clone() - displacement_planewave_image.to(wavevector_templates.device)
        analytic = p_xx_p_multiwave_vectorized(
            wavevector_templates,
            wavevector_images,
            template_ctf_angle,
            image_ctf_angle,
            cc._gamma.to(params.device) * -1.0,
            torch.stack([cc.x_displacements_expt_scale, cc.y_displacements_expt_scale]).T.to(params.device),
            polar_grid.radius_max
        )

        assert_close(res.cross_correlation_SMdw, analytic.cpu())
        torch.set_default_dtype(torch.float32)

