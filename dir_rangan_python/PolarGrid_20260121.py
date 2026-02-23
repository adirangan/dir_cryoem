from dir_matlab_macros import * ;
from torch import Tensor
# from typing import Self  # note this requires python >= 3.11;
# for python 3.10 use
from typing_extensions import Self

class PolarGrid():
    """Class implementing 2D polar-coordinate grid.

    Attributes:
        is_uniform (bool): Flag indicating whether the grid is uniform
        n_k_p_r (int): Number of frequency bands / radii
        k_p_r_ (Tensor): Frequency (k) value for each band. 1d real tensor, expected
            to have n_k_p_r elements.
        k_p_r_max (int): Highest frequency (k) value on the polar grid. Should be
            equal to the maximum value of k_p_r_.
        template_k_eq_d (float): Equatorial distance between frequencies for a
            uniform grid [??]
        n_w_ (Tensor): Number of angles per shell (will be the same value for
            every r in the case of a uniform grid). 1d real tensor of length
            equal to n_k_p_r.
        weight_2d_k_p_r_ (Tensor): Quadrature-weight for each frequency band
            on a 2d polar grid. 1D real tensor of length equal to the number
            of frequency bands, n_k_p_r. Sums to pi * k_p_r_max^2.
        weight_2d_k_p_wk_ (Tensor): Quadrature-weight for each quadrature
            point on the grid. Unrolled 1D real tensor of total length equal to
            the number of quadrature points, which is n_w_[0] * n_k_p_r for
            uniform grids. For uniform grids, the nth element of this tensor
            represents the weight at the point whose radial dimension is
            floor(n / n_w_[0]) and whose rotational dimension is n % n_w_[0].
            Sums to k_p_r_max^2 / 4pi.
        k_p_r_wk_ (Tensor): Radial dimension of each quadrature point on the
            grid, unrolled to 1D. 1D real tensor of total length equal
            to the number of quadrature points (n_w_[0] * n_k_p_r for uniform
            grids). In the case of uniform grids, the first n_w_[0] elements
            will equal k_p_r_[0], the second n_w_[0] elements will equal
            k_p_r_[1], etc.
        k_p_w_wk_ (Tensor): Angular dimension of each quadrature point on the
            grid, unrolled to 1D. 1D real tensor of total length equal to the
            number of quadrature points (n_w_[0] * n_k_p_r for uniform grids).
            In the case of uniform grids, subsequent elements of this tensor
            will increase by steps of size template_k_eq_d until they reach
            2pi, at which point they'll start over from 0; there should be
            n_k_p_r such repetitions.
        k_c_0_wk_ (Tensor): Cartesian x-coordinate of each quarature point on
            the grid, unrolled to 1D. 1D real tensor of total length equal to
            the number of quadrature points (as above). For an element n
            of this tensor, k_c_0_wk_[n] = k_p_r_wk_[n] * cos(k_p_w_wk_[n]).
        k_c_1_wk_ (Tensor): Cartesian y-coordinate of each quarature point on
            the grid, unrolled to 1D. 1D real tensor of total length equal to
            the number of quadrature points (as above). For an element n
            of this tensor, k_c_0_wk_[n] = k_p_r_wk_[n] * sin(k_p_w_wk_[n]).
    """

    is_uniform: bool
    n_k_p_r: int
    k_p_r_: Tensor
    k_p_r_max: int
    template_k_eq_d: float
    n_w_: Tensor
    weight_2d_k_p_r_: Tensor
    weight_2d_k_p_wk_: Tensor
    k_p_r_wk_: Tensor
    k_p_w_wk_: Tensor
    k_c_0_wk_: Tensor
    k_c_1_wk_: Tensor


    def __init__(self,
        is_uniform: bool,
        n_k_p_r: int,
        k_p_r_: Tensor,
        k_p_r_max: int,
        template_k_eq_d: float,
        n_w_: Tensor,
        weight_2d_k_p_r_: Tensor,
        weight_2d_k_p_wk_: Tensor,
        k_p_r_wk_: Tensor,
        k_p_w_wk_: Tensor,
        k_c_0_wk_: Tensor,
        k_c_1_wk_: Tensor,
    ) -> None:
        self.is_uniform = is_uniform
        self.n_k_p_r = n_k_p_r
        self.k_p_r_ = k_p_r_
        self.k_p_r_max = k_p_r_max
        self.template_k_eq_d = template_k_eq_d
        self.n_w_ = n_w_
        self.weight_2d_k_p_r_ = weight_2d_k_p_r_
        self.weight_2d_k_p_wk_ = weight_2d_k_p_wk_
        self.k_p_r_wk_ = k_p_r_wk_
        self.k_p_w_wk_ = k_p_w_wk_
        self.k_c_0_wk_ = k_c_0_wk_
        self.k_c_1_wk_ = k_c_1_wk_


    @classmethod
    def make_uniform_grid(
        cls,
        n_k_p_r: int = 0,
        k_p_r_: Tensor = torch.zeros(0),
        k_p_r_max: int = 0,
        template_k_eq_d: float = -1.,
        n_w_0in_: Tensor = torch.zeros(0),
        weight_3d_k_p_r_: Tensor = torch.zeros(0)
    ) -> Self:
        # NOTE there are shortcuts for this sort of long list, there's a better way to do this
        # But I don't want to introduce any more complication than I need to

        (n_w_,
        weight_2d_k_p_r_,
        weight_2d_k_p_wk_,
        k_p_r_wk_,
        k_p_w_wk_,
        k_c_0_wk_,
        k_c_1_wk_,
        ) = get_weight_2d_2(
            flag_verbose = 0,
            n_k_p_r = n_k_p_r,
            k_p_r_ = k_p_r_,
            k_p_r_max = k_p_r_max,
            template_k_eq_d = template_k_eq_d,
            n_w_0in_ = n_w_0in_,
            weight_3d_k_p_r_ = weight_3d_k_p_r_,
        )

        return cls(
            is_uniform = True,
            n_k_p_r = n_k_p_r,
            k_p_r_ = k_p_r_,
            k_p_r_max = k_p_r_max,
            template_k_eq_d = template_k_eq_d,
            n_w_ = n_w_,
            weight_2d_k_p_r_ = weight_2d_k_p_r_,
            weight_2d_k_p_wk_ = weight_2d_k_p_wk_,
            k_p_r_wk_ = k_p_r_wk_,
            k_p_w_wk_ = k_p_w_wk_,
            k_c_0_wk_ = k_c_0_wk_,
            k_c_1_wk_ = k_c_1_wk_,
        )


## NOTE: I am not super fussy about whether this code lives here or
# elsewhere, but it should be the way we initialize uniform polar grids
# everywhere. To the extent that it's designed to handle non-uniform
# grids, we should probably break that functionality out into a different
# function--as far as I can tell that's the purpose of allowing the user
# to pass in a custom tensor of angular elements and weights, but it kind
# of just complicates things.
# Alternatively, since setting the sines and cosines and radial and angular
# coordinates has to happen regardless, maybe we should break the two
# branches of "prearranged grid" and "compute grid from consistent equatorial
# distance" out into separate helper functions?

def get_weight_2d_2(
        flag_verbose: int = 0,
        n_k_p_r: int = 0,
        k_p_r_: Tensor = torch.zeros(0),
        k_p_r_max: int = 0,
        template_k_eq_d: float = -1.,
        n_w_0in_: Tensor = torch.zeros(0),
        weight_3d_k_p_r_: Tensor = torch.zeros(0)
):
    if flag_verbose > 0: print(f" %% [entering get_weight_2d_2]")

    if n_k_p_r < 1 or k_p_r_.numel() == 0 or weight_3d_k_p_r_.numel() == 0:
        # This should be an error rather than a warning, because the subsequent
        # code will actually error out if the values aren't set.
        raise ValueError(f" %% Warning, precomputation required")

    n_w_ = torch.zeros(n_k_p_r, dtype=torch.int32)

    if template_k_eq_d <= 0:
        if n_w_0in_.numel() == 0:
            if flag_verbose > 0: print(f" %% calculating n_w_0in_")
            l_max_upb = matlab_scalar_round(2 * pi * k_p_r_max)
            n_w_max = 2 * (l_max_upb + 1)
            n_w_0in_ = n_w_max * torch.ones(n_k_p_r).to(dtype=torch.int32)
        if flag_verbose > 0: print(f" %% template_k_eq_d <= 0")
        n_w_ = n_w_0in_.detach().clone()
        assert numel(n_w_) == n_k_p_r and torch.min(n_w_) > 0
    else:
        if flag_verbose > 0:
            print(f" %% template_k_eq_d > 0")
        for nk_p_r in range(n_k_p_r):
            k_p_r = k_p_r_[nk_p_r].item()
            n_equator = 3 + matlab_scalar_round(2 * pi * k_p_r / template_k_eq_d)
            n_polar_a = 3 + matlab_scalar_round(n_equator / 2)
            n_w_[nk_p_r] = 2 * n_polar_a
        # I *think* the above loop can be written as:
        # equator_factor = 2 * pi / template_k_eq_d
        # n_w_ = 3 + matlab_scalar_round(matlab_scalar_round(k_p_r_ * equator_factor)/2 + 1.5)
        ## TODO: Check the above, & see if possible to simplify the repeated rounding ops

    n_w_sum = int(torch.sum(n_w_).item())
    n_w_csum_ = cumsum_0(n_w_)

    if flag_verbose > 0:
        n_w_max = int(torch.max(n_w_).item())
        print(f" %% n_w_max {n_w_max} n_w_sum {n_w_sum}")

    weight_2d_k_p_r_ = 2 * pi * torch.reshape(weight_3d_k_p_r_, (n_k_p_r,)) \
                              / torch.maximum(torch.tensor([1e-12]),k_p_r_)
    weight_2d_k_p_r_ = weight_2d_k_p_r_.to(dtype=torch.float32)

    weight_2d_k_p_wk_ = torch.zeros(n_w_sum, dtype=torch.float32)
    for nk_p_r in range(n_k_p_r):
        tmp_index_ = int(n_w_csum_[nk_p_r].item()) + torch.arange(int(n_w_[nk_p_r].item()),dtype=torch.int32)
        weight_2d_k_p_wk_[tmp_index_] = weight_2d_k_p_r_[nk_p_r].item() / max(1,int(n_w_[nk_p_r].item())) / (2 * pi) ** 2
    # NOTE: The above loop can also probably be written as a vectorized pytorch operation

    k_c_0_wk_ = torch.zeros(n_w_sum, dtype=torch.float32)
    k_c_1_wk_ = torch.zeros(n_w_sum, dtype=torch.float32)
    k_p_r_wk_ = torch.zeros(n_w_sum, dtype=torch.float32)
    k_p_w_wk_ = torch.zeros(n_w_sum, dtype=torch.float32)

    # NOTE: There's vectorized versions of this already in the existing EMPM repo, we can
    # copy it over later
    na = 0
    for nk_p_r in range(n_k_p_r):
        n_w = int(n_w_[nk_p_r].item())
        k_p_r = k_p_r_[nk_p_r].item()
        for nw in range(n_w):
            gamma_z = (2 * pi) * nw / max(1, n_w)
            cc = np.cos(gamma_z)
            sc = np.sin(gamma_z)
            k_c_0_wk_[na] = k_p_r * cc
            k_c_1_wk_[na] = k_p_r * sc
            k_p_r_wk_[na] = k_p_r
            k_p_w_wk_[na] = gamma_z
            na += 1

    if flag_verbose > 0: print(f" %% [finished get_weight_2d_2]")

    return (
        n_w_,
        weight_2d_k_p_r_,
        weight_2d_k_p_wk_,
        k_p_r_wk_,
        k_p_w_wk_,
        k_c_0_wk_,
        k_c_1_wk_,
    )
