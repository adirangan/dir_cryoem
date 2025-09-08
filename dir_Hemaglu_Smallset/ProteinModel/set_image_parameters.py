import os
from numpy import pi

from cryolike.metadata import ImageDescriptor

verbose = True
if verbose:
    print("Setting parameters...")

n_voxels = 256
voxel_size = 1.03
precision = 'single' # 'single' or 'double'
resolution_factor = 0.5#1.0
viewing_distance = 1.0 / (4.0 * pi)
n_inplanes = 256
folder_output = 'output'

os.makedirs(folder_output, exist_ok=True)
image_parameters_filename = os.path.join(folder_output, "parameters.npz")
image_parameters = ImageDescriptor.from_individual_values(
    n_pixels = n_voxels,
    pixel_size = voxel_size,
    resolution_factor = resolution_factor,
    precision = precision,
    viewing_distance = viewing_distance,
    n_inplanes = n_inplanes,
    use_protein_residue_model = True,
    atom_shape = 'gaussian'
)
if verbose:
    image_parameters.print()
image_parameters.save(image_parameters_filename)
