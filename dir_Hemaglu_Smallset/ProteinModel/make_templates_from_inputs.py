import os
import numpy as np
from cryolike.stacks import make_templates_from_inputs

verbose = True
if verbose:
    print("Making templates...")

map_folder = "./models/"
map_file_list = ["6wxb_CA_def.pdb",  "6wxb_CA.pdb"]
list_of_inputs = [os.path.join(map_folder, map_file) for map_file in map_file_list]

folder_params = 'output/'
image_parameters_filename = os.path.join(folder_params, "parameters.npz")

folder_output = 'output/templates/'

make_templates_from_inputs(
    list_of_inputs = list_of_inputs,
    image_parameters_file=image_parameters_filename,
    folder_output = folder_output,
    verbose = verbose
)
