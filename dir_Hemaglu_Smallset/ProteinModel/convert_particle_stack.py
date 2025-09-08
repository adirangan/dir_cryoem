import os, sys
from cryolike.convert_particle_stacks import convert_particle_stacks_from_star_files

pixel_size = 1.03 
n_stack = 10
input_folder = "./particles/" 
particle_file_list = [os.path.join(input_folder, "batch_%d.mrc" % i) for i in range(n_stack)]
star_file_list = [os.path.join(input_folder, "batch_%d.star" % i) for i in range(n_stack)]
#star_file = ["./mrcs/particles.star"]


convert_particle_stacks_from_star_files(
    params_input= "./output/parameters.npz",
    folder_output = "./output/particles/",
    particle_file_list = particle_file_list,
    star_file_list = star_file_list,
    pixel_size = pixel_size,
    defocus_angle_is_degree = True,
    phase_shift_is_degree = True,
    skip_exist = False,
) 
