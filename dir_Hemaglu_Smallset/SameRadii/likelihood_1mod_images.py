import sys, os
from cryolike.run_likelihood import run_likelihood

#tag_templates = "templates_pix128_gaussian_r_3_martini_50"


i_template = int(sys.argv[1])
tag_templates = sys.argv[2]
#tag_templates = "templates_pix128_gaussian_r_3_martini_50"

if tag_templates[-1] == '/':
    tag_templates = tag_templates[:-1]

folder_particles = "./output/particles/" 
folder_output = "./likelihood_%s/" % (tag_templates)

print("folder_output = %s" % folder_output)

run_likelihood(
    params_input = "./output/parameters.npz",
    folder_templates = "./output/templates/",
    folder_particles = folder_particles, 
    folder_output = folder_output,
    i_template = i_template,
    n_stacks = 10,
    skip_exist = False,
    n_templates_per_batch = 16,
    n_images_per_batch = 128,
    search_batch_size = True,
    max_displacement_pixels = 31,
    n_displacements_x = 32,
    n_displacements_y = 32,
    return_likelihood_integrated_pose_fourier = True,
    return_likelihood_optimal_pose_physical = False,
    return_likelihood_optimal_pose_fourier = True,
    verbose = True
)
