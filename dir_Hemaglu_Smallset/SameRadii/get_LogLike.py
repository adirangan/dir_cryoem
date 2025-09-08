import os, sys
import torch
import numpy as np

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

#number of templates
numtemp = int(sys.argv[1])
tag_templates = sys.argv[2]

# Define the base directory (adjust this to your actual base directory)
base_dir = "./likelihood_%s/" % (tag_templates)

# Initialize a list to store the concatenated results for each template folder
all_templates_concat = []

# Iterate over the 80 template directories (from template0 to template79)
for x in range(0,numtemp):
    template_dir = os.path.join(base_dir, f"template{x}", "log_likelihood")

    # Initialize a list to store the 7 files within each template folder
    files_concat = []

    # Iterate over the 7 files (from log_likelihood_phys_S_stack0.pt to log_likelihood_phys_S_stack6.pt)
    for y in range(10):
#        file_path = os.path.join(template_dir, f"log_likelihood_optimal_fourier_stack_00000{y}.pt")
        file_path = os.path.join(template_dir, f"log_likelihood_integrated_fourier_stack_00000{y}.pt")
        # Load the .pt file (assuming each file contains a tensor)
        tensor = torch.load(file_path)

        # Append the tensor to the list
        files_concat.append(tensor)



    # Concatenate the 7 tensors along the first dimension
    concatenated = torch.cat(files_concat, dim=0)

    # Convert to NumPy and append the result to the final list
    concatenated_numpy = concatenated.cpu().numpy()
    all_templates_concat.append(concatenated_numpy)

# Convert the final list of arrays into a single array (concatenating along rows)
final_result_numpy = np.vstack(all_templates_concat)

# Save the logLike to a text file
np.savetxt("LogLikeMat", final_result_numpy.T)


