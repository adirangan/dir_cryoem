import pandas as pd
import mrcfile
import numpy as np
from pathlib import Path
from collections import defaultdict

# Path to your star file
star_path = Path("particles_smallset.star")

def parse_star_particles_section(star_path):
    with open(star_path, "r") as f:
        lines = f.readlines()

    # Find the start of data_particles
    particle_start = next(i for i, line in enumerate(lines) if "data_particles" in line)
    
    # Find the next loop_ after data_particles
    loop_start = next(i for i in range(particle_start, len(lines)) if lines[i].strip() == "loop_")

    # Parse headers
    headers = []
    i = loop_start + 1
    while lines[i].strip().startswith("_"):
        headers.append(lines[i].strip())
        i += 1

    # Now parse data
    data = []
    while i < len(lines) and not lines[i].startswith("data_") and lines[i].strip():
        data.append(lines[i].strip().split())
        i += 1

    df = pd.DataFrame(data, columns=headers)
    return df

df = parse_star_particles_section(star_path)

# Extract rlnImageName column
rln_image_col = [col for col in df.columns if "rlnImageName" in col][0]
entries = df[rln_image_col]

# Build mapping from mrcs path to list of zero-based indices
particles_dict = defaultdict(list)
for entry in entries:
    idx_str, path_str = entry.split("@")
    idx = int(idx_str) - 1  # Convert to 0-based index
    path = Path(path_str)
    particles_dict[path].append(idx)

# Create output directory
output_dir = Path("filtered_mrcs")
output_dir.mkdir(exist_ok=True)

# Filter and write new .mrcs files
for mrcs_path, indices in particles_dict.items():
    print(f"Processing {mrcs_path} with {len(indices)} particles")
    with mrcfile.open(mrcs_path, permissive=True) as mrc:
        data = mrc.data
        selected = data[indices]
        output_path = output_dir / mrcs_path.name
        with mrcfile.new(output_path, overwrite=True) as out_mrc:
            out_mrc.set_data(np.asarray(selected, dtype=np.float32))

print("Done. Filtered .mrcs files saved to 'filtered_mrcs/'")

