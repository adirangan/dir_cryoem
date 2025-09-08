import mrcfile
import matplotlib.pyplot as plt
import numpy as np

# Open MRC file
with mrcfile.open('batch_9.mrc', permissive=True) as mrc:
    volume = mrc.data  # shape: (Z, Y, X)

# Take 100 slices evenly spaced along the z-axis
num_slices = 100
z_indices = np.linspace(0, volume.shape[0] - 1, num_slices, dtype=int)
slices = volume[z_indices]

# Plot in 10x10 grid
fig, axes = plt.subplots(10, 10, figsize=(15, 15))
for ax, img in zip(axes.flat, slices):
    ax.imshow(img, cmap='gray')
    ax.axis('off')

plt.tight_layout()
plt.savefig('mrc_grid.png', dpi=300)
plt.close()

