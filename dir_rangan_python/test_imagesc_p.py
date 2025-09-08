import numpy as np

'''
n_r = 16;
grid_p_ = transpose([0:n_r-1]/max(1,n_r));
n_w_ = ceil(1+n_r*grid_p_);
n_A = sum(n_w_);
S_p_ = zeros(n_A,1);
ic=0;
for nr=0:n_r-1;
p = grid_p_(1+nr);
n_w = n_w_(1+nr);
for nw=0:n_w-1;
w = (2*pi*nw)/max(1,n_w);
S_p_(1+ic) = p*cos(w);
ic=ic+1;
end;%for nw=0:n_w-1;
end;%for nr=0:n_r-1;
imagesc_p(n_r,grid_p_,n_w_,n_A,S_p_,[-1,+1],colormap_beach());
'''
import matplotlib.pyplot as plt
from imagesc_p import imagesc_p
from colormap_beach import colormap_beach
from colormap_pm import colormap_pm
from colormap_80s import colormap_80s
from colormap_81s import colormap_81s

n_r = 16
grid_p_ = np.transpose(np.arange(0, n_r) / max(1, n_r))
n_w_ = np.ceil(1 + n_r * grid_p_).astype(int)
n_A = np.sum(n_w_)
S_p_ = np.zeros(n_A)
ic = 0

for nr in range(n_r):
    p = grid_p_[nr]
    n_w = n_w_[nr]
    for nw in range(n_w):
        w = (2 * np.pi * nw) / max(1, n_w)
        S_p_[ic] = p * np.cos(w)
        ic += 1



fig = plt.figure(figsize=(8, 6))
p_row = 2; p_col = 2; np=0;
for np in range(p_row*p_col):
    ax = fig.add_subplot(p_row,p_col,1+np);
    if np==0: imagesc_p(ax,n_r, grid_p_, n_w_, n_A, S_p_, [-1, +1], colormap_beach()); ax.set_title('colormap_beach')
    if np==1: imagesc_p(ax,n_r, grid_p_, n_w_, n_A, S_p_, [-1, +1], colormap_pm()); ax.set_title('colormap_pm')
    if np==2: imagesc_p(ax,n_r, grid_p_, n_w_, n_A, S_p_, [-1, +1], colormap_80s()); ax.set_title('colormap_80s')
    if np==3: imagesc_p(ax,n_r, grid_p_, n_w_, n_A, S_p_, [-1, +1], colormap_81s()); ax.set_title('colormap_81s')
    ax.set_xlim([-1, 1]); ax.set_ylim([-1, 1]); 
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xticklabels([]); ax.set_yticklabels([]);
plt.show()
