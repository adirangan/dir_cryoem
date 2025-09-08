import numpy as np

'''
n_x = 24+1; x_ = transpose(linspace(-1,+1,n_x));
n_y = 16+1; y_ = transpose(linspace(-1,+1,n_y));
[x__,y__] = ndgrid(x_,y_);
S_c__ = (x__ + y__)./2;
imagesc_c(n_x,x_,n_y,y_,S_c__,[-1,+1],colormap_beach());
axis image; axisnotick;
return;
'''
import matplotlib.pyplot as plt
from imagesc_c import imagesc_c
from colormap_beach import colormap_beach
from colormap_pm import colormap_pm
from colormap_80s import colormap_80s
from colormap_81s import colormap_81s

n_x = 24 + 1
x_ = np.linspace(-1, +1, n_x)
n_y = 16 + 1
y_ = np.linspace(-1, +1, n_y)
x__, y__ = np.meshgrid(x_, y_, indexing='ij')
S_c__ = (x__ + y__) / 2

fig = plt.figure(figsize=(8, 6))
p_row = 2; p_col = 2; np=0;
for np in range(p_row*p_col):
    ax = fig.add_subplot(p_row,p_col,1+np);
    if np==0: imagesc_c(ax,n_x,x_,n_y,y_,S_c__,[-1,+1],colormap_beach()); ax.set_title('colormap_beach')
    if np==1: imagesc_c(ax,n_x,x_,n_y,y_,S_c__,[-1,+1],colormap_pm()); ax.set_title('colormap_pm')
    if np==2: imagesc_c(ax,n_x,x_,n_y,y_,S_c__,[-1,+1],colormap_80s()); ax.set_title('colormap_80s')
    if np==3: imagesc_c(ax,n_x,x_,n_y,y_,S_c__,[-1,+1],colormap_81s()); ax.set_title('colormap_81s')
    ax.set_xlim([-1, 1]); ax.set_ylim([-1, 1]); 
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xticklabels([]); ax.set_yticklabels([]);
plt.show()
