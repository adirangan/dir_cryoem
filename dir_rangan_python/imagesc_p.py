import numpy as np
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

'''
function clim = imagesc_p(n_r,grid_p_,n_w_,n_A,S_p_,clim,cra_);
% polar imagesc;

if (nargin<1);
n_r = 8;
grid_p_ = transpose([0:n_r-1]/max(1,n_r));
n_w_ = ceil(1+8*grid_p_);
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
return;
end;%if (nargin<1);

if (nargin<6); clim = []; end;
if (nargin<7); cra_ = []; end;
if isempty(clim); 
clim = mean(S_p_) + std(S_p_)*2.5*[-1,1]; 
end;%if isempty(clim);
if isempty(cra_); 
cra_ = colormap_beach();
end;%if isempty(clim);

r0_ = zeros(1,n_A); r1_ = zeros(1,n_A);
w0_ = zeros(1,n_A); w1_ = zeros(1,n_A);
ncra = size(cra_,1);
C_ = zeros(1,n_A,3);
ic=0;
for nr=0:n_r-1;
r = grid_p_(1+nr);
if nr==0; r_pre = grid_p_(1); else r_pre = 0.5*(grid_p_(1+nr-1) + grid_p_(1+nr)); end;
if nr==n_r-1; r_pos = grid_p_(end); else r_pos = 0.5*(grid_p_(1+nr+1) + grid_p_(1+nr)); end;
if n_r==1; r_pre = grid_p_(end)/4.0; r_pos = grid_p_(end); end;
n_w = n_w_(1+nr);
dw = 2*pi/n_w;
for nw=0:n_w_(1+nr)-1;
w_pre = nw*dw - 0.5*dw;
w_pos = nw*dw + 0.5*dw;
r0_(1+ic) = r_pre; r1_(1+ic) = r_pos;
w0_(1+ic) = w_pre; w1_(1+ic) = w_pos;
nc = max(1,min(ncra,floor(ncra*(S_p_(1+ic) - min(clim))/diff(clim))));
C_(1,1+ic,1:3) = cra_(nc,:);
ic=ic+1;
end;%for nw=0:n_w_(1+nr)-1;
end;%for nr=0:n_r-1;
X_ = [r0_.*cos(w0_) ; r1_.*cos(w0_) ; r1_.*cos(w1_) ; r0_.*cos(w1_) ; r0_.*cos(w0_) ];
Y_ = [r0_.*sin(w0_) ; r1_.*sin(w0_) ; r1_.*sin(w1_) ; r0_.*sin(w1_) ; r0_.*sin(w0_) ];
p=patch(X_,Y_,C_); set(p,'EdgeColor','none');
'''
import matplotlib.pyplot as plt

def imagesc_p(ax,n_r, grid_p_, n_w_, n_A, S_p_, clim=None, cra_=None):
    if clim is None:
        clim = np.mean(S_p_) + 2.5 * np.std(S_p_) * np.array([-1, 1])
    if cra_ is None:
        cra_ = plt.cm.viridis(np.linspace(0, 1, 256))  # Default colormap

    r0_ = np.zeros(n_A)
    r1_ = np.zeros(n_A)
    w0_ = np.zeros(n_A)
    w1_ = np.zeros(n_A)
    ncra = cra_.shape[0]
    C_ = np.zeros((n_A, 3))
    ic = 0
    patches = []

    for nr in range(n_r):
        r = grid_p_[nr]
        r_pre = grid_p_[0] if nr == 0 else 0.5 * (grid_p_[nr - 1] + grid_p_[nr])
        r_pos = grid_p_[-1] if nr == n_r - 1 else 0.5 * (grid_p_[nr + 1] + grid_p_[nr])
        if n_r == 1:
            r_pre = grid_p_[-1] / 4.0
            r_pos = grid_p_[-1]
        n_w = n_w_[nr]
        dw = 2 * np.pi / n_w

        for nw in range(n_w):
            w_pre = nw * dw - 0.5 * dw
            w_pos = nw * dw + 0.5 * dw
            r0_[ic] = r_pre
            r1_[ic] = r_pos
            w0_[ic] = w_pre
            w1_[ic] = w_pos
            nc = max(0, min(ncra - 1, int(ncra * (S_p_[ic] - clim[0]) / (clim[1] - clim[0]))))
            C_[ic, :] = cra_[nc, :]
            x = [r0_[ic] * np.cos(w0_[ic]), r1_[ic] * np.cos(w0_[ic]),
                 r1_[ic] * np.cos(w1_[ic]), r0_[ic] * np.cos(w1_[ic]), r0_[ic] * np.cos(w0_[ic])]
            y = [r0_[ic] * np.sin(w0_[ic]), r1_[ic] * np.sin(w0_[ic]),
                 r1_[ic] * np.sin(w1_[ic]), r0_[ic] * np.sin(w1_[ic]), r0_[ic] * np.sin(w0_[ic])]
            polygon = Polygon(np.column_stack((x, y)), True)
            patches.append(polygon)
            ic += 1

    p = PatchCollection(patches, edgecolor='none', linewidth=0)
    p.set_facecolor(C_)
    ax.add_collection(p)
    ax.set_aspect('equal')
