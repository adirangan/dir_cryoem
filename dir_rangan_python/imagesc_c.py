import numpy as np ; pi = np.pi ; i = 1j ; import torch ; import timeit ;
from matlab_index_2d_0 import matlab_index_2d_0 ;
from matlab_index_3d_0 import matlab_index_3d_0 ;
from matlab_index_4d_0 import matlab_index_4d_0 ;
from matlab_scalar_round import matlab_scalar_round ;
from i4_torch_arange import i4_torch_arange ;
from fnorm_disp import fnorm_disp ;
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
mtr = lambda a : tuple(reversed(a)) ; #<-- matlab-arranged size (i.e., tuple(reversed(...))). ;
msr = lambda str : str[::-1] ; #<-- for einsum (i.e., string reversed (...)). ;
mts = lambda a : tuple(len(a) - x - 1 for x in a) ; #<-- for permute (i.e., tuple subtract (...)). ;
n_1 = 1; n_2 = 2; n_3 = 3;

'''
function imagesc_c(n_x,x_,n_y,y_,S_c__,clim,cra_);
% cartesian imagesc;
% assumes S_c__(nx,ny) = S_c__(nx + ny*n_x);

if (nargin<1);
n_x = 24+1; x_ = transpose(linspace(-1,+1,n_x));
n_y = 16+1; y_ = transpose(linspace(-1,+1,n_y));
[x__,y__] = ndgrid(x_,y_);
S_c__ = (x__ + y__)./2;
imagesc_c(n_x,x_,n_y,y_,S_c__,[-1,+1],colormap_beach());
axis image; axisnotick;
return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); x_=[]; end; na=na+1;
if (nargin<1+na); n_y=[]; end; na=na+1;
if (nargin<1+na); y_=[]; end; na=na+1;
if (nargin<1+na); S_c__=[]; end; na=na+1;
if (nargin<1+na); clim=[]; end; na=na+1;
if (nargin<1+na); cra_=[]; end; na=na+1;

if isempty(S_c__); S_c__ = zeros(2,2); end;
if isempty(n_x); if size(S_c__,2)>1; n_x = size(S_c__,1); else; n_x = fix(sqrt(numel(S_c__))); end; end;
if isempty(n_y); if size(S_c__,2)>1; n_y = size(S_c__,2); else; n_x = fix(numel(S_c__)/n_x); end; end;
if isempty(x_); x_ = linspace(-1,+1,n_x); end;
if isempty(y_); y_ = linspace(-1,+1,n_y); end;

if isempty(clim); 
clim = mean(S_c__(:)) + std(S_c__(:))*2.5*[-1,1]; 
end;%if isempty(clim);
if isempty(cra_); 
cra_ = colormap_beach();
end;%if isempty(clim);
dx_ = diff(x_); dy_ = diff(y_);
x0_ = zeros(1,n_x*n_y);
x1_ = zeros(1,n_x*n_y);
y0_ = zeros(1,n_x*n_y);
y1_ = zeros(1,n_x*n_y);
ncra = size(cra_,1);
C_ = zeros(1,n_x*n_y,3);
ic=0;
for ny=0:n_y-1;
if (ny==0); y_pre = y_(1)-0.5*dy_(1); else; y_pre = 0.5*(y_(1+ny-1)+y_(1+ny)); end;
if (ny==n_y-1); y_pos = y_(end)+0.5*dy_(end); else; y_pos = 0.5*(y_(1+ny+1)+y_(1+ny)); end;
for nx=0:n_x-1;
if (nx==0); x_pre = x_(1)-0.5*dx_(1); else; x_pre = 0.5*(x_(1+nx-1)+x_(1+nx)); end;
if (nx==n_x-1); x_pos = x_(end)+0.5*dx_(end); else; x_pos = 0.5*(x_(1+nx+1)+x_(1+nx)); end;
x0_(1+ic) = x_pre; x1_(1+ic) = x_pos;
y0_(1+ic) = y_pre; y1_(1+ic) = y_pos;
nc = max(1,min(ncra,floor(ncra*(S_c__(1+ic) - min(clim))/diff(clim))));
C_(1,1+ic,1:3) = cra_(nc,:);
ic = ic+1;
end;%for nx=0:n_x-1;
end;%for ny=0:n_y-1;
X_ = [x0_ ; x0_ ; x1_ ; x1_ ; x0_];
Y_ = [y0_ ; y1_ ; y1_ ; y0_ ; y0_];
p=patch(X_,Y_,C_); set(p,'EdgeColor','none');
'''
import matplotlib.pyplot as plt

def imagesc_c(
        ax,
        n_x=None,
        x_=None,
        n_y=None,
        y_=None,
        S_c__=None,
        clim_=None,
        c_3c__=None,
):
    if S_c__ is None: S_c__ = torch.zeros(mtr((2, 2))).to(dtype=torch.float32);
    if n_x is None: n_x = int(S_c__.size()[0]) ;
    if n_y is None: n_y = int(S_c__.size()[1]) ;
    if x_ is None: x_ = torch.linspace(-1, +1, n_x).to(dtype=torch.float32) ;
    if y_ is None: y_ = torch.linspace(-1, +1, n_y).to(dtype=torch.float32) ;
    if clim_ is None: clim_ = torch.mean(S_c__).item() + 2.5 * torch.std(S_c__).item() * torch.tensor([-1,+1]).to(dtype=torch.float32);
    if c_3c__ is None: c_3c__ = colormap_beach() ; # Default colormap

    dx_ = torch.diff(x_.ravel()).to(dtype=torch.float32) ;
    x_pre_ = x_ -  0.5 * torch.concatenate( ( torch.tensor([dx_[0]]) , dx_ ) , 0 );
    x_pos_ = x_ +  0.5 * torch.concatenate( ( dx_ , torch.tensor([dx_[n_x-1-1]]) ) , 0 );
    dy_ = torch.diff(y_.ravel()).to(dtype=torch.float32) ;
    y_pre_ = y_ -  0.5 * torch.concatenate( ( torch.tensor([dy_[0]]) , dy_ ) , 0 );
    y_pos_ = y_ +  0.5 * torch.concatenate( ( dy_ , torch.tensor([dy_[n_y-1-1]]) ) , 0 );
    y_prepre__,x_prepre__ = torch.meshgrid(y_pre_,x_pre_,indexing='ij');
    y_prepos__,x_prepos__ = torch.meshgrid(y_pos_,x_pre_,indexing='ij');
    y_pospre__,x_pospre__ = torch.meshgrid(y_pre_,x_pos_,indexing='ij');
    y_pospos__,x_pospos__ = torch.meshgrid(y_pos_,x_pos_,indexing='ij');
    n_c = int(c_3c__.size()[0]);
    nc_ = torch.maximum(torch.tensor(0.0),torch.minimum(torch.tensor(n_c-1),torch.floor(n_c*(S_c__.ravel()-clim_[0].item())/np.maximum(1e-12,clim_[1].item()-clim_[0].item())))).to(dtype=torch.int32);
    C_3a__ = torch.zeros(mtr((n_3,n_x*n_y))).to(dtype=torch.float32);
    patches = [];
    for na in range(n_x*n_y):
        nc = int(nc_[na].item());
        C_3a__[na,:] = c_3c__[nc,:];
        polygon = Polygon( [ [x_prepre__.ravel()[na].item(),y_prepre__.ravel()[na].item()] , [x_prepos__.ravel()[na].item(),y_prepos__.ravel()[na].item()] , [x_pospos__.ravel()[na].item(),y_pospos__.ravel()[na].item()] , [x_pospre__.ravel()[na].item(),y_pospre__.ravel()[na].item()] ], closed=True ) ;
        patches.append(polygon) ;
    #end;%for na in range(n_x*n_y):

    p = PatchCollection(patches, edgecolor='none', linewidth=0) ;
    p.set_facecolor(C_3a__.numpy()) ;
    ax.add_collection(p) ;
    ax.set_aspect('equal') ;
