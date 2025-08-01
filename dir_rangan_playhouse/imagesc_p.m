function clim = imagesc_p(n_r,grid_p_,n_w_,n_A,S_p_,clim,c__);
% polar imagesc;

if (nargin<1);
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
return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); n_r=[]; end; na=na+1;
if (nargin<1+na); grid_p_=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); n_A=[]; end; na=na+1;
if (nargin<1+na); S_p_=[]; end; na=na+1;
if (nargin<1+na); clim_=[]; end; na=na+1;
if (nargin<1+na); c__=[]; end; na=na+1;
if isempty(clim); 
clim = mean(S_p_) + std(S_p_)*2.5*[-1,1]; 
end;%if isempty(clim);
if isempty(c__); 
c__ = colormap_beach();
end;%if isempty(clim);

r0_ = zeros(1,n_A); r1_ = zeros(1,n_A);
w0_ = zeros(1,n_A); w1_ = zeros(1,n_A);
ncra = size(c__,1);
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
C_(1,1+ic,1:3) = c__(nc,:);
ic=ic+1;
end;%for nw=0:n_w_(1+nr)-1;
end;%for nr=0:n_r-1;
X_ = [r0_.*cos(w0_) ; r1_.*cos(w0_) ; r1_.*cos(w1_) ; r0_.*cos(w1_) ; r0_.*cos(w0_) ];
Y_ = [r0_.*sin(w0_) ; r1_.*sin(w0_) ; r1_.*sin(w1_) ; r0_.*sin(w1_) ; r0_.*sin(w0_) ];
p=patch(X_,Y_,C_); set(p,'EdgeColor','none');
