function imagesc_c_mask(n_x,x_,n_y,y_,S_c_,clim_,colormap__,radius_mask);
% cartesian imagesc;
% assumes S_c_(nx,ny) = S_c_(nx + ny*n_x);

na=0;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); x_=[]; end; na=na+1;
if (nargin<1+na); n_y=[]; end; na=na+1;
if (nargin<1+na); y_=[]; end; na=na+1;
if (nargin<1+na); S_c_=[]; end; na=na+1;
if (nargin<1+na); clim_=[]; end; na=na+1;
if (nargin<1+na); colormap__=[]; end; na=na+1;
if (nargin<1+na); radius_mask=[]; end; na=na+1;

if isempty(radius_mask); radius_mask = 1; end;

if isempty(clim_); 
clim_ = mean(S_c_(:)) + std(S_c_(:))*2.5*[-1,1]; 
end;%if isempty(clim_);
if isempty(colormap__); 
colormap__ = colormap_beach();
end;%if isempty(colormap__);

dx_ = diff(x_); dy_ = diff(y_);
x0_ = zeros(1,n_x*n_y);
x1_ = zeros(1,n_x*n_y);
y0_ = zeros(1,n_x*n_y);
y1_ = zeros(1,n_x*n_y);
ncra = size(colormap__,1);
C_ = zeros(1,n_x*n_y,3);
ic=0; jc=0;
for ny=0:n_y-1;
if (ny==0); y_pre = y_(1)-0.5*dy_(1); else; y_pre = 0.5*(y_(1+ny-1)+y_(1+ny)); end;
if (ny==n_y-1); y_pos = y_(end)+0.5*dy_(end); else; y_pos = 0.5*(y_(1+ny+1)+y_(1+ny)); end;
y_mid = 0.5*(y_pre + y_pos);
for nx=0:n_x-1;
if (nx==0); x_pre = x_(1)-0.5*dx_(1); else; x_pre = 0.5*(x_(1+nx-1)+x_(1+nx)); end;
if (nx==n_x-1); x_pos = x_(end)+0.5*dx_(end); else; x_pos = 0.5*(x_(1+nx+1)+x_(1+nx)); end;
x_mid = 0.5*(x_pre + x_pos);
if (sqrt(x_mid.^2 + y_mid.^2)<=radius_mask);
x0_(1+jc) = x_pre; x1_(1+jc) = x_pos;
y0_(1+jc) = y_pre; y1_(1+jc) = y_pos;
nc = max(1,min(ncra,floor(ncra*(S_c_(1+ic) - min(clim_))/diff(clim_))));
C_(1,1+jc,1:3) = colormap__(nc,:);
jc = jc+1;
end;%if (sqrt(x_mid.^2 + y_mid.^2)<=radius_mask);
ic = ic+1;
end;%for nx=0:n_x-1;
end;%for ny=0:n_y-1;
x0_ = x0_(1:jc); x1_ = x1_(1:jc);
y0_ = y0_(1:jc); y1_ = y1_(1:jc);
C_ = C_(1,1:jc,:);
X_ = [x0_ ; x0_ ; x1_ ; x1_ ; x0_];
Y_ = [y0_ ; y1_ ; y1_ ; y0_ ; y0_];
p=patch(X_,Y_,C_); set(p,'EdgeColor','none');

