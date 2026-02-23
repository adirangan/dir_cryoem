function imagesc_c(n_x_0,x_0_,n_x_1,x_1_,S_c__,clim_,c_use__);
% cartesian imagesc;
% assumes S_c__(nx_0,nx_1) = S_c__(nx_0 + nx_1*n_x_0);

if (nargin<1);
n_x_0 = 24+1; x_0_ = transpose(linspace(-1,+1,n_x_0));
n_x_1 = 16+1; x_1_ = transpose(linspace(-1,+1,n_x_1));
[x_0__,x_1__] = ndgrid(x_0_,x_1_);
S_c__ = (x_0__ + x_1__)./2;
imagesc_c(n_x_0,x_0_,n_x_1,x_1_,S_c__,[-1,+1],colormap_beach());
axis image; axisnotick;
return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); n_x_0=[]; end; na=na+1;
if (nargin<1+na); x_0_=[]; end; na=na+1;
if (nargin<1+na); n_x_1=[]; end; na=na+1;
if (nargin<1+na); x_1_=[]; end; na=na+1;
if (nargin<1+na); S_c__=[]; end; na=na+1;
if (nargin<1+na); clim_=[]; end; na=na+1;
if (nargin<1+na); c_use__=[]; end; na=na+1;

if isempty(S_c__); S_c__ = zeros(2,2); end;
if isempty(n_x_0); if size(S_c__,2)>1; n_x_0 = size(S_c__,1); else; n_x_0 = fix(sqrt(numel(S_c__))); end; end;
if isempty(n_x_1); if size(S_c__,2)>1; n_x_1 = size(S_c__,2); else; n_x_0 = fix(numel(S_c__)/n_x_0); end; end;
if isempty(x_0_); x_0_ = linspace(-1,+1,n_x_0); end;
if isempty(x_1_); x_1_ = linspace(-1,+1,n_x_1); end;

if isempty(clim_); 
clim_ = mean(S_c__(:)) + std(S_c__(:))*2.5*[-1,1]; 
end;%if isempty(clim_);
if isempty(c_use__); 
c_use__ = colormap_beach();
end;%if isempty(clim_);

n_c_use = size(c_use__,1);
x_0_ = x_0_(:); x_1_ = x_1_(:);

x_0_ext_ = [x_0_(1+0) - (x_0_(1+1)-x_0_(1+0)) ; x_0_ ; x_0_(end) + (x_0_(end)-x_0_(end-1)) ];
x_1_ext_ = [x_1_(1+0) - (x_1_(1+1)-x_1_(1+0)) ; x_1_ ; x_1_(end) + (x_1_(end)-x_1_(end-1)) ];
[x_0_ext__,x_1_ext__] = ndgrid(x_0_ext_,x_1_ext_);
x_00__ = 0.5*(x_0_ext__(1 + [0:n_x_0-1],2 + [0:n_x_1-1]) + x_0_ext__(2 + [0:n_x_0-1],2 + [0:n_x_1-1]));
x_01__ = 0.5*(x_0_ext__(2 + [0:n_x_0-1],2 + [0:n_x_1-1]) + x_0_ext__(3 + [0:n_x_0-1],2 + [0:n_x_1-1]));
x_10__ = 0.5*(x_1_ext__(2 + [0:n_x_0-1],1 + [0:n_x_1-1]) + x_1_ext__(2 + [0:n_x_0-1],2 + [0:n_x_1-1]));
x_11__ = 0.5*(x_1_ext__(2 + [0:n_x_0-1],2 + [0:n_x_1-1]) + x_1_ext__(2 + [0:n_x_0-1],3 + [0:n_x_1-1]));
x_00_ = reshape(x_00__,[1,n_x_0*n_x_1]);
x_01_ = reshape(x_01__,[1,n_x_0*n_x_1]);
x_10_ = reshape(x_10__,[1,n_x_0*n_x_1]);
x_11_ = reshape(x_11__,[1,n_x_0*n_x_1]);
nc__ = max(0,min(n_c_use-1,floor(n_c_use*(S_c__-min(clim_))/diff(clim_))));
C___ = reshape(c_use__(1+nc__(:),:),[1,n_x_0*n_x_1,3]);
X_0__ = [x_00_ ; x_00_ ; x_01_ ; x_01_ ; x_00_];
X_1__ = [x_10_ ; x_11_ ; x_11_ ; x_10_ ; x_10_];
p=patch(X_0__,X_1__,C___); set(p,'EdgeColor','none');

%{
x_00_ = zeros(1,n_x_0*n_x_1);
x_01_ = zeros(1,n_x_0*n_x_1);
x_10_ = zeros(1,n_x_0*n_x_1);
x_11_ = zeros(1,n_x_0*n_x_1);
C_ = zeros(1,n_x_0*n_x_1,3);
ic=0;
for nx_1=0:n_x_1-1;
if (nx_1==0); x_1_pre = x_1_(1)-0.5*dx_1_(1); else; x_1_pre = 0.5*(x_1_(1+nx_1-1)+x_1_(1+nx_1)); end;
if (nx_1==n_x_1-1); x_1_pos = x_1_(end)+0.5*dx_1_(end); else; x_1_pos = 0.5*(x_1_(1+nx_1+1)+x_1_(1+nx_1)); end;
for nx_0=0:n_x_0-1;
if (nx_0==0); x_0_pre = x_0_(1)-0.5*dx_0_(1); else; x_0_pre = 0.5*(x_0_(1+nx_0-1)+x_0_(1+nx_0)); end;
if (nx_0==n_x_0-1); x_0_pos = x_0_(end)+0.5*dx_0_(end); else; x_0_pos = 0.5*(x_0_(1+nx_0+1)+x_0_(1+nx_0)); end;
x_00_(1+ic) = x_0_pre; x_01_(1+ic) = x_0_pos;
x_10_(1+ic) = x_1_pre; x_11_(1+ic) = x_1_pos;
nc = max(1,min(n_c_use,floor(n_c_use*(S_c__(1+ic) - min(clim_))/diff(clim_))));
C_(1,1+ic,1:3) = c_use__(nc,:);
ic = ic+1;
end;%for nx_0=0:n_x_0-1;
end;%for nx_1=0:n_x_1-1;
X_0_ = [x_00_ ; x_00_ ; x_01_ ; x_01_ ; x_00_];
X_1_ = [x_10_ ; x_11_ ; x_11_ ; x_10_ ; x_10_];
p=patch(X_0_,X_1_,C_); set(p,'EdgeColor','none');
%}
