function isosurface_cut_f_x_u_1(f_x_u___,percent_threshold_,n_cut,n_contour,c_80s__);

na=0;
if (nargin<1+na); f_x_u___=[]; end; na=na+1;
if (nargin<1+na); percent_threshold_=[]; end; na=na+1;
if (nargin<1+na); n_cut = []; end; na=na+1;
if (nargin<1+na); n_contour = []; end; na=na+1;
if (nargin<1+na); c_80s__ = []; end; na=na+1;
if isempty(percent_threshold_); percent_threshold_=94; end;
if isempty(n_cut); n_contour = 3; end;
if isempty(n_contour); n_contour = 32; end;
if isempty(c_80s__); c_80s__ = colormap_80s; end;
if numel(percent_threshold_)==1; percent_threshold_ = [percent_threshold_,100]; end;

xlim_ = [-1,+1];
ylim_ = [-1,+1];
zlim_ = [-1,+1];

if (ndims(f_x_u___)~=3);
n_d = round(numel(f_x_u___).^(1/3));
f_x_u___ = reshape(f_x_u___,[n_d,n_d,n_d]);
end;%if (ndims(f_x_u___)~=3);
assert(ndims(f_x_u___)==3);

if (~isreal(f_x_u___));
f_x_u___ = real(f_x_u___);
end;%if (~isreal(f_x_u___));

[n_x_c,n_y_c,n_z_c] = size(f_x_u___);
X_ = transpose(linspace(min(xlim_),max(xlim_),n_x_c));
Y_ = transpose(linspace(min(ylim_),max(ylim_),n_y_c));
Z_ = transpose(linspace(min(zlim_),max(zlim_),n_z_c));
X_cut_ = transpose(linspace(min(xlim_),max(xlim_),2+n_cut));
X_cut_ = X_cut_(2:end-1);
dX_cut = 0; if (numel(X_cut_)>1); dX_cut = mean(diff(X_cut_)); end;
X_expand = 1.00;

f_cut_xyz___ = interp1(X_,f_x_u___,X_cut_);

c_80s__ = colormap_80s;
n_c_80s = size(c_80s__,1);
v_avg = mean(f_x_u___(:)); v_std = std(f_x_u___(:)); 
v_min = min(f_x_u___(:)); v_max = max(f_x_u___(:));
vlim_ = v_avg + 1.5*v_std*[-1,1];
%vlim_ = [v_min , v_max];
if ~isempty(percent_threshold_);
vlim_ = prctile(f_x_u___(:),percent_threshold_);
end;%if ~isempty(percent_threshold_);
v = min(vlim_);

for ncut=0:n_cut-1;
X_cut = X_cut_(1+ncut);
f_cut_yz__ = squeeze(f_cut_xyz___(1+ncut,:,:));
[c_cut__] = contourc(f_cut_yz__,linspace(min(vlim_),max(vlim_),n_contour));
[n_item,val_i_,length_i_,x_il__,y_il__] = cell_from_contour_0(c_cut__);
hold on;
for nitem=0:n_item-1;
val = val_i_(1+nitem);
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*(val-min(vlim_))/diff(vlim_))));
tmp_length = length_i_(1+nitem);
tmp_x_ = x_il__{1+nitem}; tmp_x_ = (tmp_x_ - n_x_c/2)*2/n_x_c;
tmp_y_ = y_il__{1+nitem}; tmp_y_ = (tmp_y_ - n_y_c/2)*2/n_y_c;
p = patch( zeros(tmp_length,1) + X_expand*dX_cut*(1+ncut) -  X_expand*dX_cut*0.125*nitem/n_item , tmp_x_ , tmp_y_ , c_80s__(1+nc_80s,:));
if (nitem>0); set(p,'LineStyle','none'); end;
end;%for nitem=0:n_item-1;
hold off;
end;%for ncut=0:n_cut-1;

xlim([0,X_expand*dX_cut*n_cut]);%xlim(xlim_);
ylim(ylim_);
zlim(zlim_);
xlabel('x');ylabel('y');zlabel('z');
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
%daspect([1,4,4]); 
view([-35,20]); 
axis vis3d;
camlight left; lighting gouraud;
figbig;
