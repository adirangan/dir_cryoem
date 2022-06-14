function isosurface_cut_f_x_u_0(f_x_u___,percent_threshold_,n_contour,c_80s__);

na=0;
if (nargin<1+na); f_x_u___=[]; end; na=na+1;
if (nargin<1+na); percent_threshold_=[]; end; na=na+1;
if (nargin<1+na); n_contour = []; end; na=na+1;
if (nargin<1+na); c_80s__ = []; end; na=na+1;
if isempty(percent_threshold_); percent_threshold_=94; end;
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
[X___,Y___,Z___] = meshgrid(X_,Y_,Z_); %<-- note peculiarites of meshgrid. ;
index_h0_ = 0:n_x_c/2-1; index_h1_ = n_x_c/2 + index_h0_;
X0___ = X___(:,1+index_h0_,:); Y0___ = Y___(:,1+index_h0_,:); Z0___ = Z___(:,1+index_h0_,:); f0___ = f_x_u___(:,1+index_h0_,:);
X1___ = X___(:,1+index_h1_,:); Y1___ = Y___(:,1+index_h1_,:); Z1___ = Z___(:,1+index_h1_,:); f1___ = f_x_u___(:,1+index_h1_,:);

f_h01__ = squeeze(mean(f_x_u___(:,n_x_c/2-1 + [0,1],:),2));
f_h0__ = squeeze(mean(f_x_u___(:,n_x_c/2-1 + [0],:),2));
f_h1__ = squeeze(mean(f_x_u___(:,n_x_c/2-1 + [1],:),2));

xgap = 1.00;

n_c_80s = size(c_80s__,1);
v_avg = mean(f_x_u___(:)); v_std = std(f_x_u___(:)); 
v_min = min(f_x_u___(:)); v_max = max(f_x_u___(:));
vlim_ = v_avg + 1.5*v_std*[-1,1];
%vlim_ = [v_min , v_max];
if ~isempty(percent_threshold_);
vlim_ = prctile(f_x_u___(:),percent_threshold_);
end;%if ~isempty(percent_threshold_);

v = min(vlim_);
%%%%%%%%;
dx_h = xgap/2;
dx_g = diff(xlim_)/max(1,n_x_c-1);
%%%%%%%%;
[c_h0__] = contourc(f_h0__,linspace(min(vlim_),max(vlim_),n_contour)); [n_item,val_i_,length_i_,x_il__,y_il__] = cell_from_contour_0(c_h0__);
hold on;
for nitem=0:n_item-1;
val = val_i_(1+nitem);
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*(val-min(vlim_))/diff(vlim_))));
tmp_length = length_i_(1+nitem);
tmp_x_ = x_il__{1+nitem}; tmp_x_ = (tmp_x_ - n_x_c/2)*2/n_x_c;
tmp_y_ = y_il__{1+nitem}; tmp_y_ = (tmp_y_ - n_y_c/2)*2/n_y_c;
p = patch(zeros(tmp_length,1) - dx_h - 0.5*dx_g + dx_g*nitem/n_item,tmp_y_,tmp_x_,c_80s__(1+nc_80s,:)); set(p,'LineStyle','none');
%p = patch(zeros(tmp_length,1) + dx_h + 0.5*dx_g - dx_g*nitem/n_item,tmp_y_,tmp_x_,c_80s__(1+nc_80s,:)); set(p,'LineStyle','none');
end;%for nitem=0:n_item-1;
hold off;
%%%%%%%%;
[c_h1__] = contourc(f_h1__,linspace(min(vlim_),max(vlim_),n_contour)); [n_item,val_i_,length_i_,x_il__,y_il__] = cell_from_contour_0(c_h1__);
hold on;
for nitem=0:n_item-1;
val = val_i_(1+nitem);
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*(val-min(vlim_))/diff(vlim_))));
tmp_length = length_i_(1+nitem);
tmp_x_ = x_il__{1+nitem}; tmp_x_ = (tmp_x_ - n_x_c/2)*2/n_x_c;
tmp_y_ = y_il__{1+nitem}; tmp_y_ = (tmp_y_ - n_y_c/2)*2/n_y_c;
%p = patch(zeros(tmp_length,1) - dx_h - 0.5*dx_g + dx_g*nitem/n_item,tmp_y_,tmp_x_,c_80s__(1+nc_80s,:)); set(p,'LineStyle','none');
p = patch(zeros(tmp_length,1) + dx_h + 0.5*dx_g - dx_g*nitem/n_item,tmp_y_,tmp_x_,c_80s__(1+nc_80s,:)); set(p,'LineStyle','none');
end;%for nitem=0:n_item-1;
hold off;
%%%%%%%%;
for nh=0:1;
if (nh==0); Xj___ = X0___; Yj___ = Y0___; Zj___ = Z0___; fj___ = f0___; dx_h = -xgap/2; end;%if (nh==0);
if (nh==1); Xj___ = X1___; Yj___ = Y1___; Zj___ = Z1___; fj___ = f1___; dx_h = +xgap/2; end;%if (nh==1);
hpatch = patch(isosurface(dx_h + Xj___,Yj___,Zj___,fj___,v));
hold on;
isonormals(dx_h + Xj___,Yj___,Zj___,fj___,hpatch);
hold off;
nc_80s = 0;
hpatch.FaceColor = c_80s__(1+nc_80s,:);
hpatch.EdgeColor = 'none';
hpatch.FaceAlpha = 1;
end;%for nh=0:1;
%%%%%%%%;
xlim(1.25*xlim_);
ylim(ylim_);
zlim(zlim_);
xlabel('x');ylabel('y');zlabel('z');
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
%daspect([1,4,4]); 
view([-35,20]); 
axis vis3d;
camlight left; lighting gouraud;
figbig;
