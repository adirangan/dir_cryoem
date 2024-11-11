function ...
[ ...
 parameter ...
,v_ ...
] = ...
isosurface_f_x_u_2( ...
 parameter ...
,f_x_u___ ...
,g_x_u___ ...
);

str_thisfunction = 'isosurface_f_x_u_2';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); f_x_u___=[]; end; na=na+1;
if (nargin<1+na); g_x_u___=[]; end; na=na+1;

flag_multicolor = ~isempty(g_x_u___);

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
tolerance_master = parameter.tolerance_master;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'flag_projection'); parameter.flag_projection = 0; end;
flag_projection = parameter.flag_projection;
if ~isfield(parameter,'flag_collapse'); parameter.flag_collapse = 1; end;
flag_collapse = parameter.flag_collapse;
if ~isfield(parameter,'flag_boxgrid'); parameter.flag_boxgrid = 1; end;
flag_boxgrid = parameter.flag_boxgrid;
if ~isfield(parameter,'x_0_lim_'); parameter.x_0_lim_ = [-1,+1]; end;
x_0_lim_ = parameter.x_0_lim_;
if ~isfield(parameter,'x_1_lim_'); parameter.x_1_lim_ = [-1,+1]; end;
x_1_lim_ = parameter.x_1_lim_;
if ~isfield(parameter,'x_2_lim_'); parameter.x_2_lim_ = [-1,+1]; end;
x_2_lim_ = parameter.x_2_lim_;
if ~isfield(parameter,'percent_threshold_'); parameter.percent_threshold_ = 98.5; end;
percent_threshold_ = parameter.percent_threshold_;
if ~isfield(parameter,'vval_'); parameter.vval_ = prctile(real(f_x_u___(:)),percent_threshold_); end;
vval_ = parameter.vval_;
if ~isfield(parameter,'vlim_'); parameter.vlim_ = prctile(real(f_x_u___(:)), [1,99]); end;
vlim_ = parameter.vlim_;
if ~isfield(parameter,'v_alpha'); parameter.v_alpha = 1; end;
v_alpha = parameter.v_alpha;
if ~isfield(parameter,'c_use__'); parameter.c_use__ = flipud(colormap('spring')); end;
c_use__ = parameter.c_use__; n_c_use = size(c_use__,1);
if ~isfield(parameter,'vlim_g_'); parameter.vlim_g_ = []; end;
vlim_g_ = parameter.vlim_g_;
if ~isfield(parameter,'c_use_g__'); parameter.c_use_g__ = flipud(colormap_pm); end;
c_use_g__ = parameter.c_use_g__; n_c_use_g = size(c_use_g__,1);

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if (ndims(f_x_u___)~=3);
n_d = round(numel(f_x_u___).^(1/3));
if (numel(f_x_u___)~=n_d^3); disp(sprintf(' %% Warning, improper dimension in %s',str_thisfunction)); end;
f_x_u___ = reshape(f_x_u___,[n_d,n_d,n_d]);
end;%if (ndims(f_x_u___)~=3);
assert(ndims(f_x_u___)==3);

if (~isreal(f_x_u___));
f_x_u___ = real(f_x_u___);
end;%if (~isreal(f_x_u___));

[n_x_0,n_x_1,n_x_2] = size(f_x_u___);
x_0_ = transpose(linspace(min(x_0_lim_),max(x_0_lim_),1+n_x_0)); x_0_ = x_0_(1:end-1); dx_0 = mean(diff(x_0_));
x_1_ = transpose(linspace(min(x_1_lim_),max(x_1_lim_),1+n_x_1)); x_1_ = x_1_(1:end-1); dx_1 = mean(diff(x_1_));
x_2_ = transpose(linspace(min(x_2_lim_),max(x_2_lim_),1+n_x_2)); x_2_ = x_2_(1:end-1); dx_2 = mean(diff(x_2_));
[x_0___,x_1___,x_2___] = meshgrid(x_0_,x_1_,x_2_); %<-- note that isosurface requires meshgrid, not ndgrid. ;

%%%%%%%%;
if flag_multicolor;
if (ndims(g_x_u___)~=3);
n_d = round(numel(g_x_u___).^(1/3));
if (numel(g_x_u___)~=n_d^3); disp(sprintf(' %% Warning, improper dimension in %s',str_thisfunction)); end;
g_x_u___ = reshape(g_x_u___,[n_d,n_d,n_d]);
end;%if (ndims(g_x_u___)~=3);
assert(ndims(g_x_u___)==3);
if (~isreal(g_x_u___));
g_x_u___ = real(g_x_u___);
end;%if (~isreal(g_x_u___));
[n_x_g_0,n_x_g_1,n_x_g_2] = size(g_x_u___);
x_g_0_ = transpose(linspace(min(x_0_lim_),max(x_0_lim_),1+n_x_g_0)); x_g_0_ = x_g_0_(1:end-1); dx_g_0 = mean(diff(x_g_0_));
x_g_1_ = transpose(linspace(min(x_1_lim_),max(x_1_lim_),1+n_x_g_1)); x_g_1_ = x_g_1_(1:end-1); dx_g_1 = mean(diff(x_g_1_));
x_g_2_ = transpose(linspace(min(x_2_lim_),max(x_2_lim_),1+n_x_g_2)); x_g_2_ = x_g_2_(1:end-1); dx_g_2 = mean(diff(x_g_2_));
[x_g_0___,x_g_1___,x_g_2___] = meshgrid(x_g_0_,x_g_1_,x_g_2_); %<-- note that isosurface requires meshgrid, not ndgrid. ;
if isempty(vlim_g_); vlim_g_ = prctile(real(g_x_u___(:)), [1,99]); end;
end;%if flag_multicolor;

if  flag_projection;
%%%%%%%%;
f_x_u_01__ = mean(f_x_u___,1+2);
nc_01__ = max(0,min(n_c_use-1,floor(n_c_use*(f_x_u_01__(:)-min(vlim_))/max(1e-12,diff(vlim_)))));
[x_0_01__,x_1_01__] = ndgrid(x_0_,x_1_);
tmp_x_0__ = transpose( [ x_0_01__(:)-0.5*dx_0 , x_0_01__(:)+0.5*dx_0 , x_0_01__(:)+0.5*dx_0 , x_0_01__(:)-0.5*dx_0 ] );
tmp_x_1__ = transpose( [ x_1_01__(:)-0.5*dx_1 , x_1_01__(:)-0.5*dx_1 , x_1_01__(:)+0.5*dx_1 , x_1_01__(:)+0.5*dx_1 ] );
tmp_x_2__ = transpose( (min(x_2_lim_)+0.0*dx_2)*ones(n_x_0*n_x_1,4) );
tmp_c___ = reshape(c_use__(1+nc_01__(:),:),[1,n_x_0*n_x_1,3]);
hold on;
patch(tmp_x_0__,tmp_x_1__,tmp_x_2__,tmp_c___,'LineStyle','none','EdgeColor','none');
hold off;
%%%%%%%%;
f_x_u_02__ = mean(f_x_u___,1+1);
nc_02__ = max(0,min(n_c_use-1,floor(n_c_use*(f_x_u_02__(:)-min(vlim_))/max(1e-12,diff(vlim_)))));
[x_0_02__,x_2_02__] = ndgrid(x_0_,x_2_);
tmp_x_0__ = transpose( [ x_0_02__(:)-0.5*dx_0 , x_0_02__(:)+0.5*dx_0 , x_0_02__(:)+0.5*dx_0 , x_0_02__(:)-0.5*dx_0 ] );
tmp_x_2__ = transpose( [ x_2_02__(:)-0.5*dx_2 , x_2_02__(:)-0.5*dx_2 , x_2_02__(:)+0.5*dx_2 , x_2_02__(:)+0.5*dx_2 ] );
tmp_x_1__ = transpose( (max(x_1_lim_)-0.0*dx_1)*ones(n_x_0*n_x_2,4) );
tmp_c___ = reshape(c_use__(1+nc_02__(:),:),[1,n_x_0*n_x_2,3]);
hold on;
patch(tmp_x_0__,tmp_x_1__,tmp_x_2__,tmp_c___,'LineStyle','none','EdgeColor','none');
hold off;
%%%%%%%%;
f_x_u_12__ = mean(f_x_u___,1+0);
nc_12__ = max(0,min(n_c_use-1,floor(n_c_use*(f_x_u_12__(:)-min(vlim_))/max(1e-12,diff(vlim_)))));
[x_1_12__,x_2_12__] = ndgrid(x_1_,x_2_);
tmp_x_1__ = transpose( [ x_1_12__(:)-0.5*dx_1 , x_1_12__(:)+0.5*dx_1 , x_1_12__(:)+0.5*dx_1 , x_1_12__(:)-0.5*dx_1 ] );
tmp_x_2__ = transpose( [ x_2_12__(:)-0.5*dx_2 , x_2_12__(:)-0.5*dx_2 , x_2_12__(:)+0.5*dx_2 , x_2_12__(:)+0.5*dx_2 ] );
tmp_x_0__ = transpose( (max(x_0_lim_)-0.0*dx_0)*ones(n_x_1*n_x_2,4) );
tmp_c___ = reshape(c_use__(1+nc_12__(:),:),[1,n_x_1*n_x_2,3]);
hold on;
patch(tmp_x_0__,tmp_x_1__,tmp_x_2__,tmp_c___,'LineStyle','none','EdgeColor','none');
hold off;
%%%%%%%%;
end;%if  flag_projection;

%%%%%%%%;
hold on;
if flag_boxgrid; plot_box_grid_0; end;% if flag_collapse
%%%%%%%%;
n_vval = numel(vval_);
for nvval=n_vval-1:-1:0;
vval = vval_(1+nvval);
[tmp_faces_,tmp_vertices_] = isosurface(x_0___,x_1___,x_2___,permute(f_x_u___,[2,1,3]),vval);
n_face = size(tmp_faces_,1);
n_vertex = size(tmp_vertices_,1);
if ~flag_multicolor;
nc_f = max(0,min(n_c_use-1,floor(n_c_use*(vval-min(vlim_))/max(1e-12,diff(vlim_)))));
tmp_cdata_f_ = repmat(c_use__(1+nc_f,:),[n_face,1]);
tmp_cdata_v_ = repmat(c_use__(1+nc_f,:),[n_vertex,1]);
end;%if ~flag_multicolor;
if  flag_multicolor;
nc_use_g = max(0,min(n_c_use_g-1,floor(n_c_use_g*(vval-min(vlim_g_))/max(1e-12,diff(vlim_g_)))));
face_center_0_v_ = mean(reshape(tmp_vertices_(tmp_faces_(:),1+0),[n_face,3]),2);
face_center_1_v_ = mean(reshape(tmp_vertices_(tmp_faces_(:),1+1),[n_face,3]),2);
face_center_2_v_ = mean(reshape(tmp_vertices_(tmp_faces_(:),1+2),[n_face,3]),2);
g_f_ = interp3(x_g_0___,x_g_1___,x_g_2___,g_x_u___,face_center_0_v_,face_center_1_v_,face_center_2_v_);
nc_f_ = max(0,min(n_c_use_g-1,floor(n_c_use_g*(g_f_-min(vlim_g_))/max(1e-12,diff(vlim_g_)))));
tmp_cdata_f_ = reshape(c_use_g__(1+nc_f_,:),[n_face,3]);
g_v_ = interp3(x_g_0___,x_g_1___,x_g_2___,g_x_u___,tmp_vertices_(:,1+0),tmp_vertices_(:,1+1),tmp_vertices_(:,1+2));
nc_v_ = max(0,min(n_c_use_g-1,floor(n_c_use_g*(g_v_-min(vlim_g_))/max(1e-12,diff(vlim_g_)))));
tmp_cdata_v_ = reshape(c_use_g__(1+nc_v_,:),[n_vertex,3]);
end;%if  flag_multicolor;
%%%%%%%%;
if flag_collapse;
tmp_vertices_0_ = bsxfun(@plus,[max(x_0_lim_),0,0],bsxfun(@times,[-1e-3,+1,+1],tmp_vertices_));
spatch = patch('Faces',tmp_faces_,'Vertices',tmp_vertices_0_,'FaceVertexCData',tmp_cdata_v_,'EdgeColor','none','FaceColor','interp');
alpha(spatch,v_alpha);
tmp_vertices_1_ = bsxfun(@plus,[0,max(x_1_lim_),0],bsxfun(@times,[+1,-1e-3,+1],tmp_vertices_));
spatch = patch('Faces',tmp_faces_,'Vertices',tmp_vertices_1_,'FaceVertexCData',tmp_cdata_v_,'EdgeColor','none','FaceColor','interp');
alpha(spatch,v_alpha);
tmp_vertices_2_ = bsxfun(@plus,[0,0,min(x_2_lim_)],bsxfun(@times,[+1,+1,-1e-3],tmp_vertices_));
spatch = patch('Faces',tmp_faces_,'Vertices',tmp_vertices_2_,'FaceVertexCData',tmp_cdata_v_,'EdgeColor','none','FaceColor','interp');
alpha(spatch,v_alpha);
end;%if flag_collapse;
%%%%%%%%;
hpatch = patch('Faces',tmp_faces_,'Vertices',tmp_vertices_,'FaceVertexCData',tmp_cdata_v_,'EdgeColor','none','FaceColor','interp');
alpha(hpatch,v_alpha);
isonormals(x_0___,x_1___,x_2___,permute(f_x_u___,[2,1,3]),hpatch);
xlim(x_0_lim_); ylim(x_1_lim_); zlim(x_2_lim_);
xlabel('x0');ylabel('x1');zlabel('x2');
set(gca,'XTick',[],'YTick',[],'ZTick',[]); grid on;
view([-65,20]); 
axis vis3d;
end;%for nvval=n_vval-1:-1:0;
if (n_vval<=1);camlight left; lighting gouraud; end;
hold off;
figbig;
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
