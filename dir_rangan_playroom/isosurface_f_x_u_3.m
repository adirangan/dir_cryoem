function ...
[ ...
 parameter ...
,d_v_ ...
] = ...
isosurface_f_x_u_3( ...
 parameter ...
,f_x_u___ ...
,g_x_u___ ...
);

str_thisfunction = 'isosurface_f_x_u_3';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); f_x_u___=[]; end; na=na+1;
if (nargin<1+na); g_x_u___=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
tolerance_master = parameter.tolerance_master;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'flag_plot'); parameter.flag_plot = 1; end;
flag_plot = parameter.flag_plot;
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
if ~isfield(parameter,'c_use__'); parameter.c_use__ = colormap('spring'); parameter.c_use__ = flipud(parameter.c_use__); end;
c_use__ = parameter.c_use__; n_c_use = size(c_use__,1);
if ~isfield(parameter,'vlim_g_'); parameter.vlim_g_ = []; end;
vlim_g_ = parameter.vlim_g_;
if ~isfield(parameter,'c_use_g__'); parameter.c_use_g__ = colormap_81s; parameter.c_use_g__ = flipud(parameter.c_use_g__); end;
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
if isempty(vlim_g_); vlim_g_ = [0,0.5]; end;

if  flag_plot &  flag_projection;
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
end;%if  flag_plot &  flag_projection;

d_v_ = [];
%%%%%%%%;
if  flag_plot;
hold on;
if flag_boxgrid; plot_box_grid_0; end;% if flag_collapse
end;%if  flag_plot;
%%%%%%%%;
n_vval = numel(vval_);
for nvval=n_vval-1:-1:0;
vval = vval_(1+nvval);
[tmp_faces_f_,tmp_vertices_f_] = isosurface(x_0___,x_1___,x_2___,permute(f_x_u___,[2,1,3]),vval);
[tmp_faces_g_,tmp_vertices_g_] = isosurface(x_g_0___,x_g_1___,x_g_2___,permute(g_x_u___,[2,1,3]),vval);
n_face_f = size(tmp_faces_f_,1);
n_vertex_f = size(tmp_vertices_f_,1);
%%%%;

max(tmp_faces_f_(:)),;
n_vertex_f,;
tmp_faces_f_(2,:),
tmp_vertices_f_(2,:),

[ij_v8__,d_v8__] = knnsearch(tmp_vertices_g_,tmp_vertices_f_,'K',8);
a0_v_ = tmp_vertices_g_(ij_v8__(:,1+0),:);
a1_v_ = tmp_vertices_g_(ij_v8__(:,1+1),:);
a2_v_ = tmp_vertices_g_(ij_v8__(:,1+2),:);
a3_v_ = tmp_vertices_g_(ij_v8__(:,1+3),:);
a4_v_ = tmp_vertices_g_(ij_v8__(:,1+4),:);
a5_v_ = tmp_vertices_g_(ij_v8__(:,1+5),:);
a6_v_ = tmp_vertices_g_(ij_v8__(:,1+6),:);
a7_v_ = tmp_vertices_g_(ij_v8__(:,1+7),:);
abar_v_ = (a0_v_ + a1_v_ + a2_v_)/3.0;
d01_v_ = a0_v_ - a1_v_;
d12_v_ = a1_v_ - a2_v_;
d20_v_ = a2_v_ - a0_v_;
dperp012_v_ = cross(d01_v_,d12_v_,2);
dperp012_v_l2_ = sqrt(sum(abs(dperp012_v_).^2,2));
dperp120_v_ = cross(d12_v_,d20_v_,2);
dperp120_v_l2_ = sqrt(sum(abs(dperp120_v_).^2,2));
dperp201_v_ = cross(d20_v_,d01_v_,2);
dperp201_v_l2_ = sqrt(sum(abs(dperp201_v_).^2,2));
%%%%;
dperp_v_ = dperp012_v_ ;
dperp_v_l2_ = dperp012_v_l2_ ;
tmp_index_ = efind(dperp120_v_l2_>dperp_v_l2_);
dperp_v_(1+tmp_index_,:) = dperp120_v_(1+tmp_index_,:) ;
dperp_v_l2_(1+tmp_index_) = dperp120_v_l2_(1+tmp_index_,:) ;
tmp_index_ = efind(dperp201_v_l2_>dperp_v_l2_);
dperp_v_(1+tmp_index_,:) = dperp201_v_(1+tmp_index_,:) ;
dperp_v_l2_(1+tmp_index_) = dperp201_v_l2_(1+tmp_index_,:) ;
%%%%;
dperp_v_ = bsxfun(@rdivide,dperp_v_,max(1e-12,dperp_v_l2_));
d_v_ = abs(sum((tmp_vertices_f_ - abar_v_).*dperp_v_,2));
subplot(1,2,1);
plot(dperp_v_l2_,'o');
subplot(1,2,2);
plot(d_v8__(:,1),d_v_,'o');
return;

%d_v_ = d_v8__(:,1);
nc_use_v_ = max(0,min(n_c_use_g-1,floor(n_c_use_g*(d_v_-min(vlim_g_))/diff(vlim_g_))));
tmp_cdata_v_ = reshape(c_use_g__(1+nc_use_v_,:),[n_vertex_f,3]);
%%%%%%%%;
if  flag_plot &  flag_collapse;
tmp_vertices_f_0_ = bsxfun(@plus,[max(x_0_lim_),0,0],bsxfun(@times,[+1e-3,-1,+1],tmp_vertices_f_));
spatch = patch('Faces',tmp_faces_f_,'Vertices',tmp_vertices_f_0_,'FaceVertexCData',tmp_cdata_v_,'EdgeColor','none','FaceColor','interp');
alpha(spatch,v_alpha);
tmp_vertices_f_1_ = bsxfun(@plus,[0,max(x_1_lim_),0],bsxfun(@times,[-1,+1e-3,+1],tmp_vertices_f_));
spatch = patch('Faces',tmp_faces_f_,'Vertices',tmp_vertices_f_1_,'FaceVertexCData',tmp_cdata_v_,'EdgeColor','none','FaceColor','interp');
alpha(spatch,v_alpha);
tmp_vertices_f_2_ = bsxfun(@plus,[0,0,min(x_0_lim_)],bsxfun(@times,[+1,+1,-1e-3],tmp_vertices_f_));
spatch = patch('Faces',tmp_faces_f_,'Vertices',tmp_vertices_f_2_,'FaceVertexCData',tmp_cdata_v_,'EdgeColor','none','FaceColor','interp');
alpha(spatch,v_alpha);
end;%if  flag_plot &  flag_collapse;
%%%%%%%%;
if  flag_plot;
hpatch = patch('Faces',tmp_faces_f_,'Vertices',tmp_vertices_f_,'FaceVertexCData',tmp_cdata_v_,'EdgeColor','none','FaceColor','interp');
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
end;%if  flag_plot;
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
