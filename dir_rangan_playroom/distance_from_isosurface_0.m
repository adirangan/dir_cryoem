function ...
[ ...
 parameter ...
,d_a_ ...
] = ...
distance_from_isosurface_0( ...
 parameter ...
,faces_b_f3__ ...
,vertices_b_v3__ ...
,faces_a_f3__ ...
,vertices_a_v3__ ...
);

str_thisfunction = 'distance_from_isosurface_0';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
rng(0);
%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
n_x_u_pack = 64;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack); x_u_0_lim_ = [-1,+1]*x_p_r_max;
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack); x_u_1_lim_ = [-1,+1]*x_p_r_max;
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack); x_u_2_lim_ = [-1,+1]*x_p_r_max;
[x_u_0___,x_u_1___,x_u_2___] = meshgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;  %<-- note that isosurface requires meshgrid, not ndgrid. ;
xxx_u_weight_ = (2*x_p_r_max/n_x_u_pack)^3;
%%%%;
n_source = 32;
source_a_3s__ = zeros(3,n_source);
source_b_3s__ = zeros(3,n_source);
for nsource=0:n_source-1;
source_a_3s__(:,1+nsource) = 0.5*randn(3,1);
source_b_3s__(:,1+nsource) = source_a_3s__(:,1+nsource) + 0.5*1e-3*randn(3,1);
end;%for nsource=0:n_source-1;
%%%%;
g___ = @(source_,sigma_x) 1/sqrt(2*pi).^3 ./ max(1e-12,sigma_x.^3) .* exp(- ( (x_u_0___ - source_(1+0)).^2 + (x_u_1___ - source_(1+1)).^2 + (x_u_2___ - source_(1+2)).^2 ) / max(1e-12,sigma_x.^2) );
%%%%;
sigma_x = diameter_x_c/8;
a_x_u___ = zeros(n_x_u_pack,n_x_u_pack,n_x_u_pack);
b_x_u___ = zeros(n_x_u_pack,n_x_u_pack,n_x_u_pack);
for nsource=0:n_source-1;
source_a_3_ = source_a_3s__(:,1+nsource);
a_x_u___ = a_x_u___ + g___(source_a_3_,sigma_x);
source_b_3_ = source_b_3s__(:,1+nsource);
b_x_u___ = b_x_u___ + g___(source_b_3_,sigma_x);
end;%for nsource=0:n_source-1;
%%%%;
vval_a = prctile(a_x_u___,75,'all');
[faces_a_f3__,vertices_a_v3__] = isosurface(x_u_0___,x_u_1___,x_u_2___,permute(a_x_u___,[2,1,3]),vval_a);
vval_b = prctile(b_x_u___,75,'all');
[faces_b_f3__,vertices_b_v3__] = isosurface(x_u_0___,x_u_1___,x_u_2___,permute(b_x_u___,[2,1,3]),vval_b);
%%%%%%%%;
%distance_from_isosurface_0(struct('flag_verbose',1,'flag_check',1),faces_b_f3__,vertices_b_v3__,faces_a_v3__,vertices_a_v3__);
%%%%%%%%;
figure(1);clf;figmed;
%%%%%%%%;
subplot(1,3,1);
hold on;
plot_box_grid_0;
%%%%;
v_alpha_a = 0.5;
cdata_a_ = bsxfun(@plus,zeros(size(vertices_a_v3__)),[1,0,1]);
hpatch_a = patch('Faces',faces_a_f3__,'Vertices',vertices_a_v3__,'FaceVertexCData',cdata_a_,'EdgeColor','none','FaceColor','interp');
alpha(hpatch_a,v_alpha_a);
isonormals(x_u_0___,x_u_1___,x_u_2___,permute(a_x_u___,[2,1,3]),hpatch_a);
%%%%;
v_alpha_b = 0.5;
cdata_b_ = bsxfun(@plus,zeros(size(vertices_b_v3__)),[1,1,0]);
hpatch_b = patch('Faces',faces_b_f3__,'Vertices',vertices_b_v3__,'FaceVertexCData',cdata_b_,'EdgeColor','none','FaceColor','interp');
alpha(hpatch_b,v_alpha_b);
isonormals(x_u_0___,x_u_1___,x_u_2___,permute(b_x_u___,[2,1,3]),hpatch_b);
%%%%;
xlim(x_u_0_lim_); ylim(x_u_1_lim_); zlim(x_u_2_lim_);
xlabel('x0');ylabel('x1');zlabel('x2');
set(gca,'XTick',[],'YTick',[],'ZTick',[]); grid on;
view([-65,20]); 
axis vis3d;
camlight left; lighting gouraud;
hold off;
%%%%%%%%;
subplot(1,3,2);
vlim_g_ = [0,1.5e-3];
isosurface_f_x_u_3(struct('vval_',[vval_a],'vlim_g_',vlim_g_),a_x_u___,b_x_u___);
xlabel(''); ylabel(''); zlabel('');
%%%%%%%%;
subplot(1,3,3);
vlim_g_ = [0,2.5e-3];
isosurface_f_x_u_3(struct('vval_',[vval_a],'vlim_g_',vlim_g_,'flag_version',1),a_x_u___,b_x_u___);
xlabel(''); ylabel(''); zlabel('');
%%%%%%%%;
figbig;
%%%%%%%%;
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_check'); parameter.flag_check=0; end;
flag_check=parameter.flag_check;
if ~isfield(parameter,'n_graph_laplacian'); parameter.n_graph_laplacian=0; end;
n_graph_laplacian=parameter.n_graph_laplacian;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_faces_b = size(faces_b_f3__,1);
n_vertices_b = size(vertices_b_v3__,1);
n_vertices_a = size(vertices_a_v3__,1);
d_a_ = zeros(n_vertices_a,1);

faces_b_3f__ = transpose(faces_b_f3__);
vertices_b_3v__ = transpose(vertices_b_v3__);
ij_faces_b_3f__ = bsxfun(@times,ones(3,1),1:n_faces_b);
ij_vertices_b_3f__ = faces_b_3f__;
ij_faces_b_from_vertices_b_fv__ = sparse(ij_faces_b_3f__(:),ij_vertices_b_3f__(:),1,n_faces_b,n_vertices_b);
% Note: each face has 3 vertices, so each row has 3 nonzero columns. ;
n_faces_per_vertex = full(max(sum(ij_faces_b_from_vertices_b_fv__,1)));
if (flag_verbose> 0); disp(sprintf(' %% n_faces_per_vertex: %d',n_faces_per_vertex)); end;

face_avg_b_3f__ = permute(reshape(mean(reshape(vertices_b_v3__(faces_b_3f__(:),:),[3,n_faces_b,3]),1),[n_faces_b,3]),[2,1]);
if flag_check;
nf=min(14,n_faces_b);
tmp_avg_3_ = mean(vertices_b_v3__(faces_b_3f__(:,1+nf),:),1);
face_avg_3_ = face_avg_b_3f__(:,1+nf);
fnorm_disp(flag_verbose,'tmp_avg_3_',tmp_avg_3_(:),'face_avg_3_',face_avg_3_(:));
end;%if flag_check;

%%%%;
n_k = n_faces_per_vertex; %<-- number of nearest faces to check. ;
[ij_vk__,d_vk__] = knnsearch(transpose(face_avg_b_3f__),vertices_a_v3__,'K',n_k);
vertex_33kv____ = zeros(3,3,n_k,n_vertices_a); %<-- (x0,x1,x2), vertex_b (i.e., triangle), neighbor index, vertex_a. ;
for nvertices_a=0:n_vertices_a-1;
for nk=0:n_k-1;
vertex_33kv____(:,:,1+nk,1+nvertices_a) = vertices_b_3v__(:,faces_b_3f__(:,ij_vk__(1+nvertices_a,1+nk)));
end;%for nk=0:n_k-1;
%vertex_33kv____(:,:,1+0,1+nvertices_a) = vertices_b_3v__(:,faces_b_3f__(:,ij_vk__(1+nvertices_a,1+0)));
%vertex_33kv____(:,:,1+1,1+nvertices_a) = vertices_b_3v__(:,faces_b_3f__(:,ij_vk__(1+nvertices_a,1+1)));
%vertex_33kv____(:,:,1+2,1+nvertices_a) = vertices_b_3v__(:,faces_b_3f__(:,ij_vk__(1+nvertices_a,1+2)));
%vertex_33kv____(:,:,1+3,1+nvertices_a) = vertices_b_3v__(:,faces_b_3f__(:,ij_vk__(1+nvertices_a,1+3)));
end;%for nvertices_a=0:n_vertices_a-1;
%%%%;
extern_33kv____ = zeros(3,3,n_k,n_vertices_a);
extern_33kv____ = repmat(reshape(transpose(vertices_a_v3__),[3,1,1,n_vertices_a]),[1,3,n_k,1]);
%%%%;
ver0_33kv____ = (1.0/1.0)*circshift(vertex_33kv____,+0,2);
ver1_33kv____ = (1.0/1.0)*circshift(vertex_33kv____,+1,2);
ver2_33kv____ = (1.0/1.0)*circshift(vertex_33kv____,+2,2);
distance_extern_ver0_13kv____ = sqrt(sum((extern_33kv____ - ver0_33kv____).^2,1));
distance_extern_ver1_13kv____ = sqrt(sum((extern_33kv____ - ver1_33kv____).^2,1));
distance_extern_ver2_13kv____ = sqrt(sum((extern_33kv____ - ver2_33kv____).^2,1));
%%%%;
edg01_33kv____ = (1.0/2.0)*(ver0_33kv____ + ver1_33kv____);
edg12_33kv____ = (1.0/2.0)*(ver1_33kv____ + ver2_33kv____);
edg20_33kv____ = (1.0/2.0)*(ver2_33kv____ + ver0_33kv____);
ver0_edg01_33kv____ = ver0_33kv____ - edg01_33kv____;
ver1_edg01_33kv____ = ver1_33kv____ - edg01_33kv____;
ver1_edg12_33kv____ = ver1_33kv____ - edg12_33kv____;
ver2_edg12_33kv____ = ver2_33kv____ - edg12_33kv____;
ver2_edg20_33kv____ = ver2_33kv____ - edg20_33kv____;
ver0_edg20_33kv____ = ver0_33kv____ - edg20_33kv____;
v_para_01_33kv____ = ver0_edg01_33kv____ - ver1_edg01_33kv____;
v_para_01_l2_13kv____ = sqrt(sum((v_para_01_33kv____).^2,1));
v_para_01_norm_33kv____ = bsxfun(@rdivide,v_para_01_33kv____,max(1e-12,v_para_01_l2_13kv____));
v_para_12_33kv____ = ver1_edg12_33kv____ - ver2_edg12_33kv____;
v_para_12_l2_13kv____ = sqrt(sum((v_para_12_33kv____).^2,1));
v_para_12_norm_33kv____ = bsxfun(@rdivide,v_para_12_33kv____,max(1e-12,v_para_12_l2_13kv____));
v_para_20_33kv____ = ver2_edg20_33kv____ - ver0_edg20_33kv____;
v_para_20_l2_13kv____ = sqrt(sum((v_para_20_33kv____).^2,1));
v_para_20_norm_33kv____ = bsxfun(@rdivide,v_para_20_33kv____,max(1e-12,v_para_20_l2_13kv____));
extern_edg01_33kv____ = extern_33kv____ - edg01_33kv____;
extern_edg12_33kv____ = extern_33kv____ - edg12_33kv____;
extern_edg20_33kv____ = extern_33kv____ - edg20_33kv____;
extern_edg01_dot_v_para_01_norm_13kv____ = sum(extern_edg01_33kv____.*v_para_01_norm_33kv____,1);
extern_edg12_dot_v_para_12_norm_13kv____ = sum(extern_edg12_33kv____.*v_para_12_norm_33kv____,1);
extern_edg20_dot_v_para_20_norm_13kv____ = sum(extern_edg20_33kv____.*v_para_20_norm_33kv____,1);
distance_extern_edg01_ext_13kv____ = sqrt(sum((extern_edg01_33kv____ - bsxfun(@times,extern_edg01_dot_v_para_01_norm_13kv____,v_para_01_norm_33kv____)).^2,1));
distance_extern_edg12_ext_13kv____ = sqrt(sum((extern_edg12_33kv____ - bsxfun(@times,extern_edg12_dot_v_para_12_norm_13kv____,v_para_12_norm_33kv____)).^2,1));
distance_extern_edg20_ext_13kv____ = sqrt(sum((extern_edg20_33kv____ - bsxfun(@times,extern_edg20_dot_v_para_20_norm_13kv____,v_para_20_norm_33kv____)).^2,1));
ver0_edg01_dot_v_para_01_norm_13kv____ = sum(ver0_edg01_33kv____.*v_para_01_norm_33kv____,1);
ver1_edg01_dot_v_para_01_norm_13kv____ = sum(ver1_edg01_33kv____.*v_para_01_norm_33kv____,1);
ver1_edg12_dot_v_para_12_norm_13kv____ = sum(ver1_edg12_33kv____.*v_para_12_norm_33kv____,1);
ver2_edg12_dot_v_para_12_norm_13kv____ = sum(ver2_edg12_33kv____.*v_para_12_norm_33kv____,1);
ver2_edg20_dot_v_para_20_norm_13kv____ = sum(ver2_edg20_33kv____.*v_para_20_norm_33kv____,1);
ver0_edg20_dot_v_para_20_norm_13kv____ = sum(ver0_edg20_33kv____.*v_para_20_norm_33kv____,1);
%%%%;
distance_extern_edg01_flag_13kv____ = (extern_edg01_dot_v_para_01_norm_13kv____>=min(ver0_edg01_dot_v_para_01_norm_13kv____,ver1_edg01_dot_v_para_01_norm_13kv____)) & (extern_edg01_dot_v_para_01_norm_13kv____<=max(ver0_edg01_dot_v_para_01_norm_13kv____,ver1_edg01_dot_v_para_01_norm_13kv____)) ;
distance_extern_edg01_use_13kv____ = (0+distance_extern_edg01_flag_13kv____).*distance_extern_edg01_ext_13kv____ + (1-distance_extern_edg01_flag_13kv____).*min(distance_extern_ver0_13kv____,distance_extern_ver1_13kv____);
distance_extern_edg12_flag_13kv____ = (extern_edg12_dot_v_para_12_norm_13kv____>=min(ver1_edg12_dot_v_para_12_norm_13kv____,ver2_edg12_dot_v_para_12_norm_13kv____)) & (extern_edg12_dot_v_para_12_norm_13kv____<=max(ver1_edg12_dot_v_para_12_norm_13kv____,ver2_edg12_dot_v_para_12_norm_13kv____)) ;
distance_extern_edg12_use_13kv____ = (0+distance_extern_edg12_flag_13kv____).*distance_extern_edg12_ext_13kv____ + (1-distance_extern_edg12_flag_13kv____).*min(distance_extern_ver1_13kv____,distance_extern_ver2_13kv____);
distance_extern_edg20_flag_13kv____ = (extern_edg20_dot_v_para_20_norm_13kv____>=min(ver2_edg20_dot_v_para_20_norm_13kv____,ver0_edg20_dot_v_para_20_norm_13kv____)) & (extern_edg20_dot_v_para_20_norm_13kv____<=max(ver2_edg20_dot_v_para_20_norm_13kv____,ver0_edg20_dot_v_para_20_norm_13kv____)) ;
distance_extern_edg20_use_13kv____ = (0+distance_extern_edg20_flag_13kv____).*distance_extern_edg20_ext_13kv____ + (1-distance_extern_edg20_flag_13kv____).*min(distance_extern_ver2_13kv____,distance_extern_ver0_13kv____);
%%%%;
distance_extern_edge_use_13kv____ = min(min(distance_extern_edg01_use_13kv____,distance_extern_edg12_use_13kv____),distance_extern_edg20_use_13kv____);
%%%%;
triavg_33kv____ = (1.0/3.0)*(ver0_33kv____ + ver1_33kv____ + ver2_33kv____);
%%;
ver0_triavg_33kv____ = ver0_33kv____ - triavg_33kv____;
ver1_triavg_33kv____ = ver1_33kv____ - triavg_33kv____;
ver2_triavg_33kv____ = ver2_33kv____ - triavg_33kv____;
extern_triavg_33kv____ = extern_33kv____ - triavg_33kv____;
%%;
v_para_01_33kv____ = ver0_triavg_33kv____ - ver1_triavg_33kv____;
v_para_01_l2_13kv____ = sqrt(sum((v_para_01_33kv____).^2,1));
v_para_01_norm_33kv____ = bsxfun(@rdivide,v_para_01_33kv____,max(1e-12,v_para_01_l2_13kv____));
v_para_12_33kv____ = ver1_triavg_33kv____ - ver2_triavg_33kv____;
v_para_12_l2_13kv____ = sqrt(sum((v_para_12_33kv____).^2,1));
v_para_12_norm_33kv____ = bsxfun(@rdivide,v_para_12_33kv____,max(1e-12,v_para_12_l2_13kv____));
v_para_20_33kv____ = ver2_triavg_33kv____ - ver0_triavg_33kv____;
v_para_20_l2_13kv____ = sqrt(sum((v_para_20_33kv____).^2,1));
v_para_20_norm_33kv____ = bsxfun(@rdivide,v_para_20_33kv____,max(1e-12,v_para_20_l2_13kv____));
%%;
ver0_perp_01_33kv____ = ver0_triavg_33kv____ - bsxfun(@times,sum(ver0_triavg_33kv____.*v_para_01_33kv____,1),v_para_01_33kv____);
ver1_perp_01_33kv____ = ver1_triavg_33kv____ - bsxfun(@times,sum(ver1_triavg_33kv____.*v_para_01_33kv____,1),v_para_01_33kv____);
v_perp_01_33kv____ = (1.0/2.0)*(ver0_perp_01_33kv____ + ver1_perp_01_33kv____);
%%;
ver1_perp_12_33kv____ = ver1_triavg_33kv____ - bsxfun(@times,sum(ver1_triavg_33kv____.*v_para_12_33kv____,1),v_para_12_33kv____);
ver2_perp_12_33kv____ = ver2_triavg_33kv____ - bsxfun(@times,sum(ver2_triavg_33kv____.*v_para_12_33kv____,1),v_para_12_33kv____);
v_perp_12_33kv____ = (1.0/2.0)*(ver1_perp_12_33kv____ + ver2_perp_12_33kv____);
%%;
ver2_perp_20_33kv____ = ver2_triavg_33kv____ - bsxfun(@times,sum(ver2_triavg_33kv____.*v_para_20_33kv____,1),v_para_20_33kv____);
ver0_perp_20_33kv____ = ver0_triavg_33kv____ - bsxfun(@times,sum(ver0_triavg_33kv____.*v_para_20_33kv____,1),v_para_20_33kv____);
v_perp_20_33kv____ = (1.0/2.0)*(ver2_perp_20_33kv____ + ver0_perp_20_33kv____);
%%;
dperp_33kv____ = ...
 cross(ver0_triavg_33kv____-ver1_triavg_33kv____,ver1_triavg_33kv____-ver2_triavg_33kv____,1) ...
+cross(ver1_triavg_33kv____-ver2_triavg_33kv____,ver2_triavg_33kv____-ver0_triavg_33kv____,1) ...
+cross(ver2_triavg_33kv____-ver0_triavg_33kv____,ver0_triavg_33kv____-ver1_triavg_33kv____,1) ...
;
dperp_l2_13kv____ = sqrt(sum((dperp_33kv____).^2,1));
dperp_norm_33kv____ = bsxfun(@rdivide,dperp_33kv____,max(1e-12,dperp_l2_13kv____));
%%;
distance_extern_triavg_ext_13kv____ = abs(sum(extern_triavg_33kv____.*dperp_norm_33kv____,1));
distance_extern_triavg_flag_13kv____ = ...
  sum(extern_triavg_33kv____.*v_perp_01_33kv____,1)<=sum(v_perp_01_33kv____.^2,1) ...
& sum(extern_triavg_33kv____.*v_perp_12_33kv____,1)<=sum(v_perp_12_33kv____.^2,1) ...
& sum(extern_triavg_33kv____.*v_perp_20_33kv____,1)<=sum(v_perp_20_33kv____.^2,1) ...
;
distance_extern_triavg_use_13kv____ = (0+distance_extern_triavg_flag_13kv____).*distance_extern_triavg_ext_13kv____ + (1-distance_extern_triavg_flag_13kv____).*distance_extern_edge_use_13kv____;
%%;
%%%%%%%%;
d_a_ = reshape(min(reshape(distance_extern_triavg_use_13kv____,[3*n_k,n_vertices_a]),[],1),[n_vertices_a,1]);

if n_graph_laplacian> 0;
n_faces_a = size(faces_a_f3__,1);
faces_a_3f__ = transpose(faces_a_f3__);
vertices_a_3v__ = transpose(vertices_a_v3__);
ij_faces_a_3f__ = bsxfun(@times,ones(3,1),1:n_faces_a);
ij_vertices_a_3f__ = faces_a_3f__;
ij_faces_a_from_vertices_a_fv__ = sparse(ij_faces_a_3f__(:),ij_vertices_a_3f__(:),1,n_faces_a,n_vertices_a);
ij_vertices_a_from_vertices_a_vv__ = transpose(ij_faces_a_from_vertices_a_fv__)*ij_faces_a_from_vertices_a_fv__;
ij_vertices_a_nrm_v_ = sum(ij_vertices_a_from_vertices_a_vv__,2);
for ngraph_laplacian=0:n_graph_laplacian-1;
d_a_ = bsxfun(@rdivide,ij_vertices_a_from_vertices_a_vv__*d_a_,max(1,ij_vertices_a_nrm_v_));
end;%for ngraph_laplacian=0:n_graph_laplacian-1;
end;%if n_graph_laplacian> 0;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
