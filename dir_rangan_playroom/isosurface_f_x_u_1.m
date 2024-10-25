function ...
[ ...
 parameter ...
,v_ ...
] = ...
isosurface_f_x_u_1( ...
 parameter ...
,f_x_u___ ...
);

str_thisfunction = 'isosurface_f_x_u_1';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); f_x_u___=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
tolerance_master = parameter.tolerance_master;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
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
if ~isfield(parameter,'c_use__'); parameter.c_use__ = flipud(colormap('spring')); end;
c_use__ = parameter.c_use__; n_c_use = size(c_use__,1);
if ~isfield(parameter,'flag_solid'); parameter.flag_solid = 1; end;
flag_solid = parameter.flag_solid;

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

hold on;
n_vval = numel(vval_);
for nvval=n_vval-1:-1:0;
vval = vval_(1+nvval);
hpatch = patch(isosurface(x_0___,x_1___,x_2___,permute(f_x_u___,[2,1,3]),vval)); 
isonormals(x_0___,x_1___,x_2___,permute(f_x_u___,[2,1,3]),hpatch);
nc_use = max(0,min(n_c_use-1,floor(n_c_use*(vval-min(vlim_))/max(1e-12,diff(vlim_)))));
hpatch.FaceColor = c_use__(1+nc_use,:); 
hpatch.EdgeColor = 'none'; 
if ~flag_solid; hpatch.FaceAlpha = max(0.0,min(1.0,(nvval/(n_vval-1)).^4 * 1.0)); end;
xlim(x_0_lim_); ylim(x_1_lim_); zlim(x_2_lim_);
xlabel('x0');ylabel('x1');zlabel('x2');
set(gca,'XTick',[],'YTick',[],'ZTick',[]); grid on;
view([-65,20]); 
axis vis3d;
end;%for nvval=n_vval-1:-1:0;
if (flag_solid | n_vval<=1);camlight left; lighting gouraud; end;
hold off;
figbig;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
