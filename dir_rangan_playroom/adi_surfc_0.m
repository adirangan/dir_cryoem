function ...
adi_surfc_0( ...
 parameter ...
, x_0_ ...
, x_1_ ...
, A__ ...
);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); x_0_=[]; end; na=na+1;
if (nargin<1+na); x_1_=[]; end; na=na+1;
if (nargin<1+na); A__=[]; end; na=na+1;

[n_row,n_col] = size(A__);
if isempty(x_0_); x_0_ = 1+transpose(0:n_row-1); end;
if isempty(x_1_); x_1_ = 1+transpose(0:n_col-1); end;

if isempty(parameter); parameter = struct('type','parameter'); end;

if ~isfield(parameter,'flag_0base'); parameter.flag_0base = 0; end;
if ~isfield(parameter,'colormap__'); parameter.colormap__ = colormap_80s; end;
if ~isfield(parameter,'pstep'); parameter.pstep = 2.5; end;
if ~isfield(parameter,'linewidth_use'); parameter.linewidth_use = 2; end;
if ~isfield(parameter,'linewidth_bot'); parameter.linewidth_bot = 3; end;
if ~isfield(parameter,'angle_az'); parameter.angle_az = -40; end;
if ~isfield(parameter,'angle_el'); parameter.angle_el = +30; end;

flag_0base = parameter.flag_0base;
colormap__ = parameter.colormap__; n_c = size(colormap__,1);
pstep = parameter.pstep;
linewidth_use = parameter.linewidth_use;
linewidth_bot = parameter.linewidth_bot;
angle_az = parameter.angle_az;
angle_el = parameter.angle_el;

n_row_end = n_row/1; n_row_mid = floor(n_row/2);
n_col_end = n_col/1; n_col_mid = floor(n_col/2);

colormap(colormap__);
hold on;
index_row_half0_ = 0:n_row_mid-1+1; index_row_half1_ = n_row_mid-1:n_row-1;
index_col_half0_ = 0:n_col_mid-1+1; index_col_half1_ = n_col_mid-1:n_col-1;
%s00=surf(x_1_(1+index_col_half0_),x_0_(1+index_row_half0_),A__(1+index_row_half0_,1+index_col_half0_)); set(s00,'LineStyle','none');
s01=surf(x_1_(1+index_col_half1_),x_0_(1+index_row_half0_),A__(1+index_row_half0_,1+index_col_half1_)); set(s01,'LineStyle','none');
s10=surf(x_1_(1+index_col_half0_),x_0_(1+index_row_half1_),A__(1+index_row_half1_,1+index_col_half0_)); set(s10,'LineStyle','none');
s11=surf(x_1_(1+index_col_half1_),x_0_(1+index_row_half1_),A__(1+index_row_half1_,1+index_col_half1_)); set(s11,'LineStyle','none');
zlim_A_ = [ min(A__,[],'all') , max(A__,[],'all') ] ; if flag_0base; zlim_A_ = [ 0 , max(A__,[],'all') ] ; end;
zmin_plot_A  = min(zlim_A_) - 0.50*diff(zlim_A_);
%%%%;
plim_ = prctile(A__,[0:pstep:100],'all'); if flag_0base; plim_ = [0 ; plim_(:)]; end;
[n_contour_A,contour_A_cell___,contour_A_level_] = contour_cell_from_matrix_0(contourc(x_1_,x_0_,A__,plim_));
contour_A_lim_ = [min(contour_A_level_) , max(contour_A_level_)];
for ncontour_A=0:n_contour_A-1;
nc = max(0,min(n_c-1,floor(n_c* ( contour_A_level_(1+ncontour_A) - min(contour_A_lim_) )/diff(contour_A_lim_))));
contour_A_x_ = contour_A_cell___{1+ncontour_A}(:,1+0);
contour_A_y_ = contour_A_cell___{1+ncontour_A}(:,1+1);
contour_A_z_ = zmin_plot_A*ones(size(contour_A_x_)) + 1e-6;
plot3(contour_A_x_,contour_A_y_,contour_A_z_,'-','LineWidth',linewidth_use,'Color',colormap__(1+nc,:));
end;%for ncontour_A=0:n_contour_A-1;
%%%%;
%line(x_1_(1+0*n_col_end)*[1,1]+1e-2,x_0_(0+1*n_row_mid)*[1,1]+0e-12,[zmin_plot_A,A__(0+1*n_row_mid,1+0*n_col_end)],'LineStyle',':','Color','k','LineWidth',linewidth_use);
%line(x_1_(0+1*n_col_end)*[1,1]-1e-2,x_0_(0+1*n_row_mid)*[1,1]+0e-12,[zmin_plot_A,A__(0+1*n_row_mid,0+1*n_col_end)],'LineStyle',':','Color','k','LineWidth',linewidth_use);
%line([x_1_(1+0*n_col_end),x_1_(0+1*n_col_end)],x_0_(0+1*n_row_mid)*[1,1]+0e-12,[zmin_plot_A,zmin_plot_A],'LineStyle',':','Color','k','LineWidth',linewidth_use);
line(x_1_(0+1*n_col_mid)*[1,1]+1e-12,x_0_(0+1*n_row_mid)*[1,1]+0e-12,[zmin_plot_A,A__(0+1*n_row_mid,0+1*n_col_mid)],'LineStyle',':','Color','k','LineWidth',linewidth_use);
line(x_1_(1+0*n_col_end)*[1,1]+1e-12,x_0_(0+1*n_row_mid)*[1,1]+0e-12,[zmin_plot_A,A__(0+1*n_row_mid,1+0*n_col_end)],'LineStyle',':','Color','k','LineWidth',linewidth_use);
line(x_1_(0+1*n_col_mid)*[1,1]-1e-12,x_0_(1+0*n_row_end)*[1,1]+0e-12,[zmin_plot_A,A__(1+0*n_row_end,0+1*n_col_mid)],'LineStyle',':','Color','k','LineWidth',linewidth_use);
line([x_1_(1+0*n_col_end),x_1_(0+1*n_col_mid)],x_0_(0+1*n_row_mid)*[1,1]+0e-12,[zmin_plot_A,zmin_plot_A],'LineStyle',':','Color','k','LineWidth',linewidth_use);
line(x_1_(0+1*n_col_mid)*[1,1]+0e-12,[x_0_(1+0*n_row_end),x_0_(0+1*n_row_mid)],[zmin_plot_A,zmin_plot_A],'LineStyle',':','Color','k','LineWidth',linewidth_use);
if flag_0base;
line(x_1_(1+0*n_col_end)*[1,1]-1e-12,[x_0_(0+1*n_row_mid),x_0_(0+1*n_row_end)],0*[1,1],'LineStyle','-','Color',colormap__(1,:),'LineWidth',linewidth_bot);
line([x_1_(0+1*n_col_mid),x_1_(0+1*n_col_end)],x_0_(1+0*n_row_end)*[1,1]-1e-12,0*[1,1],'LineStyle','-','Color',colormap__(1,:),'LineWidth',linewidth_bot);
line(x_1_(0+1*n_col_end)*[1,1]-1e-12,[x_0_(1+0*n_row_end),x_0_(0+1*n_row_end)],0*[1,1],'LineStyle','-','Color',colormap__(1,:),'LineWidth',linewidth_bot);
line([x_1_(1+0*n_col_end),x_1_(0+1*n_col_end)],x_0_(0+1*n_row_end)*[1,1]-1e-12,0*[1,1],'LineStyle','-','Color',colormap__(1,:),'LineWidth',linewidth_bot);
end;%if flag_0base;
hold off;
%xlim([min(x_0_),max(x_0_)]); %xlabel('x');
%ylim([min(x_1_),max(x_1_)]); %ylabel('y');
zlim([zmin_plot_A-1e-12,max(zlim_A_)]); %zlabel('z');
set(gca,'XTick',[]); set(gca,'YTick',[]); set(gca,'ZTick',[]);
view(angle_az,angle_el);
set(gca,'FontSize',24);
axis vis3d;
%%%%%%%%;
hold off;

