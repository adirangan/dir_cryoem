function plot_box_grid_0(parameter);
na=0;
if nargin<1+na; parameter=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'x_0_lim_'); parameter.x_0_lim_ = [-1,+1]; end;
x_0_lim_ = parameter.x_0_lim_;
if ~isfield(parameter,'x_1_lim_'); parameter.x_1_lim_ = [-1,+1]; end;
x_1_lim_ = parameter.x_1_lim_;
if ~isfield(parameter,'x_2_lim_'); parameter.x_2_lim_ = [-1,+1]; end;
x_2_lim_ = parameter.x_2_lim_;
if ~isfield(parameter,'n_line'); parameter.n_line = 1+8; end; %<-- parameter_bookmark. ;
n_line = parameter.n_line;
if ~isfield(parameter,'linecolor_use'); parameter.linecolor_use = 0.65*[1,1,1]; end; %<-- parameter_bookmark. ;
linecolor_use = parameter.linecolor_use;
if ~isfield(parameter,'linewidth_use'); parameter.linewidth_use = 0.5; end; %<-- parameter_bookmark. ;
linewidth_use = parameter.linewidth_use;
if ~isfield(parameter,'flag_solid'); parameter.flag_solid = 1; end; %<-- parameter_bookmark. ;
flag_solid = parameter.flag_solid;
if ~isfield(parameter,'patchcolor_use'); parameter.patchcolor_use = 0.95*[1,1,1]; end; %<-- parameter_bookmark. ;
patchcolor_use = parameter.patchcolor_use;

hold on;
%%%%%%%%;
if flag_solid;
x_0_ = [x_0_lim_(1);x_0_lim_(2);x_0_lim_(2);x_0_lim_(1);x_0_lim_(1)];
x_1_ = [x_1_lim_(1);x_1_lim_(1);x_1_lim_(2);x_1_lim_(2);x_1_lim_(1)];
x_2_ = [x_2_lim_(1);x_2_lim_(1);x_2_lim_(1);x_2_lim_(1);x_2_lim_(1)];
patch(x_0_,x_1_,x_2_,patchcolor_use);
x_0_ = [x_0_lim_(2);x_0_lim_(2);x_0_lim_(2);x_0_lim_(2);x_0_lim_(2)];
x_1_ = [x_1_lim_(1);x_1_lim_(2);x_1_lim_(2);x_1_lim_(1);x_1_lim_(1)];
x_2_ = [x_2_lim_(1);x_2_lim_(1);x_2_lim_(2);x_2_lim_(2);x_2_lim_(1)];
patch(x_0_,x_1_,x_2_,patchcolor_use);
x_0_ = [x_0_lim_(1);x_0_lim_(2);x_0_lim_(2);x_0_lim_(1);x_0_lim_(1)];
x_1_ = [x_1_lim_(2);x_1_lim_(2);x_1_lim_(2);x_1_lim_(2);x_1_lim_(2)];
x_2_ = [x_2_lim_(1);x_2_lim_(1);x_2_lim_(2);x_2_lim_(2);x_2_lim_(1)];
patch(x_0_,x_1_,x_2_,patchcolor_use);
end;%if flag_solid;
%%%%%%%%;
x_0_ = repmat(x_0_lim_(:),[1,n_line]);
x_1_ = repmat(linspace(min(x_1_lim_),max(x_1_lim_),n_line),[2,1]);
x_2_ = min(x_2_lim_)*ones(2,n_line)+1e-3;
line(x_0_,x_1_,x_2_,'Color',linecolor_use,'LineWidth',linewidth_use);
x_0_ = repmat(linspace(min(x_0_lim_),max(x_0_lim_),n_line),[2,1]);
x_1_ = repmat(x_1_lim_(:),[1,n_line]);
x_2_ = min(x_2_lim_)*ones(2,n_line)+1e-3;
line(x_0_,x_1_,x_2_,'Color',linecolor_use,'LineWidth',linewidth_use);
%%%%%%%%;
x_2_ = repmat(x_2_lim_(:),[1,n_line]);
x_1_ = repmat(linspace(min(x_1_lim_),max(x_1_lim_),n_line),[2,1]);
x_0_ = max(x_0_lim_)*ones(2,n_line)-1e-3;
line(x_0_,x_1_,x_2_,'Color',linecolor_use,'LineWidth',linewidth_use);
x_2_ = repmat(linspace(min(x_2_lim_),max(x_2_lim_),n_line),[2,1]);
x_1_ = repmat(x_1_lim_(:),[1,n_line]);
x_0_ = max(x_0_lim_)*ones(2,n_line)-1e-3;
line(x_0_,x_1_,x_2_,'Color',linecolor_use,'LineWidth',linewidth_use);
%%%%%%%%;
x_2_ = repmat(x_2_lim_(:),[1,n_line]);
x_0_ = repmat(linspace(min(x_0_lim_),max(x_0_lim_),n_line),[2,1]);
x_1_ = max(x_1_lim_)*ones(2,n_line)-1e-3;
line(x_0_,x_1_,x_2_,'Color',linecolor_use,'LineWidth',linewidth_use);
x_2_ = repmat(linspace(min(x_2_lim_),max(x_2_lim_),n_line),[2,1]);
x_0_ = repmat(x_0_lim_(:),[1,n_line]);
x_1_ = max(x_1_lim_)*ones(2,n_line)-1e-3;
line(x_0_,x_1_,x_2_,'Color',linecolor_use,'LineWidth',linewidth_use);
%%%%%%%%;
hold off;

%{
xlim(x_0_lim_); ylim(x_1_lim_); zlim(x_2_lim_);
xlabel('x0');ylabel('x1');zlabel('x2');
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
view([-65,20]); 
axis vis3d;
%}



