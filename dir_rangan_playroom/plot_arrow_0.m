function ...
[ ...
 parameter ...
] = ...
plot_arrow_0( ...
 parameter ...
,xy_base_ ...
,w_base ...
,r_base ...
);

str_thisfunction = 'plot_arrow_0';

if (nargin< 1);
disp(sprintf(' %% testing %s',str_thisfunction));
n_gamma_z = 16;
r_base = 0.45*2*pi*4.0/n_gamma_z;
for ngamma_z=0:n_gamma_z-1;
gamma_z = (2*pi)*ngamma_z/n_gamma_z;
plot_arrow_0([],4*[cos(gamma_z);sin(gamma_z)],gamma_z+pi/2,r_base);
end;%for ngamma_z=0:n_gamma_z-1;
xlim([-5,+5]); ylim([-5,+5]); 
set(gca,'XTick',[-5:+5]); set(gca,'YTick',[-5:+5]);
axis square; grid on;
disp('returning'); return;
end;%if (nargin< 1);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); xy_base_=[]; end; na=na+1;
if (nargin<1+na); w_base=[]; end; na=na+1;
if (nargin<1+na); r_base=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'arrow_tail_width'); parameter.arrow_tail_width = 0.125; end;
arrow_tail_width = parameter.arrow_tail_width;
if ~isfield(parameter,'arrow_tail_length'); parameter.arrow_tail_length = 1.00; end;
arrow_tail_length = parameter.arrow_tail_length;
if ~isfield(parameter,'arrow_head_width'); parameter.arrow_head_width = 0.25; end;
arrow_head_width = parameter.arrow_head_width;
if ~isfield(parameter,'arrow_head_length'); parameter.arrow_head_length = 1.00; end;
arrow_head_length = parameter.arrow_head_length;
if ~isfield(parameter,'arrow_head_start'); parameter.arrow_head_start = 0.25; end;
arrow_head_start = parameter.arrow_head_start;
if ~isfield(parameter,'arrow_head_final'); parameter.arrow_head_final = 1.00; end;
arrow_head_final = parameter.arrow_head_final;
if ~isfield(parameter,'arrow_linewidth_use'); parameter.arrow_linewidth_use = 2; end;
arrow_linewidth_use = parameter.arrow_linewidth_use;
if ~isfield(parameter,'arrow_facecolor_use'); parameter.arrow_facecolor_use = [0.85,0.85,0.85]; end;
arrow_facecolor_use = parameter.arrow_facecolor_use;
if ~isfield(parameter,'arrow_edgecolor_use'); parameter.arrow_edgecolor_use = [0.65,0.65,0.65]; end;
arrow_edgecolor_use = parameter.arrow_edgecolor_use;

if isempty(xy_base_); xy_base_ = [0.0;0.0]; end;
if isempty(w_base); w_base = 0.0; end;
if isempty(r_base); r_base = 0.5; end;
xy_base_ = xy_base_(:);

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

x0_ = [ ...
  -arrow_tail_length*r_base ...
, +(arrow_head_length-1.0)*r_base+1.0*r_base*arrow_head_start ...
, +(arrow_head_length-1.0)*r_base+1.0*r_base*arrow_head_start ...
, +(arrow_head_length-1.0)*r_base+1.0*r_base*arrow_head_final ...
 ];
y0_ = [ ...
 +r_base*arrow_tail_width ...
, +r_base*arrow_tail_width ...
, +r_base*(arrow_tail_width + arrow_head_width) ...
, 0.0 ...
 ];
x1_ = [+x0_,flip(+x0_)];
y1_ = [+y0_,flip(-y0_)];
xy1__ = [x1_ ; y1_];
R__ = [cos(w_base) , -sin(w_base) ; +sin(w_base) , cos(w_base)];
xyr__ =  R__ * xy1__;
arrow = patch(transpose(xyr__(1+0,:))+xy_base_(1+0),transpose(xyr__(1+1,:))+xy_base_(1+1),arrow_facecolor_use);
set(arrow,'LineWidth',arrow_linewidth_use);
set(arrow,'EdgeColor',arrow_edgecolor_use);

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
