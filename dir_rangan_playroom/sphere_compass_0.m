function ...
[ ...
 parameter ...
] = ...
sphere_compass_0( ...
 parameter ...
,polar_a ...
,azimu_b ...
,delta_polar_a ...
,delta_azimu_b ...
,delta_gamma_z ...
);

str_thisfunction = 'sphere_compass_0';

%%%%%%%%;
if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
parameter = struct('type','parameter');
parameter.flag_2d_vs_3d = 0;
%%%%;
figure(1);clf;figmed;
%%%%;
subplot(1,2,1);
plot_sphere_grid_0;
hold on;
sphere_compass_0(parameter);
hold off;
xlabel('x0');
ylabel('x1');
zlabel('x2');
axis equal; axis vis3d;
%%%%;
parameter.flag_2d_vs_3d = 1;
subplot(1,2,2);
hold on;
sphere_compass_0(parameter);
hold off;
xlim([0,2*pi]); ylim([0,1*pi]);
xlabel('azimu_b','Interpreter','none');
ylabel('polar_a','Interpreter','none');
axisnotick;
%%%%;
disp('returning'); return;
end;%if nargin<1;
%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); polar_a=[]; end; na=na+1;
if (nargin<1+na); azimu_b=[]; end; na=na+1;
if (nargin<1+na); delta_polar_a=[]; end; na=na+1;
if (nargin<1+na); delta_azimu_b=[]; end; na=na+1;
if (nargin<1+na); delta_gamma_z=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'k_p_r_max'); parameter.k_p_r_max=1.0; end;
k_p_r_max=parameter.k_p_r_max;
if ~isfield(parameter,'c_use__'); parameter.c_use__=colormap_beach(); end;
c_use__=parameter.c_use__; n_c_use = size(c_use__,1);
if ~isfield(parameter,'flag_2d_vs_3d'); parameter.flag_2d_vs_3d=0; end;
flag_2d_vs_3d=parameter.flag_2d_vs_3d;
if ~isfield(parameter,'compass_r_base'); parameter.compass_r_base=1.0/(2*pi)/2; end;
compass_r_base=parameter.compass_r_base;
if ~isfield(parameter,'compass_n_side_base'); parameter.compass_n_side_base=64; end;
compass_n_side_base=parameter.compass_n_side_base;
if ~isfield(parameter,'compass_base_linecolor'); parameter.compass_base_linecolor=[0.6,0.6,0.6]; end;
compass_base_linecolor=parameter.compass_base_linecolor;
if ~isfield(parameter,'compass_base_linewidth'); parameter.compass_base_linewidth=0.5; end;
compass_base_linewidth=parameter.compass_base_linewidth;
if ~isfield(parameter,'compass_pointer_linecolor'); parameter.compass_pointer_linecolor=[0.2,0.2,0.2]; end;
compass_pointer_linecolor=parameter.compass_pointer_linecolor;
if ~isfield(parameter,'compass_pointer_linewidth'); parameter.compass_pointer_linewidth=0.5; end;
compass_pointer_linewidth=parameter.compass_pointer_linewidth;
if ~isfield(parameter,'compass_pointer_patchcolor'); parameter.compass_pointer_patchcolor=[0.8,0.2,1.0]; end;
compass_pointer_patchcolor=parameter.compass_pointer_patchcolor;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(polar_a); polar_a = +1.0*pi/3.0; end;
if isempty(azimu_b); azimu_b = -2.0*pi/7.0; end;
if isempty(delta_polar_a); delta_polar_a = -0.15; end;
if isempty(delta_azimu_b); delta_azimu_b = +0.20; end;
if isempty(delta_gamma_z); delta_gamma_z = -0.50; end;

pointer_template_s2__ = [ ...
  -0.85 , +0.00 ...
; -0.85 , +0.25 ...
; +0.05 , +0.25 ...
; -0.25 , +0.65 ...
; +1.00 , +0.00 ...
; -0.25 , -0.65 ...
; +0.05 , -0.25 ...
; -0.85 , -0.25 ...
];
compass_n_side_pointer = size(pointer_template_s2__,1);
pointer_template_2s__ = transpose(pointer_template_s2__);

compass_c_base_ = c_use__(round(n_c_use/2),:);
ca = cos(polar_a); sa = sin(polar_a);
cb = cos(azimu_b); sb = sin(azimu_b);
Rz__ = [ +cb , -sb , 0 ; +sb , +cb , 0 ; 0 , 0 , 1];
Ry__ = [ +ca , 0 , +sa ; 0 , 1 , 0 ; -sa , 0 , +ca ];
R__ = Rz__*Ry__ ;
delta_x0 = +cb*ca*delta_polar_a - sb*sa*delta_azimu_b; %<-- (dxda + dxdb)*(cbsa) ;
delta_x1 = +sb*ca*delta_polar_a + cb*sa*delta_azimu_b; %<-- (dxda + dxdb)*(sbsa) ;
delta_x2 = -   sa*delta_polar_a +     0*delta_azimu_b; %<-- (dxda + dxdb)*(ca) ;
delta_xx_ = inv(R__)*[delta_x0;delta_x1;delta_x2]; delta_xx_omega = atan2(delta_xx_(1+1),delta_xx_(1+0));
Rx__ = [ cos(delta_xx_omega) , -sin(delta_xx_omega) ; +sin(delta_xx_omega) , cos(delta_xx_omega) ];
compass_base_3s__ = zeros(3,compass_n_side_base);
for compass_nside_base=0:compass_n_side_base-1+1;
compass_w = (2*pi*compass_nside_base)/max(1,compass_n_side_base);
compass_b_ = [ compass_r_base*cos(compass_w) ; compass_r_base*sin(compass_w) ; sqrt(1-compass_r_base.^2) * 1.0 ] ;
compass_base_3s__(:,1+compass_nside_base) = compass_b_;
end;%for compass_nside_base=0:compass_n_side_base-1+1;
%%%%;
compass_pointer_2s__ = compass_r_base * Rx__ * pointer_template_2s__;
compass_pointer_3s__ = [ compass_pointer_2s__ ; sqrt(1-sum(compass_pointer_2s__.^2,1)) ];
%%%%;

%%%%;
compass_base_3s__ = pagemtimes(R__,compass_base_3s__); %<-- list of vertices in compass base. ;
compass_pointer_3s__ = pagemtimes(R__,compass_pointer_3s__); %<-- list of vertices in compass pointer. ;
%%%%;


if flag_2d_vs_3d==0;
%%%%;
l = line( ...
 transpose(compass_base_3s__(1+0,:)) ...
,transpose(compass_base_3s__(1+1,:)) ...
,transpose(compass_base_3s__(1+2,:)) ...
);
set(l,'LineStyle','-','Color',compass_base_linecolor,'LineWidth',compass_base_linewidth);
%%%%;
p = patch( ...
 transpose(compass_pointer_3s__(1+0,:)) ...
,transpose(compass_pointer_3s__(1+1,:)) ...
,transpose(compass_pointer_3s__(1+2,:)) ...
,compass_pointer_patchcolor ...
);
set(p,'LineStyle','-','EdgeColor',compass_pointer_linecolor,'LineWidth',compass_pointer_linewidth);
%%%%;
end;%if flag_2d_vs_3d==0;

if flag_2d_vs_3d==1;
%%%%;
compass_base_xs_ = reshape(compass_base_3s__(1+0,:),[1+compass_n_side_base,1]);
compass_base_ys_ = reshape(compass_base_3s__(1+1,:),[1+compass_n_side_base,1]);
compass_base_zs_ = reshape(compass_base_3s__(1+2,:),[1+compass_n_side_base,1]);
compass_base_rs_ = reshape(sum(compass_base_3s__.^2,1),[1+compass_n_side_base,1]);
compass_base_ps_ = reshape(sum(compass_base_3s__(1:2,:).^2,1),[1+compass_n_side_base,1]);
compass_base_bs_ = atan2(compass_base_ys_,compass_base_xs_);
compass_base_as_ = atan2(compass_base_ps_,compass_base_zs_);
[compass_base_as_,compass_base_bs_] = periodize_polar_a_azimu_b_0(compass_base_as_,compass_base_bs_);
%%%%;
compass_pointer_xs_ = reshape(compass_pointer_3s__(1+0,:),[compass_n_side_pointer,1]);
compass_pointer_ys_ = reshape(compass_pointer_3s__(1+1,:),[compass_n_side_pointer,1]);
compass_pointer_zs_ = reshape(compass_pointer_3s__(1+2,:),[compass_n_side_pointer,1]);
compass_pointer_rs_ = reshape(sum(compass_pointer_3s__.^2,1),[compass_n_side_pointer,1]);
compass_pointer_ps_ = reshape(sum(compass_pointer_3s__(1:2,:).^2,1),[compass_n_side_pointer,1]);
compass_pointer_bs_ = atan2(compass_pointer_ys_,compass_pointer_xs_);
compass_pointer_as_ = atan2(compass_pointer_ps_,compass_pointer_zs_);
[compass_pointer_as_,compass_pointer_bs_] = periodize_polar_a_azimu_b_0(compass_pointer_as_,compass_pointer_bs_);
%%%%;
flag_border = 0;
flag_border = flag_border | (max(compass_base_bs_)-min(compass_base_bs_))>=pi/2 ;
flag_border = flag_border | (max(compass_base_as_)-min(compass_base_as_))>=pi/2 ;
flag_border = flag_border | (max(compass_pointer_bs_)-min(compass_pointer_bs_))>=pi/2 ;
flag_border = flag_border | (max(compass_pointer_as_)-min(compass_pointer_as_))>=pi/2 ;
%%%%;
if ~flag_border;
l = line( ...
 compass_base_bs_ ...
,compass_base_as_ ...
);
set(l,'LineStyle','-','Color',compass_base_linecolor,'LineWidth',compass_base_linewidth);
p = patch( ...
 compass_pointer_bs_ ...
,compass_pointer_as_ ...
,compass_pointer_patchcolor ...
);
set(p,'LineStyle','-','EdgeColor',compass_pointer_linecolor,'LineWidth',compass_pointer_linewidth);
end;%if ~flag_border;
%%%%;
end;%if flag_2d_vs_3d==1;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;


