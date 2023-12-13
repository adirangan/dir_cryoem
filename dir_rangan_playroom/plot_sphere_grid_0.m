function plot_sphere_grid_0(parameter);
na=0;
if nargin<1+na; parameter=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'n_polar_a'); parameter.n_polar_a = 7; end; %<-- parameter_bookmark. ;
n_polar_a = parameter.n_polar_a;
if ~isfield(parameter,'n_azimu_b'); parameter.n_azimu_b = 12; end; %<-- parameter_bookmark. ;
n_azimu_b = parameter.n_azimu_b;
if ~isfield(parameter,'n_gamma_z'); parameter.n_gamma_z = 128; end; %<-- parameter_bookmark. ;
n_gamma_z = parameter.n_gamma_z;
if ~isfield(parameter,'k_max'); parameter.k_max = 1.0; end; %<-- parameter_bookmark. ;
k_max = parameter.k_max;
if ~isfield(parameter,'k_mid'); parameter.k_mid = 0.875*parameter.k_max; end; %<-- parameter_bookmark. ;
k_mid = parameter.k_mid;
if ~isfield(parameter,'facealpha_use'); parameter.facealpha_use = 1.00; end; %<-- parameter_bookmark. ;
facealpha_use = parameter.facealpha_use;
if ~isfield(parameter,'facecolor_use'); parameter.facecolor_use = 0.85*[1,1,1]; end; %<-- parameter_bookmark. ;
facecolor_use = parameter.facecolor_use;
if ~isfield(parameter,'linecolor_a'); parameter.linecolor_a = 0.75*[1,1,1]; end; %<-- parameter_bookmark. ;
linecolor_a = parameter.linecolor_a;
if ~isfield(parameter,'linecolor_b'); parameter.linecolor_b = 0.95*[1,1,1]; end; %<-- parameter_bookmark. ;
linecolor_b = parameter.linecolor_b;
if ~isfield(parameter,'linewidth_use'); parameter.linewidth_use = 2; end; %<-- parameter_bookmark. ;
linewidth_use = parameter.linewidth_use;
if ~isfield(parameter,'flag_solid'); parameter.flag_solid = 0; end; %<-- parameter_bookmark. ;
flag_solid = parameter.flag_solid;

r_max = k_max*(1-1e-2);
r_mid = k_mid*(1-1e-2);

gamma_z_ = linspace(0,2*pi,n_gamma_z+1); gamma_z_ = transpose(gamma_z_(1:n_gamma_z+1));
cc_ = cos(gamma_z_);
cc_pre_ = 0.5*cc_ + 0.5*circshift(cc_,+1);
cc_pos_ = 0.5*cc_ + 0.5*circshift(cc_,-1);
sc_ = sin(gamma_z_);
sc_pre_ = 0.5*sc_ + 0.5*circshift(sc_,+1);
sc_pos_ = 0.5*sc_ + 0.5*circshift(sc_,-1);

Rz = @(azimu_b) ...
[ +cos(azimu_b) -sin(azimu_b) 0 ; ...
  +sin(azimu_b) +cos(azimu_b) 0 ; ...
   0             0            1 ; ...
] ;

Ry = @(polar_a) ...
[ +cos(polar_a) 0 +sin(polar_a) ; ...
   0            1  0            ; ...
  -sin(polar_a) 0 +cos(polar_a) ; ...
];

hold on;
%%%%%%%%;
if flag_solid;
ca_ = zeros(n_polar_a,1);
sa_ = zeros(n_polar_a,1);
for npolar_a=0:n_polar_a-1;
polar_a = npolar_a*pi/max(1,n_polar_a-1);
ca_(1+npolar_a) = cos(polar_a);
sa_(1+npolar_a) = sin(polar_a);
end;%for npolar_a=0:n_polar_a-1;
tmp_x__ = zeros(2*n_polar_a,n_azimu_b);
tmp_y__ = zeros(2*n_polar_a,n_azimu_b);
tmp_z__ = zeros(2*n_polar_a,n_azimu_b);
for nazimu_b=0:n_azimu_b-1;
azimu_b_pre = (nazimu_b+0)*2*pi/max(1,n_azimu_b);
azimu_b_pos = (nazimu_b+1)*2*pi/max(1,n_azimu_b);
cb_pre = cos(azimu_b_pre); sb_pre = sin(azimu_b_pre);
cb_pos = cos(azimu_b_pos); sb_pos = sin(azimu_b_pos);
tmp_x__(:,1+nazimu_b) = [ cb_pre * sa_ ; cb_pos * flipud(sa_) ];
tmp_y__(:,1+nazimu_b) = [ sb_pre * sa_ ; sb_pos * flipud(sa_) ];
tmp_z__(:,1+nazimu_b) = [ ca_ ; flipud(ca_) ];
end;%for nazimu_b=0:n_azimu_b-1;
tmp_c_ = repmat(reshape(facecolor_use,[1,1,3]),[1,n_azimu_b,1]);
patch(tmp_x__,tmp_y__,tmp_z__,tmp_c_,'LineStyle','none','FaceAlpha',facealpha_use);
end;%if flag_solid;
%%%%%%%%;
for npolar_a=0:n_polar_a-1;
polar_a = npolar_a*pi/max(1,n_polar_a-1);
ca = cos(polar_a); sa = sin(polar_a);
plot3(r_max*sa*cc_,r_max*sa*sc_,r_max*ca*ones(n_gamma_z+1,1),'Color',linecolor_a,'LineWidth',linewidth_use);
plot3(r_mid*sa*cc_,r_mid*sa*sc_,r_max*ca*ones(n_gamma_z+1,1),'Color',linecolor_a,'LineWidth',linewidth_use);
tmp_x_ = [ ...
  reshape(r_mid*sa*cc_pre_,[1,n_gamma_z+1]) ...
; reshape(r_mid*sa*cc_pos_,[1,n_gamma_z+1]) ...
; reshape(r_max*sa*cc_pos_,[1,n_gamma_z+1]) ...
; reshape(r_max*sa*cc_pre_,[1,n_gamma_z+1]) ...
];
tmp_y_ = [ ...
  reshape(r_mid*sa*sc_pre_,[1,n_gamma_z+1]) ...
; reshape(r_mid*sa*sc_pos_,[1,n_gamma_z+1]) ...
; reshape(r_max*sa*sc_pos_,[1,n_gamma_z+1]) ...
; reshape(r_max*sa*sc_pre_,[1,n_gamma_z+1]) ...
];
tmp_z_ = [ ...
  r_max*ca*ones(1,n_gamma_z+1) ...
; r_max*ca*ones(1,n_gamma_z+1) ...
; r_max*ca*ones(1,n_gamma_z+1) ...
; r_max*ca*ones(1,n_gamma_z+1) ...
];
tmp_c_ = repmat(reshape(facecolor_use,[1,1,3]),[1,n_gamma_z+1,1]);
patch(tmp_x_,tmp_y_,tmp_z_,tmp_c_,'LineStyle','none','FaceAlpha',facealpha_use);
end;%for npolar_a=0:n_polar_a-1;
%%%%%%%%;
for nazimu_b=0:n_azimu_b-1;
azimu_b = nazimu_b*2*pi/max(1,n_azimu_b);
tmp_ = Rz(azimu_b)*[transpose(cc_);zeros(1,n_gamma_z+1);transpose(sc_)];
plot3(r_max*tmp_(1+0,:),r_max*tmp_(1+1,:),r_max*tmp_(1+2,:),'Color',linecolor_b,'LineWidth',linewidth_use);
clear tmp_;
end;%for nazimu_b=0:n_azimu_b-1;
%%%%%%%%;
hold off;


