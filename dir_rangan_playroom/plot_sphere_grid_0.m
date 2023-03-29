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
if ~isfield(parameter,'linecolor_use'); parameter.linecolor_use = 0.85*[1,1,1]; end; %<-- parameter_bookmark. ;
linecolor_use = parameter.linecolor_use;
if ~isfield(parameter,'linewidth_use'); parameter.linewidth_use = 2; end; %<-- parameter_bookmark. ;
linewidth_use = parameter.linewidth_use;

r_max = k_max*(1-1e-2);

gamma_z_ = linspace(0,2*pi,n_gamma_z+1); gamma_z_ = transpose(gamma_z_(1:n_gamma_z+1));
cc_ = cos(gamma_z_);
sc_ = sin(gamma_z_);

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
for npolar_a=0:n_polar_a-1;
polar_a = npolar_a*pi/max(1,n_polar_a-1);
ca = cos(polar_a); sa = sin(polar_a);
plot3(r_max*sa*cc_,r_max*sa*sc_,r_max*ca*ones(n_gamma_z+1,1),'Color',linecolor_use,'LineWidth',linewidth_use);
end;%for npolar_a=0:n_polar_a-1;
%%%%%%%%;
for nazimu_b=0:n_azimu_b-1;
azimu_b = nazimu_b*2*pi/max(1,n_azimu_b);
tmp_ = Rz(azimu_b)*[transpose(cc_);zeros(1,n_gamma_z+1);transpose(sc_)];
plot3(r_max*tmp_(1+0,:),r_max*tmp_(1+1,:),r_max*tmp_(1+2,:),'Color',linecolor_use,'LineWidth',linewidth_use);
clear tmp_;
end;%for nazimu_b=0:n_azimu_b-1;
%%%%%%%%;
hold off;


