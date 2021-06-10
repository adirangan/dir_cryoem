function P1__ = mp_init_surface_repulsion_wrap_0(P0__,n_neighbor,n_iteration);
% wrapper for mp_init_surface_repulsion. ;
if (nargin<2); n_neighbor = []; end;
if (nargin<3); n_iteration = []; end;
if isempty(n_neighbor); n_neighbor = 1; end;
if isempty(n_iteration); n_iteration = 1024; end;

n_M = size(P0__,1);
P1__ = P0__./repmat(sqrt(sum(P0__.^2,2)),[1,3]); %<-- normalize. ;
dP_avg = sqrt(4*pi^2 / n_M); %<-- average distance between uniformly spaced points on the sphere. ;
dt = dP_avg/16; %<-- timestep for evolution. ;
for niteration=0:n_iteration-1;
[dP1__,dP2_l2_min] = mp_init_surface_repulsion_0(n_M,P1__,n_neighbor);
dt = dP2_l2_min/8; %<-- timestep. ;
P1__ = P1__ + dt*dP1__;
end;%for niteration=0:n_iteration-1;

