function ...
[ ...
 parameter ...
,a_x_u___ ...
,n_Atom ...
,Atom_dilated_X_ ...
,Atom_dilated_Y_ ...
,Atom_dilated_Z_ ...
,Atom_occupancy_ ...
] = ...
a_x_u___from_pdb_0( ...
 parameter ...
,fname_pdb ...
,n_x_u ...
,Pixel_Spacing ...
);

str_thisfunction = 'a_x_u___from_pdb_0';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
platform = 'rusty';%platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data/rangan'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home/rangan'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home/rangan'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home/rangan'; end;
if (strcmp(platform,'ceph')); setup_ceph; string_root = 'mnt/home/rangan/ceph'; end;
fname_pdb = sprintf('/%s/dir_cryoem/dir_RNA0/models/open.pdb',string_root);
parameter = struct('type','parameter');
parameter.flag_disp = 1;
parameter.flag_center = 1;
n_x_u = 128;
Pixel_Spacing = 2.676;
[ ...
 parameter ...
,a_x_u___ ...
] = ...
a_x_u___from_pdb_0( ...
 parameter ...
,fname_pdb ...
,n_x_u ...
,Pixel_Spacing ...
);
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); fname_pdb=[]; end; na=na+1;
if (nargin<1+na); n_x_u=[]; end; na=na+1;
if (nargin<1+na); Pixel_Spacing=[]; end; na=na+1;


if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp=0; end;
flag_disp=parameter.flag_disp; nf=0;
if ~isfield(parameter,'flag_center'); parameter.flag_center=1; end;
flag_center=parameter.flag_center; nf=0;
if ~isfield(parameter,'sigma_angstrom'); parameter.sigma_angstrom=8.0; end;
sigma_angstrom=parameter.sigma_angstrom;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(Pixel_Spacing); Pixel_Spacing = 1.0; end;

half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
half_diameter_angstrom = n_x_u/2 * Pixel_Spacing;
x_p_r_max = 1.0;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u);
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u^3;
xxx_u_weight_ = (2*x_p_r_max/n_x_u)^3;
sigma_gaussian = sigma_angstrom/max(1e-12,half_diameter_angstrom)*x_p_r_max ;
g = @(x0,x1,x2) 1/sqrt(2*pi)^3 / max(1e-12,sigma_gaussian^3) * exp( -( (x_u_0___ - x0).^2 + (x_u_1___ - x1).^2 + (x_u_2___ - x2).^2 )/max(1e-12,2*sigma_gaussian^2));
%%%%;
pdb = pdbread(fname_pdb);
Atom_ = pdb.Model.Atom;
n_Atom = numel(Atom_);
Atom_X_ = zeros(n_Atom,1);
Atom_Y_ = zeros(n_Atom,1);
Atom_Z_ = zeros(n_Atom,1);
Atom_occupancy_ = zeros(n_Atom,1);
for nAtom=0:n_Atom-1;
Atom_X_(1+nAtom) = Atom_(1+nAtom).X;
Atom_Y_(1+nAtom) = Atom_(1+nAtom).Y;
Atom_Z_(1+nAtom) = Atom_(1+nAtom).Z;
Atom_occupancy_(1+nAtom) = Atom_(1+nAtom).occupancy;
end;%for nAtom=0:n_Atom-1;
%%%%;
Atom_centered_X_ = Atom_X_;
Atom_centered_Y_ = Atom_Y_;
Atom_centered_Z_ = Atom_Z_;
if flag_center;
Atom_avg_X_ = mean(Atom_centered_X_.*Atom_occupancy_);
Atom_avg_Y_ = mean(Atom_centered_Y_.*Atom_occupancy_);
Atom_avg_Z_ = mean(Atom_centered_Z_.*Atom_occupancy_);
Atom_centered_X_ = Atom_centered_X_ - Atom_avg_X_;
Atom_centered_Y_ = Atom_centered_Y_ - Atom_avg_Y_;
Atom_centered_Z_ = Atom_centered_Z_ - Atom_avg_Z_;
end;%if flag_center;
Atom_dilated_X_ = Atom_centered_X_ / max(1e-12,half_diameter_angstrom) ;
Atom_dilated_Y_ = Atom_centered_Y_ / max(1e-12,half_diameter_angstrom) ;
Atom_dilated_Z_ = Atom_centered_Z_ / max(1e-12,half_diameter_angstrom) ;
%%%%;
a_x_u___ = zeros(n_x_u,n_x_u,n_x_u);
for nAtom=0:n_Atom-1;
x = Atom_dilated_X_(1+nAtom);
y = Atom_dilated_Y_(1+nAtom);
z = Atom_dilated_Z_(1+nAtom);
occupancy = Atom_occupancy_(1+nAtom);
if (flag_verbose>1); disp(sprintf(' %% nAtom %d/%d: (%0.6f,%0.6f,%0.6f)',nAtom,n_Atom,x,y,z)); end;
a_x_u___ = a_x_u___ + occupancy*g(x,y,z);
end;%for nAtom=0:n_Atom-1;
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;figmed;
markersize_use = 16;
subplot(1,2,1);
plot3(Atom_dilated_X_,Atom_dilated_Y_,Atom_dilated_Z_,'.','MarkerSize',markersize_use,'MarkerFaceColor','c');
xlabel('x'); ylabel('y'); zlabel('z');
xlim([-1,+1]);
ylim([-1,+1]);
zlim([-1,+1]);
axisnotick3d; axis vis3d; grid on;
subplot(1,2,2);
isosurface_f_x_u_1([],a_x_u___);
end;%if flag_disp;
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
