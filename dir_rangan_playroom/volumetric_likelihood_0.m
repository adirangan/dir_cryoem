function ...
[ ...
 parameter ...
] = ...
volumetric_likelihood_0( ...
 parameter ...
,n_x_u_0in ...
,a_x_u_load_0in_ ...
,n_x_u_pack_0in ...
,a_x_u_pack_0in_ ...
,k_p_r_max_0in ...
,k_eq_d_0in ...
,n_M_0in ...
,M_x_u_0in_xxM___ ...
,euler_polar_a_0in_M_ ...
,euler_azimu_b_0in_M_ ...
,euler_gamma_z_0in_M_ ...
);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_x_u_0in=[]; end; na=na+1;
if (nargin<1+na); a_x_u_load_0in_=[]; end; na=na+1;
if (nargin<1+na); n_x_u_pack_0in=[]; end; na=na+1;
if (nargin<1+na); a_x_u_pack_0in_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max_0in=[]; end; na=na+1;
if (nargin<1+na); k_eq_d_0in=[]; end; na=na+1;
if (nargin<1+na); n_M_0in=[]; end; na=na+1;
if (nargin<1+na); M_x_u_0in_xxM___=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_0in_M_=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_0in_M_=[]; end; na=na+1;
if (nargin<1+na); euler_gamma_z_0in_M_=[]; end; na=na+1;

% parameter=[]; n_x_u_0in=[]; a_x_u_load_0in_=[]; n_x_u_pack_0in=[]; a_x_u_pack_0in_=[]; k_p_r_max_0in=[]; k_eq_d_0in=[]; n_M_0in=[]; M_x_u_0in_xxM___=[]; euler_polar_a_0in_M_=[]; euler_azimu_b_0in_M_=[]; euler_gamma_z_0in_M_=[];

str_thisfunction = 'volumetric_likelihood_0';

if isempty(n_x_u_0in); n_x_u_0in = 256; end; %<-- this is a typical resolution for a volume loaded from the protein database. ;
if isempty(a_x_u_load_0in_); a_x_u_load_0in_ = []; end; %<-- this is the loaded volume in x-space cartesian coordinates. ;
if isempty(n_x_u_pack_0in); n_x_u_pack_0in = 64; end; %<-- this is a lower resolution for a downsampled volume we will work with. ;
if isempty(a_x_u_pack_0in_); a_x_u_pack_0in_ = []; end; %<-- this is the downsampled volume in x-space cartesian coordinates. ;
if isempty(k_p_r_max_0in); k_p_r_max_0in = 48.0/(2*pi); end; %<-- k_p_r_max is the maximum frequency. ;
if isempty(k_eq_d_0in); k_eq_d_0in = 1.0/(2*pi); end; %<-- k_eq_d is the distance (along the equator) between gridpoints. ;
if isempty(n_M_0in); n_M_0in = []; end; %<-- this is the number of images. ;
if isempty(M_x_u_0in_xxM___); M_x_u_0in_xxM___ = []; end; %<-- this is the image-stack, assumed to be of dimension [n_x_u_pack,n_x_u_pack,n_M]. ;
if isempty(euler_polar_a_0in_M_); euler_polar_a_0in_M_ = []; end; %<-- this are the euler_polar_a angles for each image. ;
if isempty(euler_azimu_b_0in_M_); euler_azimu_b_0in_M_ = []; end; %<-- this are the euler_azimu_b angles for each image. ;
if isempty(euler_gamma_z_0in_M_); euler_gamma_z_0in_M_ = []; end; %<-- this are the euler_gamma_z angles for each image. ;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master=1e-2; end;
tolerance_master=parameter.tolerance_master;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=1; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_setpath'); parameter.flag_setpath=1; end;
flag_setpath=parameter.flag_setpath;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp=1; end;
flag_disp=parameter.flag_disp; nf=0;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if flag_setpath;
%%%%%%%%;
% Define paths. ;
% You will have to write your own (short) 'platform.type' file.
% Then you may have to write your own 'setup_local.m' and define your own string_root. ;
%%%%%%%%;
tmp_t = tic();
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (flag_verbose> 0); disp(sprintf(' %% defining paths for %s',platform)); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
tmp_t = toc(tmp_t); if (flag_verbose> 0); disp(sprintf(' %% defining paths: %0.2fs',tmp_t)); end;
%%%%%%%%;
end;%if flag_setpath;

%%%%%%%%;
% Define grids for x_c_ (cartesian), ;
% sometimes called x_u_ (for uniform). ;
% This is in contrast to x_p_ (polar) (which is rare). ;
% Note: the x_c_ grids used in this test file are only used for visualization, ;
% and are not exactly the x_c_ grids that are 'pixel-aligned' to the images. ;
% The pixel-aligned x_c_ grids used when processing the images are described immediately below: ;
%%%%%%%%;
% Each image is assumed to be a torus of side-length 2 lying in (periodic) [-1,+1]^{2}. ;
% half_diameter_x_c is typically 1.0. ;
% n_x_M_u: integer number of pixels on the side of an image. ;
% typically n_x_M_u is 256 or 360 in a typical mrc or mrcs image-stack from EMPIAR. ;
% x_c_0_ and x_c_1_ are the first and second coordinates used to describe the image-pixels. ;
% Note that these are left-aligned to each (square) pixel: ;
% x_c_0_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c; %<-- records the array of left-sides of each pixel. ;
% x_c_1_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c; %<-- records the array of bottom-sides of each pixel. ;
% so, for example, x_c_0_ ranges from -1.0 all the way up to (+1.0 - diameter_x_c/n_x_M_u). ;
% This is so that we can use the type-1 and type-2 nufft, where the first omitted spatial-point (i.e., +1.0) ;
% is assumed to be the same as the first included spatial-point (i.e. -1.0). ;
% With the above notation, a real-space image M_x_c__ ;
% (read from an image-stack using rlnImageName_from_star_1 or rlnImageName_from_star_2) ;
% has the following organization: ;
% M_x_c__(1+nx0,1+nx1) is the value in pixel (1+nx0,1+nx1), ;
% where nx0 and nx1 are zero-based. ;
% To convert from M_x_c__ (2d array in real-space) to M_k_p_ (unrolled array in fourier-space polar coordinates): ;
% dx = diameter_x_c/n_x_M_u; %<-- pixel-width. ;
% M_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c__,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2 ;
% n_k_p_r: integer number of radii (i.e., rings) in the polar-grid. ;
% k_p_r_: double array of size n_k_p_r; k_p_r_(1+nk_p_r) is the radius of ring indexed by nk_p_r (zero-based). ;
% n_w_: integer array of size n_k_p_r; n_w_(1+nk_p_r) is the number of equispaced angular points on the ring of radius k_p_r_(1+nk_p_r). ;
%         angles (implicitly defined) are assumed to take the form of: ;
%         gamma_z_: double array of size n_w = n_w_(1+nk_p_r). ; Note that this is different for each ring. ;
%                 gamma_z_(1+nw) = (2*pi)*cast(nw,'double')/max(1.0,cast(n_w,'double')); %<-- where nw is zero-based. ;
%                 Thus, gamma_z_ (for a particular image-ring) ranges from 0 all the way up to (2*pi - dgamma), ;
%                 where dgamma is (2*pi)/n_w. ;
% The normalization constant at the end (i.e., sqrt(n_x_M_u^2)*dx^2) is just a convention. ;
%%%%%%%%;
n_x_u = n_x_u_0in; %<-- this is a typical resolution for a volume loaded from the protein database. ;
half_diameter_x_c = 1.0d0; %<-- assume the volume sits in a box of side-length 2. ;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
n_x_u_pack = n_x_u_pack_0in; %<-- this is a lower resolution for a downsampled volume we will work with. ;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack); %<-- downsampled grid for x_u_0_. ;
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack); %<-- downsampled grid for x_u_1_. ;
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack); %<-- downsampled grid for x_u_2_. ;
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;
xxx_u_weight_ = (2*x_p_r_max/n_x_u_pack)^3;
%%%%%%%%;
n_pack = n_x_u/n_x_u_pack;
pack_row_ij_ = zeros(n_x_u_pack,1);
pack_col_ij_ = zeros(n_x_u_pack,1);
pack_val_ij_ = zeros(n_x_u_pack,1);
na=0;
for nx_u=0:n_x_u-1;
pack_row_ij_(1+na) = 1+nx_u;
pack_col_ij_(1+na) = 1+floor(nx_u/n_pack);
pack_val_ij_(1+na) = 1/n_pack;
na=na+1;
end;%for nx_u=0:n_x_u-1;
x_u_pack_ = sparse(pack_row_ij_,pack_col_ij_,pack_val_ij_,n_x_u,n_x_u_pack);
%%%%%%%%;
% If there were a volume a_x_u_load_ provided on a uniform grid of size n_x_u^3, ;
% this volume could be converted to a volume a_x_u_pack_ on the grid of size n_x_u_pack^3 via: ;
%%%%%%%%;
if ~isempty(a_x_u_pack_0in_); a_x_u_pack_ = a_x_u_pack_0in_; end;
if ( ~isempty(a_x_u_load_0in_) &  isempty(a_x_u_pack_0in_) );
a_x_u_pack_ = reshape(a_x_u_load_0in_,[n_x_u*n_x_u,n_x_u])*x_u_pack_;
a_x_u_pack_ = reshape(permute(reshape(a_x_u_pack_,[n_x_u,n_x_u,n_x_u_pack]),[3,1,2]),[n_x_u*n_x_u_pack,n_x_u])*x_u_pack_;
a_x_u_pack_ = reshape(permute(reshape(a_x_u_pack_,[n_x_u_pack,n_x_u,n_x_u_pack]),[3,1,2]),[n_x_u_pack*n_x_u_pack,n_x_u])*x_u_pack_;
a_x_u_pack_ = permute(reshape(a_x_u_pack_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),[3,1,2]);
end;%if ( ~isempty(a_x_u_load_0in_) &  isempty(a_x_u_pack_0in_) );
if (  isempty(a_x_u_load_0in_) &  isempty(a_x_u_pack_0in_) );
rng(0);
a_x_u_pack_ = zeros(n_x_u_pack,n_x_u_pack,n_x_u_pack);
n_source = 16; tmp_sigma = x_p_r_max/8;
for nsource=0:n_source-1;
tmp_delta_ = x_p_r_max*randn(3,1); tmp_delta_f = fnorm(tmp_delta_); if tmp_delta_f>x_p_r_max/sqrt(2); tmp_delta_ = tmp_delta_/tmp_delta_f*x_p_r_max/sqrt(2); end;
tmp_x_u_r___ = sqrt( (x_u_0___-tmp_delta_(1+0)).^2 + (x_u_1___-tmp_delta_(1+1)).^2 + (x_u_2___-tmp_delta_(1+2)).^2 );
tmp_x_u_g___ = 1/sqrt(2*pi)^3 / tmp_sigma^3 * exp(-tmp_x_u_r___.^2 / (2*tmp_sigma^2));
a_x_u_pack_ = a_x_u_pack_ + tmp_x_u_g___;
clear tmp_delta_ tmp_x_u_r___ tmp_x_u_g___ ;
end;%for nsource=0:n_source-1;
end;%if (  isempty(a_x_u_load_0in_) &  isempty(a_x_u_pack_0in_) );
%%%%%%%%;
% Here are grids for k_c_ (cartesian), sometimes referred to as k_u_ (uniform). ;
%%%%%%%%;
k_u_0_ = periodize(0:n_x_u_pack-1,-n_x_u_pack/2,+n_x_u_pack/2)/2; %<-- box has diameter 2. ;
k_u_1_ = periodize(0:n_x_u_pack-1,-n_x_u_pack/2,+n_x_u_pack/2)/2; %<-- box has diameter 2. ;
k_u_2_ = periodize(0:n_x_u_pack-1,-n_x_u_pack/2,+n_x_u_pack/2)/2; %<-- box has diameter 2. ;
[k_u_0___,k_u_1___,k_u_2___] = ndgrid(k_u_0_,k_u_1_,k_u_2_); n_kkk_u = n_x_u_pack^3;
%%%%%%%%;
if (flag_disp> 1);
%%%%%%%%;
% visualize x_c_ grid. ;
% This does not matter too much since we only really use this for visualization. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
markersize_use = 8;
plot3(x_u_0___(:),x_u_1___(:),x_u_2___(:),'ko','MarkerSize',markersize_use,'MarkerFaceColor','c');
axis(half_diameter_x_c*[-1,+1,-1,+1,-1,+1]); axis vis3d;
xlabel('x0'); ylabel('x1'); zlabel('x2');
title('x_c_','Interpreter','none');
end;%if (flag_disp> 0);
%%%%%%%%%;
if (flag_disp> 0);
%%%%%%%%;
% visualize a_x_u_pack_. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
isosurface_f_x_u_1([],a_x_u_pack_); axisnotick; axis vis3d;
xlabel('x0'); ylabel('x1'); zlabel('x2');
title('a_x_u_pack_','Interpreter','none');
end;%if (flag_disp> 0);
%%%%%%%%%;

%%%%%%%%;
% Now set up grids for (3d) k_p_ (polar). ;
% For historical reasons, the suffix '_all_' refers to the unrolled gridpoints associated with the sphere in k_p_ (in 3d). ;
% This can be confusing later on, since the suffix '_all_' is also used to refer to ;
% the unrolled gridpoints associated with a disk in k_p_ (in 2d), ;
% as well as the unrolled list of viewing-angles associated with a spherical shell (for templates). ;
%%%%%%%%;
% Meanwhile, the suffix '_k_' refers to values associated with the list of radii k_p_r_. ;
% Unfortunately, I misnamed the variable 'weight_shell_k_' long ago, and it should really be called 'weight_shell_all_'. ;
% Fortunately, this variable is not used frequently. ;
%%%%%%%%;
% The k_p_ grid itself is defined in terms of the following: ;
% n_k_p_r: integer number of shells for the spherical grid. ; %<-- strongly assume this is the same as the n_k_p_r used for the (2d) k_p_. ;
% k_p_r_: double array of size n_k_p_r; k_p_r_(1+nk_p_r) is the radius of the shell with index nk_p_r (zero-based). ;
%         strongly assume that this k_p_r_ is the same as the one used for the (2d) k_p_. ;
% n_polar_a_k_: integer array of size n_k_p_r; n_polar_a_k_(1+nk_p_r) is the number of latitudes (i.e., polar_a values) ;
%         associated with shell nk_p_r. ;
% polar_a_ka__: cell array of size n_k_p_r. polar_a_ka__{1+nk_p_r} is the array of polar_a values for shell nk_p_r. ;
%         For shell nk_p_r, ;
%         n_polar_a = n_polar_a_k_(1+nk_p_r) is an integer number of latitudes for that shell. ;
%         polar_a_a_ := polar_a_ka__{1+nk_p_r} is an array of size n_polar_a. ; 
%                 polar_a_a_(1+npolar_a) is a double storing the polar_a for latitude-line npolar_a. ;
% n_azimu_b_ka__: cell array of size n_k_p_r. ;
%         For shell nk_p_r, ;
%         n_polar_a = n_polar_a_k_(1+nk_p_r) is an integer number of latitudes for that shell. ;
%         n_azimu_b_a_ := n_azimu_b_ka__{1+nk_p_r} is an array of size n_polar_a. ; 
%                 n_azimu_b_a_(1+npolar_a) is an integer storing the number of azimu_b values for latitude-line npolar_a. ;
% weight_3d_k_all_: double array of size n_k_all, storing quadrature (integration) weights for each point in the unrolled array. ;
%%%%%%%%;
tmp_t = tic;
k_p_r_max = k_p_r_max_0in; %<-- k_p_r_max is the maximum frequency. ;
% Note that we technically use the 'wavenumber' rather than the 'frequency' (related by 2*pi), so the exponentials involved in our calculations will look like exp(+i*2*pi*k_p_r*x_p_r). ;
k_eq_d = k_eq_d_0in; %<-- k_eq_d is the distance (along the equator) between gridpoints. ;
% Note that (due to the support of a_x_u_) we expect the bandlimit of a_k_p_ to be 1.0 in terms of frequency (or 1.0/(2*pi) in terms of wavenumber). ;
[ ...
 n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,~ ...
,~ ...
,~ ...
,~ ...
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
] = ...
sample_sphere_7( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,'L' ...
) ; %<-- sum(weight_3d_k_p_r_)*(4*pi) = (4/3)*pi*k_p_r_max^3 --> sum(weight_3d_k_p_r_) = (1/3)*k_p_r_max^3 ;
tmp_t = toc(tmp_t); if (flag_verbose> 0); disp(sprintf(' %% sample_sphere_7: time %0.2fs',tmp_t)); end;
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_quad_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_x_u_pack_(:).*xxx_u_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); if (flag_verbose> 0); disp(sprintf(' %% xxnufft3d3: a_k_p_quad_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_u_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_quad_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); if (flag_verbose> 0); disp(sprintf(' %% xxnufft3d3: a_x_u_reco_ time %0.2fs',tmp_t)); end;
disp(sprintf(' %% xxnufft3d3: a_x_u_reco error: %0.16f',fnorm(a_x_u_pack_(:)-a_x_u_reco_)/fnorm(a_x_u_pack_(:))));
%%%%%%%%;
if (flag_disp> 1);
figure(1+nf);nf=nf+1;clf;figsml;
makersize_use = 4;
nk_p_r = n_k_p_r-1;
tmp_index_ = n_k_all_csum_(1+nk_p_r):n_k_all_csum_(1+nk_p_r+1)-1;
plot3(k_c_0_all_(1+tmp_index_),k_c_1_all_(1+tmp_index_),k_c_2_all_(1+tmp_index_),'k.','MarkerSize',markersize_use);
axis(k_p_r_max*[-1,+1,-1,+1,-1,+1]); axis vis3d;
xlabel('k0'); ylabel('k1'); zlabel('k2');
title(sprintf(' shell nk_p_r=%d k_p_ plotted in k_c_',nk_p_r),'Interpreter','none');
end;%if flag_disp;
%%%%%%%%;
if (flag_disp> 0);
%%%%%%%%;
% visualize a_x_u_reco_. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
isosurface_f_x_u_1([],a_x_u_reco_); axisnotick; axis vis3d;
xlabel('x0'); ylabel('x1'); zlabel('x2');
title('a_x_u_reco_','Interpreter','none');
end;%if (flag_disp> 0);
%%%%%%%%%;

%%%%%%%%;
% Now set up grids for k_Y_ (spherical_harmonic). ;
%%%%%%%%;
l_max_upb = round(2*pi*k_p_r_max); %<-- typically sufficient for 2-3 digits of precision. ;
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
l_max_(1+nk_p_r) = max(0,min(l_max_upb,1+ceil(2*pi*k_p_r_(1+nk_p_r))));
end;%for nk_p_r=0:n_k_p_r-1;
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
Y_l_max_val_ = zeros(n_lm_sum,1);
Y_l_val_ = zeros(n_lm_sum,1);
Y_m_val_ = zeros(n_lm_sum,1);
Y_k_val_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
tmp_l_val_ = zeros(n_lm_(1+nk_p_r),1);
tmp_m_val_ = zeros(n_lm_(1+nk_p_r),1);
na=0; 
for l_val=0:l_max;
for m_val=-l_val:+l_val;
tmp_l_val_(1+na) = l_val;
tmp_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
Y_l_max_val_(1+tmp_index_) = l_max;
Y_l_val_(1+tmp_index_) = tmp_l_val_;
Y_m_val_(1+tmp_index_) = tmp_m_val_;
Y_k_val_(1+tmp_index_) = k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
weight_Y_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
weight_Y_(1+tmp_index_) = weight_3d_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_t = tic;
[a_k_Y_quad_yk_] = convert_k_p_to_spharm_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_quad_);
tmp_t = toc(tmp_t); if (flag_verbose> 0); disp(sprintf(' %% a_k_Y_quad_yk_ time %0.2fs',tmp_t)); end;
disp(sprintf(' %% Here we should ensure that a_k_Y_ on the outer shells has decayed to the desired precision.'));
disp(sprintf(' %% Moreover, we should ensure that a_k_Y_(l,m) has decayed for large l,m at each shell.'));
tmp_t = tic;
[a_k_p_reco_] = convert_spharm_to_k_p_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_quad_yk_);
tmp_t = toc(tmp_t); if (flag_verbose> 0); disp(sprintf(' %% a_k_Y_quad_yk_ --> a_k_p_reco_ time %0.2fs',tmp_t)); end;
disp(sprintf(' %% xxnufft3d3: a_k_p_reco error: %0.16f',fnorm(a_k_p_quad_-a_k_p_reco_)/fnorm(a_k_p_quad_))); %<-- this should be 2-3 digits. ;
%%%%%%%%;
% Here a_k_Y_quad_yk_ is a one-dimensional array, as indicated by the final '_'. ;
% I.e., a_k_Y_quad_yk_ stores the unrolled array of spherical-harmonic coefficients ;
% (with spherical-harmonic index y varying quickly, and the radial index k varying slowly) ;
% for each shell one after another. ;
% We can instead create an array a_k_Y_quad_yk__ which is of dimension 2, ;
% (again with spherical-harmonic index y varying quickly, and the radial index k varying slowly). ;
% Often this conversion can be done using matlab 'reshape'. ;
% However this only works for rectangular arrays. ;
% In this case we actually have fewer total spherical-harmonic coefficients for the smaller shells. ;
% (i.e., we determine l_max_(1+nk_p_r) adaptively). ;
% So the linear a_k_p_quad_yk_ is smaller than full a_k_p_quad_yk__, ;
% the latter of which is rectangular (2d) and padded with zeros. ;
%%%%%%%%;
a_k_Y_quad_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_quad_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_quad_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;

%%%%%%%%;
% Use a_k_Y_quad_yk_ to generate templates S_k_p_wkS__ on a 2d k_p_ grid with uniform n_w_. ;
% The array n_w_ is shorthand for n_gamma_k_ (or n_omega_k_); ;
% i.e., the number of gridpoints of (2d) k_p_ on the ring of radius k_p_r_(1+nk_p_r) (in 2d). ;
%%%%%%%%;
% Here 'template' means 'projection of volume at some viewing-angle'. ;
% Each template S_k_p_wk_ (or S_k_p_ for short) corresponds to the 1-dimensonal array of ;
% coefficients, each associated with a gridpoint in (2d) k_p_ space. ;
% For this 1d-array, the w-index (corresponding to polar-angle in the image) varies quickly, ;
% and the k-index (corresponding to the radius of the image-ring) varies slowly. ;
% The same convention is used for images (labeled M_k_p_wk_ or M_k_p_ later on). ;
% There are other representations of each S_k_p_wk_, such as: ;
% S_k_p_wk__ <-- 2d rectangular array with polar-angle varying quickly and k varying slowly. ;
% S_k_q_wk__ or S_k_q_ <-- 1d array storing fourier-bessel coefficients. ;
% S_x_p_wx_ <-- 1d array storing template in real-space polar-coordinates (rare). ;
% S_x_c__ or S_x_u__ <-- 2d array storing template in real-space cartesian-coordinates. ;
% S_x_c_ or S_x_u_ <-- 1d unrolled version of S_x_c__ or S_x_u__. ; 
%%%%%%%%;
% While an adaptive n_w_ might be tempting, ;
% blas3 acceleration of subsequent computations benefits from a uniform n_w_. ;
% For this reason I assume n_w_ is simply n_w_max*ones(n_k_p_r,1). ;
%%%%%%%%;
% Unfortunately, for historical reasons the suffix '_all_' here refers to the list of viewing-angles. ;
% (associated with the stack of templates S_k_p_). ;
% I have tried to disambiguate these variables from those associated with the spherical grid above ;
% by using the prefix 'viewing' or 'template_viewing'. ;
% So here n_S := n_viewing_all is the number of templates generates. ;
% The viewing-angle for template index nS is stored in: ;
% viewing_azimu_b_all_(1+nS) and viewing_polar_a_all_(1+nS). ;
%%%%%%%%;
n_w_max = 2*(l_max_max+1);
n_w_0in_ = n_w_max*ones(n_k_p_r,1);
template_k_eq_d = -1;
viewing_k_eq_d = 1.0; %<-- subsample just a few templates for visualization. ;
tmp_t = tic();
[ ...
 n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,~ ...
,~ ...
,~ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
] = ...
sample_shell_5( ...
 k_p_r_max ...
,viewing_k_eq_d ...
,'L' ...
) ; %<-- obtain viewing angles on outer shell. ;
n_viewing_azimu_b_sum = sum(n_viewing_azimu_b_);
n_viewing_azimu_b_csum_ = cumsum([0;n_viewing_azimu_b_]);
if (flag_verbose); disp(sprintf(' %% n_viewing_all %d n_viewing_polar_a %d n_viewing_azimu_b_sum %d',n_viewing_all,n_viewing_polar_a,n_viewing_azimu_b_sum)); end;
%%%%%%%%;
% Here we generate weights for k_p_ (in 2-dimensions). ;
% Unfortunately (again due to historical reasons), the suffix '_all_' refers to the unrolled list of k_p_ (in 2d). ;
%%%%%%%%;
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
] = ...
get_weight_2d_1( ...
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
if (flag_verbose); disp(sprintf(' %% n_w_max %d n_w_sum %d',n_w_max,n_w_sum)); end;
tmp_t = tic();
[ ...
 S_k_p_wkS__ ...
,~ ...
,n_S ...
,~ ...
,~ ...
] = ...
pm_template_2( ...
 flag_verbose ...
,l_max ...
,n_k_p_r ...
,reshape(a_k_Y_quad_yk__,[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
);
tmp_t = toc(tmp_t); if (flag_verbose> 0); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t)); end;
S_k_p_wkS__ = reshape(S_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
n_S = n_viewing_all; n_w_max = max(n_w_);
%%%%%%%%;

if ~isempty(n_M_0in); n_M = n_M_0in; end;
if  isempty(n_M_0in); n_M = n_S; end;
if ~isempty(M_x_u_0in_xxM___);
M_x_u_xxM___ = M_x_u_0in_xxM___;
end;%if ~isempty(M_x_u_0in_xxM___);
if  isempty(M_x_u_0in_xxM___);
n_M = n_S;
M_x_u_xxM___ = zeros(n_x_u_pack,n_x_u_pack,n_M);
for nM=0:n_M-1;
S_k_p_wk_ = S_k_p_wkS__(:,1+nM);
S_x_c_xx__ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u_pack ...
,diameter_x_c ...
,n_x_u_pack ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,S_k_p_wk_ ...
);
M_x_u_xxM___(:,:,1+nM) = reshape(S_x_c_xx__,[n_x_u_pack,n_x_u_pack]);
clear S_k_p_wk_ S_x_c_xx__;
end;%for nM=0:n_M-1;
end;%if  isempty(M_x_u_0in_xxM___);

if flag_disp;
%%%%%%%%;
% View various images. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 3; p_col = 5; np=0;
for np=0:p_row*p_col-1;
subplot(p_row,p_col,1+np);
nM = max(0,min(n_M,floor(n_M*np/(p_row*p_col))));
M_x_u_xx__ = M_x_u_xxM___(:,:,1+nM);
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,real(M_x_u_xx__),[],colormap_beach);
title(sprintf('nM %d real(M_x_u_xx__)',nM),'Interpreter','none'); axis image; axisnotick;
end;%for np=0:p_row*p_col-1;
end;%if flag_disp;

%%%%%%%%;
% Now we construct the stack of images in k-space polar coordinates. ;
%%%%%%%%;
tmp_t = tic();
M_k_p_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_x_u_xx__ = M_x_u_xxM___(:,:,1+nM);
M_k_p_wk_ = ...
interp_x_c_to_k_p_xxnufft( ...
 n_x_u_pack ...
,diameter_x_c ...
,n_x_u_pack ...
,diameter_x_c ...
,M_x_u_xx__ ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
) ;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose> 0); disp(sprintf(' %% interp_x_c_to_k_p_xxnufft: %0.2fs',tmp_t)); end;

if ~isempty(euler_polar_a_0in_M_); euler_polar_a_M_ = euler_polar_a_0in_M_; end;
if  isempty(euler_polar_a_0in_M_); euler_polar_a_M_ = viewing_polar_a_all_; end;
if ~isempty(euler_azimu_b_0in_M_); euler_azimu_b_M_ = euler_azimu_b_0in_M_; end;
if  isempty(euler_azimu_b_0in_M_); euler_azimu_b_M_ = viewing_azimu_b_all_; end;
if ~isempty(euler_gamma_z_0in_M_); euler_gamma_z_M_ = euler_gamma_z_0in_M_; end;
if  isempty(euler_gamma_z_0in_M_); euler_gamma_z_M_ = zeros(n_M,1); end;

%%%%%%%%;
% Currently: ;
% weight_3d_k_p_r_ is a double-array of size n_k_p_r. ;
%         This stores integration-weights for 3d-integration. ;
%         weight_3d_k_p_r_(1+nk_p_r) is the 3d-integration weight ;
%         associated with the shell of radius k_p_r_(1+nk_p_r). ;
% weight_3d_k_all_ is a double-array of size n_k_all. ;
%         This stores integration-weights for 3d-integration. ;
%         weight_3d_k_all_(1+nk_all) is the 3d-integration weight ;
%         associated with the point in (3d) k_p_ associated with index nk_all. ;
%         That is to say, the point associated with the k_c_ position: ;
%         [ k_c_0_all_(1+nk_all), k_c_1_all_(1+nk_all), k_c_2_all_(1+nk_all) ], ;
%         or with the spherical-coordinates: ;
%         [ k_p_r_all_(1+nk_all), k_p_azimu_b_all_(1+nk_all), k_p_polar_a_all_(1+nk_all) ]. ;
% weight_2d_k_p_r_ is a double-array of size n_k_p_r. ;
%         This stores integration-weights for 2d-integration. ;
%         weight_2d_k_p_r_(1+nk_p_r) is the 2d-integration weight ;
%         associated with the ring of radius k_p_r_(1+nk_p_r). ;
% weight_2d_k_all_ is a double-array of size n_w_sum. ;
%         This stores integration-weights for 2d-integration. ;
%         weight_2d_k_all_(1+n_w_csum_(1+nk_p_r) + nw) is the 2d-integration weight ;
%         associated with the point in (2d) k_p_ associated with radius k_p_r_(1+nk_p_r) ;
%         and polar-angle 2*pi*nw/n_w_(1+nk_p_r). ;
% viewing_weight_all_ is a double-array of size n_S. ;
%         This stores integration-weights for spherical-shell integration ;
%         associated with the templates S_k_p_wkS__. ;
%         viewing_weight_all_(1+nS) is the spherical-shell-integration weight ;
%         associated with template S_k_p_wkS__(:,1+nS). ;
%%%%%%%%;
% For future reference we construct weight_Y_2d_. ;
% weight_Y_2d_ is a double-array of size n_lm_sum. ;
%         This stores the values of weight_2d_k_p_r_, ;
%         referenced using the ordering of spherical-harmonic coefficients. ;
%         Thus, weight_Y_2d_(1+nk_all) will be equal to weight_2d_k_p_r_(1+nk_p_r), ;
%         where k_p_r_(1+nk_p_r) is equal to k_p_r_all_(1+nk_all). ;
%%%%%%%%;
weight_Y_2d_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
weight_Y_2d_(1+tmp_index_) = weight_2d_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
% For future reference we construct weight_2d_k_p_wkS__. ;
% weight_2d_k_p_wkS__ is a double-array of size (n_w_sum,n_S). ;
%         This stores a tensor-grid of 2d-integration weights. ;
%         weight_2d_k_p_wkS__(1+n_w_csum_(1+nk_p_r)+nw,1+nS) is proportional to ;
%         the 2d-integration weight weight_2d_k_all_(1+n_w_csum_(1+nk_p_r)+nw) ;
%         multiplied by the spherical-shell-integration weight viewing_weight_all_(1+nS). ;
%%%%%%%%;
weight_2d_k_p_wkS__ = reshape(weight_2d_k_all_*4*pi^2,[n_w_sum,1])*reshape(viewing_weight_all_/(k_p_r_max^2),[1,n_S]);
%%%%;
E_3d_k_Y = sum(weight_Y_);
E_3d_k_p = sum(weight_3d_k_all_);
E_2d_k_p = sum(weight_2d_k_p_wkS__, 'all' );
E_2d_k_Y = sum(weight_Y_2d_);
%%%%;
L2_3d_k_Y = sum(abs(a_k_Y_quad_yk_).^2 .* weight_Y_);
flag_check=0;
if flag_check;
L2_3d_k_p = sum(abs(a_k_p_quad_).^2 .* weight_3d_k_all_);
end;%if flag_check;
%%%%;
L2_2d_k_p = sum(abs(S_k_p_wkS__).^2 .* weight_2d_k_p_wkS__,'all');
L2_2d_k_Y = sum(abs(a_k_Y_quad_yk_).^2 .* weight_Y_2d_);
if (flag_verbose); disp(sprintf(' %% L2_2d_k_p %0.6f vs L2_2d_k_Y %0.6f: %0.16f',L2_2d_k_p,L2_2d_k_Y,fnorm(L2_2d_k_p - L2_2d_k_Y)/fnorm(L2_2d_k_Y))); end;
%%%%%%%%;

%%%%%%%%;
% Now reconstruct volumetric terms. ;
% Because of the way qbp_6 is structured, ;
% these calculations benefit from the (nearly) uniform-distribution of viewing-angles. ;
% A highly nonuniform-distribution of viewing-angles will introduce errors (due to qbp_6), ;
% in which case one should use a different algorithm for reconstruction (e.g., cg_lsq_6). ;
% Note that this ignores the CTF for now (see empty entries below). ;
%%%%%%%%;
% Note that the calculations below can be sped up by choosing a uniform l_max_. ;
% That is, qbp_6 is quite a bit slower when l_max_ is nonuniform (as is the case here). ;
% The same holds for convert_spharm_to_k_p. ;
%%%%%%%%;
tmp_t = tic();
qbp_eps = 1e-3;
[ ...
 M_avg_k_Y_yk_ ...
] = ...
qbp_6( ...
 qbp_eps ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,[] ...
,[] ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose> 0); disp(sprintf(' %% qbp_6: %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
tmp_parameter = struct('type','parameter');
tmp_parameter.flag_loading_svd_vs_iterate = 1;
tmp_parameter.flag_loading_skip_loading = 1;
[ ...
 tmp_parameter ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,R_k_p_wkM__ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
] = ...
get_loading_qbp_2( ...
 tmp_parameter ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,l_max_max*ones(n_k_p_r,1) ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,[] ...
,[] ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose> 0); disp(sprintf(' %% get_loading_qbp_2 time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
qbp_eps = 1e-3;
[ ...
 M_std_k_Y_yk_ ...
] = ...
qbp_6( ...
 qbp_eps ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,sqrt(abs(R_k_p_wkM__).^2) ...
,[] ...
,[] ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose> 0); disp(sprintf(' %% qbp_6: %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
qbp_eps = 1e-3;
[ ...
 M_var_k_Y_yk_ ...
] = ...
qbp_6( ...
 qbp_eps ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,abs(R_k_p_wkM__).^2 ...
,[] ...
,[] ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose> 0); disp(sprintf(' %% qbp_6: %0.2fs',tmp_t)); end;
%%%%%%%%;
% Now convert the k_Y representation to k_p. ;
% This can be sped up by using convert_spharm_to_k_p_4 (or one of the more recent versions) ;
% and performing the necessary precomputations. ;
%%%%%%%%;
tmp_t = tic;
[M_avg_k_p_] = convert_spharm_to_k_p_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,M_avg_k_Y_yk_);
tmp_t = toc(tmp_t); disp(sprintf(' %% M_avg_k_Y_yk_ --> M_avg_k_p_ time %0.2fs',tmp_t));
tmp_t = tic;
[M_std_k_p_] = convert_spharm_to_k_p_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,M_std_k_Y_yk_);
tmp_t = toc(tmp_t); disp(sprintf(' %% M_std_k_Y_yk_ --> M_std_k_p_ time %0.2fs',tmp_t));
tmp_t = tic;
[M_var_k_p_] = convert_spharm_to_k_p_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,M_var_k_Y_yk_);
tmp_t = toc(tmp_t); disp(sprintf(' %% M_var_k_Y_yk_ --> M_var_k_p_ time %0.2fs',tmp_t));
%%%%%%%%;

%%%%%%%%;
% Now we reconstruct viewing-angle distribution on the surface of the sphere, and store the results in jl_ab_all_. ;
% jl_ab_all_ is a double-array of size n_S. ;
%         jl_ab_all_(1+nS) stores the probability of observing the viewing-angle ;
%         associated with [ viewing_azimu_b_(1+nS) , viewing_polar_a_(1+nS) ]. ;
% We construct jl_ab_all_ with a simple nearest-neighbor algorithm. ;
%%%%%%%%;
jl_ab_all_ = zeros(n_S,1); 
viewing_k_c_0_all_ = cos(viewing_azimu_b_all_).*sin(viewing_polar_a_all_);
viewing_k_c_1_all_ = sin(viewing_azimu_b_all_).*sin(viewing_polar_a_all_);
viewing_k_c_2_all_ = cos(viewing_polar_a_all_);
viewing_k_c_3S__ = [reshape(viewing_k_c_0_all_,[1,n_S]) ; reshape(viewing_k_c_1_all_,[1,n_S]) ; reshape(viewing_k_c_2_all_,[1,n_S]) ];
euler_k_c_0_M_ = cos(euler_azimu_b_M_).*sin(euler_polar_a_M_);
euler_k_c_1_M_ = sin(euler_azimu_b_M_).*sin(euler_polar_a_M_);
euler_k_c_2_M_ = cos(euler_polar_a_M_);
euler_k_c_3M__ = [reshape(euler_k_c_0_M_,[1,n_M]) ; reshape(euler_k_c_1_M_,[1,n_M]) ; reshape(euler_k_c_2_M_,[1,n_M]) ];
tmp_sigma = 2*k_eq_d;
for nS=0:n_S-1; 
viewing_k_c_ = viewing_k_c_3S__(:,1+nS);
viewing_distance_M_ = sqrt(sum(abs(bsxfun(@minus,euler_k_c_3M__,viewing_k_c_)).^2,1));
jl_ab_all_(1+nS) = numel(efind(viewing_distance_M_<=tmp_sigma));
end;%for nS=0:n_S-1; 
jl_ab_all_ = jl_ab_all_/max(1e-12,sum(jl_ab_all_.*viewing_weight_all_)/(k_p_r_max^2));
if (flag_verbose); disp(sprintf(' %% sum(jl_ab_all_.*viewing_weight_all_)/(k_p_r_max^2): %0.16f (should be 1)',sum(jl_ab_all_.*viewing_weight_all_)/(k_p_r_max^2))); end;
%%%%%%%%;
% Now construct shell-complement of viewing-angle distribution. ;
% m_val_ is a double-array of size (1+2*l_max) storing the array of spherical-harmonic orders. ;
% Note that l_max_ is often chosen adaptively, so not all orders will be used on all shells. ;
% jl_hk_Y_ is a double-complex-array of size n_lm_max. ;
%         jl_hk_Y_ stores the volumetric-description of the function ;
%         \tilde{J}(\hk) from the overleaf notes. ;
%         Briefly, the function \tilde{J}(\hk) refers to the probability of observing
%         a particular viewing-angle on the sphere (indicated by \hk). ;
%         Note that \tilde{J}(\hk) is defined on a spherical-shell. ;
% jl_hk_lm__ is a double-complex-array of size(1+l_max,numel(m_val_)). ;
%         jl_hk_lm__ stores the values of jl_hk_Y_ in a rectangular (2d) array. ;
%%%%%%%%;
m_val_ = -l_max:+l_max;
%%%%%%%%;
qbp_eps = 1e-3;
[ ...
 jl_hk_Y_ ...
,n_quad_from_data_q_ ...
] = ...
qbp_6( ...
 qbp_eps ...
,1 ...
,1 ...
,l_max ...
,n_w_max ...
,n_S ...
,ones(n_w_max,1)*transpose(jl_ab_all_) ...
,[] ...
,[] ...
,viewing_polar_a_all_ ...
,viewing_azimu_b_all_ ...
,zeros(n_S,1) ...
,[] ...
,[] ...
,[]...
);
jl_hk_lm__=zeros(1+l_max,n_m_max);
for l_val=0:l_max;
index_0in_ = l_val^2 + [0:1+2*l_val-1];
index_out_ = l_max + [-l_val:+l_val];
jl_hk_lm__(1+l_val,1+index_out_) = jl_hk_Y_(1+index_0in_);
end;%for l_val=0:l_max;
if (flag_verbose); disp(sprintf(' %% jl_hk_lm__(1,1+l_max)*sqrt(4*pi): %0.16f (should be within 2-3 digits of 1)',jl_hk_lm__(1,1+l_max)*sqrt(4*pi))); end;
%%%%%%%%;
% jl_hk_k_all_ is a double-array of size n_S. ;
%         jl_hk_k_all_(1+nS) stores the probability of observing viewing-angle for nS. ;
%%%%%%%%;
jl_hk_k_all_ = ...
convert_spharm_to_k_p_2( ...
 0 ...
,n_viewing_all ...
,[0;n_viewing_all] ...
,ones(n_viewing_all,1) ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,[] ...
,viewing_weight_all_ ...
,1 ...
,1 ...
,1 ...
,l_max ...
,jl_hk_Y_ ...
);
if (flag_verbose); disp(sprintf(' %% sum(real(jl_hk_k_all_).*viewing_weight_all_/k_p_r_max^2): %0.16f (should be within 2-3 digits of 1)',sum(real(jl_hk_k_all_).*viewing_weight_all_/k_p_r_max^2))); end;
%%%%%%%%;
% jl_hk_yk_ is a double-array of size n_lm_sum, ;
% referencing the above distribution on a spherical-harmonic grid (including all shells). ;
%%%%%%%%;
jl_hk_yk_ = zeros(n_lm_sum,1);
for nlm_sum=0:n_lm_sum-1;
jl_hk_yk_(1+nlm_sum) = jl_hk_lm__(1+Y_l_val_(1+nlm_sum),1+l_max+Y_m_val_(1+nlm_sum));
end;%for nlm_sum=0:n_lm_sum-1;
%%%%%%%%;
if flag_disp;
%%%%%%%%;
% Demonstrate the nonuniform viewing-angle distribution (JL). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
imagesc_polar_a_azimu_b_0(viewing_polar_a_all_,viewing_azimu_b_all_,real(jl_ab_all_),[0,0.1],[],0);
xlabel('k0'); ylabel('k1'); zlabel('k2');
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
title('jl_ab_all_','Interpreter','none');
end;%if flag_disp;
%%%%%%%%;
% Now convert k_Y to k_p. ;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% jl_hk_k_p_')); end;
tmp_t = tic();
[jl_hk_k_p_] = ...
convert_spharm_to_k_p_2( ...
 0*flag_verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,jl_hk_yk_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% jl_hk_yk_ --> jl_hk_k_p_ time %0.2fs',tmp_t));
%%%%%%%%;

%%%%%%%%;
% We construct the riesz integration-weights on the sphere. ;
% These are associated with the riesz-potential 1/k^2.5, ;
% or a weighting-function (for the squared-L2-norm) of 1/k. ;
%%%%%%%%;
weight_3d_riesz_k_all_ = weight_3d_k_all_;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
weight_3d_k_p_r = weight_3d_k_p_r_(1+nk_p_r);
weight_2d_k_p_r = weight_2d_k_p_r_(1+nk_p_r);
tmp_index_ = n_k_all_csum_(1+nk_p_r):n_k_all_csum_(1+nk_p_r+1)-1;
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(weight_3d_k_all_(1+tmp_index_))/(weight_3d_k_p_r * 4*pi): %0.16f',k_p_r,sum(weight_3d_k_all_(1+tmp_index_))/(4*pi*weight_3d_k_p_r))); end;
weight_3d_riesz_k_all_(1+tmp_index_) = weight_3d_k_all_(1+tmp_index_) * weight_2d_k_p_r / max(1e-16,weight_3d_k_p_r);
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(weight_3d_riesz_k_all_(1+tmp_index_))/(weight_2d_k_p_r * 4*pi): %0.16f',k_p_r,sum(weight_3d_riesz_k_all_(1+tmp_index_))/(4*pi*weight_2d_k_p_r))); end;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;

%%%%%%%%;
% Now we performt the likelihood calculation. ;
%%%%%%%%;
L2_2d_order1_k_p = sum( bsxfun(@times,bsxfun(@times,abs(bsxfun(@minus,M_avg_k_p_,a_k_p_quad_)).^2,reshape(jl_hk_k_p_,[n_k_all,1])),weight_3d_riesz_k_all_) , 'all' );
L2_2d_order2_k_p = sum( bsxfun(@times,bsxfun(@times,abs(M_std_k_p_).^2,reshape(jl_hk_k_p_,[n_k_all,1])),weight_3d_riesz_k_all_) , 'all' );
L2_2d = (L2_2d_order1_k_p*n_M + L2_2d_order2_k_p);

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
