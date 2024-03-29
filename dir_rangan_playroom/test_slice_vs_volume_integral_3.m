global_parameter=struct('type','parameter');
global_parameter.tolerance_master = 1e-6;
global_parameter.flag_verbose = 1;
flag_verbose = global_parameter.flag_verbose;
str_thisfunction = 'test_slice_vs_volume_integral_3';

if (flag_verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
% Define paths. ;
% You will have to write your own (short) 'platform.type' file.
% Then you may have to write your own 'setup_local.m' and define your own string_root. ;
%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

%%%%%%%%;
% Here are default values for the 'global_parameter'. ;
%%%%%%%%;
if isempty(global_parameter); global_parameter = struct('type','parameter'); end;
if (~isfield(global_parameter,'tolerance_master')); global_parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_verbose')); global_parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_recalc')); global_parameter.flag_recalc = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_replot')); global_parameter.flag_replot = 0; end; %<-- parameter_bookmark. ;
tolerance_master = global_parameter.tolerance_master;
flag_verbose = global_parameter.flag_verbose;
flag_recalc = global_parameter.flag_recalc;
flag_replot = global_parameter.flag_replot;
nf=0;

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
n_x_u = 256; %<-- this is a typical resolution for a volume loaded from the protein database. ;
half_diameter_x_c = 1.0d0; %<-- assume the volume sits in a box of side-length 2. ;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
n_x_u_pack = 64; %<-- this is a lower resolution for a downsampled volume we will work with. ;
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
% a_x_u_pack_ = reshape(a_x_u_load_,[n_x_u*n_x_u,n_x_u])*x_u_pack_;
% a_x_u_pack_ = reshape(permute(reshape(a_x_u_pack_,[n_x_u,n_x_u,n_x_u_pack]),[3,1,2]),[n_x_u*n_x_u_pack,n_x_u])*x_u_pack_;
% a_x_u_pack_ = reshape(permute(reshape(a_x_u_pack_,[n_x_u_pack,n_x_u,n_x_u_pack]),[3,1,2]),[n_x_u_pack*n_x_u_pack,n_x_u])*x_u_pack_;
% a_x_u_pack_ = permute(reshape(a_x_u_pack_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),[3,1,2]);
%%%%%%%%;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack); %<-- downsampled grid for x_u_0_. ;
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack); %<-- downsampled grid for x_u_1_. ;
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack); %<-- downsampled grid for x_u_2_. ;
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;
xxx_u_weight_ = (2*x_p_r_max/n_x_u_pack)^3;
%%%%%%%%;
% Here are grids for k_c_ (cartesian), sometimes referred to as k_u_ (uniform). ;
%%%%%%%%;
k_u_0_ = periodize(0:n_x_u_pack-1,-n_x_u_pack/2,+n_x_u_pack/2)/2; %<-- box has diameter 2. ;
k_u_1_ = periodize(0:n_x_u_pack-1,-n_x_u_pack/2,+n_x_u_pack/2)/2; %<-- box has diameter 2. ;
k_u_2_ = periodize(0:n_x_u_pack-1,-n_x_u_pack/2,+n_x_u_pack/2)/2; %<-- box has diameter 2. ;
[k_u_0___,k_u_1___,k_u_2___] = ndgrid(k_u_0_,k_u_1_,k_u_2_); n_kkk_u = n_x_u_pack^3;
%%%%%%%%;
flag_disp=1;
if flag_disp;
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
end;%if flag_disp;
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
k_p_r_max = 48.0/(2*pi); %<-- k_p_r_max is the maximum frequency. ;
% Note that we technically use the 'wavenumber' rather than the 'frequency' (related by 2*pi), so the exponentials involved in our calculations will look like exp(+i*2*pi*k_p_r*x_p_r). ;
k_eq_d = 1.0/(2*pi); %<-- k_eq_d is the distance (along the equator) between gridpoints. ;
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
,J_node_ ...
,J_weight_ ...
,J_chebfun_ ...
,J_polyval_ ...
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
tmp_t = toc(tmp_t); disp(sprintf(' %% sample_sphere_7: time %0.2fs',tmp_t));
%%%%%%%%;
flag_disp=1;
if flag_disp;
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
% Now define a random volume a_k_Y_base_yk_. ;
%%%%%%%%;
rng(0);
a_k_Y_base_yk_ = ...
 ( randn(n_lm_sum,1) + i*randn(n_lm_sum,1) ) ...
.* exp(-(Y_l_val_ + abs(Y_m_val_)).^2./(2*(Y_l_max_val_/4).^2)) ...
.* exp(-(Y_k_val_-k_p_r_max/2).^2./(2*(k_p_r_max/8)).^2) ...
;
flag_check=0;
if flag_check;
%%%%%%%%;
% Check to see if the conversion between k_p_ and k_Y_ is accurate. ;
% Note that this uses convert_spharm_to_k_p_1 and convert_k_p_to_spharm_1, which are quite slow. ;
% This is because these functions recalculate the associated Ylm__ for each nk_p_r. ;
% There are faster versions (see convert_spharm_to_k_p_4, etc.) which precompute the necessary Ylm__. ;
%%%%%%%%;
tmp_t = tic;
[a_k_p_quad_] = convert_spharm_to_k_p_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_base_yk_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_base_yk_ --> a_k_p_quad_ time %0.2fs',tmp_t));
tmp_t = tic;
[a_k_Y_quad_yk_] = convert_k_p_to_spharm_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_p_quad --> a_k_Y_quad_yk_ time %0.2fs',tmp_t));
disp(sprintf(' %% a_k_Y_base_yk_ vs a_k_Y_quad_yk_: %0.16f',fnorm(a_k_Y_base_yk_-a_k_Y_quad_yk_)/fnorm(a_k_Y_quad_yk_))); %<-- this should be at least 2-3 digits. ;
end;%if flag_check;
%%%%%%%%;
% Because the 'recovered' a_k_Y_quad_yk_ is so close to the original a_k_Y_base_yk_, ;
% we simply set the two to be equal for now. ;
%%%%%%%%;
a_k_Y_quad_yk_ = a_k_Y_base_yk_; %<-- approximation a_k_Y_quad_yk_ is quite close to a_k_Y_base_yk_. ;
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
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t));
S_k_p_wkS__ = reshape(S_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
n_S = n_viewing_all; n_w_max = max(n_w_);
%%%%%%%%;
flag_check=0;
if flag_check;
%%%%%%%%;
% Test pm_template_2 vs get_template_1. ;
% The latter allows for adaptive grids n_w_, but is significantly slower. ;
%%%%%%%%;
tmp_t = tic();
[ ...
 tmp_S_k_p_wkS__ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,template_k_c_0_wkS__ ...
,template_k_c_1_wkS__ ...
,template_k_c_2_wkS__ ...
] = ...
get_template_1( ...
 0*flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_quad_yk_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_max*ones(n_k_p_r,1) ...
);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% get_template_1: %0.2fs',tmp_t)); end;
if (flag_verbose); disp(sprintf(' %% n_viewing_all %d n_viewing_polar_a %d n_w_max %d',n_viewing_all,n_viewing_polar_a,max(n_w_))); end;
disp(sprintf(' %% S_k_p_wkS__ vs tmp_S_k_p_wkS__: %0.16f',fnorm(S_k_p_wkS__-tmp_S_k_p_wkS__)/fnorm(S_k_p_wkS__)));
end;%if flag_check;
%%%%%%%%;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% View one of the templates. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
nS=0; %<-- template index to view. ;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
S_x_c_ = ...
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
subplot(2,2,1); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wk_),[],colormap_80s); title(sprintf('nS %d real(S_k_p_wk_)',nS),'Interpreter','none'); axis image; axisnotick;
subplot(2,2,2); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(S_k_p_wk_),[],colormap_80s); title(sprintf('nS %d imag(S_k_p_wk_)',nS),'Interpreter','none'); axis image; axisnotick;
subplot(2,2,3); imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,real(S_x_c_),[],colormap_beach); title(sprintf('nS %d real(S_x_c_)',nS),'Interpreter','none'); axis image; axisnotick;
subplot(2,2,4); imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,imag(S_x_c_),[],colormap_beach); title(sprintf('nS %d imag(S_x_c_)',nS),'Interpreter','none'); axis image; axisnotick;
end;%if flag_disp;

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

E_3d_k_Y = sum(weight_Y_);
E_3d_k_p = sum(weight_3d_k_all_);
E_2d_k_p = sum(weight_2d_k_p_wkS__, 'all' );
E_2d_k_Y = sum(weight_Y_2d_);

L2_3d_k_Y = sum(abs(a_k_Y_base_yk_).^2 .* weight_Y_);
flag_check=0;
if flag_check;
L2_3d_k_p = sum(abs(a_k_p_quad_).^2 .* weight_3d_k_all_);
end;%if flag_check;

L2_2d_k_p = sum(abs(S_k_p_wkS__).^2 .* weight_2d_k_p_wkS__,'all');
L2_2d_k_Y = sum(abs(a_k_Y_base_yk_).^2 .* weight_Y_2d_);
if (flag_verbose); disp(sprintf(' %% L2_2d_k_p %0.6f vs L2_2d_k_Y %0.6f: %0.16f',L2_2d_k_p,L2_2d_k_Y,fnorm(L2_2d_k_p - L2_2d_k_Y)/fnorm(L2_2d_k_Y))); end;

%%%%%%%%;
% Now we add several sets of images. ;
% Typically n_M refers to the number of images, and nM is the image index. ;
% Temporarily (i.e., in the code below) we will use n_M to refer to the number of image-sets, ;
% with nM referring to the index of the image-set. ;
%%%%;
% For this particular script, we are interested in checking the equivalence of image- vs volumetric likelihood calculations, ;
% and so we construct the images in such a way that we immediately know the back-propagated volume (1st-order) and volumetric-variance (2nd-order) ;
% associated with the images. ;
% To do this, we actually generate the images in 'sets' or batches, where each image-set corresponds to a volume ;
% which itself generates a full set of templates (i.e., one for each viewing-angle). ;
% Thus, each image-set will have n_S images in it. ;
% Note also that, due to this particular construction, each image-set ;
% will have a distribution of viewing-angles that is close to uniform ;
% (because the templates were generated with a uniform distribution of viewing-angles). ;
%%%%%%%%;
% Thus, if we have n_M image-sets, each corresponding to a particular volume, ;
% the back-propagated volume constructed from the full image-pool will simply be the set-wise average volume (i.e., 1st-order). ;
% Moreover, the back-propagated volumetric variance will simply be the set-wise variance of those volumes (i.e., 2nd-order). ;
% We structure the image pool this way so that there is a very small error associated with back-propagation. ;
% We use n_M> 1 so that the variance is nonzero. ;
%%%%%%%%;
n_M = 5; %<-- number of individual 'independent' image-sets. ;
b_k_Y_quad_ykM__ = zeros(n_lm_sum,n_M); %<-- volume in k_Y_ coordinates for each image-set. ;
b_k_Y_quad_ykM___ = zeros(n_lm_max,n_k_p_r,n_M);
for nM=0:n_M-1;
rng(1+nM);
b_k_Y_quad_ykM__(:,1+nM) = ...
 ( randn(n_lm_sum,1) + i*randn(n_lm_sum,1) ) ...
.* exp(-(Y_l_val_ + abs(Y_m_val_)).^2./(2*(Y_l_max_val_/4).^2)) ...
.* exp(-(Y_k_val_-k_p_r_max/2).^2./(2*(k_p_r_max/8)).^2) ...
;
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
b_k_Y_quad_ykM___(1:n_lm_(1+nk_p_r),1+nk_p_r,1+nM) = b_k_Y_quad_ykM__(1+tmp_index_,1+nM);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nM=0:n_M-1;
%%%%%%%%;
b_avg_k_Y_quad_yk_ = mean(b_k_Y_quad_ykM__,2); %<-- average volume (averaged over image-sets). ;
b_avg_k_Y_quad_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
b_avg_k_Y_quad_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = b_avg_k_Y_quad_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
% generate templates M_k_p_wkSM___ on k_p_ grid with uniform n_w_. ;
% Note that each image-set has n_S images in it (with image-index mirroring the template-index). ;
% Thus, image nS from image-set nM is stored in M_k_p_wkSM___(:,1+nS,1+nM). ;
%%%%%%%%;
M_k_p_wkSM___ = zeros(n_w_sum,n_S,n_M);
for nM=0:n_M-1;
if (flag_verbose); disp(sprintf(' %% nM %d/%d',nM,n_M)); end;
tmp_t = tic();
[ ...
 M_k_p_wkS__ ...
] = ...
pm_template_2( ...
 flag_verbose ...
,l_max ...
,n_k_p_r ...
,reshape(b_k_Y_quad_ykM___(:,:,1+nM),[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
);
M_k_p_wkS__ = reshape(M_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
M_k_p_wkSM___(:,:,1+nM) = M_k_p_wkS__; clear M_k_p_wkS__;
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t));
end;%for nM=0:n_M-1;
%%%%%%%%;
% We also construct the set-averaged image-set M_avg_k_p_wkS__. ;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% M_avg')); end;
tmp_t = tic();
[ ...
 M_avg_k_p_wkS__ ...
] = ...
pm_template_2( ...
 flag_verbose ...
,l_max ...
,n_k_p_r ...
,reshape(b_avg_k_Y_quad_yk__,[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
);
M_avg_k_p_wkS__ = reshape(M_avg_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
%%%%%%%%;
% Note that this set-averaged image-set could also have been generated ;
% by taking the set-wise average of M_k_p_wkSM___. ;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% fnorm(M_avg_k_p_wkS__ - mean(M_k_p_wkSM___,3)): %0.16f',fnorm(M_avg_k_p_wkS__ - mean(M_k_p_wkSM___,3)))); end;
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t));

%%%%%%%%;
% Now we check the volumetric formulation, assuming a uniform viewing-angle-distribution. ;
% This compares formul 'image-wise calculation' with 'volumetric-calculation' in the section entitled: ;
% '[simple volumetric formulation with uniformly distributed viewing-angles]' in the overleaf. ;
% The 'set-volume' associated with each image-set is b_k_Y_quad_ykM__(:,1+nM), ;
% while the 'set-avg-volume' is b_avg_k_Y_quad_yk_. ;
% The 'true-volume' associated with the templates is a_k_Y_base_yk_. ;
%%%%%%%%;
% Note that, because of the way this calculation is structured, ;
% we have very little 'error' in the volume-formation for each set. ;
% The volumes are known ahead of time, ;
% and the images and templates are produced from those volumes. ;
% Later on we should repeat this calculation while reconstructing the volume from the images ;
% (as one might do when dealing with real data). ;
% We actually do something along these lines (i.e., reconstructing the 2nd-order term from noisy images) below. ;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% Assuming uniform distribution of viewing-angles: ')); end;
if (flag_verbose); disp(sprintf(' %% (multi volume to ensure accurate variance calculation): ')); end;
%%%%%%%%;
% In each calculation below 'L2' refers to a squared-L2-norm. ;
% L2_2d_alM_k_p contrasts the images and the templates in (2d) k_p_ coordinates. ;
%%%%%%%%;
L2_2d_alM_k_p = sum( bsxfun(@times,abs(bsxfun(@minus,M_k_p_wkSM___,S_k_p_wkS__)).^2,weight_2d_k_p_wkS__) , 'all' ); %<-- image-wise calculation. ;
%%%%%%%%;
% L2_2d_alM_k_Y contrasts the set-volumes against the true-volume. ;
%%%%%%%%;
L2_2d_alM_k_Y = sum( bsxfun(@times,abs(bsxfun(@minus,b_k_Y_quad_ykM__,a_k_Y_base_yk_)).^2,weight_Y_2d_) , 'all' ); %<-- image-wise calculation using spherical harmonics. This kind of calculation is not typically accessible with arbitrary images. ;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% L2_2d_alM_k_p %0.6f vs L2_2d_alM_k_Y %0.6f: %0.16f',L2_2d_alM_k_p,L2_2d_alM_k_Y,fnorm(L2_2d_alM_k_p - L2_2d_alM_k_Y)/fnorm(L2_2d_alM_k_Y))); end;
%%%%%%%%;
% L2_2d_avg_k_p contrasts the set-averaged images against the templates. ;
%%%%%%%%;
L2_2d_avg_k_p = sum( bsxfun(@times,abs(bsxfun(@minus,M_avg_k_p_wkS__,S_k_p_wkS__)).^2,weight_2d_k_p_wkS__) , 'all' );
%%%%%%%%;
% L2_2d_avg_k_Y contrasts the set-avg-volume against the true-volume. ;
%%%%%%%%;
L2_2d_avg_k_Y = sum( bsxfun(@times,abs(bsxfun(@minus,b_avg_k_Y_quad_yk_,a_k_Y_base_yk_)).^2,weight_Y_2d_) , 'all' );
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% L2_2d_avg_k_p %0.6f vs L2_2d_avg_k_Y %0.6f: %0.16f',L2_2d_avg_k_p,L2_2d_avg_k_Y,fnorm(L2_2d_avg_k_p - L2_2d_avg_k_Y)/fnorm(L2_2d_avg_k_Y))); end;
%%%%%%%%;
% L2_2d_var_k_p calculates the second-order volume-term from the images. ;
%%%%%%%%;
L2_2d_var_k_p = sum( bsxfun(@times,abs(bsxfun(@minus,M_k_p_wkSM___,M_avg_k_p_wkS__)).^2,weight_2d_k_p_wkS__) , 'all' );
%%%%%%%%;
% L2_2d_var_k_Y calculates the second-order volume-term from the set-volumes. ;
%%%%%%%%;
L2_2d_var_k_Y = sum( bsxfun(@times,abs(bsxfun(@minus,b_k_Y_quad_ykM__,b_avg_k_Y_quad_yk_)).^2,weight_Y_2d_) , 'all' );
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% L2_2d_var_k_p %0.6f vs L2_2d_var_k_Y %0.6f: %0.16f',L2_2d_var_k_p,L2_2d_var_k_Y,fnorm(L2_2d_var_k_p - L2_2d_var_k_Y)/fnorm(L2_2d_var_k_Y))); end;
%%%%%%%%;
% Now we compare: ;
% L2_2d_alM_k_Y (i.e., the log-likelihood calculated using images) to ;
% (L2_2d_avg_k_Y*n_M + L2_2d_var_k_Y) (i.e., the the log-likelihood calculated using volumes). ;
% This is the 'volumetric-calculation' in the overleaf. ;
%%%%%%%%;
disp(sprintf(' %% L2_2d_alM_k_Y vs (L2_2d_avg_k_Y*n_M + L2_2d_var_k_Y): %0.16f',fnorm(L2_2d_alM_k_Y-(L2_2d_avg_k_Y*n_M + L2_2d_var_k_Y))/fnorm(L2_2d_alM_k_Y)));
%%%%%%%%;

%%%%%%%%;
% Checking uniform distribution again, but this time using only a single image-set. ;
% In this calculation we will construct the second-order volume-term from the images. ;
%%%%%%%%;
c_k_Y_quad_yk_ = b_k_Y_quad_ykM__(:,1); %<-- pick one of the set-volumes. ;
N_k_p_wkS__ = M_k_p_wkSM___(:,:,1); %<-- and the associated collection of n_S images. ;
n_N = 7; %<-- here n_N is the multiplier used when creating our 'full' image-set N_k_p_wkSN___. ;
N_k_p_wkSN___ = repmat(N_k_p_wkS__,[1,1,n_N]); %<-- Initialize the full image-set. ;
rng(0);
Ndif_wkSN___ = 0.5*randn(size(N_k_p_wkSN___)); %<-- Noise added to the images. ;
Ndif_wkSN___ = bsxfun(@minus,Ndif_wkSN___,mean(Ndif_wkSN___,3)); %<-- Noise is mean-0. ;
N_k_p_wkSN___ = N_k_p_wkSN___ + Ndif_wkSN___;
N_avg_k_p_wkS__ = mean(N_k_p_wkSN___,3); %<-- average-image for each viewing-angle. ;
%%%%%%%%;
% here nN refers to the image-index in our expanded image-set of size n_S*n_N. ;
%%%%%%%%;
euler_polar_a_N_ = reshape(repmat(viewing_polar_a_all_,[1,n_N]),[n_S*n_N,1]);
euler_azimu_b_N_ = reshape(repmat(viewing_azimu_b_all_,[1,n_N]),[n_S*n_N,1]);
euler_gamma_z_N_ = zeros(n_S*n_N,1);
%%%%%%%%;
% Now we reconstruct the second-order volume-terms. ;
% Because of the way qbp_6 is structured, ;
% these calculations benefit from the (nearly) uniform-distribution of viewing-angles. ;
% A highly nonuniform-distribution of viewing-angles will introduce errors (due to qbp_6), ;
% in which case one should use a different algorithm for reconstruction. ;
%%%%%%%%;
qbp_eps = 1e-3;
[ ...
 N_std_k_Y_yk_ ...
] = ...
qbp_6( ...
 qbp_eps ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_N*n_S ...
,reshape(sqrt(abs(Ndif_wkSN___).^2),[n_w_sum,n_S*n_N]) ...
,[] ...
,[] ...
,euler_polar_a_N_ ...
,euler_azimu_b_N_ ...
,euler_gamma_z_N_ ...
);
%%%%%%%%;
qbp_eps = 1e-3;
[ ...
 N_var_k_Y_yk_ ...
] = ...
qbp_6( ...
 qbp_eps ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_N*n_S ...
,reshape(abs(Ndif_wkSN___).^2,[n_w_sum,n_S*n_N]) ...
,[] ...
,[] ...
,euler_polar_a_N_ ...
,euler_azimu_b_N_ ...
,euler_gamma_z_N_ ...
);
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% Assuming uniform distribution of viewing-angles (single volume): ')); end;
%%%%;
L2_2d_alM_k_p = sum( bsxfun(@times,abs(bsxfun(@minus,N_k_p_wkSN___,S_k_p_wkS__)).^2,weight_2d_k_p_wkS__) , 'all' );
%%%%;
L2_2d_avg_k_p = sum( bsxfun(@times,abs(bsxfun(@minus,N_avg_k_p_wkS__,S_k_p_wkS__)).^2,weight_2d_k_p_wkS__) , 'all' );
L2_2d_avg_k_Y = sum( bsxfun(@times,abs(bsxfun(@minus,c_k_Y_quad_yk_,a_k_Y_base_yk_)).^2,weight_Y_2d_) , 'all' );
if (flag_verbose); disp(sprintf(' %% L2_2d_avg_k_p %0.6f vs L2_2d_avg_k_Y %0.6f: %0.16f',L2_2d_avg_k_p,L2_2d_avg_k_Y,fnorm(L2_2d_avg_k_p - L2_2d_avg_k_Y)/fnorm(L2_2d_avg_k_Y))); end;
%%%%;
L2_2d_var_k_p = sum( bsxfun(@times,abs(bsxfun(@minus,N_k_p_wkSN___,N_avg_k_p_wkS__)).^2,weight_2d_k_p_wkS__) , 'all' );
L2_2d_var_k_Y = sum( bsxfun(@times,abs(N_std_k_Y_yk_).^2,weight_Y_2d_*(n_N*pi/2)) , 'all' );
%L2_2d_var_k_Y = sum( bsxfun(@times,abs(N_var_k_Y_yk_).^1,weight_Y_2d_*n_N) , 'all' );
if (flag_verbose); disp(sprintf(' %% L2_2d_var_k_p %0.6f vs L2_2d_var_k_Y %0.6f: %0.16f',L2_2d_var_k_p,L2_2d_var_k_Y,fnorm(L2_2d_var_k_p - L2_2d_var_k_Y)/fnorm(L2_2d_var_k_Y))); end;
%%%%;
disp(sprintf(' %% L2_2d_alM_k_p vs (L2_2d_avg_k_p*n_M + L2_2d_var_k_p): %0.16f',fnorm(L2_2d_alM_k_p-(L2_2d_avg_k_p*n_N + L2_2d_var_k_p))/fnorm(L2_2d_alM_k_p)));
%%%%;

%%%%%%%%;
% Now we will check the image-wise and volumetric-calculations in the section: ;
% 'More general volume formulation' in the overleaf. ;
% We set up a nonuniform distribution of viewing-angles. ;
% This is defined using the function JL_AB_, with values stored in jl_ab_all_. ;
% jl_ab_all_ is a double-array of size n_S. ;
%         jl_ab_all_(1+nS) stores the probability of observing the viewing-angle ;
%         associated with [ viewing_azimu_b_(1+nS) , viewing_polar_a_(1+nS) ]. ;
%%%%%%%%;
% m_val_ is a double-array of size (1+2*l_max) storing the array of spherical-harmonic orders. ;
% Note that l_max_ is often chosen adaptively, so not all orders will be used on all shells. ;
%%%%%%%%;
m_val_ = -l_max:+l_max;
%%%%%%%%;
JL_AB_const_a = 3; JL_AB_const_b = 4;
JL_AB_ = @(a,b) (JL_AB_const_b + cos(b)).*(JL_AB_const_a + sin(a)) / (JL_AB_const_b * 2*pi * (JL_AB_const_a*2 + pi/2 ));
jl_ab_all_ = JL_AB_(viewing_polar_a_all_,viewing_azimu_b_all_);
if (flag_verbose); disp(sprintf(' %% sum(jl_ab_all_.*viewing_weight_all_)/(k_p_r_max^2): %0.16f (should be 1)',sum(jl_ab_all_.*viewing_weight_all_)/(k_p_r_max^2))); end;
%%%%%%%%;
% j1_ab_all_ stores the uniform distribution. ;
%%%%%%%%;
j1_ab_all_ = ones(n_S,1)/(4*pi);
%%%%%%%%;
% Now construct shell-complement of viewing-angle distribution. ;
% jl_hk_Y_ is a double-complex-array of size n_lm_max. ;
%         jl_hk_Y_ stores the volumetric-description of the function ;
%         \tilde{J}(\hk) from the overleaf notes. ;
%         Briefly, the function \tilde{J}(\hk) refers to the probability of observing
%         a particular viewing-angle on the sphere (indicated by \hk). ;
%         Note that \tilde{J}(\hk) is defined on a spherical-shell. ;
% jl_hk_lm__ is a double-complex-array of size(1+l_max,numel(m_val_)). ;
%         jl_hk_lm__ stores the values of jl_hk_Y_ in a rectangular (2d) array. ;
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
% j1_hk_Y_ and j1_hk_lm__ are the analogous arrays for the uniform distribution. ;
%%%%%%%%;
qbp_eps = 1e-3;
[ ...
 j1_hk_Y_ ...
,n_quad_from_data_q_ ...
] = ...
qbp_6( ...
 qbp_eps ...
,1 ...
,1 ...
,l_max ...
,n_w_max ...
,n_S ...
,ones(n_w_max,1)*transpose(j1_ab_all_) ...
,[] ...
,[] ...
,viewing_polar_a_all_ ...
,viewing_azimu_b_all_ ...
,zeros(n_S,1) ...
,[] ...
,[] ...
,[]...
);
j1_hk_lm__=zeros(1+l_max,n_m_max);
for l_val=0:l_max;
index_0in_ = l_val^2 + [0:1+2*l_val-1];
index_out_ = l_max + [-l_val:+l_val];
j1_hk_lm__(1+l_val,1+index_out_) = j1_hk_Y_(1+index_0in_);
end;%for l_val=0:l_max;
if (flag_verbose); disp(sprintf(' %% j1_hk_lm__(1,1+l_max)*sqrt(4*pi): %0.16f (should be 1)',j1_hk_lm__(1,1+l_max)*sqrt(4*pi))); end;
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
% j1_hk_k_all_ is the analogous array for the uniform-distribution. ;
%%%%%%%%;
j1_hk_k_all_ = ...
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
,j1_hk_Y_ ...
);
if (flag_verbose); disp(sprintf(' %% sum(real(j1_hk_k_all_).*viewing_weight_all_/k_p_r_max^2): %0.16f (should be 1)',sum(real(j1_hk_k_all_).*viewing_weight_all_/k_p_r_max^2))); end;
%%%%%%%%;
% jl_hk_yk_ and j1_hk_yk_ are double-arrays of size n_lm_sum, ;
% referencing the above distributions on the spherical-harmonic grid (including all shells). ;
%%%%%%%%;
jl_hk_yk_ = zeros(n_lm_sum,1);
for nlm_sum=0:n_lm_sum-1;
jl_hk_yk_(1+nlm_sum) = jl_hk_lm__(1+Y_l_val_(1+nlm_sum),1+l_max+Y_m_val_(1+nlm_sum));
end;%for nlm_sum=0:n_lm_sum-1;
j1_hk_yk_ = zeros(n_lm_sum,1);
for nlm_sum=0:n_lm_sum-1;
j1_hk_yk_(1+nlm_sum) = j1_hk_lm__(1+Y_l_val_(1+nlm_sum),1+l_max+Y_m_val_(1+nlm_sum));
end;%for nlm_sum=0:n_lm_sum-1;
%%%%%%%%;
flag_disp=1;
if flag_disp;
%%%%%%%%;
% Demonstrate the nonuniform viewing-angle distribution (JL). ;
% Compare with the uniform  viewing-angle distribution (J1). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
imagesc_polar_a_azimu_b_0(viewing_polar_a_all_,viewing_azimu_b_all_,real(jl_ab_all_),[0,0.1],[],0);
xlabel('k0'); ylabel('k1'); zlabel('k2');
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
title('jl_ab_all_','Interpreter','none');
subplot(1,2,2);
imagesc_polar_a_azimu_b_0(viewing_polar_a_all_,viewing_azimu_b_all_,real(j1_ab_all_),[0,0.1],[],0);
xlabel('k0'); ylabel('k1'); zlabel('k2');
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
title('j1_ab_all_','Interpreter','none');
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% Now we construct representations of the distribution jl and j1 in (3d) k_p_ coordinates. ;
% We also contruct representations of the various set-volumes in (3d) k_p_ coordinates. ;
% Since the code convert_spharm_to_k_p_2 is slow (see more recent versions), ;
% we store the results. ;
% Note that loading jl_hk_k_all_ and j1_hk_k_all_ (which are arrays of size n_lm_sum) ;
% will overwrite the above definitions (which were arrays of size n_S). ;
%%%%%%%%;
dir_base = sprintf('/%s/rangan/dir_cryoem/dir_rangan_playroom',string_root);
dir_pm_mat = sprintf('%s/dir_pm_mat',dir_base);
fname_mat = sprintf('%s/test_slice_vs_volume_integral_k_p_1.mat',dir_pm_mat);
if (~exist(fname_mat,'file'));
if (flag_verbose); disp(sprintf(' %% %s not found, creating',fname_mat)); end;
%%%%%%%%;
% Now convert all k_Y to k_p. ;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% jl_hk_k_all_')); end;
[jl_hk_k_all_] = ...
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
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% j1_hk_k_all_')); end;
[j1_hk_k_all_] = ...
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
,j1_hk_yk_ ...
);
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% a_k_p_quad_k_all_')); end;
[a_k_p_base_k_all_] = ...
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
,a_k_Y_base_yk_ ...
);
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% c_k_p_quad_k_all_')); end;
[c_k_p_quad_k_all_] = ...
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
,b_k_Y_quad_ykM__(:,1) ...
);
%%%%%%%%;
b_k_p_quad_k_all_kM__ = zeros(n_k_all,n_M);
for nM=0:n_M-1;
if (flag_verbose); disp(sprintf(' %% nM %d/%d b_k_p_quad_k_all_kM__',nM,n_M)); end;
[b_k_p_quad_k_all_kM__(:,1+nM)] = ...
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
,b_k_Y_quad_ykM__(:,1+nM) ...
);
end;%for nM=0:n_M-1;
save(fname_mat ...
,'n_k_all' ...
,'n_k_all_csum_' ...
,'k_p_r_all_' ...
,'k_p_azimu_b_all_' ...
,'k_p_polar_a_all_' ...
,'weight_3d_k_all_' ...
,'weight_shell_k_' ...
,'n_k_p_r' ...
,'k_p_r_' ...
,'weight_3d_k_p_r_' ...
,'l_max_' ...
,'jl_hk_k_all_' ...
,'j1_hk_k_all_' ...
,'a_k_p_base_k_all_' ...
,'c_k_p_quad_k_all_' ...
,'b_k_p_quad_k_all_kM__' ...
);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
if (flag_verbose); disp(sprintf(' %% %s found, not creating',fname_mat)); end;
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

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
% First we constrast the image-calculation and the volumetric-calculation ;
% when assuming a (nearly) uniform-distribution of viewing-angles. ;
%%%%%%%%;
L2_2d_MvS_k_p = sum( bsxfun(@times,abs(bsxfun(@minus,M_k_p_wkSM___,S_k_p_wkS__)).^2,weight_2d_k_p_wkS__) , 'all' );
L2_2d_bva_k_p = sum( bsxfun(@times,abs(bsxfun(@minus,b_k_p_quad_k_all_kM__,a_k_p_base_k_all_)).^2,weight_3d_riesz_k_all_) , 'all' );
if (flag_verbose); disp(sprintf(' %% L2_2d_MvS_k_p %0.6f vs L2_2d_bva_k_p %0.6f: %0.16f',L2_2d_MvS_k_p,L2_2d_bva_k_p,fnorm(L2_2d_MvS_k_p - L2_2d_bva_k_p)/fnorm(L2_2d_MvS_k_p))); end;
%%%%%%%%;
% Now we repeat these calculations explicitly using the j1 values (i.e. uniform viewing-angle). ;
%%%%%%%%;
L2_2d_MvS_k_p = sum( bsxfun(@times,bsxfun(@times,abs(bsxfun(@minus,M_k_p_wkSM___,S_k_p_wkS__)).^2,reshape(j1_ab_all_,[1,n_S,1])),weight_2d_k_p_wkS__) , 'all' );
L2_2d_bva_k_p = sum( bsxfun(@times,bsxfun(@times,abs(bsxfun(@minus,b_k_p_quad_k_all_kM__,a_k_p_base_k_all_)).^2,reshape(j1_hk_k_all_,[n_k_all,1])),weight_3d_riesz_k_all_) , 'all' );
if (flag_verbose); disp(sprintf(' %% (uniform weight) L2_2d_MvS_k_p %0.6f vs L2_2d_bva_k_p %0.6f: %0.16f',L2_2d_MvS_k_p,L2_2d_bva_k_p,fnorm(L2_2d_MvS_k_p - L2_2d_bva_k_p)/fnorm(L2_2d_MvS_k_p))); end;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% Assuming non-uniform distribution of viewing-angles (multi volume): ')); end;
if (flag_verbose); disp(sprintf(' %% (multi volume to ensure accurate variance calculation): ')); end;
%%%%%%%%;
% Now we repeat the calculations again, this time using the non-uniform jl values. ;
% Note that this calculation requires both a first-order volumetric-term (L2_2d_bva_avg_k_p) ;
% as well as a second-order volumetric-term (L2_2d_bva_var_k_p). ;
%%%%%%%%;
L2_2d_MvS_k_p = sum( bsxfun(@times,bsxfun(@times,abs(bsxfun(@minus,M_k_p_wkSM___,S_k_p_wkS__)).^2,reshape(jl_ab_all_,[1,n_S,1])),weight_2d_k_p_wkS__) , 'all' );
L2_2d_bva_k_p = sum( bsxfun(@times,bsxfun(@times,abs(bsxfun(@minus,b_k_p_quad_k_all_kM__,a_k_p_base_k_all_)).^2,reshape(jl_hk_k_all_,[n_k_all,1])),weight_3d_riesz_k_all_) , 'all' );
if (flag_verbose); disp(sprintf(' %% (nonuniform weight) L2_2d_MvS_k_p %0.6f vs L2_2d_bva_k_p %0.6f: %0.16f (should be good to 2-3 digits)',L2_2d_MvS_k_p,L2_2d_bva_k_p,fnorm(L2_2d_MvS_k_p - L2_2d_bva_k_p)/fnorm(L2_2d_MvS_k_p))); end;
%%%%;
b_avg_k_p_quad_k_all_ = mean(b_k_p_quad_k_all_kM__,2);
L2_2d_MvS_avg_k_p = sum( bsxfun(@times,bsxfun(@times,abs(bsxfun(@minus,M_avg_k_p_wkS__,S_k_p_wkS__)).^2,reshape(jl_ab_all_,[1,n_S,1])),weight_2d_k_p_wkS__) , 'all' );
L2_2d_bva_avg_k_p = sum( bsxfun(@times,bsxfun(@times,abs(bsxfun(@minus,b_avg_k_p_quad_k_all_,a_k_p_base_k_all_)).^2,reshape(jl_hk_k_all_,[n_k_all,1])),weight_3d_riesz_k_all_) , 'all' );
if (flag_verbose); disp(sprintf(' %% (nonuniform weight) L2_2d_MvS_avg_k_p %0.6f vs L2_2d_bva_avg_k_p %0.6f: %0.16f (should be good to 2-3 digits)',L2_2d_MvS_avg_k_p,L2_2d_bva_avg_k_p,fnorm(L2_2d_MvS_avg_k_p - L2_2d_bva_avg_k_p)/fnorm(L2_2d_MvS_avg_k_p))); end;
%%%%;
L2_2d_MvS_var_k_p = sum( bsxfun(@times,bsxfun(@times,abs(bsxfun(@minus,M_k_p_wkSM___,M_avg_k_p_wkS__)).^2,reshape(jl_ab_all_,[1,n_S,1])),weight_2d_k_p_wkS__) , 'all' );
L2_2d_bva_var_k_p = sum( bsxfun(@times,bsxfun(@times,abs(bsxfun(@minus,b_k_p_quad_k_all_kM__,b_avg_k_p_quad_k_all_)).^2,reshape(jl_hk_k_all_,[n_k_all,1])),weight_3d_riesz_k_all_) , 'all' );
if (flag_verbose); disp(sprintf(' %% L2_2d_MvS_var_k_p %0.6f vs L2_2d_bva_var_k_p %0.6f: %0.16f (should be good to 2-3 digits)',L2_2d_MvS_var_k_p,L2_2d_bva_var_k_p,fnorm(L2_2d_MvS_var_k_p - L2_2d_bva_var_k_p)/fnorm(L2_2d_bva_var_k_p))); end;
%%%%%%%%;
% And finally we compare the image- and volumetric-calculations of the likelihood. ;
%%%%%%%%;
disp(sprintf(' %% L2_2d_MvS_k_p vs (L2_2d_bva_avg_k_p*n_M + L2_2d_bva_var_k_p): %0.16f (should be good to 2-3 digits)',fnorm(L2_2d_MvS_k_p-(L2_2d_bva_avg_k_p*n_M + L2_2d_bva_var_k_p))/fnorm(L2_2d_MvS_k_p)));
%%%%%%%%;
