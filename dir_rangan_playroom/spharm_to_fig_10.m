function [ FA_x_c_ , X_0_ , X_1_ , X_2_ ]  = spharm_to_fig_10(n_k_p_r,k_p_r_,l_max_,a_k_Y_,res,n_iso,angle_,delta_,transf_type);
% Plots the spherical harmonic representation passed in by modsph. ;
% Displays n_iso isosurfaces (default 3). ;
% Allows for angle and delta transformation (in either order). ;
% ;
% Calculates real-space function via nufft3d3. ;
% (Can be used to test functions such as convert_spharm_to_k_p_0.m and convert_k_p_to_spharm_0.m.) ;
% 
% ;
% n_k_p_r = integer maximum k_p_r ;
% k_p_r_ = integer array of length n_k_p_r; k_p_r_(nk_p_r) = k_p_r_value for shell nk_p_r ;
% l_max_ = integer array of length n_k_p_r; l_max_(nk_p_r) = spherical harmonic order on shell nk_p_r; l_max_(nk_p_r) corresponds to n_lm_(nk_p_r) = (l_max_(nk_p_r)+1)^2 coefficients ;
% a_k_Y_ = complex array of length \sum_{nk_p_r} (n_lm_(nk_p_r)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% res = resolution to plot ;
% n_iso = number of isosurfaces to plot (default 5). ;
% angle_ = real vector of euler-angles to apply for rotation [ alpha , beta , gamma ] (default [0,0,0]);
% delta_ = real vector of displacements to apply [delta_x , delta_y , delta_z ] (default (0,0,0));
% transf_type = either 'angle first' or 'delta first' for transformation order (default 'angle first'). ;
% test with: ;
%{
  spharm_to_fig_10();
  %}

path(path,'/data/rangan/dir_cryoem/dir_nufftall-1.33/')

verbose=0;

%verbose=0; transf_type = 'angle first'; delta_ = [0,0,0]; angle_ = [0,0,0]; n_iso = 5; res = 64;
if nargin<9; transf_type = 'angle first'; end;
if nargin<8; delta_ = 1*0.125*[1,2,3]; end;
if nargin<7; angle_ = 1*[1*pi/4,1*pi/4,1*pi/6]; end;
if nargin<6; n_iso = 5; end;
if nargin<5; res = 64; end;
if nargin<4;
isph_start_ = MDA_read_i4('./dir_mda6/isph_start_.mda');
nterms_sph_ = MDA_read_i4('./dir_mda6/nterms_sph_.mda');
modsph_A_ori_ = MDA_read_c16('./dir_mda6/modsph_A_ori_.mda');
n_k_p_r = length(isph_start_);
k_p_r_ = 1:n_k_p_r;
l_max_ = nterms_sph_;
a_k_Y_ = modsph_A_ori_;
end;%if nargin<4;

n_lm_ = (l_max_+1).^2;
sample_d = 0.5;

if strcmp(transf_type,'angle first');
a_k_Y_ = rotate_spharm_to_spharm_0(verbose,n_k_p_r,k_p_r_,l_max_,a_k_Y_,angle_);
% slower, since this involves two unnecessary spherical-harmonic-transforms. ;
%[a_k_Y_] = transf_spharm_to_spharm_0(verbose,n_k_p_r,k_p_r_,l_max_,a_k_Y_,delta_,sample_d);
%[n_all,n_sub_,k_p_r_all_,azimu_b_all_,polar_a_all_,weight_all_,a_all_,kx_all_,ky_all_,kz_all_] = convert_spharm_to_k_p_0(verbose,n_k_p_r,k_p_r_,l_max_,a_k_Y_,sample_d); 
% faster, since we only transform to k_p_r_p once. ;
[n_all,n_sub_,k_p_r_all_,azimu_b_all_,polar_a_all_,weight_all_,a_all_,kx_all_,ky_all_,kz_all_] = convert_spharm_to_k_p_0(verbose,n_k_p_r,k_p_r_,l_max_,a_k_Y_,sample_d); 
a_all_ = transf_k_p_to_k_p_0(verbose,n_all,k_p_r_all_,azimu_b_all_,polar_a_all_,a_all_,delta_);
end;%if strcmp(transf_type,'angle first'); end;

if strcmp(transf_type,'delta first');
% not terribly efficient. ;
a_k_Y_ = transf_spharm_to_spharm_0(verbose,n_k_p_r,k_p_r_,l_max_,a_k_Y_,delta_,sample_d);
a_k_Y_ = rotate_spharm_to_spharm_0(verbose,n_k_p_r,k_p_r_,l_max_,a_k_Y_,angle_);
[n_all,n_sub_,k_p_r_all_,azimu_b_all_,polar_a_all_,weight_all_,a_all_,kx_all_,ky_all_,kz_all_] = convert_spharm_to_k_p_0(verbose,n_k_p_r,k_p_r_,l_max_,a_k_Y_,sample_d); 
end;%if strcmp(transf_type,'delta first');

k_p_r_max = k_p_r_(n_k_p_r);
kx_all_resampled_ = pi*kx_all_/k_p_r_max ;
ky_all_resampled_ = pi*ky_all_/k_p_r_max ;
kz_all_resampled_ = pi*kz_all_/k_p_r_max ;
n_m = res.^3;
m_x_ = linspace(-k_p_r_max/2,+k_p_r_max/2,res+1); m_x_ = m_x_(1:res);
m_y_ = linspace(-k_p_r_max/2,+k_p_r_max/2,res+1); m_y_ = m_y_(1:res);
m_z_ = linspace(-k_p_r_max/2,+k_p_r_max/2,res+1); m_z_ = m_z_(1:res);
[M_X_,M_Y_,M_Z_] = meshgrid(m_x_,m_y_,m_z_);
[FA_x_c_nu3d3_,ier] = nufft3d3(n_all,kx_all_resampled_,ky_all_resampled_,kz_all_resampled_,a_all_.*weight_all_,-1,1e-12,n_m,M_X_(:),M_Y_(:),M_Z_(:));
FA_x_c_nu3d3_ = real(reshape(FA_x_c_nu3d3_,res,res,res));

x_max = +1.0;
x_0_ = linspace(-x_max,+x_max,res); x_1_ = linspace(-x_max,+x_max,res); x_2_ = linspace(-x_max,+x_max,res);
[X_0_,X_1_,X_2_] = meshgrid(x_0_,x_1_,x_2_);
weight_k_p_r_ = ones(n_k_p_r,1);
FB_x_c_direct_ = real(convert_spharm_to_x_c_0(verbose,n_k_p_r,k_p_r_/(2*pi),ones(n_k_p_r,1),l_max_,a_k_Y_,X_0_,X_1_,X_2_));

plot_flag=1;
if plot_flag;
% plot molecule_A in real space ;
subplot(1,2,1);
cra = colormap('spring'); cra = cra(end:-1:1,:);  ncra = size(cra,1);
F_tmp_ = FA_x_c_nu3d3_;
v_avg = mean(F_tmp_(:)); v_std = std(F_tmp_(:)); v_min = min(F_tmp_(:)); v_max = max(F_tmp_(:));
vlim = v_avg + 2.5*v_std*[-1,1];
if (n_iso>1); v_ = linspace(vlim(1),vlim(2),n_iso); end;
if (n_iso==1); v_ = [v_avg]; end;
disp(num2str(v_));
for niso=n_iso:-1:1;
v = v_(niso);
hpatch = patch(isosurface(X_0_,X_1_,X_2_,FA_x_c_nu3d3_,v)); isonormals(X_0_,X_1_,X_2_,FA_x_c_nu3d3_,hpatch);
nc = max(1,min(ncra,floor(ncra*niso/n_iso)));
hpatch.FaceColor = cra(nc,:);
hpatch.EdgeColor = 'none';
hpatch.FaceAlpha = (niso/n_iso).^4 * 1.0;
xlim(1*[-1,1]);ylim(1*[-1,1]);zlim(1*[-1,1]);
xlabel('x');ylabel('y');zlabel('z');
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for niso=1:n_iso;
title(sprintf('FA_x nu3d3, r_max %d',k_p_r_max));
%%%%%%%%;
subplot(1,2,2);
%cra = colormap('spring'); cra = cra(end:-1:1,:);  ncra = size(cra,1);

cra = colormap(colormap_beach()); ncra = size(cra,1);
F_tmp_ = FB_x_c_direct_;
v_avg = mean(F_tmp_(:)); v_std = std(F_tmp_(:)); v_min = min(F_tmp_(:)); v_max = max(F_tmp_(:));
vlim = v_avg + 2.5*v_std*[-1,1];
if (n_iso>1); v_ = linspace(vlim(1),vlim(2),n_iso); end;
if (n_iso==1); v_ = [v_avg]; end;
disp(num2str(v_));
for niso=n_iso:-1:1;
v = v_(niso);
hpatch = patch(isosurface(X_0_,X_1_,X_2_,FB_x_c_direct_,v)); isonormals(X_0_,X_1_,X_2_,FB_x_c_direct_,hpatch);
nc = max(1,min(ncra,floor(ncra*niso/n_iso)));
hpatch.FaceColor = cra(nc,:);
hpatch.EdgeColor = 'none';
hpatch.FaceAlpha = (niso/n_iso).^4 * 1.0;
xlim(1*[-1,1]);ylim(1*[-1,1]);zlim(1*[-1,1]);
xlabel('x');ylabel('y');zlabel('z');
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for niso=1:n_iso;
title(sprintf('FB_x direct, r_max %d',k_p_r_max));
figbig;

end;%if plot_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
