function [ FA_x_c , X_1_ , X_2_ , X_3_ ]  = spharm_to_fig_3(n_k,k_,n_l_,a_,res,alpha_,delta_);
% plots the spherical harmonic representation passed in by modsph ;
% calculates real-space function directly, as well as via nufft ;
% ;
% n_k = integer maximum k ;
% k_ = integer array of length n_k; k_(nk) = k_value for shell nk ;
% n_l_ = integer array of length n_k; n_l_(nk) = spherical harmonic order on shell nk; n_l_(nk) corresponds to n_lm_(nk) = (n_l_(nk)+1)^2 coefficients ;
% a_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% res = resolution to plot ;
% alpha_ = real vector of euler-angles to apply for rotation [ alpha , beta , gamma ];
% delta_ = real vector of displacements to apply [delta_x , delta_y , delta_z ];
% test with: ;
%{
  spharm_to_fig_3();
  %}

path(path,'/data/rangan/dir_cryoem/nufftall-1.33/')

verbose=1;

if nargin<7; delta_ = 0*0.125*[1,2,3]; end;
if nargin<6; alpha_ = [0*pi/4,0*pi/4,0*pi/6]; end;
if nargin<5; res = 40; end;
if nargin<4;
isph_start_ = MDA_read_i4('./dir_mda6/isph_start_.mda');
nterms_sph_ = MDA_read_i4('./dir_mda6/nterms_sph_.mda');
modsph_A_ori_ = MDA_read_c16('./dir_mda6/modsph_A_ori_.mda');
n_k = length(isph_start_);
k_ = 1:n_k;
n_l_ = nterms_sph_;
a_ = modsph_A_ori_;
end;%if nargin<4;

n_lm_ = (n_l_+1).^2;

if (norm(alpha_)>0);
if (verbose); disp(sprintf(' %% rotating molecule by: [%0.2f %0.2f %0.2f]',+alpha_)); end;
if (verbose); disp(sprintf(' %% rotating coordinate_frame by: [%0.2f %0.2f %0.2f]',-alpha_(3),-alpha_(2),-alpha_(1))); end;
b_ = zeros(size(a_));
for nk=1:n_k;
n_l = n_l_(nk); n_lm = n_lm_(nk); ix_base = sum(n_lm_(1:nk-1));
a_k_ = a_(ix_base + (1:n_lm)); b_k_ = zeros(size(a_k_));
W_beta_ = wignerd_b(n_l,-alpha_(2));
for nl=0:n_l;
W_alpha = diag(exp(+i*[-nl:+nl]*-alpha_(3)));
W_gamma = diag(exp(+i*[-nl:+nl]*-alpha_(1)));
a_k_tmp = a_k_(1+nl*(nl+1) + (-nl:+nl));
a_k_tmp = reshape(a_k_tmp,2*nl+1,1);
b_k_tmp = W_gamma*W_beta_{1+nl}*W_alpha*a_k_tmp;
b_k_(1+nl*(nl+1) + (-nl:+nl)) = b_k_tmp;
end;%for nl=0:n_l;
b_(ix_base + (1:n_lm)) = b_k_;
end;%for nk=1:n_k;
a_ = b_;
end;%if (norm(alpha_)>0);

k_max = k_(n_k);
k_1_ = linspace(-k_max,+k_max,res); k_2_ = linspace(-k_max,+k_max,res); k_3_ = linspace(-k_max,+k_max,res);
[K_1_,K_2_,K_3_] = meshgrid(k_1_,k_2_,k_3_);

K_RAD_min = 1; K_RAD_max = n_k;
K_THETA_ = atan2(K_2_,K_1_);
K_RXY_ = max(1e-6,sqrt(K_2_.^2 + K_1_.^2));
%K_PHI_ = atan2(K_RXY_,K_3_);
K_RADIUS_ = min(K_RAD_max+1,max(1e-6,sqrt(K_3_.^2 + K_2_.^2 + K_1_.^2)));
K_PHI_ = acos(K_3_./K_RADIUS_);
K_RAD_PRE_ = min(K_RAD_max+1,max(K_RAD_min,floor(K_RADIUS_)));
K_RAD_POS_ = min(K_RAD_max+1,max(K_RAD_min,ceil(K_RADIUS_)));
DK_PRE_ = abs(K_RAD_PRE_ - K_RADIUS_);
DK_POS_ = abs(K_RAD_POS_ - K_RADIUS_);
DK_A_ = DK_PRE_;
DK_B_ = DK_POS_;
DK_C_ = DK_A_ + DK_B_;
DK_ALPHA_ = zeros(size(DK_A_)); DK_BETA_ = ones(size(DK_B_));
ij_tmp = find(DK_C_>0);
DK_ALPHA_(ij_tmp) = DK_A_(ij_tmp)./DK_C_(ij_tmp);
DK_BETA_(ij_tmp) = DK_B_(ij_tmp)./DK_C_(ij_tmp);
MOL_A_K_ = zeros(res,res,res);

for nk=1:n_k;
k = k_(nk); 
K_RAD_PRE_IJ = find(K_RAD_PRE_==k); K_RAD_POS_IJ = find(K_RAD_POS_==k);
length_pre = length(K_RAD_PRE_IJ); length_pos = length(K_RAD_POS_IJ);
K_PHI_PRE_ = K_PHI_(K_RAD_PRE_IJ); K_PHI_POS_ = K_PHI_(K_RAD_POS_IJ);
K_THETA_PRE_ = transpose(K_THETA_(K_RAD_PRE_IJ)); K_THETA_POS_ = transpose(K_THETA_(K_RAD_POS_IJ));
DK_BETA_PRE_ = transpose(DK_BETA_(K_RAD_PRE_IJ)); DK_BETA_POS_ = transpose(DK_BETA_(K_RAD_POS_IJ));
DK_ALPHA_PRE_ = transpose(DK_ALPHA_(K_RAD_PRE_IJ)); DK_ALPHA_POS_ = transpose(DK_ALPHA_(K_RAD_POS_IJ));
n_l = n_l_(nk); n_lm = n_lm_(nk); ix_base = sum(n_lm_(1:nk-1));
a_k_ = a_(ix_base + (1:n_lm)); 
l_ = []; m_ = []; for nl=0:n_l; l_ = [l_ , nl*ones(1,2*nl+1) ]; m_ = [m_ , [-nl:+nl] ]; end;%for nl=0:n_l;
if (verbose>0); disp(sprintf(' %% nk %d k %d/%d: length_pre %d length_pos %d; n_l %d n_lm %d ix_base %d',nk,k,k_max,length_pre,length_pos,n_l,n_lm,ix_base)); end;
for nl=0:n_l;
l_val = nl;
if (verbose>2); disp(sprintf(' %% nk %d k %d/%d: nl %d l_val %d',nk,k,k_max,nl,l_val)); end;
if (length_pre>0); Llm_pre__=legendre(l_val,cos(K_PHI_PRE_),'unnorm'); end;
if (length_pos>0); Llm_pos__=legendre(l_val,cos(K_PHI_POS_),'unnorm'); end;
A_pre_ = zeros(1,length_pre); A_pos_ = zeros(1,length_pos);
for m_val = -l_val:+l_val;
ix = 1+l_val*(l_val+1)+m_val;
m_abs = abs(m_val);
if (length_pre>0); if (l_val>0); Llm_pre_ = squeeze(Llm_pre__(1+m_abs,:,:)); end; if (l_val==0); Llm_pre_ = Llm_pre__(:,:); end; end;
if (length_pos>0); if (l_val>0); Llm_pos_ = squeeze(Llm_pos__(1+m_abs,:,:)); end; if (l_val==0); Llm_pos_ = Llm_pos__(:,:); end; end;
a1=((2*l_val+1)/(4*pi));
a2=exp(lfactorial(l_val-m_abs) - lfactorial(l_val+m_abs));
c=sqrt(a1*a2);
s=(-1)^((m_val<0)*m_val); % needed to preserve condon-shortley phase. ;
%s=1; % original phase ;
if (length_pre>0); 
Ylm_pre_ = s*c*Llm_pre_.*exp(+i*m_val*K_THETA_PRE_); 
A_pre_ = A_pre_ + s*a_k_(ix)*Ylm_pre_.*DK_BETA_PRE_;
end;%if (length_pre>0); 
if (length_pos>0); 
Ylm_pos_ = s*c*Llm_pos_.*exp(+i*m_val*K_THETA_POS_); 
A_pos_ = A_pos_ + s*a_k_(ix)*Ylm_pos_.*DK_ALPHA_POS_;
end;%if (length_pos>0); 
end;%for m_val = -l_val:+l_val;
if (length_pre>0); 
MOL_A_K_(K_RAD_PRE_IJ) = MOL_A_K_(K_RAD_PRE_IJ) + transpose(A_pre_);
end;%if (length_pre>0); 
if (length_pos>0);
MOL_A_K_(K_RAD_POS_IJ) = MOL_A_K_(K_RAD_POS_IJ) + transpose(A_pos_);
end;%if (length_pos>0);
end;%for nl=0:n_l;
end;%for nk=1:n_k;

if (norm(delta_)>0);
if (verbose); disp(sprintf(' %% translating by [%0.2f %0.2f %0.2f]',delta_)); end;
for nk1=1:length(k_1_);
k_1 = k_1_(nk1);
for nk2=1:length(k_2_);
k_2 = k_2_(nk2);
for nk3=1:length(k_3_);
k_3 = k_3_(nk3);
kd = k_1*delta_(1) + k_2*delta_(2) + k_3*delta_(3);
MOL_A_K_(nk1,nk2,nk3) = MOL_A_K_(nk1,nk2,nk3) * exp(i*kd);
end;end;end;% for nk3,nk2,nk1;
end;%if (norm(delta_)>0);

FA_k_c = MOL_A_K_; FA_l_k_c = abs(FA_k_c); FA_w_k_c = angle(FA_k_c);

plot_flag=1;
if plot_flag;
% plot molecule_A in fourier space;
subplot(2,2,1);
cra = colormap('spring'); cra = cra(end:-1:1,:);  ncra = size(cra,1);
F_tmp = abs(FA_k_c);
v_avg = mean(F_tmp(:)); v_std = std(F_tmp(:)); v_min = min(F_tmp(:)); v_max = max(F_tmp(:));
vlim = v_avg + 2.5*v_std*[-1,1];
v_ = linspace(vlim(1),vlim(2),7);
for nv=length(v_):-1:1;
v = v_(nv);
hpatch = patch(isosurface(K_1_,K_2_,K_3_,F_tmp,v)); isonormals(K_1_,K_2_,K_3_,F_tmp,hpatch);
nc = max(1,min(ncra,floor(ncra*nv/length(v_))));
hpatch.FaceColor = cra(nc,:);
hpatch.EdgeColor = 'none';
hpatch.FaceAlpha = (nv/length(v_)).^4 * 1.0;
xlim(k_max*[-1,1]);ylim(k_max*[-1,1]);zlim(k_max*[-1,1]);
xlabel('x');ylabel('y');zlabel('z');
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for nv=1:length(v_);
title(sprintf('FA_k, r_max %d res %d',k_max,res));
% end plot;
end;%if plot_flag;

sample_rate = res / k_max;
k_max = k_(n_k);
k_1_p_ = linspace(-sample_rate*k_max,+sample_rate*k_max,sample_rate*res); k_2_p_ = linspace(-sample_rate*k_max,+sample_rate*k_max,sample_rate*res); k_3_p_ = linspace(-sample_rate*k_max,+sample_rate*k_max,sample_rate*res);
[K_1_p_,K_2_p_,K_3_p_] = meshgrid(k_1_p_,k_2_p_,k_3_p_);
ij_tmp = (1:res) - floor(res/2) + floor(sample_rate*res/2);
FA_k_c_p = zeros(sample_rate*res,sample_rate*res,sample_rate*res);
FA_k_c_p(ij_tmp,ij_tmp,ij_tmp) = FA_k_c;

x_max = +1.0;
x_1_ = linspace(-x_max,+x_max,res); x_2_ = linspace(-x_max,+x_max,res); x_3_ = linspace(-x_max,+x_max,res);
[X_1_,X_2_,X_3_] = meshgrid(x_1_,x_2_,x_3_);

FA_x_c = real(recenter3(fftn(recenter3(FA_k_c_p)))); FA_x_c = FA_x_c(ij_tmp,ij_tmp,ij_tmp);

plot_flag=1;
if plot_flag;
% plot molecule_A in real space ;
subplot(2,2,2);
cra = colormap('spring'); cra = cra(end:-1:1,:);  ncra = size(cra,1);
F_tmp = FA_x_c;
v_avg = mean(F_tmp(:)); v_std = std(F_tmp(:)); v_min = min(F_tmp(:)); v_max = max(F_tmp(:));
vlim = v_avg + 2.5*v_std*[-1,1];
v_ = linspace(vlim(1),vlim(2),7);
for nv=length(v_):-1:1;
v = v_(nv);
hpatch = patch(isosurface(X_1_,X_2_,X_3_,FA_x_c,v)); isonormals(X_1_,X_2_,X_3_,FA_x_c,hpatch);
nc = max(1,min(ncra,floor(ncra*nv/length(v_))));
hpatch.FaceColor = cra(nc,:);
hpatch.EdgeColor = 'none';
hpatch.FaceAlpha = (nv/length(v_)).^4 * 1.0;
xlim(1*[-1,1]);ylim(1*[-1,1]);zlim(1*[-1,1]);
xlabel('x');ylabel('y');zlabel('z');
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for nv=1:length(v_);
title(sprintf('FA_x, r_max %d',k_max));
% end plot;
end;%if plot_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

n_all_ = zeros(n_k,1);
for nk=1:n_k;
k = k_(nk); 
n_all_(nk) =  sample_shell(k,1);
end%for nk=1:n_k;
n_all = sum(n_all_);
k_all_ = zeros(n_all,1); theta_all_ = zeros(n_all,1); phi_all_ = zeros(n_all,1);
MOL_A_K_P_ = zeros(n_all,1);

n_sub = 0;
for nk=1:n_k;
k = k_(nk); 
[length_sub,K_THETA_,K_PHI_] = sample_shell(k,1);
n_l = n_l_(nk); n_lm = n_lm_(nk); ix_base = sum(n_lm_(1:nk-1));
a_k_ = a_(ix_base + (1:n_lm)); 
l_ = []; m_ = []; for nl=0:n_l; l_ = [l_ , nl*ones(1,2*nl+1) ]; m_ = [m_ , [-nl:+nl] ]; end;%for nl=0:n_l;
if (verbose>0); disp(sprintf(' %% nk %d k %d/%d: n_l %d n_lm %d ix_base %d',nk,k,k_max,n_l,n_lm,ix_base)); end;
for nl=0:n_l;
l_val = nl;
if (verbose>2); disp(sprintf(' %% nk %d k %d/%d: nl %d l_val %d',nk,k,k_max,nl,l_val)); end;
Llm__=legendre(l_val,cos(K_PHI_),'unnorm');
A_ = zeros(1,length_sub);
for m_val = -l_val:+l_val;
ix = 1+l_val*(l_val+1)+m_val;
m_abs = abs(m_val);
if (length_sub>0); if (l_val>0); Llm_ = squeeze(Llm__(1+m_abs,:,:)); end; if (l_val==0); Llm_ = Llm__(:,:); end; end;
a1=((2*l_val+1)/(4*pi));
a2=exp(lfactorial(l_val-m_abs) - lfactorial(l_val+m_abs));
c=sqrt(a1*a2);
s=(-1)^((m_val<0)*m_val); % needed to preserve condon-shortley phase. ;
%s=1; % original phase ;
if (length_sub>0); 
Ylm_ = s*c*Llm_.*exp(+i*m_val*transpose(K_THETA_)); 
A_ = A_ + s*a_k_(ix)*Ylm_;
end;%if (length_sub>0); 
end;%for m_val = -l_val:+l_val;
if (length_sub>0); 
MOL_A_K_P_(1 + n_sub + (0:length_sub-1)) = MOL_A_K_P_(1 + n_sub + (0:length_sub-1)) + transpose(A_);
end;%if (length_sub>0); 
end;%for nl=0:n_l;
k_all_(1 + n_sub + (0:length_sub-1)) = k*ones(length_sub,1);
theta_all_(1 + n_sub + (0:length_sub-1)) = K_THETA_;
phi_all_(1 + n_sub + (0:length_sub-1)) = K_PHI_;
n_sub = n_sub + length_sub;
end;%for nk=1:n_k;
kx_all_ = k_all_ .* cos(theta_all_) .* sin(phi_all_);
ky_all_ = k_all_ .* sin(theta_all_) .* sin(phi_all_);
kz_all_ = k_all_ .* cos(phi_all_);

if (norm(delta_)>0);
if (verbose); disp(sprintf(' %% translating by [%0.2f %0.2f %0.2f]',delta_)); end;
for nall=1:n_all;
kd = kx_all_(nall)*delta_(1) + ky_all_(nall)*delta_(2) + kz_all_(nall)*delta_(3);
MOL_A_K_P_(nall) = MOL_A_K_P_(nall) * exp(i*kd);
end;%for nall=1:n_all;
end;%if (norm(delta_)>0);

FA_k_p = MOL_A_K_P_;

ij_tmp = (1:res) - floor(res/2) + floor(sample_rate*res/2);
FA_x_c = real(recenter3(fftn(recenter3(FA_k_c_p)))); FA_x_c = FA_x_c(ij_tmp,ij_tmp,ij_tmp);
sample_rate = res / k_max; mx1 = sample_rate*res;
[FA_x_c_nuff1,ier] = nufft3d1(numel(FA_k_c_p),pi*K_1_p_(:)/res,pi*K_2_p_(:)/res,pi*K_3_p_(:)/res,FA_k_c_p(:),-1,1e-12,mx1,mx1,mx1);
FA_x_c_nuff1 = real(permute(reshape(FA_x_c_nuff1,mx1,mx1,mx1),[2,1,3]))*numel(FA_k_c_p); FA_x_c_nuff1 = FA_x_c_nuff1(ij_tmp,ij_tmp,ij_tmp);
kx_all_resampled_ = pi*kx_all_/k_max ;
ky_all_resampled_ = pi*ky_all_/k_max ;
kz_all_resampled_ = pi*kz_all_/k_max ;
n_m = res.^3;
m_x_ = linspace(-k_max/2,+k_max/2,res+1); m_x_ = m_x_(1:res);
m_y_ = linspace(-k_max/2,+k_max/2,res+1); m_y_ = m_y_(1:res);
m_z_ = linspace(-k_max/2,+k_max/2,res+1); m_z_ = m_z_(1:res);
[M_X_,M_Y_,M_Z_] = meshgrid(m_x_,m_y_,m_z_);
[FA_x_c_nuff2,ier] = nufft3d3(n_all,kx_all_resampled_,ky_all_resampled_,kz_all_resampled_,FA_k_p,-1,1e-12,n_m,M_X_(:),M_Y_(:),M_Z_(:));
FA_x_c_nuff2 = real(reshape(FA_x_c_nuff2,res,res,res))*n_all;
z_ = round(linspace(1/4*res,3/4*res,6));
for nz=1:length(z_);
z = z_(nz);
subplot(3,length(z_),nz+0*length(z_)); imagesc(real(squeeze(FA_x_c(:,:,z)))); colorbar;
subplot(3,length(z_),nz+1*length(z_)); imagesc(real(squeeze(FA_x_c_nuff1(:,:,z)))); colorbar;
subplot(3,length(z_),nz+2*length(z_)); imagesc(real(squeeze(FA_x_c_nuff2(:,:,z)))); colorbar;
end;%for nz=1:length(z_);

x_max = +1.0;
x_1_ = linspace(-x_max,+x_max,res); x_2_ = linspace(-x_max,+x_max,res); x_3_ = linspace(-x_max,+x_max,res);
[X_1_,X_2_,X_3_] = meshgrid(x_1_,x_2_,x_3_);

plot_flag=1;
if plot_flag;
% plot molecule_A in real space ;
subplot(2,2,3);
cra = colormap('spring'); cra = cra(end:-1:1,:);  ncra = size(cra,1);
F_tmp = FA_x_c_nuff1;
v_avg = mean(F_tmp(:)); v_std = std(F_tmp(:)); v_min = min(F_tmp(:)); v_max = max(F_tmp(:));
vlim = v_avg + 2.5*v_std*[-1,1];
v_ = linspace(vlim(1),vlim(2),7);
for nv=length(v_):-1:1;
v = v_(nv);
hpatch = patch(isosurface(X_1_,X_2_,X_3_,FA_x_c_nuff1,v)); isonormals(X_1_,X_2_,X_3_,FA_x_c_nuff1,hpatch);
nc = max(1,min(ncra,floor(ncra*nv/length(v_))));
hpatch.FaceColor = cra(nc,:);
hpatch.EdgeColor = 'none';
hpatch.FaceAlpha = (nv/length(v_)).^4 * 1.0;
xlim(1*[-1,1]);ylim(1*[-1,1]);zlim(1*[-1,1]);
xlabel('x');ylabel('y');zlabel('z');
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for nv=1:length(v_);
title(sprintf('FA_x nuff1, r_max %d',k_max));
% end plot;
end;%if plot_flag;

plot_flag=1;
if plot_flag;
% plot molecule_A in real space ;
subplot(2,2,4);
cra = colormap('spring'); cra = cra(end:-1:1,:);  ncra = size(cra,1);
F_tmp = FA_x_c_nuff2;
v_avg = mean(F_tmp(:)); v_std = std(F_tmp(:)); v_min = min(F_tmp(:)); v_max = max(F_tmp(:));
vlim = v_avg + 2.5*v_std*[-1,1];
v_ = linspace(vlim(1),vlim(2),7);
for nv=length(v_):-1:1;
v = v_(nv);
hpatch = patch(isosurface(X_1_,X_2_,X_3_,FA_x_c_nuff2,v)); isonormals(X_1_,X_2_,X_3_,FA_x_c_nuff2,hpatch);
nc = max(1,min(ncra,floor(ncra*nv/length(v_))));
hpatch.FaceColor = cra(nc,:);
hpatch.EdgeColor = 'none';
hpatch.FaceAlpha = (nv/length(v_)).^4 * 1.0;
xlim(1*[-1,1]);ylim(1*[-1,1]);zlim(1*[-1,1]);
xlabel('x');ylabel('y');zlabel('z');
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for nv=1:length(v_);
title(sprintf('FA_x nuff2, r_max %d',k_max));
% end plot;
end;%if plot_flag;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

function [n_all,theta_all_,phi_all_] = sample_shell(r,d) ;
% generate theta_ and phi_ arrays sampled on shell of radius r at equatorial_distance d ;
% test with: ;
%{
  r=3.0; d=1.0;
  [n_all,theta_all_,phi_all_] = sample_shell(r,d) ;
  disp(sprintf(' %% n_all %d %d %d',n_all,length(theta_all_),length(phi_all_)));
  x_ = cos(theta_all_).*sin(phi_all_);
  y_ = sin(theta_all_).*sin(phi_all_);
  z_ = cos(phi_all_);
  plot3(x_,y_,z_,'.');axis vis3d; 
  xlabel('x');
  ylabel('y');
  zlabel('z');
 %}
n_equator = 1+round(2*pi*r/d);
n_phi = 1+round(n_equator/2);
cphi_ = linspace(-1,1,n_phi+2); cphi_ = cphi_(2:end-1);
phi_ = acos(cphi_);
n_all = 0 ;
for nphi = 1:n_phi;
phi = phi_(nphi);
sphi = sin(phi);
n_theta = 1+round(2*pi*sphi*r/d);
n_all = n_all + n_theta;
end;%for nphi = 1:n_phi;
theta_all_ = zeros(n_all,1); phi_all_ = zeros(n_all,1); ix=0;
for nphi = 1:n_phi;
phi = phi_(nphi);
sphi = sin(phi);
n_theta = 1+round(2*pi*sphi*r/d);
theta_ = linspace(0,2*pi,n_theta+1); theta_ = theta_(1:end-1);
theta_all_(1+ix+(0:n_theta-1)) = transpose(theta_);
phi_all_(1+ix+(0:n_theta-1)) = phi*ones(n_theta,1);
ix = ix+n_theta;
end;%for nphi = 1:n_phi;
