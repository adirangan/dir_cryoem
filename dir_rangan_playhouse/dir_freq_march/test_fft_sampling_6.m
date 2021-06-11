function test_fft_sampling_6(k_r,n_polar_a,model_type);
%{

  model_type = 'helix'; n_polar_a = 32; for k_r = 8:4:48; test_fft_sampling_6(k_r,n_polar_a,model_type); end;
  model_type = 'lissajous'; n_polar_a = 32; for k_r = 8:4:48; test_fft_sampling_6(k_r,n_polar_a,model_type); end;
  %%%%%%%%;
  model_type = 'random'; n_polar_a = 32; for k_r = 8:4:48; test_fft_sampling_6(k_r,n_polar_a,model_type); end;
  model_type = 'lissajous'; n_polar_a = 32; for k_r = 8:4:48; test_fft_sampling_6(k_r,n_polar_a,model_type); end;
  model_type = 'helix'; n_polar_a = 32; for k_r = 8:4:48; test_fft_sampling_6(k_r,n_polar_a,model_type); end;

  %}

if nargin<1; k_r = 24; end;
if nargin<2; n_polar_a = 32; end;
if nargin<3; model_type = 'helix'; end;

%%%%%%%%;
% generating inital volume as sum of gaussians. ;
%%%%%%%%;
%%%%%%%%;
if (strcmp(model_type,'random'));
rng(1);
n_source = 32;
source_x_ = zeros(n_source,1); source_y_ = zeros(n_source,1); source_z_ = zeros(n_source,1);
source_r_ = zeros(n_source,1); source_polar_a_ = zeros(n_source,1); source_azimu_b_ = zeros(n_source,1);
for nsource=1:n_source;
tmp_r = 1;
tmp_w = 2*pi*nsource/n_source;
tmp_x = 0.75*(2*rand()-1);
tmp_y = 0.75*(2*rand()-1);
tmp_z = 0.75*(2*rand()-1);
source_x_(nsource) = tmp_x;
source_y_(nsource) = tmp_y;
source_z_(nsource) = tmp_z;
source_r_(nsource) = sqrt(tmp_x.^2 + tmp_y.^2 + tmp_z.^2); 
source_polar_a_(nsource) = atan2(sqrt(tmp_x.^2 + tmp_y.^2),tmp_z);
source_azimu_b_(nsource) = atan2(tmp_y,tmp_x);
end;%for nsource=1:n_source;
end;%if (strcmp(model_type,'random'));
%%%%%%%%;
if (strcmp(model_type,'helix'));
n_source = 64;
source_x_ = zeros(n_source,1); source_y_ = zeros(n_source,1); source_z_ = zeros(n_source,1);
source_r_ = zeros(n_source,1); source_polar_a_ = zeros(n_source,1); source_azimu_b_ = zeros(n_source,1);
for nsource=1:n_source;
tmp_r = nsource/n_source;
tmp_w = 2*pi*nsource/n_source;
tmp_x = 0.75*tmp_r*cos(3*tmp_w);
tmp_y = 0.75*tmp_r*sin(3*tmp_w);
tmp_z = 0.75*(2*tmp_r-1);
source_x_(nsource) = tmp_x;
source_y_(nsource) = tmp_y;
source_z_(nsource) = tmp_z;
source_r_(nsource) = sqrt(tmp_x.^2 + tmp_y.^2 + tmp_z.^2); 
source_polar_a_(nsource) = atan2(sqrt(tmp_x.^2 + tmp_y.^2),tmp_z);
source_azimu_b_(nsource) = atan2(tmp_y,tmp_x);
end;%for nsource=1:n_source;
end;%if (strcmp(model_type,'helix'));
%%%%%%%%;
if (strcmp(model_type,'lissajous'));
n_source = 128;
source_x_ = zeros(n_source,1); source_y_ = zeros(n_source,1); source_z_ = zeros(n_source,1);
source_r_ = zeros(n_source,1); source_polar_a_ = zeros(n_source,1); source_azimu_b_ = zeros(n_source,1);
for nsource=1:n_source;
tmp_r = 1;
tmp_w = 2*pi*nsource/n_source;
tmp_x = 0.75*tmp_r*cos(7*tmp_w);
tmp_y = 0.75*tmp_r*sin(5*tmp_w);
tmp_z = 0.75*tmp_r*cos(3*tmp_w);
source_x_(nsource) = tmp_x;
source_y_(nsource) = tmp_y;
source_z_(nsource) = tmp_z;
source_r_(nsource) = sqrt(tmp_x.^2 + tmp_y.^2 + tmp_z.^2); 
source_polar_a_(nsource) = atan2(sqrt(tmp_x.^2 + tmp_y.^2),tmp_z);
source_azimu_b_(nsource) = atan2(tmp_y,tmp_x);
end;%for nsource=1:n_source;
end;%if (strcmp(model_type,'lissajous'));
%%%%%%%%;
flag_plot=0;
if flag_plot;
hold on;
plot3(source_x_,source_y_,source_z_,'o'); 
plot3(source_r_.*cos(source_azimu_b_).*sin(source_polar_a_),source_r_.*sin(source_azimu_b_).*sin(source_polar_a_),source_r_.*cos(source_polar_a_),'.'); 
hold off;
axis vis3d; view(10,25);
set(gcf,'Position',1+[0,0,512,512]);
fname_tmp = sprintf('%s/test_fft_sampling_6_%s',pwd,model_type);
disp(sprintf(' %% writing %s.jpg',fname_tmp));
print('-djpeg',sprintf('%s.jpg',fname_tmp));
disp('returning'); return;
end;%if flag_plot;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now creating fourier values on shell of radius k_r.; 
%%%%%%%%;
sigma = 0.1;
Fk = (1/sqrt(2*pi).^3) / sigma.^3 * sigma.^3 * exp(-0.5 * k_r.^2 ./ (1 ./ sigma).^2 );
%n_polar_a = 16;
polar_a_k_ = linspace(0,pi,n_polar_a);
n_azimu_b_ = zeros(n_polar_a,1);
azimu_b_k__ = cell(n_polar_a,1);
n_equator = 128*2;
gamma_ = linspace(0,2*pi,n_equator+1);
n_all=0;
for npolar_a=1:n_polar_a;
n_azimu_b_(npolar_a) = max(1,round(n_polar_a*2*sin(polar_a_k_(npolar_a))));
tmp_azimu_b_ = linspace(0,2*pi,n_azimu_b_(npolar_a)+1); tmp_azimu_b_ = tmp_azimu_b_(1:end-1);
azimu_b_k__{npolar_a} = tmp_azimu_b_;
n_all = n_all + n_azimu_b_(npolar_a);
end;%for npolar_a=1:n_polar_a;
k__ = zeros(n_all,3);
F_ = zeros(n_all,1);
equator_F__ = zeros(n_all,n_equator+1);
equator_x__ = zeros(n_all,n_equator+1);
equator_y__ = zeros(n_all,n_equator+1);
equator_z__ = zeros(n_all,n_equator+1);
equator_a__ = zeros(n_all,n_equator+1);
equator_b__ = zeros(n_all,n_equator+1);
polar_a_v__ = zeros(n_all,1);
azimu_b_v__ = zeros(n_all,1);
na = 0;
for npolar_a=1:n_polar_a;
polar_a_k = polar_a_k_(npolar_a);
for nazimu_b=1:n_azimu_b_(npolar_a);
azimu_b_k = azimu_b_k__{npolar_a}(nazimu_b);
polar_a_v__(1+na) = polar_a_k;
azimu_b_v__(1+na) = azimu_b_k;
k_x = k_r*cos(azimu_b_k)*sin(polar_a_k);
k_y = k_r*sin(azimu_b_k)*sin(polar_a_k);
k_z = k_r*cos(polar_a_k);
k__(1+na,:) = [k_x,k_y,k_z];
tmp_w_ = [k_x,k_y,k_z]; tmp_w_ = tmp_w_/norm(tmp_w_);
tmp_u_ = cross(tmp_w_,[1,0,0]); tmp_u_ = tmp_u_/norm(tmp_u_);
tmp_v_ = cross(tmp_w_,tmp_u_); tmp_v_ = tmp_v_/norm(tmp_v_);
equator_x__(1+na,:) = k_r * (cos(gamma_)*tmp_u_(1) + sin(gamma_)*tmp_v_(1));
equator_y__(1+na,:) = k_r * (cos(gamma_)*tmp_u_(2) + sin(gamma_)*tmp_v_(2));
equator_z__(1+na,:) = k_r * (cos(gamma_)*tmp_u_(3) + sin(gamma_)*tmp_v_(3));
equator_a__(1+na,:) = real(acos(equator_z__(1+na,:)/k_r));
equator_b__(1+na,:) = real(atan2(equator_y__(1+na,:)/k_r,equator_x__(1+na,:)/k_r));
for nsource=1:n_source;
F_(1+na) = F_(1+na) + Fk * exp( i* ( k_x*source_x_(nsource) + k_y*source_y_(nsource) + k_z*source_z_(nsource) ) ) ;
equator_F__(1+na,:) = equator_F__(1+na,:) + Fk * exp( i* ( equator_x__(1+na,:)*source_x_(nsource) + equator_y__(1+na,:)*source_y_(nsource) + equator_z__(1+na,:)*source_z_(nsource) ) ) ;
end;%for nsource=1:n_source;
na = na+1;
end;%for nazimu_b=1:n_azimu_b_(npolar_a);
end;%for npolar_a=1:n_polar_a;
color_p_v__ = colorpsphere_tetra(polar_a_v__,azimu_b_v__);
color_s_v__ = colorsphere_tetra(polar_a_v__,azimu_b_v__);
%%%%%%%%;

%%%%%%%%;
% Calculate correlations X__ = <A1,A2>/|A1||A2|. ;
% as well as squared distances Y__ = |A1 - A2|^2 . 
%%%%%%%%;
equator_F_q__ = fft(equator_F__(:,1:n_equator),[],2);
AA_ = zeros(n_all,1);
for na1=1:n_all;
AA_(na1) = sum(conj(equator_F__(na1,1:n_equator)).*equator_F__(na1,1:n_equator)) * (2*pi) / n_equator;
end;%for na1=1:n_all;
X__ = zeros(n_all,n_all);
Y__ = zeros(n_all,n_all);
for na1=1:n_all;
for na2=na1:n_all;
tmp_x = max(real(innerproduct_q(equator_F_q__(na1,:),equator_F_q__(na2,:))));
Y__(na1,na2) = AA_(na1) - 2*tmp_x + AA_(na2);
Y__(na2,na1) = Y__(na1,na2);
X__(na1,na2) = tmp_x / sqrt((AA_(na1)*AA_(na2)));
X__(na2,na1) = X__(na1,na2);
end;%for na2=na1:n_all;
end;%for na1=1:n_all;
flag_plot=1;
if flag_plot;
clf;
imagesc(X__,[-1,1]); colormap(colormap_beach()); axis image; %colorbar;
figbig;
fname_tmp = sprintf('%s/test_fft_sampling_6_%s_X_k%da%d',pwd,model_type,round(k_r),n_polar_a);
disp(sprintf(' %% writing %s.jpg',fname_tmp));
print('-djpeg',sprintf('%s.jpg',fname_tmp));
end;%if flag_plot;

%%%%%%%%;
flag_plot=0;
if flag_plot;
%%%%%%%%;
sigma_ = 0.01:0.01:0.24;
for ns=1:length(sigma_);
sigma = sigma_(ns);
tmp_x__ = exp(-(1-X__).^2/(2*sigma^2));
tmp_x__ = tmp_x__./repmat(sum(tmp_x__,2),1,n_all);
tmp_k__ = tmp_x__*k__;
subplot(4,6,ns);
hold on; scatter3(tmp_k__(:,1),tmp_k__(:,2),tmp_k__(:,3),15,color_p_v__,'filled'); hold off;
axis vis3d; view(10,25);
title(sprintf('\\sigma %0.2f',sigma));
end;%for ns=1:length(sigma_);
%%%%%%%%;
set(gcf,'Position',1+[0,0,1024*2,1024]);
fname_tmp = sprintf('%s/test_fft_sampling_6_%s_scatterplot_k%da%d',pwd,model_type,round(k_r),n_polar_a);
disp(sprintf(' %% writing %s.jpg',fname_tmp));
print('-djpeg',sprintf('%s.jpg',fname_tmp));
end;%if flag_plot;
%%%%%%%%;

%%%%%%%%;
% plot mesh grid? ;
%%%%%%%%;
flag_plot=1;
if flag_plot;
clf;
n_K = 5; [IK_,ID_] = knnsearch(k__,k__,'K',n_K);
n_p = 2;%n_p = n_all*n_K*2;
tmp_x_ = zeros(1,n_p); tmp_y_ = zeros(1,n_p); tmp_z_ = zeros(1,n_p); tmp_c_ = zeros(1,n_p,3);
sigma_ = linspace(0.01,0.24,6);
for ns=1:length(sigma_);
sigma = sigma_(ns);
tmp_x__ = exp(-(1-X__).^2/(2*sigma^2));
tmp_x__ = tmp_x__./repmat(sum(tmp_x__,2),1,n_all);
tmp_k__ = tmp_x__*k__;
subplot(2,3,ns);
%%%%%%%%;
hold on;
for na=1:n_all;
for nk=1:n_K;
if (IK_(na,nk)~=na);
nak=0;
nak = nak+1;
nb = na;         tmp_x_(1,nak) = tmp_k__(nb,1); tmp_y_(1,nak) = tmp_k__(nb,2); tmp_z_(1,nak) = tmp_k__(nb,3); tmp_c_(1,nak,:) = color_p_v__(nb,:); nak = nak+1;
nb = IK_(na,nk); tmp_x_(1,nak) = tmp_k__(nb,1); tmp_y_(1,nak) = tmp_k__(nb,2); tmp_z_(1,nak) = tmp_k__(nb,3); tmp_c_(1,nak,:) = color_p_v__(nb,:); nak = nak+1;
surface([tmp_x_;tmp_x_],[tmp_y_;tmp_y_],[tmp_z_;tmp_z_],[tmp_c_;tmp_c_],'facecol','no','edgecol','interp','linew',2.0);
end;%if (IK_(na,nk)~=na);
end;%for nk=1:n_K;
end;%for na=1:n_all;
hold off;
%%%%%%%%;
axis vis3d; view(10,25);
title(sprintf('\\sigma %0.2f',sigma));
end;%for ns=1:length(sigma_);
figbig;
fname_tmp = sprintf('%s/test_fft_sampling_6_%s_mesh_k%da%d',pwd,model_type,round(k_r),n_polar_a);
disp(sprintf(' %% writing %s.jpg',fname_tmp));
print('-djpeg',sprintf('%s.jpg',fname_tmp));
end;%if flag_plot;
