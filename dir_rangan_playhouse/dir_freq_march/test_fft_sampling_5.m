%%%%%%%%;
% generating inital volume as sum of gaussians. ;
%%%%%%%%;
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
%%%%%%%%;
flag_plot=0;
if flag_plot;
hold on;
plot3(source_x_,source_y_,source_z_,'o'); 
plot3(source_r_.*cos(source_azimu_b_).*sin(source_polar_a_),source_r_.*sin(source_azimu_b_).*sin(source_polar_a_),source_r_.*cos(source_polar_a_),'.'); 
hold off;
axis vis3d;
disp('returning'); return;
end;%if flag_plot;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now creating fourier values on shell of radius k_r.; 
%%%%%%%%;
sigma = 0.1;
k_r = 4;
Fk = (1/sqrt(2*pi).^3) / sigma.^3 * sigma.^3 * exp(-0.5 * k_r.^2 ./ (1 ./ sigma).^2 );
%n_polar_a = 32;
n_polar_a = 16;
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
na = 0;
for npolar_a=1:n_polar_a;
polar_a_k = polar_a_k_(npolar_a);
for nazimu_b=1:n_azimu_b_(npolar_a);
azimu_b_k = azimu_b_k__{npolar_a}(nazimu_b);
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
%%%%%%%%;

%%%%%%%%;
% calculate correlations. ;
%%%%%%%%;
equator_F_q__ = fft(equator_F__(:,1:n_equator),[],2);
A_ = zeros(n_all,1);
for na1=1:n_all;
A_(na1) = sum(conj(equator_F__(na1,1:n_equator)).*equator_F__(na1,1:n_equator)) * (2*pi) / n_equator;
end;%for na1=1:n_all;
X__ = zeros(n_all,n_all);
for na1=1:n_all;
for na2=na1:n_all;
tmp_x = max(real(innerproduct_q(equator_F_q__(na1,:),equator_F_q__(na2,:))))/sqrt((A_(na1)*A_(na2)));
X__(na1,na2) = tmp_x;
X__(na2,na1) = tmp_x;
end;%for na2=na1:n_all;
end;%for na1=1:n_all;
imagesc(X__,[0,1]); colormap(colormap_beach()); axis image; colorbar;

%%%%%%%%;
flag_plot=1;
if flag_plot;
c_hsv_ = colormap('hsv');
%%%%%%%%;
subplot(2,2,3);
hold on;
na = 0;
for npolar_a=1:n_polar_a;
polar_a_k = polar_a_k_(npolar_a);
for nazimu_b=1:n_azimu_b_(npolar_a);
azimu_b_k = azimu_b_k__{npolar_a}(nazimu_b);
k_x = k_r*cos(azimu_b_k)*sin(polar_a_k);
k_y = k_r*sin(azimu_b_k)*sin(polar_a_k);
k_z = k_r*cos(polar_a_k);
tmp_x_ = exp(-(1-X__(1+na,:)).^2/0.05); tmp_x_ = tmp_x_/sum(tmp_x_);
tmp_k_ = tmp_x_*k__;
%plot3(k_x,k_y,k_z,'.','Color',colorsphere_trio(polar_a_k,azimu_b_k,c_hsv_),'MarkerSize',25);
plot3(tmp_k_(1),tmp_k_(2),tmp_k_(3),'.','Color',colorsphere_tetra(polar_a_k,azimu_b_k),'MarkerSize',25);
na = na+1;
end;%for nazimu_b=1:n_azimu_b_(npolar_a);
end;%for npolar_a=1:n_polar_a;
hold off;
axis vis3d; view(10,25);
title('viewing angles distorted by similarity profile');
end;%if flag_plot

%%%%%%%%;
flag_plot=1;
if flag_plot;
c_hsv_ = colormap('hsv');
%%%%%%%%;
subplot(2,2,1);
hold on;
na = 0;
for npolar_a=1:n_polar_a;
polar_a_k = polar_a_k_(npolar_a);
for nazimu_b=1:n_azimu_b_(npolar_a);
azimu_b_k = azimu_b_k__{npolar_a}(nazimu_b);
k_x = k_r*cos(azimu_b_k)*sin(polar_a_k);
k_y = k_r*sin(azimu_b_k)*sin(polar_a_k);
k_z = k_r*cos(polar_a_k);
%plot3(k_x,k_y,k_z,'.','Color',colorsphere_trio(polar_a_k,azimu_b_k,c_hsv_),'MarkerSize',25);
plot3(k_x,k_y,k_z,'.','Color',colorsphere_tetra(polar_a_k,azimu_b_k),'MarkerSize',25);
na = na+1;
end;%for nazimu_b=1:n_azimu_b_(npolar_a);
end;%for npolar_a=1:n_polar_a;
hold off;
axis vis3d; view(10,25);
title('viewing angle color code');
end;%if flag_plot

%%%%%%%%;
flag_plot=1;
if flag_plot;
c_hsv_ = colormap('hsv');
subplot(2,2,2);
hold on;
na = 0;
for npolar_a=1:n_polar_a;
polar_a_k = polar_a_k_(npolar_a);
for nazimu_b=1:n_azimu_b_(npolar_a);
azimu_b_k = azimu_b_k__{npolar_a}(nazimu_b);
%tmp_c_ = colorsphere_trio(polar_a_k,azimu_b_k,c_hsv_);
tmp_c_ = colorsphere_tetra(polar_a_k,azimu_b_k);
plot(real(equator_F__(1+na,:)),imag(equator_F__(1+na,:)),'-','Color',tmp_c_,'LineWidth',0.5);
na = na+1;
end;%for nazimu_b=1:n_azimu_b_(npolar_a);
end;%for npolar_a=1:n_polar_a;
hold off;
title('color by viewing angle (normal to great circle)');
end;%if flag_plot;

%%%%%%%%;
flag_plot=1;
if flag_plot;
c_hsv_ = colormap('hsv');
subplot(2,2,4);
hold on;
na = 0;
for npolar_a=1:n_polar_a;
polar_a_k = polar_a_k_(npolar_a);
for nazimu_b=1:n_azimu_b_(npolar_a);
azimu_b_k = azimu_b_k__{npolar_a}(nazimu_b);
%tmp_c_ = colorsphere_trio(equator_a__(1+na,:),equator_b__(1+na,:),c_hsv_);
tmp_c_ = colorsphere_tetra(equator_a__(1+na,:),equator_b__(1+na,:));
tmp_c_ = real(reshape(tmp_c_,1,n_equator+1,3));
tmp_x_ = real(equator_F__(1+na,:));
tmp_y_ = imag(equator_F__(1+na,:));
surface([tmp_x_;tmp_x_],[tmp_y_;tmp_y_],zeros(2,n_equator+1),[tmp_c_;tmp_c_],'facecol','no','edgecol','interp','linew',0.5);
na = na+1;
end;%for nazimu_b=1:n_azimu_b_(npolar_a);
end;%for npolar_a=1:n_polar_a;
hold off;
title('color by point position on sphere');
end;%if flag_plot;

if flag_plot;
set(gcf,'Position',1+[0,0,1024*2,1024]);
end;%if flag_plot;

fname_tmp = sprintf('%s/test_fft_sampling_5_k%d',pwd,round(k_r));
disp(sprintf(' %% writing %s.jpg',fname_tmp));
print('-djpeg',sprintf('%s.jpg',fname_tmp));

