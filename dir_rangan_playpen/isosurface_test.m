[X,Y,Z] = meshgrid(linspace(-1,1,64)); 
R = sqrt(X.^2+Y.^2+Z.^2); A = atan2(Y,X); B = atan2(Z,sqrt(X.^2+Y.^2));
n_source = 32;
sigma = 1/32;
sigma_x_ = sigma*ones(1,n_source); sigma_y_ = sigma*ones(1,n_source); sigma_z_ = sigma*ones(1,n_source);
center_r_ = 0.5*(0:n_source-1)/max(1,n_source-1);
center_x_ = center_r_.*cos(2*pi*2*(0:n_source-1)/(n_source));
center_y_ = center_r_.*sin(2*pi*2*(0:n_source-1)/(n_source));
center_z_ = 2.0*center_r_ - 0.5 ;
tmp_gaus = 0 * (1/sqrt(pi)/sigma)^3 * exp(-R.^2/2/sigma.^2);
for nsigma=1:length(sigma_x_);
sigma_x = sigma_x_(nsigma); sigma_y = sigma_y_(nsigma); sigma_z = sigma_z_(nsigma);
center_x = center_x_(nsigma); center_y = center_y_(nsigma); center_z = center_z_(nsigma);
tmp_gaus = tmp_gaus + (1/sqrt(pi)/sigma_x)*(1/sqrt(pi)/sigma_y)*(1/sqrt(pi)/sigma_z) * exp(-(X-center_x).^2/2/sigma_x.^2 -(Y-center_y).^2/2/sigma_y.^2 -(Z-center_z).^2/2/sigma_z.^2);
end;%for nsigma=1:length(sigma_x_);
param_A = 8; param_B = 8;
%tmp_crop = 0.5*(1+erf(-16.0*(R - 0.25*(1 + 0.5*cos(param_A*A)) - 0.25*(1 + 0.5*cos(param_B*B)))));
tmp_crop = 0.5*(1+erf(-16.0*(R - 0.125*(1 + 0.5*cos(param_A*(A-1*pi/8))) - 0.125*(1 + 0.5*cos(param_B*B)))));
F_x_c = tmp_gaus ;
%F_x_c = tmp_crop ;
%F_x_c = tmp_crop .* tmp_gaus ;
F_k_c = recenter3(fftn(recenter3(F_x_c))); Fl_k_c = abs(F_k_c); Fw_k_c = angle(F_k_c);

subplot(1,2,1);
cra = colormap('jet'); ncra = size(cra,1);
v_avg = mean(F_x_c(:)); v_std = std(F_x_c(:)); v_min = min(F_x_c(:)); v_max = max(F_x_c(:));
vlim = v_avg + 2.5*v_std*[-1,1];
%vlim = [v_min , v_max];
v_ = linspace(vlim(1),vlim(2),5);
for nv=length(v_):-1:1
v = v_(nv);
hpatch = patch(isosurface(X,Y,Z,F_x_c,v)); isonormals(X,Y,Z,F_x_c,hpatch);
nc = max(1,min(ncra,floor(ncra*nv/length(v_))));
%hpatch.FaceColor = 'red'; 
hpatch.FaceColor = cra(nc,:);
hpatch.EdgeColor = 'none';
%hpatch.EdgeColor = cra(nc,:);
hpatch.FaceAlpha = (nv/length(v_)).^2 * 0.5;
%hpatch.EdgeAlpha = (nv/length(v_)).^2 * 1;
%title(sprintf('v %0.2f',v));
xlim(1*[-1,1]);ylim(1*[-1,1]);zlim(1*[-1,1]);
%daspect([1,4,4]); 
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for nv=1:length(v_);
title(sprintf('F(x)'));

subplot(1,2,2);
cra = colormap('jet'); ncra = size(cra,1);
v_avg = mean(Fl_k_c(:)); v_std = std(Fl_k_c(:)); v_min = min(Fl_k_c(:)); v_max = max(Fl_k_c(:));
v25 = prctile(Fl_k_c(:),25);
v50 = prctile(Fl_k_c(:),50);
v75 = prctile(Fl_k_c(:),75);
v85 = prctile(Fl_k_c(:),85);
v95 = prctile(Fl_k_c(:),95);
v97 = prctile(Fl_k_c(:),97.5);
%vlim = [v95,v97];
vlim = v_avg + 2.5*v_std*[-1,1];
%vlim = [v_min , v_max];
v_ = linspace(vlim(1),vlim(2),5);
for nv=length(v_):-1:1
v = v_(nv);
hpatch = patch(isosurface(X,Y,Z,Fl_k_c,v)); isonormals(X,Y,Z,Fl_k_c,hpatch);
nc = max(1,min(ncra,floor(ncra*nv/length(v_))));
%hpatch.FaceColor = 'red'; 
hpatch.FaceColor = cra(nc,:);
hpatch.EdgeColor = 'none';
%hpatch.EdgeColor = cra(nc,:);
hpatch.FaceAlpha = (nv/length(v_)).^2 * 0.5;
%hpatch.EdgeAlpha = (nv/length(v_)).^2 * 1;
%title(sprintf('v %0.2f',v));
xlim(1*[-1,1]);ylim(1*[-1,1]);zlim(1*[-1,1]);
%daspect([1,4,4]); 
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for nv=1:length(v_);
title(sprintf('(abs(F(k)))'));

%{
subplot(1,3,3);
cra = colormap('hsv'); ncra = size(cra,1);
vlim = [-pi,+pi];
v_ = linspace(vlim(1),vlim(2),2+1); v_ = v_(1:end-1);
for nv=length(v_):-1:1
v = v_(nv);
hpatch = patch(isosurface(X,Y,Z,Fw_k_c,v)); isonormals(X,Y,Z,Fw_k_c,hpatch);
nc = max(1,min(ncra,floor(ncra*nv/length(v_))));
%hpatch.FaceColor = 'red'; 
hpatch.FaceColor = cra(nc,:);
hpatch.EdgeColor = 'none';
%hpatch.EdgeColor = cra(nc,:);
hpatch.FaceAlpha = (nv/length(v_)).^2 * 0.5;
%hpatch.EdgeAlpha = (nv/length(v_)).^2 * 1;
%title(sprintf('v %0.2f',v));
xlim(1*[-1,1]);ylim(1*[-1,1]);zlim(1*[-1,1]);
%daspect([1,4,4]); 
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for nv=1:length(v_);
title(sprintf('abs(F(k))'));
 %}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
g_ = linspace(-1,1,125); g_ = g_(1:end-1);
[X,Y] = meshgrid(g_); 
R = sqrt(X.^2+Y.^2); A = atan2(Y,X);
%n_source = 12;
%sigma_x_ = 0.25*ones(1,n_source); sigma_y_ = 0.25*ones(1,n_source);
%center_x_ = 0.25*zeros(1,n_source); center_y_ = 0.25*zeros(1,n_source);
sigma = 0.0625;
tmp_gaus = (1/sqrt(pi)/sigma)^2 * exp(-R.^2/2/sigma.^2);
%for nsigma=1:length(sigma_x_);
%sigma_x = sigma_x_(nsigma); sigma_y = sigma_y_(nsigma);
%center_x = center_x_(nsigma); center_y = center_y_(nsigma);
%tmp_gaus = tmp_gaus + (1/sqrt(pi)/sigma_x)*(1/sqrt(pi)/sigma_y) * exp(-X.^2/2/sigma_x.^2 -Y.^2/2/sigma_y.^2);
%end;%for nsigma=1:length(sigma_x_);
param_A = 8;
%tmp_crop = 0.5*(1+erf(-16.0*(R - 0.25*(1 + 0.5*cos(param_A*A)))));
%F_x_c = tmp_crop .* tmp_gaus ;
F_x_c = tmp_gaus ;
F_k_c = recenter(fft2(recenter(F_x_c))); Fl_k_c = (max(1e-6,abs(F_k_c))); Fw_k_c = angle(F_k_c);
%F_k_c = fft2(F_x_c); Fl_k_c = log10(max(1e-6,abs(F_k_c))); Fw_k_c = angle(F_k_c);
subplot(1,2,1); imagesc(F_x_c);
subplot(1,2,2); imagesc(Fl_k_c);
return;
 %}
