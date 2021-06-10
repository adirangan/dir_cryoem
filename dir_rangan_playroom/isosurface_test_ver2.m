figure;

[X_m,Y_m,Z_m] = meshgrid(linspace(-1,1,64)); 
R_m = sqrt(X_m.^2+Y_m.^2+Z_m.^2); A_m = atan2(Y_m,X_m); B = atan2(Z_m,sqrt(X_m.^2+Y_m.^2));
n_source = 32;
sigma = 1/32;
sigma_x_ = sigma*ones(1,n_source); sigma_y_ = sigma*ones(1,n_source); sigma_z_ = sigma*ones(1,n_source);

center_A_r_ = 0.5*(0:n_source-1)/max(1,n_source-1);
center_A_w_ = 2*pi*2*(0:n_source-1)/(n_source);
center_A_x_ = center_A_r_.*cos(center_A_w_);
center_A_y_ = center_A_r_.*sin(center_A_w_);
center_A_z_ = 2.0*center_A_r_ - 0.5 ;
tmp_gaus = 0 * R_m;
for nsigma=1:length(sigma_x_);
sigma_x = sigma_x_(nsigma); sigma_y = sigma_y_(nsigma); sigma_z = sigma_z_(nsigma);
center_x = center_A_x_(nsigma); center_y = center_A_y_(nsigma); center_z = center_A_z_(nsigma);
tmp_gaus = tmp_gaus + (1/sqrt(pi)/sigma_x)*(1/sqrt(pi)/sigma_y)*(1/sqrt(pi)/sigma_z) * exp(-(X_m-center_x).^2/2/sigma_x.^2 -(Y_m-center_y).^2/2/sigma_y.^2 -(Z_m-center_z).^2/2/sigma_z.^2);
end;%for nsigma=1:length(sigma_x_);
F_A_x_c = tmp_gaus ;

tmp_gaus = 0 * R_m;
for nsigma=30:32;
sigma_x = sigma_x_(nsigma); sigma_y = sigma_y_(nsigma); sigma_z = sigma_z_(nsigma);
center_x = center_A_x_(nsigma); center_y = center_A_y_(nsigma); center_z = center_A_z_(nsigma);
tmp_gaus = tmp_gaus + (1/sqrt(pi)/sigma_x)*(1/sqrt(pi)/sigma_y)*(1/sqrt(pi)/sigma_z) * exp(-(X_m-center_x).^2/2/sigma_x.^2 -(Y_m-center_y).^2/2/sigma_y.^2 -(Z_m-center_z).^2/2/sigma_z.^2);
end;%for nsigma=1:length(sigma_x_);
for nsigma=27:29;
sigma_x = sigma_x_(nsigma); sigma_y = sigma_y_(nsigma); sigma_z = sigma_z_(nsigma);
center_x = center_A_x_(nsigma); center_y = center_A_y_(nsigma); center_z = center_A_z_(nsigma);
tmp_gaus = tmp_gaus - (1/sqrt(pi)/sigma_x)*(1/sqrt(pi)/sigma_y)*(1/sqrt(pi)/sigma_z) * exp(-(X_m-center_x).^2/2/sigma_x.^2 -(Y_m-center_y).^2/2/sigma_y.^2 -(Z_m-center_z).^2/2/sigma_z.^2);
end;%for nsigma=1:length(sigma_x_);
F_A_x_c = tmp_gaus ;

n_cut = floor(0.825*n_source);
c_A = 0:n_cut-1; c_B = n_cut:n_source-1;
center_B_r_ = [ 0.5*(c_A)/max(1,n_source-1) , 0.5*(1.0*(c_A(end)+1) + 3.5*(c_B-c_B(1)))/max(1,n_source-1) ];
center_B_w_ = [ 2*pi*2*(c_A)/(n_source) , 2*pi*2*(1.0*(c_A(end)+1) - 0.15*(c_B-c_B(1)).*(c_B-c_B(end)))/(n_source) ];
center_B_x_ = center_B_r_.*cos(center_B_w_);
center_B_y_ = center_B_r_.*sin(center_B_w_);
center_B_z_ = [ 1.0*(c_A)/max(1,n_source-1) , 1.0*(1.0*(c_A(end)+1) - 0.5*(c_B-c_B(1)).*(c_B-c_B(end)))/max(1,n_source-1) ] - 0.5;
tmp_gaus = 0 * R_m;
for nsigma=1:length(sigma_x_);
sigma_x = sigma_x_(nsigma); sigma_y = sigma_y_(nsigma); sigma_z = sigma_z_(nsigma);
center_x = center_B_x_(nsigma); center_y = center_B_y_(nsigma); center_z = center_B_z_(nsigma);
tmp_gaus = tmp_gaus + (1/sqrt(pi)/sigma_x)*(1/sqrt(pi)/sigma_y)*(1/sqrt(pi)/sigma_z) * exp(-(X_m-center_x).^2/2/sigma_x.^2 -(Y_m-center_y).^2/2/sigma_y.^2 -(Z_m-center_z).^2/2/sigma_z.^2);
end;%for nsigma=1:length(sigma_x_);
F_B_x_c = tmp_gaus ;

nF_ = 0:0;
for nF=nF_;

if (nF==0);
F_x_c = F_A_x_c ;
X = X_m; Y = Y_m; Z = Z_m;
end;% if;
if (nF==1);
F_x_c = F_B_x_c ;
X = X_m; Y = Y_m; Z = Z_m;
end;% if;
if (nF==2);
F_C_x_c = (4/7)*F_A_x_c + (3/7)*F_B_x_c ;
F_x_c = F_A_x_c - F_C_x_c ;
X = X_m; Y = Y_m; Z = Z_m;
end;% if;

F_k_c = recenter3(fftn(recenter3(F_x_c))); Fl_k_c = abs(F_k_c); 
%Fw_k_c = exp(-angle(F_k_c).^2/2/(pi/16)^2);
Fw_k_c = angle(F_k_c);

subplot(length(nF_),4,1 + 4*nF);
%cra = colormap('jet'); ncra = size(cra,1);
cra = colormap('spring'); cra = cra(end:-1:1,:); ncra = size(cra,1);
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
%hpatch.FaceAlpha = (nv/length(v_)).^2 * 0.5;
hpatch.FaceAlpha = (nv/length(v_)).^4 * 1.0;
%hpatch.EdgeAlpha = (nv/length(v_)).^2 * 1;
%title(sprintf('v %0.2f',v));
xlim(1*[-1,1]);ylim(1*[-1,1]);zlim(1*[-1,1]);
xlabel('x');ylabel('y');zlabel('z');
%daspect([1,4,4]); 
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for nv=1:length(v_);
title(sprintf('F_x(x)'));

subplot(length(nF_),4,2 + 4*nF);
%cra = colormap('jet'); ncra = size(cra,1);
cra = colormap('spring'); cra = cra(end:-1:1,:); ncra = size(cra,1);
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
%hpatch.FaceAlpha = (nv/length(v_)).^2 * 0.5;
hpatch.FaceAlpha = (nv/length(v_)).^4 * 1.0;
%hpatch.EdgeAlpha = (nv/length(v_)).^2 * 1;
%title(sprintf('v %0.2f',v));
xlim(1*[-1,1]);ylim(1*[-1,1]);zlim(1*[-1,1]);
xlabel('x');ylabel('y');zlabel('z');
%daspect([1,4,4]); 
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for nv=1:length(v_);
title(sprintf('Fl_k(k)'));

subplot(length(nF_),4,3 + 4*nF);
cra = colormap('cool'); cra = cra(end:-1:1,:); ncra = size(cra,1);
imagesc(abs(squeeze(F_k_c(:,:,end/2))));
title(sprintf('abs F_k_c at kz=0'));
axis square;

subplot(length(nF_),4,4 + 4*nF);
cra = colormap('hsv'); cra = cra(end:-1:1,:); ncra = size(cra,1);
imagesc(angle(squeeze(F_k_c(:,:,end/2))),[-pi,pi]);
title(sprintf('angle F_k_c at kz=0'));
axis square;

end;% for nF;

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
