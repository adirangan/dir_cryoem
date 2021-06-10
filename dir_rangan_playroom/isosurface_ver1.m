for nF=0:1;

if (nF==0);
F_x_c = MDA_read_c16('./S_pre_.mda'); F_x_c = real(permute(F_x_c,[2,1,3]));
X = MDA_read_r8('./x_pre_.mda'); X = permute(X,[2,1,3]);
Y = MDA_read_r8('./y_pre_.mda'); Y = permute(Y,[2,1,3]);
Z = MDA_read_r8('./z_pre_.mda'); Z = permute(Z,[2,1,3]);
end;% if;
if (nF==1);
F_x_c = MDA_read_c16('./S_pos_.mda'); F_x_c = real(permute(F_x_c,[2,1,3]));
X = MDA_read_r8('./x_pos_.mda'); X = permute(X,[2,1,3]);
Y = MDA_read_r8('./y_pos_.mda'); Y = permute(Y,[2,1,3]);
Z = MDA_read_r8('./z_pos_.mda'); Z = permute(Z,[2,1,3]);
end;% if;

F_k_c = recenter3(fftn(recenter3(F_x_c))); Fl_k_c = abs(F_k_c); Fw_k_c = angle(F_k_c);

subplot(2,2,1 + 2*nF);
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
xlabel('x');ylabel('y');zlabel('z');
%daspect([1,4,4]); 
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for nv=1:length(v_);
title(sprintf('F_x(x)'));

subplot(2,2,2 + 2*nF);
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
xlabel('x');ylabel('y');zlabel('z');
%daspect([1,4,4]); 
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for nv=1:length(v_);
title(sprintf('Fl_k(k)'));

end;% for nF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
[X_m,Y_m,Z_m] = meshgrid(linspace(-1,1,64)); 
R_m = sqrt(X_m.^2+Y_m.^2+Z_m.^2); A_m = atan2(Y_m,X_m); B = atan2(Z_m,sqrt(X_m.^2+Y_m.^2));
n_source = 128;
sigma = 1/32;
sigma_x_ = sigma*ones(1,n_source); sigma_y_ = sigma*ones(1,n_source); sigma_z_ = sigma*ones(1,n_source);
center_r_ = 2*pi*(0:n_source-1)/max(1,n_source-1);
center_x_ = cos(2*center_r_).*(3 + cos(3*center_r_))/8;
center_y_ = sin(2*center_r_).*(3 + cos(3*center_r_))/8;
center_z_ = 2.0*sin(3*center_r_)/8;
tmp_gaus = 0 * (1/sqrt(pi)/sigma)^3 * exp(-R_m.^2/2/sigma.^2);
for nsigma=1:length(sigma_x_);
sigma_x = sigma_x_(nsigma); sigma_y = sigma_y_(nsigma); sigma_z = sigma_z_(nsigma);
center_x = center_x_(nsigma); center_y = center_y_(nsigma); center_z = center_z_(nsigma);
tmp_gaus = tmp_gaus + (1/sqrt(pi)/sigma_x)*(1/sqrt(pi)/sigma_y)*(1/sqrt(pi)/sigma_z) * exp(-(X_m-center_x).^2/2/sigma_x.^2 -(Y_m-center_y).^2/2/sigma_y.^2 -(Z_m-center_z).^2/2/sigma_z.^2);
end;%for nsigma=1:length(sigma_x_);
F_x_c = tmp_gaus;
X = X_m;Y = Y_m;Z = Z_m;
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
xlabel('x');ylabel('y');zlabel('z');
%daspect([1,4,4]); 
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for nv=1:length(v_);
title(sprintf('F_x(x)'));
 %}
