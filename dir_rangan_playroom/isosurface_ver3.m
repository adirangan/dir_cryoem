function isosurface_ver3(n_,dirname,tname)
figure;

%nlow = 4;
%nmax = 22;
%n_ = [0,[nlow:2:nmax]];
prows = 2; pcols = ceil(length(n_)/prows);
for nF=1:length(n_);
ncur = n_(nF);

if (ncur==0);
F_x_c = MDA_read_c16(sprintf('%s/S_ori_.mda',dirname)); F_x_c = real(permute(F_x_c,[2,1,3]));
X = MDA_read_r8(sprintf('%s/x_pre_.mda',dirname)); X = permute(X,[2,1,3]);
Y = MDA_read_r8(sprintf('%s/y_pre_.mda',dirname)); Y = permute(Y,[2,1,3]);
Z = MDA_read_r8(sprintf('%s/z_pre_.mda',dirname)); Z = permute(Z,[2,1,3]);
elseif (ncur>0);
F_x_c = MDA_read_c16(sprintf('%s/S_%s_%d_.mda',dirname,tname,ncur)); F_x_c = real(permute(F_x_c,[2,1,3]));
X = MDA_read_r8(sprintf('%s/x_pre_.mda',dirname)); X = permute(X,[2,1,3]);
Y = MDA_read_r8(sprintf('%s/y_pre_.mda',dirname)); Y = permute(Y,[2,1,3]);
Z = MDA_read_r8(sprintf('%s/z_pre_.mda',dirname)); Z = permute(Z,[2,1,3]);
end;% if;

F_k_c = recenter3(fftn(recenter3(F_x_c))); Fl_k_c = abs(F_k_c); Fw_k_c = angle(F_k_c);
subplot(prows,pcols,nF);
cra = colormap('spring'); cra = cra(end:-1:1,:); 
%cra = colormap('cool'); 
ncra = size(cra,1);
v_avg = mean(F_x_c(:)); v_std = std(F_x_c(:)); v_min = min(F_x_c(:)); v_max = max(F_x_c(:));
vlim = v_avg + 2.5*v_std*[-1,1];
%vlim = [v_min , v_max];
v_ = linspace(vlim(1),vlim(2),7);
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
title(sprintf('F_x(x) %d',ncur));

end;% for nF;

set(gcf,'Position',[1,1,1024,768]);
print('-djpeg',sprintf('%s/S_%s_%d_%d_%d.jpg',dirname,tname,n_(2),round(mean(diff(n_(2:end)))),n_(end)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
