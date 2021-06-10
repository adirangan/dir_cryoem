function isosurface_ver8(n_,dirname,infix) ;
% for use with test_ver15.f, which has model_A,B,C, as well as Y_tru_ and Y_est_ ;
% to see model__A,model__B,model__C and blend, try: ;
% isosurface_ver8([-1,-2,-3,0],'./dir_k16_dF_s1_nF_cT_tF_hF_rF_M32_S16_s0_MS_r1','');
% to see model__A and iteration 16(tru), try: ;
% isosurface_ver8([-1,16],'./dir_k16_dF_s2_nF_cT_tF_hF_rF_M0_S0_s0_SMa0_r1','tru');
figure;

if length(n_)>1; prows = 2; pcols = ceil(length(n_)/prows); end;
if length(n_)<=1; prows = 1; pcols = 1; end;

for nF=1:length(n_);
ncur = n_(nF);

if (0);
elseif (ncur==-1);
F_x_c = MDA_read_c16(sprintf('%s/f_x_c_A___.mda',dirname)); F_x_c = real(permute(F_x_c,[2,1,3]));
X = MDA_read_r8(sprintf('%s/grid_x1_c___.mda',dirname)); X = permute(X,[2,1,3]);
Y = MDA_read_r8(sprintf('%s/grid_x2_c___.mda',dirname)); Y = permute(Y,[2,1,3]);
Z = MDA_read_r8(sprintf('%s/grid_x3_c___.mda',dirname)); Z = permute(Z,[2,1,3]);
elseif (ncur==-2);
F_x_c = MDA_read_c16(sprintf('%s/f_x_c_B___.mda',dirname)); F_x_c = real(permute(F_x_c,[2,1,3]));
X = MDA_read_r8(sprintf('%s/grid_x1_c___.mda',dirname)); X = permute(X,[2,1,3]);
Y = MDA_read_r8(sprintf('%s/grid_x2_c___.mda',dirname)); Y = permute(Y,[2,1,3]);
Z = MDA_read_r8(sprintf('%s/grid_x3_c___.mda',dirname)); Z = permute(Z,[2,1,3]);
elseif (ncur==-3);
F_x_c = MDA_read_c16(sprintf('%s/f_x_c_C___.mda',dirname)); F_x_c = real(permute(F_x_c,[2,1,3]));
X = MDA_read_r8(sprintf('%s/grid_x1_c___.mda',dirname)); X = permute(X,[2,1,3]);
Y = MDA_read_r8(sprintf('%s/grid_x2_c___.mda',dirname)); Y = permute(Y,[2,1,3]);
Z = MDA_read_r8(sprintf('%s/grid_x3_c___.mda',dirname)); Z = permute(Z,[2,1,3]);
elseif (ncur==0);
F_A_x_c = MDA_read_c16(sprintf('%s/f_x_c_A___.mda',dirname)); 
F_B_x_c = MDA_read_c16(sprintf('%s/f_x_c_B___.mda',dirname)); 
F_C_x_c = MDA_read_c16(sprintf('%s/f_x_c_C___.mda',dirname)); 
F_x_c = (3/7)*real(permute(F_A_x_c,[2,1,3])) + (2/7)*real(permute(F_B_x_c,[2,1,3])) + (2/7)*real(permute(F_C_x_c,[2,1,3]));
X = MDA_read_r8(sprintf('%s/grid_x1_c___.mda',dirname)); X = permute(X,[2,1,3]);
Y = MDA_read_r8(sprintf('%s/grid_x2_c___.mda',dirname)); Y = permute(Y,[2,1,3]);
Z = MDA_read_r8(sprintf('%s/grid_x3_c___.mda',dirname)); Z = permute(Z,[2,1,3]);
elseif (ncur>2);
F_x_c = MDA_read_c16(sprintf('%s/f_x_c_%s___%d_.mda',dirname,infix,ncur)); F_x_c = real(permute(F_x_c,[2,1,3]));
X = MDA_read_r8(sprintf('%s/grid_x1_c___.mda',dirname)); X = permute(X,[2,1,3]);
Y = MDA_read_r8(sprintf('%s/grid_x2_c___.mda',dirname)); Y = permute(Y,[2,1,3]);
Z = MDA_read_r8(sprintf('%s/grid_x3_c___.mda',dirname)); Z = permute(Z,[2,1,3]);
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
title(sprintf('%s: F_x(x) %d',infix,ncur));

end;% for nF;

if length(n_)>2; fname = sprintf('%s/f_x_c_%s___%d_%d_%d.jpg',dirname,infix,n_(2),round(mean(diff(n_(2:end)))),n_(end)); end;
if length(n_)==2; fname = sprintf('%s/f_x_c_%s___%d_%d.jpg',dirname,infix,n_(1),n_(2)); end;
if length(n_)<=1; fname = sprintf('%s/f_x_c_%s___%d.jpg',dirname,infix,n_(1)); end;

set(gcf,'Position',[1,1,1024,768]);
print('-djpeg',fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
