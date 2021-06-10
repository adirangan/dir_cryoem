function [F__] = gen_Jsvd_Jmn_0(N_pixel,eps_target);
% returns jacobi-polynomial coefficients of J_{l}. ;
% F__(m,n,l) = \int \int J_{l}(\delta k) P_{m}(\delta) P_{n}(k) \deltad\delta kdk

if nargin<1; N_pixel = 3.0*sqrt(2); end;
if nargin<2; eps_target = 1e-14; end;

%N_pixel = 3.0*sqrt(2); eps_target = 1e-14;
setup; verbose=2; K_max = 1.0; 
if (N_pixel<=5); l_max = 64;a_K = 48;b_K=48; end;
if (N_pixel> 5); l_max = 64;a_K = 48;b_K=48; end;

%%%%%%%%;
K_target = K_max;
z_target = N_pixel*pi*sqrt(2);
D_target = z_target/(2*pi*K_target);
%%%%%%%%;
r_max = 2*pi*K_target;
d_max = D_target;
a_m = r_max/2; a_r = a_m;
b_m = d_max/2; b_r = b_m;
%%%%%%%%;
[a_jx,a_jw] = jacpts(a_K,0,1); a_jw=transpose(a_jw);
a_Jv_ = cell(1+a_K,1+a_K); a_Jx = zeros(a_K,a_K);
for nkA=0:a_K;
aj = jacpoly(nkA,0,1)*sqrt(nkA+1)/sqrt(2);
if (nkA<a_K); a_Jx(1+nkA,:) = aj(a_jx); end;
a_Jv_{1+nkA} = aj;
end;%for nkA=1:a_K;
a_jt = a_jx*a_r + a_m;
%%%%%%%%;
[b_jx,b_jw] = jacpts(b_K,0,1); b_jw=transpose(b_jw);
b_Jv_ = cell(1+b_K,1+b_K); b_Jx = zeros(b_K,b_K);
for nkB=0:b_K;
bj = jacpoly(nkB,0,1)*sqrt(nkB+1)/sqrt(2);
if (nkB<b_K); b_Jx(1+nkB,:) = bj(b_jx); end;
b_Jv_{1+nkB} = bj;
end;%for nkB=1:b_K;
b_jt = b_jx*b_r + b_m;
%%%%%%%%;
[A_jt_,B_jt_] = meshgrid(a_jt,b_jt);
%%%%%%%%;
F__ = zeros(b_K,a_K,2*l_max+1);
for tmp_l = -l_max:+l_max;
F = @(a,b) besselj(tmp_l,a.*b);
F_jt_ = F(A_jt_,B_jt_);
F_ = zeros(b_K,a_K);
jw_ = b_jw*transpose(a_jw);
for nkA=0:a_K-1;for nkB=0:b_K-1;
J_tmp_ = transpose(b_Jx(1+nkB,:))*a_Jx(1+nkA,:);
S_tmp = F_jt_.*J_tmp_.*jw_;
F_(1+nkB,1+nkA) = sum(S_tmp(:));
end;end;%for nkA=0:a_K-1;for nkB=0:b_K-1;
F__(:,:,1+l_max+tmp_l) = F_;
end;%for tmp_l = -l_max:+l_max;
%%%%%%%%;

flag_plot=0;
if flag_plot;
figure();clf;
llim_ = [log10(eps_target),0];
prows = 5; pcols = 9;
tmp_l_ = round(linspace(-l_max,+l_max,prows*pcols));
colormap(colormap_beach());
for nl=1:length(tmp_l_);
tmp_l = tmp_l_(nl);
subplot(prows,pcols,nl);
imagesc(log10(abs(F__(:,:,1+l_max+tmp_l))),llim_);
%imagesc(((F__(:,:,1+l_max+tmp_l))),10.^llim_);
set(gca,'Xtick',[],'Ytick',[]);
end;%for nl=1:length(tmp_l_);
end;%if flag_plot;

flag_plot=1;
if flag_plot;
G__ = F__(:,:,1+l_max:1+2*l_max); llim_ = [log10(eps_target),0];
c_ = colormap(colormap_beach()); n_c = size(c_,1);
figure();clf;
%%%%%%%%;
subplot(1,3,1);
hold on; 
eps_contour_ = 10.^[-2,-4,-6,-8];
for ne=1:length(eps_contour_);
eps_contour = eps_contour_(ne);
nx = (log10(eps_contour)-min(llim_))/diff(llim_);
tmp_X = min(llim_) + nx*diff(llim_);
disp(sprintf(' %% contour at level %0.4f',tmp_X));
hpatch = patch(isosurface(0:b_K-1,0:a_K-1,0:+l_max,log10(abs(G__)),tmp_X));
isonormals(0:b_K-1,0:a_K-1,0:+l_max,log10(abs(G__)),hpatch);
nc = max(1,min(n_c,floor(n_c*nx)));
hpatch.FaceColor = c_(nc,:);
hpatch.EdgeColor = 'none';
hpatch.FaceAlpha = 0.45;
end;%for ne=1:length(eps_contour_);
xlim([1,b_K]-1);ylim([1,a_K]-1);zlim([0,+l_max]);
xlabel('m');ylabel('n');zlabel('l');
set(gca,'Xtick',[1,b_K]-1,'Ytick',[1,a_K]-1,'Ztick',[0,+l_max]);
tmp_l = line([0;b_K-1],[0;0],[0;0]); set(tmp_l,'Color',0.85*[1,1,1],'LineWidth',2);
tmp_l = line([0;0],[0;a_K-1],[0;0]); set(tmp_l,'Color',0.85*[1,1,1],'LineWidth',2);
tmp_l = line([0;0],[0;0],[0;+l_max]); set(tmp_l,'Color',0.85*[1,1,1],'LineWidth',2);
for na=0:4:a_K-1;
tmp_l = line([0;b_K-1],[na;na],[0;0]); set(tmp_l,'Color',0.85*[1,1,1],'LineWidth',0.5);
tmp_l = line([0;0],[na;na],[0;+l_max]); set(tmp_l,'Color',0.85*[1,1,1],'LineWidth',0.5);
end;%for na=0:4:a_K-1;
for nb=0:4:b_K-1;
tmp_l = line([nb;nb],[0;a_K-1],[0;0]); set(tmp_l,'Color',0.85*[1,1,1],'LineWidth',0.5);
tmp_l = line([nb;nb],[0;0],[0;+l_max]); set(tmp_l,'Color',0.85*[1,1,1],'LineWidth',0.5);
end;%for nb=0:4:b_K-1;
for nl=0:4:l_max;
tmp_l = line([0;b_K-1],[0;0],[nl;nl]); set(tmp_l,'Color',0.85*[1,1,1],'LineWidth',0.5);
tmp_l = line([0;0],[0;a_K-1],[nl;nl]); set(tmp_l,'Color',0.85*[1,1,1],'LineWidth',0.5);
end;%for nl=0:4:l_max;
hold off;
view([-45+90+90,20]); 
axis vis3d; camlight right; lighting gouraud;
title(sprintf('W %0.2f',N_pixel/sqrt(2)));
%%%%%%%%;
subplot(1,3,2);
G_ = squeeze(G__(:,1,:));
colormap(colormap_beach());
imagesc(log10(abs(transpose(G_))),llim_); colorbar;
set(gca,'YDir','normal');
axis image;
xlabel('m'); ylabel('l');
title(sprintf('W %0.2f',N_pixel/sqrt(2)));
%%%%%%%%;
subplot(1,3,3);
G_ = squeeze(G__(:,:,1));
colormap(colormap_beach());
imagesc(log10(abs(G_)),llim_); colorbar;
set(gca,'YDir','normal');
xlabel('n'); ylabel('m');
axis image;
title(sprintf('W %0.2f',N_pixel/sqrt(2)));
%%%%%%%%;
set(gcf,'Position',1+[0,0,512*3,512*1]);
fname_fig = sprintf('%s/dir_svd/gen_Jsvd_Jmn_N%3d',pwd,round(100*N_pixel));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if flag_plot;

%{

  llim_ = [log(eps_target),0];
  G_ = squeeze(G__(:,1,:));
  H_ = log(abs(G_));
  figure;clf;
  subplot(1,3,1);
  colormap(colormap_beach());
  imagesc(transpose(H_),llim_); colorbar;
  set(gca,'YDir','normal');
  axis image;
  xlabel('m'); ylabel('l');
  title(sprintf('W %0.2f',N_pixel/sqrt(2)));
  subplot(1,3,2);
  tmp_l_ = repmat(0:l_max,a_K,1);
  tmp_m_ = repmat(transpose(0:a_K-1),1,l_max+1);
  tmp_a = 1.25; tmp_b = 1.25; 
  tmp_c = exp(tmp_a)*(3/2 + 2*cosh(tmp_b) + 0.5*cosh(2*tmp_b));
  %tmp_c = 2;
  tmp_g_ = tmp_c - tmp_a*tmp_l_ - tmp_b*tmp_m_;
  colormap(colormap_beach());
  imagesc(transpose(tmp_g_),llim_); colorbar;
  set(gca,'YDir','normal');
  axis image;
  xlabel('m'); ylabel('l');
  title(sprintf('W %0.2f',N_pixel/sqrt(2)));
  subplot(1,3,3);
  colormap(colormap_beach());
  imagesc(transpose(tmp_g_-H_),[-7,+7]); colorbar;
  set(gca,'YDir','normal');
  axis image;
  xlabel('m'); ylabel('l');
  title(sprintf('W %0.2f',N_pixel/sqrt(2)));

  %}
