function test_spharm_0(l_,m_,a_,res_);
% l_ = array of l-values ;
% m_ = array of m-values ;
% a_ = array of coefficients ;
% note that l_(ns) should be higher than m_(ns) ;
% res_ = [ resolution_theta , resolution_phi ] ;
% test with: ;
%{
  test_spharm_0();
  % or something like: ;
  l_ = [+0,+1,+1,+1,+2,+2,+2,+2,+2];
  m_ = [+0,-1,+0,+1,-2,-1,+0,+1,+2];
  %a_ = [+0,+0,+0,+0,+0,+0,+0,+1,+0];
  a_ = linspace(-1,1,9);
  res_ = [64,64];
  test_spharm_0(l_,m_,a_,res_);
  %}

if nargin<4; res_ = [64,64]; end;
if nargin<3; a_ = linspace(-1,1,9); end;
if nargin<2; m_ = [+0,-1,+0,+1,-2,-1,+0,+1,+2]; end;
if nargin<1; l_ = [+0,+1,+1,+1,+2,+2,+2,+2,+2]; end;
if nargin<1;
isph_start_ = MDA_read_i4('./dir_mda6/isph_start_.mda');
nterms_sph_ = MDA_read_i4('./dir_mda6/nterms_sph_.mda');
modsph_A_ori_ = MDA_read_c16('./dir_mda6/modsph_A_ori_.mda');
ngridr = length(isph_start_);
ncur = 10;
n_l = nterms_sph_(ncur);
n_lm = (n_l+1).^2;
l_ = []; m_ = [];
for nl=0:n_l;
l_ = [l_ , nl*ones(1,2*nl+1) ];
m_ = [m_ , [-nl:nl] ];
end;%for nl=0:n_l;
a_ = modsph_A_ori_(isph_start_(ncur) + (1:n_lm));
end;%if nargin<1;

theta_ = linspace( 0 , 2*pi , res_(1) );  % Azimuthal/Longitude/Circumferential
phi_   = linspace( 0 ,   pi , res_(2) );  % Altitude /Latitude /Elevation
[THETA_,PHI_]=meshgrid(theta_,phi_);
  
n_t = length(l_);
Lmn_ = cell(n_t,1);
Ymn_ = cell(n_t,1);
 
for nt=1:n_t;
l_val = l_(nt);
m_val = m_(nt);
m_abs = abs(m_val);
Lmn_{nt}=legendre(l_val,cos(PHI_),'unnorm');
if l_val~=0;
Lmn_{nt}=squeeze(Lmn_{nt}(m_abs+1,:,:));
end;%if l_val~=0;
a1=((2*l_val+1)/(4*pi));
a2=exp(lfactorial(l_val-m_abs) - lfactorial(l_val+m_abs));
c=sqrt(a1*a2);
Ymn_{nt}=c*Lmn_{nt}.*exp(i*m_val*THETA_);
end;%for nt=1:n_t;

Ymn = zeros(size(Ymn_{1}));
for nt=1:n_t;
Ymn = Ymn + a_(nt) * Ymn_{nt};
end;%for nt=1:n_t;

[Xm,Ym,Zm]=sph2cart(THETA_,PHI_-pi/2,abs(Ymn).^2);
[Xr,Yr,Zr]=sph2cart(THETA_,PHI_-pi/2,real(Ymn).^2);
[Xi,Yi,Zi]=sph2cart(THETA_,PHI_-pi/2,imag(Ymn).^2);
% [Xp,Yp,Zp]=sph2cart(THETA_,PHI_-pi/2,angle(Ymn).^2);
plot_flag=1;
if plot_flag;
f=figure; axis off; hold on;
axes('position',[0.0500 0 0.2666 1]); 
surf(Xm,Ym,Zm,'linestyle','none'); 
axis equal off; %rot3d;
axis vis3d;
light; lighting phong; camzoom(1.3);
axes('position',[0.3666 0 0.2666 1]); 
surf(Xr,Yr,Zr,'linestyle','none'); 
axis equal off; %rot3d;
axis vis3d;
light; lighting phong; camzoom(1.3);
axes('position',[0.6833 0 0.2666 1]); 
surf(Xi,Yi,Zi,'linestyle','none'); 
%    surf(Xp,Yp,Zp,'linestyle','none'); 
axis equal off; %rot3d;
axis vis3d;
light; lighting phong; camzoom(1.3);
axes('position',[0 0.9 1 0.1]); axis off;
t(1)=text(0.50,0.25,'Spherical Harmonics','HorizontalAlignment','Center');
axes('position',[0 0 1 0.1]); axis off;
t(2)=text(0.20,0.25,['|Y|^2'],'HorizontalAlignment','Center');
t(3)=text(0.50,0.25,['Real(Y)^2'],'HorizontalAlignment','Center');
t(4)=text(0.80,0.25,['Imag(Y)^2'],'HorizontalAlignment','Center');
%setfig(gcf,10,5,12);
end;%if plot_flag;

[X1,Y1,Z1]=sph2cart(THETA_,PHI_-pi/2,1.0);
cra = colormap('hsv'); ncra = size(cra,1);
y_min = 0; y_max = prctile(abs(Ymn(:)),85);
C1 = zeros(res_(1),res_(2),3);
for nr=1:res_(1); for nc=1:res_(2);
y_tmp = Ymn(nr,nc);
y_ang = angle(y_tmp); y_nrm = abs(y_tmp);
nba = max(1,min(ncra,floor(ncra*(y_ang+pi)/(2*pi))));
dbn = max(0,min(1,(y_nrm-y_min)/(y_max-y_min)));
C1(nr,nc,:) = cra(nba,:)*dbn;
end;end;%for nr=1:res_(1); for nc=1:res_(2);
plot_flag=1;
if plot_flag;
figure;
s=surf(X1,Y1,Z1,C1);set(s,'EdgeColor','none');
axis equal off; %rot3d;
axis vis3d;
%light; lighting phong;
end;%if plot_flag;

return
