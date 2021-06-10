function [n_all,theta_all_,phi_all_,weight_all_] = sample_shell_2(r,d,TorL) ;
% generate theta_ and phi_ arrays sampled on shell of radius r at equatorial_distance d ;
% Note: theta_ in [0,2*pi], phi_ in [0,1*pi];
% Note that: ;
% \int (sphere of radius K) exp(+i*k*delta) = I1 = I2 = I3 = I4, where: ;
% I1 = dblquad(@(k,phi) 2*pi*cos(delta.*k.*cos(phi)).*k.^2.*sin(phi),0,Kmax,0,pi);
% I2 = 4*pi*(1/delta^1)*quad(@(k) k.^2.*sin(delta*k)./k,0,Kmax);
% I3 = 4*pi*(1/delta^3)*quad(@(k) k.*sin(k),0,Kmax*delta);
% I4 = 4*pi*(1/delta^3)*(sin(Kd) - Kd*cos(Kd)); 
% Test quadrature in sphere with: ;
%{
addpath('~/chebfun');
K_max = 2; delta = 1.2; sample_d = 0.25;
 b_K = 4; [b_jx_,b_jw_] = jacpts(b_K,0,2);
k_ = (b_jx_+1.0)*K_max/2; k_w_ = b_jw_*K_max/2;
Iq_ = zeros(b_K,1); Ix_ = zeros(b_K,1);
[n_all,theta_all_,phi_all_,weight_all_] = sample_shell_2(1,sample_d,'L') ;
for nk=1:b_K;
k = k_(nk);
Iq_(nk) = sum(cos(k*delta*cos(phi_all_)).*weight_all_);
Ix_(nk) = 4*pi*(1/delta)*sin(delta*k).*k./k.^2;
end;% for nk=1:b_K;
Iq = k_w_*Iq_;
Ix = 4*pi*(1/delta^3)*( sin(K_max*delta) - (K_max*delta)*cos(K_max*delta) );
disp(sprintf(' %% (Ix-Iq)/Ix = %0.16f/%0.16f = %0.16f',Ix-Iq,Ix,(Ix-Iq)/Ix));
  %}
%{
  % Test quadrature with r^2 weight via:
K_max = 2; delta = 1.2;
Ix = 4*pi*(1/delta^3)*( sin(K_max*delta) - (K_max*delta)*cos(K_max*delta) );
 %%%%%%%%;
b_K = 6; 
 %%%%%%%%;
[b_jx_,b_jw_] = jacpts(b_K,0,2); k_jx_ = (b_jx_+1.0)*K_max/2; k_jw_ = b_jw_*K_max/2;
Ix_ = zeros(b_K,1); for nk=1:b_K; k = k_jx_(nk); Ix_(nk) = 4*pi*(1/delta)*sin(delta*k).*k./k.^2; end;% for nk=1:b_K;
Ij = k_jw_*Ix_; 
disp(sprintf(' %% (Ix-Ij)/Ix = %0.16f/%0.16f = %0.16f',Ix-Ij,Ix,(Ix-Ij)/Ix));
 %%%%%%%%;
[b_lx_,b_lw_] = legpts(b_K);     k_lx_ = (b_lx_+1.0)*K_max/2; k_lw_ = b_lw_*K_max/2;
Ix_ = zeros(b_K,1); for nk=1:b_K; k = k_lx_(nk); Ix_(nk) = 4*pi*(1/delta)*sin(delta*k).*k; end;% for nk=1:b_K;
Il = k_lw_*Ix_;
disp(sprintf(' %% (Ix-Il)/Ix = %0.16f/%0.16f = %0.16f',Ix-Il,Ix,(Ix-Il)/Ix));
 %}     
% Test quadrature on sphere with: ;
%{
  sample_shell_2();
 %}

if (nargin<1);

d_ = 2.^[-0.5:-0.5:-2.0];
prows = length(d_); pcols = 3;
for nd=1:length(d_);
d = d_(nd);
disp(sprintf(' %% sampling at equatorial_distance d %0.3f',d));
r=3.0;

subplot(prows,pcols,1 + (nd-1)*pcols);
n_theta = 128; n_phi = 64;
theta_ = linspace(0,2*pi,n_theta); phi_ = linspace(0,pi,n_phi);
[THETA_,PHI_]=meshgrid(theta_,phi_);
f_ = feval(@f_test,PHI_); g_ = feval(@g_test,THETA_); h_ = f_.*g_;
[X_,Y_,Z_]=sph2cart(THETA_,PHI_-pi/2,1 + 0.2*h_);
surf(X_,Y_,Z_,'linestyle','none'); 
xlabel('x'); ylabel('y'); zlabel('z'); title('original function');
view([-65,20]);
axis equal; %rot3d;
axis vis3d;
light; lighting phong; camzoom(1.3);

[n_all,theta_all_,phi_all_,weight_all_] = sample_shell_2(r,d,'T') ;
weight_avg = mean(weight_all_); weight_std = std(weight_all_);
disp(sprintf(' %% n_all %d %d %d',n_all,length(theta_all_),length(phi_all_)));
x_ = cos(theta_all_).*sin(phi_all_);
y_ = sin(theta_all_).*sin(phi_all_);
z_ = cos(phi_all_);
subplot(prows,pcols,2 + (nd-1)*pcols);
cra = colormap('jet'); ncra = size(cra,1); clim = weight_avg + 1.5*weight_std*[-1,1];
hold on;
for nall = 1:n_all;
x = x_(nall);
y = y_(nall);
z = z_(nall);
nb = max(1,min(ncra,floor(ncra*(weight_all_(nall) - clim(1))/diff(clim))));
plot3(x,y,z,'.','Color',cra(nb,:),'MarkerSize',15);
end;%for nall = 1:n_all;
hold off;
view([-65,20]);
xlabel('x'); ylabel('y'); zlabel('z'); title('T nodes colored by weight');
axis equal;
axis vis3d;
f_all_ = feval(@f_test,phi_all_); g_all_ = feval(@g_test,theta_all_); h_all_ = f_all_.*g_all_;
It = dot(weight_all_,h_all_);

[n_all,theta_all_,phi_all_,weight_all_] = sample_shell_2(r,d,'L') ;
weight_avg = mean(weight_all_); weight_std = std(weight_all_);
disp(sprintf(' %% n_all %d %d %d',n_all,length(theta_all_),length(phi_all_)));
x_ = cos(theta_all_).*sin(phi_all_);
y_ = sin(theta_all_).*sin(phi_all_);
z_ = cos(phi_all_);
subplot(prows,pcols,3 + (nd-1)*pcols);
cra = colormap('jet'); ncra = size(cra,1); clim = weight_avg + 1.5*weight_std*[-1,1];
hold on;
for nall = 1:n_all;
x = x_(nall);
y = y_(nall);
z = z_(nall);
nb = max(1,min(ncra,floor(ncra*(weight_all_(nall) - clim(1))/diff(clim))));
plot3(x,y,z,'.','Color',cra(nb,:),'MarkerSize',15);
end;%for nall = 1:n_all;
hold off;
view([-65,20]);
xlabel('x'); ylabel('y'); zlabel('z'); title('L nodes colored by weight');
axis equal;
axis vis3d;
f_all_ = feval(@f_test,phi_all_); g_all_ = feval(@g_test,theta_all_); h_all_ = f_all_.*g_all_;
Il = dot(weight_all_,h_all_);

Ix = r.^2 .* (F_test(pi)-F_test(0)) .* (G_test(2*pi)-G_test(0));
disp(sprintf(' %% Ix %0.3f It %0.3f error %0.16f Il %0.3f error %0.16f',Ix,It,norm(It-Ix),Il,norm(Il-Ix)));
end;%for nd=1:length(d_);

return;
end;%if (nargin<1);

n_equator = 1+round(2*pi*r/d);
n_phi = 1+round(n_equator/2);

if TorL == 'T';
[tsch_node_,tsch_weight_] = tsch_node_weight_0(n_phi,-1,1);
phi_ = acos(tsch_node_); dphi = mean(diff(sort(phi_)));
weight_ = tsch_weight_;
end;%if TorL == 'T';

if TorL == 'L';
[lgnd_node_,lgnd_weight_] = lgnd_node_weight_0(n_phi,-1,1);
phi_ = acos(lgnd_node_); dphi = mean(diff(sort(phi_)));
weight_ = lgnd_weight_;
end;%if TorL == 'L';

n_all = 0 ;
for nphi = 1:n_phi;
phi = phi_(nphi);
sphi = sin(phi);
n_theta = 1+round(2*pi*sphi*r/d);
n_all = n_all + n_theta;
end;%for nphi = 1:n_phi;
theta_all_ = zeros(n_all,1); phi_all_ = zeros(n_all,1); weight_all_ = zeros(n_all,1); ix=0;
for nphi = 1:n_phi;
phi = phi_(nphi);
sphi = sin(phi);
n_theta = 1+round(2*pi*sphi*r/d);
theta_ = linspace(0,2*pi,n_theta+1); theta_ = theta_(1:end-1); dtheta = mean(diff(theta_));
theta_all_(1+ix+(0:n_theta-1)) = transpose(theta_);
phi_all_(1+ix+(0:n_theta-1)) = phi*ones(n_theta,1);
weight_all_(1+ix+(0:n_theta-1)) = r.^2 .* weight_(nphi) .* dtheta .* ones(n_theta,1);
%weight_all_(1+ix+(0:n_theta-1)) = r.^2 .* dphi .* sin(phi)*dtheta .* ones(n_theta,1);
ix = ix+n_theta;
end;%for nphi = 1:n_phi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

function output = f_test(phi) ;
m=6;
output = exp(-(-m*cos(phi)).^2)*m;
%output = ones(size(phi));

function output = F_test(phi) ;
m=6;
output = sqrt(pi)/2*erf(-m*cos(phi));
%output = -cos(phi);

function output = g_test(theta) ;
w1 = 3; w2 = 5;
output = (cos(w1*theta) + 1) + (cos(w2*theta) + 1);
%output = ones(size(theta));

function output= G_test(theta) ;
w1 = 3; w2 = 5;
output = (sin(w1*theta)/w1 + theta) + (sin(w2*theta)/w2 + theta);
%output = theta;
