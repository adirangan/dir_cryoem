function [n_all,theta_all_,phi_all_,weight_all_] = sample_shell_3_nojvm(r,d,TorL) ;
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

setup;
 sample_d_ = 0.5.^[-2:0.125:4]; n_sample_d = length(sample_d_);
 E_ = zeros(n_sample_d,1);
 for nsample_d = 1:n_sample_d;
 sample_d = sample_d_(nsample_d);
 K_max = 7.5; delta = 2.32; 
 b_K = 1+ceil(K_max/2/sample_d); [b_jx_,b_jw_] = jacpts(b_K,0,2);
 k_ = (b_jx_+1.0)*K_max/2; k_w_ = b_jw_*(K_max/2)^3;
 Iq_ = zeros(b_K,1); Ix_ = zeros(b_K,1);
 [n_all,theta_all_,phi_all_,weight_all_] = sample_shell_3_nojvm(1,sample_d,'L') ;
 for nk=1:b_K;
 k = k_(nk);
 Iq_(nk) = sum(cos(k*delta*cos(phi_all_)).*weight_all_);
 Ix_(nk) = 4*pi*(1/delta)*sin(delta*k).*k./k.^2;
 end;% for nk=1:b_K;
 Iq = k_w_*Iq_;
 Ix = 4*pi*(1/delta^3)*( sin(K_max*delta) - (K_max*delta)*cos(K_max*delta) );
 disp(sprintf(' %% sample_d = %0.6f n_all %d b_K %d (Ix-Iq)/Ix = %0.16f/%0.16f = %0.16f',sample_d,n_all,b_K,Ix-Iq,Ix,(Ix-Iq)/Ix));
 E_(nsample_d) = abs((Ix-Iq)/Ix);
 end;%for nsample_d = 1:n_sample_d;
% plot(-log2(sample_d_),log10(E_),'o-'); xlabel('-log2(sample_d)'); ylabel('log10(E)');

  %}
%{
  % Test quadrature with r^2 weight via:
K_max = 7; delta = 1.2;
Ix = 4*pi*(1/delta^3)*( sin(K_max*delta) - (K_max*delta)*cos(K_max*delta) );
 %%%%%%%%;
b_K = 6; 
 %%%%%%%%;
 [b_jx_,b_jw_] = jacpts(b_K,0,2); k_jx_ = (b_jx_+1.0)*K_max/2; k_jw_ = b_jw_*(K_max/2)^3;
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
  sample_shell_3_nojvm();
 %}

if (nargin<1);

d_ = 2.^[-0.5:-0.25:-2.0];
%prows = length(d_); pcols = 3;
for nd=1:length(d_);
r=3.0;
d = d_(nd); n_equator = 3+round(2*pi*r/d); n_phi = 3+round(n_equator/2);
disp(sprintf(' %% sampling at equatorial_distance d %0.3f (n_phi %d)',d,n_phi));

[n_all,theta_all_,phi_all_,weight_all_] = sample_shell_3_nojvm(r,d,'T') ;
f_all_ = feval(@f_test,phi_all_); g_all_ = feval(@g_test,theta_all_); h_all_ = f_all_.*g_all_;
It = dot(weight_all_,h_all_);
[n_all,theta_all_,phi_all_,weight_all_] = sample_shell_3_nojvm(r,d,'C') ;
f_all_ = feval(@f_test,phi_all_); g_all_ = feval(@g_test,theta_all_); h_all_ = f_all_.*g_all_;
Ic = dot(weight_all_,h_all_);
[n_all,theta_all_,phi_all_,weight_all_] = sample_shell_3_nojvm(r,d,'L') ;
f_all_ = feval(@f_test,phi_all_); g_all_ = feval(@g_test,theta_all_); h_all_ = f_all_.*g_all_;
Il = dot(weight_all_,h_all_);

Ix = r.^2 .* (F_test(pi)-F_test(0)) .* (G_test(2*pi)-G_test(0));
et = norm(It-Ix); let = -log10(et);
ec = norm(Ic-Ix); lec = -log10(ec);
el = norm(Il-Ix); lel = -log10(el);
disp(sprintf(' %% Ix %0.3f It %0.3f error %0.16f (%0.2f) Ic %0.3f error %0.16f (%0.2f) Il %0.3f error %0.16f (%0.2f)',Ix,It,et,let,Ic,ec,lec,Il,el,lel));
end;%for nd=1:length(d_);

return;
end;%if (nargin<1);

n_equator = 3+round(2*pi*r/d);
n_phi = 3+round(n_equator/2);

if TorL == 'T';
%[tsch_node_,tsch_weight_] = tsch_node_weight_0(n_phi,-1,1);
[tsch_node_,tsch_weight_] = chebpts(n_phi,1);
phi_ = acos(tsch_node_); dphi = mean(diff(sort(phi_)));
weight_ = tsch_weight_;
end;%if TorL == 'T';

if TorL == 'C';
[tsch_node_,tsch_weight_] = chebpts(n_phi,2);
phi_ = acos(tsch_node_); dphi = mean(diff(sort(phi_)));
weight_ = tsch_weight_;
end;%if TorL == 'T';

if TorL == 'L';
%[lgnd_node_,lgnd_weight_] = lgnd_node_weight_0(n_phi,-1,1);
[lgnd_node_,lgnd_weight_] = legpts(n_phi);
phi_ = acos(lgnd_node_); dphi = mean(diff(sort(phi_)));
weight_ = lgnd_weight_;
end;%if TorL == 'L';

n_all = 0 ;
for nphi = 1:n_phi;
phi = phi_(nphi);
sphi = sin(phi);
n_theta = 3+round(2*pi*sphi*r/d);
n_all = n_all + n_theta;
end;%for nphi = 1:n_phi;
theta_all_ = zeros(n_all,1); phi_all_ = zeros(n_all,1); weight_all_ = zeros(n_all,1); ix=0;
for nphi = 1:n_phi;
phi = phi_(nphi);
sphi = sin(phi);
n_theta = 3+round(2*pi*sphi*r/d);
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
