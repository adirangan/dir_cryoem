function ...
[ ...
 n_all ...
,azimu_b_all_ ...
,polar_a_all_ ...
,weight_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,n_polar_a ...
,polar_a_ ...
,n_azimu_b_ ...
] = ...
sample_shell_6( ...
 r ...
,d ...
,str_T_vs_L ...
,flag_uniform_over_polar_a ...
) ;
% generate azimu_b_ and polar_a_ arrays sampled on shell of radius r at equatorial_distance d ;
% Note: azimu_b_ in [0,2*pi], polar_a_ in [0,1*pi];
% Note that: ;
% \int (sphere of radius K) exp(+i*k*delta) = I1 = I2 = I3 = I4, where: ;
% I1 = dblquad(@(k,polar_a) 2*pi*cos(delta.*k.*cos(polar_a)).*k.^2.*sin(polar_a),0,Kmax,0,pi);
% I2 = 4*pi*(1/delta^1)*quad(@(k) k.^2.*sin(delta*k)./k,0,Kmax);
% I3 = 4*pi*(1/delta^3)*quad(@(k) k.*sin(k),0,Kmax*delta);
% I4 = 4*pi*(1/delta^3)*(sin(Kd) - Kd*cos(Kd)); 

str_thisfunction = 'sample_shell_6';

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
disp(sprintf(' %% testing quadrature with r^2 weight: '));
K_max = 7; delta = 1.2;
Ix = 4*pi*(1/delta^3)*( sin(K_max*delta) - (K_max*delta)*cos(K_max*delta) );
b_K = 6; [b_jx_,b_jw_] = jacpts(b_K,0,2); k_jx_ = (b_jx_+1.0)*K_max/2; k_jw_ = b_jw_*(K_max/2)^3;
Ix_ = zeros(b_K,1); for nk=1:b_K; k = k_jx_(nk); Ix_(nk) = 4*pi*(1/delta)*sin(delta*k).*k./k.^2; end;% for nk=1:b_K;
Ij = k_jw_*Ix_; disp(sprintf(' %% (Ix-Ij)/Ix = %0.16f/%0.16f = %0.16f',Ix-Ij,Ix,(Ix-Ij)/Ix));
[b_lx_,b_lw_] = legpts(b_K); k_lx_ = (b_lx_+1.0)*K_max/2; k_lw_ = b_lw_*K_max/2;
Ix_ = zeros(b_K,1); for nk=1:b_K; k = k_lx_(nk); Ix_(nk) = 4*pi*(1/delta)*sin(delta*k).*k; end;% for nk=1:b_K;
Il = k_lw_*Ix_; disp(sprintf(' %% (Ix-Il)/Ix = %0.16f/%0.16f = %0.16f',Ix-Il,Ix,(Ix-Il)/Ix));
%%%%%%%%;
disp(sprintf(' %% testing basic quadrature on sphere: '));
sample_d_ = 0.5.^[-2:0.125:4]; n_sample_d = length(sample_d_);
for flag_uniform_over_polar_a=[0,1];
E_ = zeros(n_sample_d,1);
for nsample_d = 1:n_sample_d;
sample_d = sample_d_(nsample_d);
K_max = 7.5; delta = 2.32; 
b_K = 1+ceil(K_max/2/sample_d); [b_jx_,b_jw_] = jacpts(b_K,0,2);
k_ = (b_jx_+1.0)*K_max/2; k_w_ = b_jw_*(K_max/2)^3;
Iq_ = zeros(b_K,1); Ix_ = zeros(b_K,1);
[n_all,azimu_b_all_,polar_a_all_,weight_all_] = sample_shell_6(1,sample_d,'L',flag_uniform_over_polar_a) ;
for nk=1:b_K;
k = k_(nk);
Iq_(nk) = sum(cos(k*delta*cos(polar_a_all_)).*weight_all_);
Ix_(nk) = 4*pi*(1/delta)*sin(delta*k).*k./k.^2;
end;% for nk=1:b_K;
Iq = k_w_*Iq_;
Ix = 4*pi*(1/delta^3)*( sin(K_max*delta) - (K_max*delta)*cos(K_max*delta) );
disp(sprintf(' %% flag_uniform_over_polar_a %d sample_d = %0.6f n_all %.4d b_K %.2d (Ix-Iq)/Ix = %+0.16f/%+0.16f = %+0.16f',flag_uniform_over_polar_a,sample_d,n_all,b_K,Ix-Iq,Ix,(Ix-Iq)/Ix));
E_(nsample_d) = abs((Ix-Iq)/Ix);
end;%for nsample_d = 1:n_sample_d;
end;%for flag_uniform_over_polar_a=[0,1];
%%%%%%%%;
disp(sprintf(' %% testing different equatorial distances: '));
d_ = 2.^[-0.5:-0.25:-2.0];
for flag_uniform_over_polar_a=[0,1];
%prows = length(d_); pcols = 3;
for nd=1:length(d_);
r=3.0;
d = d_(nd); n_equator = 3+round(2*pi*r/d); n_polar_a = 3+round(n_equator/2);
disp(sprintf(' %% sampling at equatorial_distance d %0.3f (n_polar_a %d)',d,n_polar_a));
%%%%;
[n_all,azimu_b_all_,polar_a_all_,weight_all_] = sample_shell_6(r,d,'T',flag_uniform_over_polar_a) ;
f_all_ = feval(@f_test,polar_a_all_); g_all_ = feval(@g_test,azimu_b_all_); h_all_ = f_all_.*g_all_;
It = dot(weight_all_,h_all_);
[n_all,azimu_b_all_,polar_a_all_,weight_all_] = sample_shell_6(r,d,'C',flag_uniform_over_polar_a) ;
f_all_ = feval(@f_test,polar_a_all_); g_all_ = feval(@g_test,azimu_b_all_); h_all_ = f_all_.*g_all_;
Ic = dot(weight_all_,h_all_);
[n_all,azimu_b_all_,polar_a_all_,weight_all_] = sample_shell_6(r,d,'L',flag_uniform_over_polar_a) ;
f_all_ = feval(@f_test,polar_a_all_); g_all_ = feval(@g_test,azimu_b_all_); h_all_ = f_all_.*g_all_;
Il = dot(weight_all_,h_all_);
%%%%;
Ix = r.^2 .* (F_test(pi)-F_test(0)) .* (G_test(2*pi)-G_test(0));
et = norm(It-Ix); let = -log10(et);
ec = norm(Ic-Ix); lec = -log10(ec);
el = norm(Il-Ix); lel = -log10(el);
disp(sprintf(' %% flag_uniform_over_polar_a %d n_all %d Ix %0.3f It %0.3f error %0.16f (%0.2f) Ic %0.3f error %0.16f (%0.2f) Il %0.3f error %0.16f (%0.2f)',flag_uniform_over_polar_a,n_all,Ix,It,et,let,Ic,ec,lec,Il,el,lel));
end;%for nd=1:length(d_);
end;%for flag_uniform_over_polar_a=[0,1];
%%%%%%%%;
disp(sprintf(' %% returning')); return;
end;%if (nargin<1);

if (nargin<2); d = []; end;
if (nargin<3); str_T_vs_L = []; end;
if (nargin<4); flag_uniform_over_polar_a = []; end;
if isempty(d); d = 1/(2*pi); end;
if isempty(str_T_vs_L); str_T_vs_L = 'L'; end;
if isempty(flag_uniform_over_polar_a); flag_uniform_over_polar_a = 0; end;

n_equator = 3+round(2*pi*r/max(1e-12,d));
n_polar_a = 3+round(n_equator/2);

if str_T_vs_L == 'T';
%[tsch_node_,tsch_weight_] = tsch_node_weight_0(n_polar_a,-1,1);
[tsch_node_,tsch_weight_] = chebpts(n_polar_a,1);
polar_a_ = acos(tsch_node_); dpolar_a = mean(diff(sort(polar_a_)));
weight_ = tsch_weight_;
end;%if str_T_vs_L == 'T';

if str_T_vs_L == 'C' | str_T_vs_L == 'C2';
[tsch_node_,tsch_weight_] = chebpts(n_polar_a,2);
polar_a_ = acos(tsch_node_); dpolar_a = mean(diff(sort(polar_a_)));
weight_ = tsch_weight_;
end;%if str_T_vs_L == 'C' | str_T_vs_L == 'C2';

if str_T_vs_L == 'L';
%[lgnd_node_,lgnd_weight_] = lgnd_node_weight_0(n_polar_a,-1,1);
[lgnd_node_,lgnd_weight_] = legpts(n_polar_a);
polar_a_ = acos(lgnd_node_); dpolar_a = mean(diff(sort(polar_a_)));
weight_ = lgnd_weight_;
end;%if str_T_vs_L == 'L';

n_azimu_b_max = 3+round(2*pi*r/max(1e-12,d));
n_azimu_b_ = zeros(n_polar_a,1);
n_all = 0 ;
for npolar_a = 1:n_polar_a;
polar_a = polar_a_(npolar_a);
spolar_a = sin(polar_a);
if flag_uniform_over_polar_a==0; n_azimu_b = 3+round(2*pi*spolar_a*r/max(1e-12,d)); end;
if flag_uniform_over_polar_a==1; n_azimu_b = n_azimu_b_max; end;
if (str_T_vs_L == 'C2'); n_azimu_b = n_azimu_b + mod(n_azimu_b,2); end;
if (str_T_vs_L == 'C2') & (1-abs(cos(polar_a))<1e-16); n_azimu_b = 1; end;
n_azimu_b_(npolar_a) = n_azimu_b;
n_all = n_all + n_azimu_b;
end;%for npolar_a = 1:n_polar_a;
azimu_b_all_ = zeros(n_all,1); polar_a_all_ = zeros(n_all,1); weight_all_ = zeros(n_all,1); ix=0;
for npolar_a = 1:n_polar_a;
polar_a = polar_a_(npolar_a);
spolar_a = sin(polar_a);
n_azimu_b = n_azimu_b_(npolar_a);
azimu_b_ = linspace(0,2*pi,n_azimu_b+1); azimu_b_ = azimu_b_(1:end-1); dazimu_b = mean(diff(azimu_b_));
if ~isfinite(dazimu_b); dazimu_b = 2*pi; end;
if n_azimu_b==1; dazimu_b = 2*pi; end;
azimu_b_all_(1+ix+(0:n_azimu_b-1)) = transpose(azimu_b_);
polar_a_all_(1+ix+(0:n_azimu_b-1)) = polar_a*ones(n_azimu_b,1);
weight_all_(1+ix+(0:n_azimu_b-1)) = r.^2 .* weight_(npolar_a) .* dazimu_b .* ones(n_azimu_b,1);
%weight_all_(1+ix+(0:n_azimu_b-1)) = r.^2 .* dpolar_a .* sin(polar_a)*dazimu_b .* ones(n_azimu_b,1);
ix = ix+n_azimu_b;
end;%for npolar_a = 1:n_polar_a;

if (nargout>4);
k_c_0_all_ = r.*cos(azimu_b_all_).*sin(polar_a_all_);
k_c_1_all_ = r.*sin(azimu_b_all_).*sin(polar_a_all_);
k_c_2_all_ = r.*cos(polar_a_all_);
end;%if (nargout>4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

function output = f_test(polar_a) ;
m=6;
output = exp(-(-m*cos(polar_a)).^2)*m;
%output = ones(size(polar_a));

function output = F_test(polar_a) ;
m=6;
output = sqrt(pi)/2*erf(-m*cos(polar_a));
%output = -cos(polar_a);

function output = g_test(azimu_b) ;
w1 = 3; w2 = 5;
output = (cos(w1*azimu_b) + 1) + (cos(w2*azimu_b) + 1);
%output = ones(size(azimu_b));

function output= G_test(azimu_b) ;
w1 = 3; w2 = 5;
output = (sin(w1*azimu_b)/w1 + azimu_b) + (sin(w2*azimu_b)/w2 + azimu_b);
%output = azimu_b;
