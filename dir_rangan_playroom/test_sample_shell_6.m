function ...
[ ...
 parameter ...
] = ...
test_sample_shell_6( ...
 parameter ...
);

str_thisfunction = 'test_sample_shell_6';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
disp(sprintf(' %% testing quadrature with r^2 weight: '));
k2p_max = 7; %<-- plays the role of a maximum 2*pi*k (i.e., bandlimit k). ;
delta = 1.2; %<-- plays the role of a displacement-vector magnitude for plane-wave of the form: exp(+i*2*pi*dot(k_,delta_)). ;
Ix = 4*pi*(1/delta^3)*( sin(k2p_max*delta) - (k2p_max*delta)*cos(k2p_max*delta) ); %<-- integral of plane-wave under bandlimit k, with domain normalized to have volume 4*pi/3. ;
b_K2P = 6; [b_jx_,b_jw_] = jacpts(b_K2P,0,2); k2p_jx_ = (b_jx_+1.0)*k2p_max/2; k2p_jw_ = b_jw_*(k2p_max/2)^3;
Ix_ = zeros(b_K2P,1); for nk2p=0:b_K2P-1; k2p = k2p_jx_(1+nk2p); Ix_(1+nk2p) = 4*pi*(1/delta)*sin(delta*k2p).*k2p./k2p.^2; end;% for nk2p=0:b_K2P-1; %<-- function sampled on gauss_jacobi nodes. ;
Ij = k2p_jw_*Ix_; disp(sprintf(' %% (Ix-Ij)/Ix = %0.16f/%0.16f = %0.16f',Ix-Ij,Ix,(Ix-Ij)/Ix)); %<-- gauss_jacobi errrel. ;
[b_lx_,b_lw_] = legpts(b_K2P); k2p_lx_ = (b_lx_+1.0)*k2p_max/2; k2p_lw_ = b_lw_*k2p_max/2;
Ix_ = zeros(b_K2P,1); for nk2p=0:b_K2P-1; k2p = k2p_lx_(1+nk2p); Ix_(1+nk2p) = 4*pi*(1/delta)*sin(delta*k2p).*k2p; end;% for nk2p=0:b_K2P-1; %<-- function sampled on gauss_legendre nodes. ;
Il = k2p_lw_*Ix_; disp(sprintf(' %% (Ix-Il)/Ix = %0.16f/%0.16f = %0.16f',Ix-Il,Ix,(Ix-Il)/Ix)); %<-- gauss_legendre errrel. ;
%%%%%%%%;
disp(sprintf(' %% testing basic quadrature on sphere: '));
k_eq_d_ = 0.5.^[-2:0.125:4]; n_k_eq_d = numel(k_eq_d_);
for flag_uniform_over_polar_a=[0,1];
E_df__ = zeros(n_k_eq_d,2);
for nk_eq_d = 0:n_k_eq_d-1;
k_eq_d = k_eq_d_(1+nk_eq_d);
k2p_max = 7.5; delta = 2.32; 
b_K2P = 1+ceil(k2p_max/2/k_eq_d); [b_jx_,b_jw_] = jacpts(b_K2P,0,2);
k2p_ = (b_jx_+1.0)*k2p_max/2; k2p_w_ = b_jw_*(k2p_max/2)^3;
Iq_ = zeros(b_K2P,1); Ix_ = zeros(b_K2P,1);
[n_q,azimu_b_q_,polar_a_q_,weight_shell_q_] = sample_shell_6(1,k_eq_d,'L',flag_uniform_over_polar_a) ;
for nk2p=0:b_K2P-1;
k2p = k2p_(1+nk2p);
Iq_(1+nk2p) = sum(cos(k2p*delta*cos(polar_a_q_)).*weight_shell_q_);
Ix_(1+nk2p) = 4*pi*(1/delta)*sin(delta*k2p).*k2p./k2p.^2;
end;% for nk2p=0:b_K2P-1;
Iq = k2p_w_*Iq_;
Ix = 4*pi*(1/delta^3)*( sin(k2p_max*delta) - (k2p_max*delta)*cos(k2p_max*delta) );
disp(sprintf(' %% flag_uniform_over_polar_a %d k_eq_d = %0.6f n_q %.4d b_K2P %.2d (Ix-Iq)/Ix = %+0.16f/%+0.16f = %+0.16f',flag_uniform_over_polar_a,k_eq_d,n_q,b_K2P,Ix-Iq,Ix,(Ix-Iq)/Ix));
E_df__(1+nk_eq_d,1+flag_uniform_over_polar_a) = abs((Ix-Iq)/max(1e-12,Ix));
end;%for nk_eq_d = 0:n_k_eq_d-1;
end;%for flag_uniform_over_polar_a=[0,1];
%%%%%%%%;
disp(sprintf(' %% testing different equatorial distances: '));
k_eq_d_ = 2.^[-0.5:-0.25:-2.0]; n_k_eq_d = numel(k_eq_d_);
for flag_uniform_over_polar_a=[0,1];
for nk_eq_d=0:n_k_eq_d-1;
k_p_r_max=3.0;
k_eq_d = k_eq_d_(1+nk_eq_d); n_equator = 3+round(2*pi*k_p_r_max/k_eq_d); n_polar_a = 3+round(n_equator/2);
disp(sprintf(' %% sampling at equatorial_distance k_eq_d %0.3f (n_polar_a %d)',k_eq_d,n_polar_a));
%%%%;
[n_q,azimu_b_q_,polar_a_q_,weight_shell_q_] = sample_shell_6(k_p_r_max,k_eq_d,'T',flag_uniform_over_polar_a) ;
f_q_ = feval(@f_test,polar_a_q_); g_q_ = feval(@g_test,azimu_b_q_); h_q_ = f_q_.*g_q_;
It = dot(weight_shell_q_,h_q_);
[n_q,azimu_b_q_,polar_a_q_,weight_shell_q_] = sample_shell_6(k_p_r_max,k_eq_d,'C',flag_uniform_over_polar_a) ;
f_q_ = feval(@f_test,polar_a_q_); g_q_ = feval(@g_test,azimu_b_q_); h_q_ = f_q_.*g_q_;
Ic = dot(weight_shell_q_,h_q_);
[n_q,azimu_b_q_,polar_a_q_,weight_shell_q_] = sample_shell_6(k_p_r_max,k_eq_d,'C2',flag_uniform_over_polar_a) ;
f_q_ = feval(@f_test,polar_a_q_); g_q_ = feval(@g_test,azimu_b_q_); h_q_ = f_q_.*g_q_;
Id = dot(weight_shell_q_,h_q_);
[n_q,azimu_b_q_,polar_a_q_,weight_shell_q_] = sample_shell_6(k_p_r_max,k_eq_d,'L',flag_uniform_over_polar_a) ;
f_q_ = feval(@f_test,polar_a_q_); g_q_ = feval(@g_test,azimu_b_q_); h_q_ = f_q_.*g_q_;
Il = dot(weight_shell_q_,h_q_);
%%%%;
Ix = k_p_r_max.^2 .* (F_test(pi)-F_test(0)) .* (G_test(2*pi)-G_test(0));
et = norm(It-Ix); let = -log10(et);
ec = norm(Ic-Ix); lec = -log10(ec);
ed = norm(Id-Ix); led = -log10(ed);
el = norm(Il-Ix); lel = -log10(el);
disp(sprintf(' %% flag_uniform_over_polar_a %d n_q %d Ix %0.3f It %0.3f error %0.16f (%0.2f) Ic %0.3f error %0.16f (%0.2f) Id %0.3f error %0.16f (%0.2f) Il %0.3f error %0.16f (%0.2f)',flag_uniform_over_polar_a,n_q,Ix,It,et,let,Ic,ec,lec,Id,ed,led,Il,el,lel));
end;%for nk_eq_d=0:n_k_eq_d-1;
end;%for flag_uniform_over_polar_a=[0,1];
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

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
