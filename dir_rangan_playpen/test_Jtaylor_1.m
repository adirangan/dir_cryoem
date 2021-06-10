% Testing taylor-function expansion T_ of F_ = exp(i*z*cos(2*pi*j/N - w)) ; 
% Also testing bessel-function expansion J_ of F_ (note that this does not naturally factorize) ;
% Determining what order (and number of terms) is required to reach accuracy eps out to some fixed z=z_max ;
% Typical values of eps are 0.1 and 0.01 ;
% Typical values of z_max are 2*pi (1 pixel) and 4*pi (2 pixels) ;

% First generate taylor-expansion ;
l_max = 32; % maximum order
n_T = round((l_max+1)*(l_max+2)/2);
f_ = zeros(n_T,1);
p_ = zeros(n_T,1);
b_ = zeros(n_T,1);
B_ = zeros(n_T,1);
l_ = zeros(n_T,1);
j_ = zeros(n_T,1);
%l_factorial = 1;
ci = -i/2; f = 1; ic=1;
for l=0:l_max;
%l_start = round(l*(l+1)/2);
%l_factorial = l_factorial * max(1,l);
if (l>0); f = f*ci/max(1,l); end;
j_factorial = 1;
lmj_factorial = 1;
for j=0:l;
j_factorial = j_factorial * max(1,j);
%f_(1 + l_start + j) = (-i/2)^l / l_factorial;
%f_(1 + l_start + j) = f;
%p_(1 + l_start + j) = -l + 2*j;
%b_(1 + l_start + j) = lmj_factorial / j_factorial;
%B_(1 + l_start + j) = nchoosek(l,j);
f_(ic) = f;
p_(ic) = -l + 2*j;
b_(ic) = lmj_factorial / j_factorial;
B_(ic) = nchoosek(l,j);
lmj_factorial = lmj_factorial * max(1,l-j);
%l_(1 + l_start + j) = l;
%j_(1 + l_start + j) = j;
l_(ic) = l;
j_(ic) = j;
ic = ic+1;
end;%for j=0:l;
end;%for l=0:l_max;

n_z = 4; z_ = [1.0 , pi , 2*pi , 4*pi];
%n_eps = 4; eps_ = 0.1.^[0.5 , 1.0 , 1.5 , 2.0];
n_eps = 4; eps_ = 0.1.^[1 , 2 , 3 , 4];
T_req_ = zeros(n_z,n_eps);
J_req_ = zeros(n_z,n_eps);

n_phi = 1024;
phi_ = 2*pi*(0:n_phi-1)/n_phi;

for nz=1:n_z;
z = z_(nz);
for neps = 1:n_eps;
eps = eps_(neps);

% Calculate function F_ ;
F_ = exp(-i*z*cos(phi_));
%F_ = exp(+i*z*cos(phi_));

% Calculate bessel expansion J_ ; 
J_ = zeros(size(phi_));
EJ_ = zeros(1+l_max,1);
jc=0; l=0; continue_flag=1;
while (l<=l_max & continue_flag);
if (l==0); J_ = J_ + besselj(0,z); jc = jc+1; end;
%if (l>0); J_ = J_ + exp(+i*l*pi/2)*(-1)^l*besselj(-l,z)*exp(+i*l*phi_) + exp(-i*l*pi/2)*(-1)^l*besselj(+l,z)*exp(-i*l*phi_); jc = jc+2; end;
%if (l>0); J_ = J_ + exp(+i*l*pi/2)*besselj(-l,z)*exp(+i*l*phi_) + exp(-i*l*pi/2)*besselj(+l,z)*exp(-i*l*phi_); jc = jc+2; end;
if (l>0); J_ = J_ + exp(+i*l*pi/2)*besselj(-l,z)*exp(-i*l*phi_) + exp(-i*l*pi/2)*besselj(+l,z)*exp(+i*l*phi_); jc = jc+2; end;
EJ_(1+l) = mean(abs(F_-J_).^2);
%disp(sprintf(' %% z %0.2f l %d EJ %0.4f',z,l,EJ_(l+1)));
if (EJ_(l+1)<eps); continue_flag=0; J_req_(nz,neps) = l; disp(sprintf('z_(%d) = %0.2f eps_(%d) = %0.3f J_req = %d(%d)',nz,z,neps,eps,J_req_(nz,neps),jc)); end;
l = l+1;
end;%while;

% Calculate taylor expansion T_ ;
T_ = zeros(size(phi_));
ET_ = zeros(1+l_max,1);
jc=0; l=0; continue_flag=1;
while (l<=l_max & continue_flag)
l_start = round(l*(l+1)/2);
for j=0:l;
f = f_(1 + l_start + j);
p = p_(1 + l_start + j);
b = b_(1 + l_start + j);
T_ = T_ + f*b*exp(i*p*phi_)*z^l;
jc = jc+1;
end;%for j=0:l;
ET_(1+l) = mean(abs(F_-T_).^2);
%disp(sprintf(' %% z %0.2f l %d ET %0.4f',z,l,ET_(l+1)));
if (ET_(l+1)<eps); continue_flag=0; T_req_(nz,neps) = l; disp(sprintf('z_(%d) = %0.2f eps_(%d) = %0.3f T_req = %d(%d)',nz,z,neps,eps,T_req_(nz,neps),jc)); end;
l = l+1;
end;%while;

if neps==n_eps;
subplot(2,n_z,0*n_z + nz); hold on;
plot(phi_,real(F_),'k-');
plot(phi_,real(J_),'r-');
plot(phi_,real(T_),'m-');
xlim([min(phi_),max(phi_)]);
subplot(2,n_z,1*n_z + nz); hold on;
plot(phi_,imag(F_),'k-');
plot(phi_,imag(J_),'b-');
plot(phi_,imag(T_),'g-');
xlim([min(phi_),max(phi_)]);
end;%if neps==n_eps;

end;%for neps = 1:n_eps;
end;%for nz=1:n_z;
