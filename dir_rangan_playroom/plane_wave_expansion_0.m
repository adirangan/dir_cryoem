function a_Y_ = plane_wave_expansion_0(k_,x_,l_max);
% plane-wave expansion in terms of spherical-bessel-functions and spherical-harmonics. ;
% we assume plane-wave is of the form exp( +i * < 2*pi*k_ , x_ > ). ;

na=0;
if (nargin<1+na); k_ = []; end; na=na+1;
if (nargin<1+na); x_ = []; end; na=na+1;
if (nargin<1+na); l_max = []; end; na=na+1;

if (isempty(k_)); k_ = [+1.2;+0.3;-5.6]; end; %<- wavenumber. ;
if (isempty(x_)); x_ = [-0.8;+1.3;+2.3]/3; end;  %<- space. ;
if (isempty(l_max)); l_max = 128; end; %<-- order. ;

verbose=0;

k = fnorm(k_); k_01 = fnorm(k_(1:2)); azimu_b_k = atan2(k_(2),k_(1)); polar_a_k = atan2(k_01,k_(3));
x = fnorm(x_); x_01 = fnorm(x_(1:2)); azimu_b_x = atan2(x_(2),x_(1)); polar_a_x = atan2(x_01,x_(3));
t = 2*pi*k*x;
Ix_ = exp( +i * dot(2*pi*k_,x_));
Ip_ = 0;
Ylm_k__ = get_Ylm__(1+l_max,0:l_max,1,azimu_b_k,polar_a_k,0);
Ylm_x__ = get_Ylm__(1+l_max,0:l_max,1,azimu_b_x,polar_a_x,0);
n_lm = (1+l_max)^2;
a_Y_ = zeros(n_lm,1);
na=0;
for l_val=0:l_max;
Ylm_d_k_ = Ylm_k__{1+l_val}(:);
Ylm_d_x_ = Ylm_x__{1+l_val}(:);
jl = besselj(l_val+0.5,t)*sqrt(pi/(2*t));
%Ip_ = Ip_ + sum(4*pi*(i^l_val)*jl*((Ylm_d_k_).*conj(Ylm_d_x_)));
Ip_ = Ip_ + sum(4*pi*(i^l_val)*jl*(conj(Ylm_d_k_).*(Ylm_d_x_)));
a_Y_(1+na+[0:1+2*l_val - 1]) = 4*pi*(i^l_val)*jl*conj(Ylm_d_x_);
na=na+1+2*l_val;
end;%for l_val=0:l_max;
assert(na==n_lm);
if (verbose); disp(sprintf('k: %f, x: %f, Ix: %f+%fi, l_val %d; Ip: %f+%fi, error: %0.16f',k,x,real(Ix_),imag(Ix_),l_val,real(Ip_),imag(Ip_),fnorm(Ix_-Ip_))); end;




