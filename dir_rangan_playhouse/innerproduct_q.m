function x_ = innerproduct_q(a_q_,b_q_);
% calculates the innerproduct: ;
% x_(omega) = \int_{gamma=0}^{gamma=2*pi} conj(a_p_(gamma)) * b_p_(gamma+omega) dgamma ;
% assuming that a_p_ and b_p_ are given on equispaced points on a ring. ;
% omega is assumed to be equispaced on the ring as well. ;
% The values of omega are given by: ;
% omega_ = linspace(0,2*pi,length(a_)+1) ; omega_ = omega_(1:end-1); 

if nargin<2;
n_g = 28; rng(1);
a_p_ = randn(n_g,1) + i*randn(n_g,1);
b_p_ = randn(n_g,1) + i*randn(n_g,1);
x0_ = zeros(n_g,1);
bb_p_ = [b_p_;b_p_];
for ng=0:n_g-1;
x0_(1+ng) = sum(conj(a_p_(1:n_g)).*bb_p_([1:n_g] + ng))*(2*pi) / n_g;
end;%for ng=0:n_g-1;
a_q_ = fft(a_p_); b_q_ = fft(b_p_);
x1_ = innerproduct_q(a_q_,b_q_);
disp(x0_); disp(x1_); disp(x0_./x1_);
disp('returning');return;
end;%if nargin<2;

c_q_ = conj(a_q_).*b_q_;
x_ = ifft(c_q_)*2*pi/length(a_q_);
