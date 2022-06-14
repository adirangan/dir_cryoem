function [grad_,A_x_approx_] = grad_from_1d_periodic_0(x_0_,A_x_,eps);
%%%%%%%%;
% calculates gradient from 1d array. ;
%%%%%%%%;

if (nargin<1);
n_row = 164;
x_0_ = transpose(linspace(0,1,n_row));
A_x_ = sin(2*pi*x_0_);
B_x_ = cos(2*pi*x_0_)*2*pi;
[D_x_,A_x_approx_] = grad_from_1d_periodic_0(x_0_,A_x_);
D_x_ = real(D_x_);
A_x_approx_ = real(A_x_approx_);
figure(1);clf;figmed;
subplot(1,2,1);plot(x_0_,A_x_,'k.-',x_0_,A_x_approx_,'ko-',x_0_,B_x_,'r.-');
xlabel('x');ylabel('f');legend({'A','A approx','dA'});
subplot(1,2,2);plot(x_0_,B_x_,'r.-',x_0_,D_x_,'ro-');
xlabel('x');ylabel('f');legend({'dA','dA approx'});
disp('returning'); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); x_0_=[]; end; na=na+1;
if (nargin<1+na); A_x_=[]; end; na=na+1;
if (nargin<1+na); eps=[]; end; na=na+1;

[n_row] = numel(A_x_);
if isempty(x_0_); x_0_ = 1:n_row; end;
if isempty(eps); eps = 1e-12; end;

x_min = min(x_0_);
x_max = max(x_0_);
x_lim_ = [x_min,x_max];
x_scale_ = 2*pi*(x_0_ - mean(x_lim_))/diff(x_lim_);

n_K = n_row;
n_K_half = floor(n_K/2);
k_scale_ = -n_K_half + transpose([0:n_K-1]);

A_k_ = xxnufft1d1(n_row,x_scale_,A_x_,-1,eps,n_K);
if (nargout>1); A_x_approx_ = xxnufft1d2(n_row,x_scale_,+1,eps,n_K,A_k_); end;
B_k_ = A_k_.*i.*k_scale_;
C_x_ = xxnufft1d2(n_row,x_scale_,+1,eps,n_K,B_k_);
D_x_ = C_x_*2*pi/diff(x_lim_);
grad_ = D_x_;

