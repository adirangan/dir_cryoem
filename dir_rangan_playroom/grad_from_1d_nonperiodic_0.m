function [grad_] = grad_from_1d_nonperiodic_0(x_0_,A_);
%%%%%%%%;
% calculates gradient from 1d array. ;
% assumes that the input array x_0_ is sorted. ;
% if A_ is an array, acts on each column of A_. ;
%%%%%%%%;

if (nargin<1);
n_row = 38;
x_0_ = -pi/2 + 0.4*transpose(linspace(0,1,n_row));
A_ = [ sin(2*pi*x_0_) , sin(2*pi*2*x_0_) , sin(2*pi*3*x_0_) ];
B_ = 2*pi*[ cos(2*pi*x_0_) , 2*cos(2*pi*2*x_0_) , 3*cos(2*pi*3*x_0_) ];
D_ = grad_from_1d_nonperiodic_0(x_0_,A_);
figure(1);clf;figmed;
subplot(1,2,1);plot(x_0_,A_,'k.-',x_0_,B_,'r.-');
xlabel('x');ylabel('f');legend({'A','dA'});
subplot(1,2,2);plot(x_0_,B_,'r.-',x_0_,D_,'ro-');
xlabel('x');ylabel('f');legend({'dA','dA approx'});
disp('returning'); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); x_0_=[]; end; na=na+1;
if (nargin<1+na); A_=[]; end; na=na+1;
if (nargin<1+na); eps=[]; end; na=na+1;

[n_row] = size(A_,1);
if isempty(x_0_); x_0_ = 1:n_row; end;

x_0_ = reshape(x_0_,n_row,1);
grad_0in_ = bsxfun(@rdivide,A_(3:end,:) - A_(1:end-2,:),x_0_(3:end) - x_0_(1:end-2));
grad_start = (A_(2,:) - A_(1,:))/(x_0_(2) - x_0_(1));
grad_final = (A_(end,:) - A_(end-1,:))/(x_0_(end) - x_0_(end-1));
grad_ = [grad_start;grad_0in_;grad_final];
