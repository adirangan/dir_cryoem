function c_ = colormap_gaussian_2d(x0_,x1_,sigma,gamma);

na=0;
if (nargin<1+na); x0_=[]; end; na=na+1;
if (nargin<1+na); x1_=[]; end; na=na+1;
if (nargin<1+na); sigma=[]; end; na=na+1;
if (nargin<1+na); gamma=[]; end; na=na+1;

if nargin<2;
n_x = 129; dx = 2/(n_x-1); dx2 = dx/2;
x0_ = linspace(-1,+1,n_x);
x1_ = linspace(-1,+1,n_x);
[X0__,X1__] = ndgrid(x0_,x1_);
gamma_ = 0.05:0.05:1.0; n_gamma = numel(gamma_);
figure(1);clf;figbig;
for ngamma=0:n_gamma-1;
gamma = gamma_(1+ngamma);
subplot(4,5,1+ngamma);
c_ = colormap_gaussian_2d(X0__(:),X1__(:),[],gamma);
c__ = reshape(c_,[1,n_x^2,3]);
x0__ = transpose(repmat(X0__(:),[1,4]) + dx2*repmat([-1,+1,+1,-1],[n_x^2,1]));
x1__ = transpose(repmat(X1__(:),[1,4]) + dx2*repmat([-1,-1,+1,+1],[n_x^2,1]));
patch(x0__,x1__,c__,'EdgeColor','none');
title(sprintf('$\\gamma=%0.2f$',gamma),'Interpreter','latex');
axisnotick;
axis image;
end;%for ngamma=0:n_gamma-1;
disp('returning');return;
end;%if nargin<2;

if isempty(sigma); sigma = sqrt(2)*max(std(x1_),std(x1_,1)); end;
if isempty(gamma); gamma = 0.35; end;

x0_ = x0_(:);
x1_ = x1_(:);
W_ = atan2(x1_,x0_);
R_ = sqrt(x0_.^2 + x1_.^2);
azimu_b_ = W_;
polar_a_ = 2*atan2(R_/sigma,1);
c_ = colorpsphere_8points_gamma(polar_a_,azimu_b_,gamma);

  
