function grad_2rc___ = grad_from_2d_nonperiodic_0(x_0_,x_1_,A__);
%%%%%%%%;
% calculates gradient from 2d array. ;
%%%%%%%%;

if (nargin<1);
n_row = 38;
n_col = 43;
x_0_ = -0.8 + 0.6*transpose(linspace(0,1,n_row));
x_1_ = -0.3 + 0.3*linspace(0,1,n_col);
[x_0__,x_1__] = ndgrid(x_0_,x_1_);
A__ = [ sin(2*pi*x_0__) .* sin(2*pi*2*x_1__) ];
B0__ = [ 2*pi * cos(2*pi*x_0__) .* sin(2*pi*2*x_1__) ];
B1__ = [ 2*2*pi * sin(2*pi*x_0__) .* cos(2*pi*2*x_1__) ];
D___ = grad_from_2d_nonperiodic_0(x_0_,x_1_,A__);
figure(1);clf;figmed;np=0;
subplot(2,3,1+np);np=np+1; imagesc(B0__); colorbar;
subplot(2,3,1+np);np=np+1; imagesc(squeeze(D___(1+0,:,:))); colorbar;
subplot(2,3,1+np);np=np+1; imagesc(B0__ - squeeze(D___(1+0,:,:))); colorbar;
subplot(2,3,1+np);np=np+1; imagesc(B1__); colorbar;
subplot(2,3,1+np);np=np+1; imagesc(squeeze(D___(1+1,:,:))); colorbar;
subplot(2,3,1+np);np=np+1; imagesc(B1__ - squeeze(D___(1+1,:,:))); colorbar;
disp('returning'); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); x_0_=[]; end; na=na+1;
if (nargin<1+na); x_1_=[]; end; na=na+1;
if (nargin<1+na); A__=[]; end; na=na+1;

[n_row,n_col] = size(A__);
if isempty(x_0_); x_0_ = transpose(1:n_row); end;
if isempty(x_1_); x_1_ = 1:n_col; end;

grad_2rc___ = zeros(2,n_row,n_col);
grad_2rc___(1,:,:) = grad_from_1d_nonperiodic_0(x_0_,A__);
grad_2rc___(2,:,:) = transpose(grad_from_1d_nonperiodic_0(x_1_,transpose(A__)));


