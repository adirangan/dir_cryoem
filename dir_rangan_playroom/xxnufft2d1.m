function [fk,ier] = xxnufft2d1(nj,xj,yj,cj,iflag,eps,ms,mt);

if (nargin<1);
rng(0);
nj = 59*5;
xj = pi*(2*rand(nj,1)-1);
yj = pi*(2*rand(nj,1)-1);
cj = randn(nj,1) + i*randn(nj,1);
iflag = +1;
eps = 1e-12;
ms = 17*7;
mt = 19*11;
fk0_ = finufft2d1(xj,yj,cj,iflag,eps,ms,mt) / nj ;
fk1_ = nufft2d1(nj,xj,yj,cj,iflag,eps,ms,mt) ;
disp(sprintf(' %% fk0_ vs fk1_: %0.16f',fnorm(fk0_-fk1_)/fnorm(fk0_)));
disp('returning'); return;
end;%if (nargin<1);

try;
ier=0; fk = finufft2d1(xj,yj,cj,iflag,eps,ms,mt) / nj ;
catch;
[fk,ier] = nufft2d1(nj,xj,yj,cj,iflag,eps,ms,mt) ;
end;%try;
