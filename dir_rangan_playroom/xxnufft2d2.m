function [fj,ier] = xxnufft2d2(nj,xj,yj,iflag,eps,ms,mt,fk);

if (nargin<1);
rng(0);
nj = 59*5;
xj = pi*(2*rand(nj,1)-1);
yj = pi*(2*rand(nj,1)-1);
iflag = +1;
eps = 1e-12;
ms = 17*13;
mt = 19*11;
fk = rand(ms,mt);
fj0_ = finufft2d2(xj,yj,iflag,eps,fk);
fj1_ = nufft2d2(nj,xj,yj,iflag,eps,ms,mt,fk);
disp(sprintf(' %% fj0_ vs fj1_: %0.16f',fnorm(fj0_-fj1_)/fnorm(fj0_)));
disp('returning'); return;
end;%if (nargin<1);

try;
ier=0; fj = finufft2d2(xj,yj,iflag,eps,fk);
catch;
[fj,ier] = nufft2d2(nj,xj,yj,iflag,eps,ms,mt,fk);
end;%try;
