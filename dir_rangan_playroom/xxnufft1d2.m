function [fj,ier] = xxnufft1d2(nj,xj,iflag,eps,ms,fk);

if (nargin<1);
rng(0);
nj = 59*5;
xj = pi*(2*rand(nj,1)-1);
iflag = +1;
eps = 1e-12;
ms = 17*13;
fk = rand(ms,1);
fj0_ = finufft1d2(xj,iflag,eps,fk);
fj1_ = nufft1d2(nj,xj,iflag,eps,ms,fk);
disp(sprintf(' %% fj0_ vs fj1_: %0.16f',fnorm(fj0_-fj1_)/fnorm(fj0_)));
disp('returning'); return;
end;%if (nargin<1);

try;
ier=0; fj = finufft1d2(xj,iflag,eps,fk);
catch;
[fj,ier] = nufft1d2(nj,xj,iflag,eps,ms,fk);
end;%try;
