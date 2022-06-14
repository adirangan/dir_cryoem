function [fk,ier] = xxnufft1d1(nj,xj,cj,iflag,eps,ms);

if (nargin<1);
rng(0);
nj = 59*5;
xj = pi*(2*rand(nj,1)-1);
cj = randn(nj,1) + i*randn(nj,1);
iflag = +1;
eps = 1e-12;
ms = 17*7;
fk0_ = finufft1d1(xj,cj,iflag,eps,ms) / nj ;
fk1_ = nufft1d1(nj,xj,cj,iflag,eps,ms) ;
disp(sprintf(' %% fk0_ vs fk1_: %0.16f',fnorm(fk0_-fk1_)/fnorm(fk0_)));
disp('returning'); return;
end;%if (nargin<1);

try;
ier=0; fk = finufft1d1(xj,cj,iflag,eps,ms) / nj ;
catch;
[fk,ier] = nufft1d1(nj,xj,cj,iflag,eps,ms) ;
end;%try;
