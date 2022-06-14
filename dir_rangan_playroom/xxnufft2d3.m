function [fk,ier] = xxnufft2d3(nj,xj,yj,cj,iflag,eps,nk,sk,tk);

if (nargin<1);
rng(0);
nj = 59*11;
xj = pi*(2*rand(nj,1)-1);
yj = pi*(2*rand(nj,1)-1);
cj = randn(nj,1) + i*randn(nj,1);
iflag = +1;
eps = 1e-12;
nk = 31*13;
sk = pi*(2*rand(nk,1)-1);
tk = pi*(2*rand(nk,1)-1);
fk0_ = finufft2d3(xj,yj,cj,iflag,eps,sk,tk) ;
fk1_ = nufft2d3(nj,xj,yj,cj,iflag,eps,nk,sk,tk) ;
disp(sprintf(' %% fk0_ vs fk1_: %0.16f',fnorm(fk0_-fk1_)/fnorm(fk0_)));
disp('returning'); return;
end;%if (nargin<1);

try;
ier=0; fk = finufft2d3(xj,yj,cj,iflag,eps,sk,tk) ;
catch;
[fk,ier] = nufft2d3(nj,xj,yj,cj,iflag,eps,nk,sk,tk) ;
end;%try;
