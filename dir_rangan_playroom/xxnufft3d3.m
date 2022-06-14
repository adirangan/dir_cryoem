function [fk,ier] = xxnufft3d3(nj,xj,yj,zj,cj,iflag,eps,nk,sk,tk,uk);

if (nargin<1);
rng(0);
nj = 59;
xj = pi*(2*rand(nj,1)-1);
yj = pi*(2*rand(nj,1)-1);
zj = pi*(2*rand(nj,1)-1);
cj = randn(nj,1) + i*randn(nj,1);
iflag = +1;
eps = 1e-12;
nk = 17;
sk = pi*(2*rand(nk,1)-1);
tk = pi*(2*rand(nk,1)-1);
uk = pi*(2*rand(nk,1)-1);
fk0_ = finufft3d3(xj,yj,zj,cj,iflag,eps,sk,tk,uk) ;
fk1_ = nufft3d3(nj,xj,yj,zj,cj,iflag,eps,nk,sk,tk,uk) ;
disp(sprintf(' %% fk0_ vs fk1_: %0.16f',fnorm(fk0_-fk1_)/fnorm(fk0_)));
disp('returning'); return;
end;%if (nargin<1);

try;
ier=0; fk = finufft3d3(xj,yj,zj,cj,iflag,eps,sk,tk,uk) ;
catch;
[fk,ier] = nufft3d3(nj,xj,yj,zj,cj,iflag,eps,nk,sk,tk,uk) ;
end;%try;
