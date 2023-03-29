function [fk,ier] = xxnufft1d3(nj,xj,fc,iflag,eps,nk,sk);

if (nargin<1);
rng(0);
nj = 59*5;
xj = pi*(2*rand(nj,1)-1);
fc = randn(nj,1) + i*randn(nj,1);
iflag = +1;
eps = 1e-12;
nk = 63*4;
sk = rand(nk,1);
fk0_ = finufft1d3(xj,fc,iflag,eps,sk);
fk1_ = nufft1d3(nj,xj,fc,iflag,eps,nk,sk);
disp(sprintf(' %% fk0_ vs fk1_: %0.16f',fnorm(fk0_-fk1_)/fnorm(fk0_)));
disp('returning'); return;
end;%if (nargin<1);

try;
ier=0; fk = finufft1d3(xj,fc,iflag,eps,sk);
catch;
[fk,ier] = nufft1d2(nj,xj,fc,iflag,eps,nk,sk);
end;%try;
