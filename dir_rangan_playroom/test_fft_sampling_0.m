N=16;
% set up imaginary part to obey constraint. ;
f_=0.5*randn(N,N) + i*0.5*randn(N,N);
kA0_=zeros(N,N); kA1_=zeros(N,N);
kB0_=zeros(N,N); kB1_=zeros(N,N);
kC0_=zeros(N,N); kC1_=zeros(N,N);
kD0_=zeros(N,N); kD1_=zeros(N,N);
flag_conj_=zeros(N,N);
flag_real_=zeros(N,N);
for kA0=0:N-1; for kA1=0:N-1;
kB0=kA0; if kB0>=N/2; kB0=kB0-N; end;
kB1=kA1; if kB1>=N/2; kB1=kB1-N; end;
kC0=-kB0;
kC1=-kB1;
kD0=kC0; if kD0< 0; kD0=kD0+N; end;
kD1=kC1; if kD1< 0; kD1=kD1+N; end;
kA0_(1+kA0,1+kA1) = kA0;
kA1_(1+kA0,1+kA1) = kA1;
kB0_(1+kA0,1+kA1) = kB0;
kB1_(1+kA0,1+kA1) = kB1;
kC0_(1+kA0,1+kA1) = kC0;
kC1_(1+kA0,1+kA1) = kC1;
kD0_(1+kA0,1+kA1) = kD0;
kD1_(1+kA0,1+kA1) = kD1;
if kB0<=0; %<-- lower half-plane. ;
flag_conj_(1+kD0,1+kD1) = -1;
flag_conj_(1+kA0,1+kA1) = +1;
f_(1+kD0,1+kD1) = conj(f_(1+kA0,1+kA1));
end;%if kB0<=0;
if ( (kD0==kA0) & (kD1==kA1) );
flag_real_(1+kD0,1+kD1) = 1;
f_(1+kD0,1+kD1) = 1.0*randn(1,1);
end;%if ( (kD0==kA0) & (kD1==kA1) );
end;end;%for kA0=0:N-1; for kA1=0:N-1;
% now take fourier-transform. ;
x_=fft2(f_)/N;
disp(sprintf(' %% fnorm(imag(x_)): %0.16f',fnorm(imag(x_))));

%{
d1 = 5; d2 =7;
rng(1);
C = reshape(ceil(8*rand(d1*d2,1)) + i*ceil(8*rand(d1*d2,1)),d1,d2);
B = fft2(C);
k1_ = linspace(0,2*pi,d1+1); k1_ = k1_(1:end-1);
k2_ = linspace(0,2*pi,d2+1); k2_ = k2_(1:end-1);
[K2_,K1_] = meshgrid(k2_,k1_);
A = nufft2d2(d1*d2,K1_(:),K2_(:),-1,1e-12,d1,d2,decenter2(C));
A = reshape(A,d1,d2);
%}

%{
res = 64; x_ = 0:res-1; k_ = 0:res-1;
sigma = 2.5;
f_x_ = exp(-(x_).^2/2/sigma^2) + exp(-(res-x_).^2/2/sigma^2);
subplot(2,3,1); plot(x_,f_x_,'o-'); xlim([0,res-1]);
f_k_ = fft(f_x_);
subplot(2,3,2); plot(k_,real(f_k_),'ro-',k_,imag(f_k_),'bo-'); xlim([0,res-1]);
k2_ = 0:0.25:res-1;
f_k2_ = interp1(k_,f_k_,k2_);
subplot(2,3,3); plot(k2_,real(f_k2_),'ro-',k2_,imag(f_k2_),'bo-'); xlim([0,res-1]);
f_x2_ = ifft(f_k2_);
x2_ = 0:4*(res-1); 
subplot(2,3,4); plot(x2_,real(f_x2_),'o-'); xlim([0,res-1]);
nres2 = length(k2_); res2 = 5*res;
f_k3_ = recenter1([zeros(1,nres2*2),recenter1(f_k2_),zeros(1,nres2*2)]);
nres3 = length(f_k3_);
k3_ = linspace(0,5*(res-1),nres3);
subplot(2,3,5); plot(k3_,real(f_k3_),'ro-',k3_,imag(f_k3_),'bo-'); xlim([0,5*(res-1)]);
f_x3_ = ifft(f_k3_);
%x3_ = (0:max(k3_)/mean(diff(k3_)))*mean(diff(k3_));
x3_ = linspace(0,max(k3_),1+ceil(max(k3_)/mean(diff(k3_))));
subplot(2,3,6); plot(x3_,real(f_x3_),'o-'); xlim([0,max(k3_)]);
 %}


%{
res = 512*4;
sigma = 1/32;
x_ = linspace(-1,1,res);
f_x_ = 1/sqrt(2*pi)/sigma*exp(-x_.^2/2/sigma^2);
subplot(2,4,1);plot(x_,f_x_,'.-');
f_k_ = fft(recenter1(f_x_))/res;
subplot(2,4,2);
k_ = pi*((recenter1(0:res-1))-res/2);
plot(k_,real(f_k_),'r.-',k_,imag(f_k_),'b.-'); xlim([min(k_),max(k_)]);
subplot(2,4,3);
g_k_ = 1/sqrt(2*pi)*sigma*exp(-k_.^2/2*sigma^2);
plot(k_,real(g_k_),'r.-',k_,imag(g_k_),'b.-'); xlim([min(k_),max(k_)]);
subplot(2,4,4);
g_x_ = recenter1(fft(g_k_));
plot(x_,real(g_x_),'r.-',x_,imag(g_x_),'b.-'); xlim([min(x_),max(x_)]);
subplot(2,4,5);
res2 = 1024*4;
k2_ = pi*linspace(-200,+200,res2);
g_k2_ = 1/sqrt(2*pi)*sigma*exp(-k2_.^2/2*sigma^2);
plot(k2_,real(g_k2_),'r.-',k2_,imag(g_k2_),'b.-'); xlim([min(k2_),max(k2_)]);
subplot(2,4,6);
g_x2_ = recenter1(fft(recenter1(g_k2_)));
x2_ = linspace(-1,1,res2);
plot(x2_,real(g_x2_),'r.-',x2_,imag(g_x2_),'b.-'); xlim([min(x2_),max(x2_)]);
 %}

