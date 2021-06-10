function [Et_,El_,Er_] = test_quad_3(a,c);

%K_ = 1:64;
K_ = round(2.^[0:0.125:10]);
Et_ = zeros(length(K_),1);
El_ = zeros(length(K_),1);
Er_ = zeros(length(K_),1);

for nK=1:length(K_);
K = K_(nK);
[wt,xt] = tschintmat(K,a,c);
ft = zeros(K,1); for index=1:length(xt); ft(index)=feval(@f,xt(index)); end;
It = wt*ft;
[xl,wl] = lgwt(K,a,c); wl = transpose(wl);
fl = zeros(K,1); for index=1:length(xl); fl(index)=feval(@f,xl(index)); end;
Il = wl*fl;
wr = ones(1,K)/K * (c-a) ; xr = linspace(a,c,K);
fr = zeros(K,1); for index=1:length(xr); fr(index)=feval(@f,xr(index)); end;
Ir = wr*fr;
Ix = feval(@F,c) - feval(@F,a);
Et_(nK) = abs(Ix-It); El_(nK) = abs(Ix-Il); Er_(nK) = abs(Ix-Ir);
end;%for nK=1:length(K_);

x_ = linspace(a,c,1024*2); y_ = feval(@f,x_);
subplot(1,2,1); plot(x_,y_); xlabel('x'); ylabel('f'); xlim([a,c]);
subplot(1,2,2);
plot(log2(K_),log(Et_),'b.-',log2(K_),log(El_),'r.-',log2(K_),log(Er_),'k.-'); 
xlabel('log2(K)'); ylabel('log(E)'); legend('T','L');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

function output = f(x);
m = 1.2; w = 3;
output = exp(-(m*x+m*sin(w*x)).^2).*(m + m*w*cos(w*x));

function output = F(x) ;
% integral of f(x);
m = 1.2; w = 3;
output = sqrt(pi)/2*erf(m*x+m*sin(w*x));

