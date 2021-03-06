%function test_LSQ_vs_Sampling_0();
% example for LSQ_vs_sampling. ;

i = sqrt(-1);
ni=0;
rng(ni);
K_lo = 8;
K_hi = 16;
K_al_ = -K_hi:+K_hi;
K_lo_ = -K_lo:+K_lo;
K_hi_ = [-K_hi:-K_lo-1 , +K_lo+1:+K_hi];
K_lo_al_ = K_hi - K_lo + (1:2*K_lo+1);
K_hi_al_ = [1:K_hi-K_lo , 1+K_hi+K_lo+1:2*K_hi+1];

%x__s = linspace(-pi,pi,16*K_hi+1); x__s = x__s(2:end);
x__s = linspace(-pi,pi,16*K_hi);
F_al = zeros(length(x__s),2*K_hi+1);
F_al(:,1+K_hi) = ones(length(x__s),1);
tmp_e = exp(i*1*transpose(x__s));
for nk=1+K_hi+1:2*K_hi+1;
F_al(:,nk) = F_al(:,nk-1).*tmp_e;
end;%for nk=1+K_hi+1:2*K_hi+1;
for nk=1+K_hi-1:-1:1;
F_al(:,nk) = F_al(:,nk+1)./tmp_e;
end;%for nk=1+K_hi-1:-1:1;
F_al = F_al/sqrt(2*pi);
disp_flag=0;
if disp_flag;
subplot(1,3,1);imagesc(real(F_al));
subplot(1,3,2);imagesc(imag(F_al));
subplot(1,3,3);imagesc(real(ctranspose(F_al)*F_al));
end;%if disp_flag;
G_al = (F_al - repmat(exp(-i*pi*K_al_)/sqrt(2*pi),length(x__s),1))*diag(1./(i*K_al_));
G_al(:,1+K_hi) = (x__s + pi)/sqrt(2*pi);
disp_flag=0;
if disp_flag;
subplot(1,3,1);imagesc(real(G_al));
subplot(1,3,2);imagesc(imag(G_al));
subplot(1,3,3);imagesc(real(ctranspose(G_al)*G_al));
end;%if disp_flag;

r__k = zeros(2*K_hi+1,1);
tmp = 1/(2*K_hi*sqrt(2*pi));
for nk=1:2*K_hi+1;
k_val = K_al_(nk);
%r__k(nk) = tmp*(randn(1) + i*randn(1)) / max(1,abs(k_val)^2) ;
r__k(nk) = tmp*(randn(1) + i*randn(1)) ;
end;%for nk=1:2*K_hi+1;
for nk=1:2*K_hi+1;
r__k(nk) = conj(r__k(end-(nk-1)));
end;%for nk=1:2*K_hi+1;
r__k(1+K_hi) = 1/sqrt(2*pi);
disp_flag=0;
if disp_flag;
subplot(2,2,1);plot(x__s,real(F_al*r__k),'r-',x__s,imag(F_al*r__k),'b-');
xlabel('x__s'); ylabel('r'); xlim([-pi,pi]); title('rho');
subplot(2,2,2);plot(x__s,real(G_al*r__k),'r-',x__s,imag(G_al*r__k),'b-');
xlabel('x__s'); ylabel('r'); xlim([-pi,pi]); title('int rho');
subplot(2,2,[3,4]);
tmp_z = rand(512*K_hi,1); tmp_x = interp1(real(G_al*r__k),x__s,tmp_z);
x_h = linspace(-pi,pi,128);
tmp_h = hist(tmp_x,x_h); tmp_h = tmp_h/sum(tmp_h)/mean(diff(x_h));
hold on;
plot(x__s,real(F_al*r__k),'r-',x__s,imag(F_al*r__k),'b-');
stairs(x_h,tmp_h);
hold off;
xlabel('x'); ylabel('r'); xlim([-pi,pi]); title('sampled rho');
end;%if disp_flag;

y__k = zeros(2*K_hi+1,1);
for nk=1:2*K_hi+1;
k_val = K_al_(nk);
%y__k(nk) = randn(1) + i*randn(1) / max(1,abs(k_val)^2) ;
y__k(nk) = randn(1) + i*randn(1) ;
end;%for nk=1:2*K_hi+1;
disp_flag=0;
if disp_flag;
plot(x__s,real(F_al*y__k),'r-',x__s,imag(F_al*y__k),'b-');
xlabel('x'); ylabel('y'); xlim([-pi,pi]); title('y');
end;%if disp_flag;

sigma = 0.0; %sigma = 1.0; %sigma = 0.125;
J = 1*512*K_hi;
tmp_z = rand(J,1); x__j = interp1(real(G_al*r__k),x__s,tmp_z);
disp_flag=0;
if disp_flag;
x_h = linspace(-pi,pi,128);
tmp_h = hist(x__j,x_h); tmp_h = tmp_h/sum(tmp_h)/mean(diff(x_h));
hold on;
plot(x__s,real(F_al*r__k),'r-',x__s,imag(F_al*r__k),'b-');
stairs(x_h,tmp_h);
hold off;
xlabel('x'); ylabel('r'); xlim([-pi,pi]); title('sampled rho');
end;%if disp_flag;

F__j = zeros(length(x__j),2*K_hi+1);
F__j(:,1+K_hi) = ones(length(x__j),1);
tmp_e = exp(i*1*x__j);
for nk=1+K_hi+1:2*K_hi+1;
F__j(:,nk) = F__j(:,nk-1).*tmp_e;
end;%for nk=1+K_hi+1:2*K_hi+1;
for nk=1+K_hi-1:-1:1;
F__j(:,nk) = F__j(:,nk+1)./tmp_e;
end;%for nk=1+K_hi-1:-1:1;
F__j = F__j/sqrt(2*pi);
disp_flag=0;
if disp_flag;
[~,tmp_ij] = sort(x__j); 
subplot(1,2,1);imagesc(real(F__j(tmp_ij,:)))
subplot(1,2,2);imagesc(imag(F__j(tmp_ij,:)))
end;%if disp_flag;
y__j = F__j*y__k;
e__j = sigma*(randn(size(y__j)) + i*randn(size(y__j)));
disp_flag=0;
if disp_flag;
hold on;
plot(x__s,real(F_al*y__k),'r-',x__s,imag(F_al*y__k),'b-');
plot(x__j,real(y__j+e__j),'ro',x__j,imag(y__j+e__j),'bo');
hold off;
xlabel('x'); ylabel('y'); xlim([-pi,pi]); title('sampled y+e');
end;%if disp_flag;

F__j_lo = F__j(:,K_lo_al_);
F_al_lo = F_al(:,K_lo_al_);
c_ls = F__j_lo\y__j;

disp_flag=0;
if disp_flag;
hold on;
plot(K_al_,real(ctranspose(F__j)*y__j)*2*pi/J,'rx-',K_al_,imag(ctranspose(F__j)*y__j)*2*pi/J,'bx-');
plot(K_al_,real(R_al*y__k*sqrt(2*pi)),'r.:',K_al_,imag(R_al*y__k*sqrt(2*pi)),'b.:');
hold off;
xlim([-K_hi,+K_hi]); xlabel('K'); ylabel('y_k'); title('F_j^* y_j vs R y_k');
end;%if disp_flag;

R_al = zeros(2*K_hi+1,2*K_hi+1);
for kv1=-K_hi:+K_hi;
for kv2=-K_hi:+K_hi;
dk = kv1-kv2;
if (abs(dk)<=K_hi);
R_al(1+K_hi+kv1,1+K_hi+kv2) = r__k(1+K_hi+dk);
end;%if (abs(dk)<=K_hi);
end;%for kv2=-K_hi:+K_hi;
end;%for kv1=-K_hi:+K_hi;

R_lo_lo = R_al(K_lo_al_,K_lo_al_);
R_lo_hi = R_al(K_lo_al_,K_hi_al_);
y_lo = y__k(K_lo_al_);
y_hi = y__k(K_hi_al_);
b_ls = y_lo + inv(R_lo_lo)*R_lo_hi*y_hi;

eta = 2*pi/(2*K_lo+1)/2;
x__n = 2*eta*[-K_lo:+K_lo];
F__n = zeros(length(x__n),2*K_hi+1);
F__n(:,1+K_hi) = ones(length(x__n),1);
tmp_e = exp(i*1*transpose(x__n));
for nk=1+K_hi+1:2*K_hi+1;
F__n(:,nk) = F__n(:,nk-1).*tmp_e;
end;%for nk=1+K_hi+1:2*K_hi+1;
for nk=1+K_hi-1:-1:1;
F__n(:,nk) = F__n(:,nk+1)./tmp_e;
end;%for nk=1+K_hi-1:-1:1;
F__n = F__n/sqrt(2*pi);
F__n_lo = F__n(:,K_lo_al_);
F__n_hi = F__n(:,K_hi_al_);

disp_flag=0;
if disp_flag;
tmp_k1 = repmat(transpose(K_lo_),1,numel(K_hi_));
tmp_k2 = repmat(K_hi_,numel(K_lo_),1);
tmp_dk = mod(tmp_k1-tmp_k2,2*K_lo+1)==0;
subplot(1,2,1); imagesc(tmp_dk,[0,1]);
subplot(1,2,2); imagesc(real(ctranspose(F__n_lo)*F__n_hi),[0,1]);
end;%if disp_flag;

x__b = max(1,min(2*K_lo+1,1+floor((2*K_lo+1)*(x__j+pi)/(2*pi))));
y__b = zeros(2*K_lo+1,1);
for nk=1:2*K_lo+1;
y__b(nk) = mean(y__j(find(x__b==nk)));
end;%for nk=1:2*K_lo+1;
c_i0 = 2*eta*ctranspose(F__n_lo)*y__b;

F1_n_lo = F__n_lo*diag(i*[-K_lo:K_lo]);
F2_n_lo = F1_n_lo*diag(i*[-K_lo:K_lo]);
F3_n_lo = F2_n_lo*diag(i*[-K_lo:K_lo]);
p=1;
y__b = zeros((2*K_lo+1),1);
y1_b = zeros((2*K_lo+1),1);
for nk=1:2*K_lo+1;
tmp_ij = find(x__b==nk);
[p_x_] = test_LSQ_vs_Sampling_excerpt_0(p,x__j(tmp_ij),y__j(tmp_ij));
y__b(nk) = polyval(p_x_,x__n(nk));
y1_b(nk) = polyval(polyder(p_x_),x__n(nk));
end;%for nk=1:2*K_lo+1;
c_i1 = [F__n_lo ; F1_n_lo] \ [y__b ; y1_b];

p=2;
y__b = zeros((2*K_lo+1),1);
y1_b = zeros((2*K_lo+1),1);
y2_b = zeros((2*K_lo+1),1);
for nk=1:2*K_lo+1;
tmp_ij = find(x__b==nk);
[p_x_] = test_LSQ_vs_Sampling_excerpt_0(p,x__j(tmp_ij),y__j(tmp_ij));
y__b(nk) = polyval(p_x_,x__n(nk));
y1_b(nk) = polyval(polyder(p_x_),x__n(nk));
y2_b(nk) = polyval(polyder(polyder(p_x_)),x__n(nk));
end;%for nk=1:2*K_lo+1;
c_i2 = [F__n_lo ; F1_n_lo ; F2_n_lo] \ [y__b ; y1_b ; y2_b];

p=3;
y__b = zeros((2*K_lo+1),1);
y1_b = zeros((2*K_lo+1),1);
y2_b = zeros((2*K_lo+1),1);
y3_b = zeros((2*K_lo+1),1);
for nk=1:2*K_lo+1;
tmp_ij = find(x__b==nk);
[p_x_] = test_LSQ_vs_Sampling_excerpt_0(p,x__j(tmp_ij),y__j(tmp_ij));
y__b(nk) = polyval(p_x_,x__n(nk));
y1_b(nk) = polyval(polyder(p_x_),x__n(nk));
y2_b(nk) = polyval(polyder(polyder(p_x_)),x__n(nk));
y3_b(nk) = polyval(polyder(polyder(polyder(p_x_))),x__n(nk));
end;%for nk=1:2*K_lo+1;
c_i3 = [F__n_lo ; F1_n_lo ; F2_n_lo ; F3_n_lo] \ [y__b ; y1_b ; y2_b ; y3_b];

D_al = zeros(2*K_hi+1,1);
for kv=-K_hi:+K_hi;
D_al(1+K_hi+kv) = (2/kv)*sin(kv*eta);
end;%for kv=-K_hi:+K_hi;
D_al(1+K_hi) = 2*eta;
D_lo = D_al(K_lo_al_);
D_hi = D_al(K_hi_al_);
r_lo = r__k(K_lo_al_);
r_hi = r__k(K_hi_al_);
R_hi_lo = R_al(K_hi_al_,K_lo_al_);
R_hi_hi = R_al(K_hi_al_,K_hi_al_);

z_vec = F__n_lo*diag(D_lo)*r_lo + F__n_hi*diag(D_hi)*r_hi;
s_vec = (F__n_lo*diag(D_lo)*(R_lo_lo*y_lo + R_lo_hi*y_hi) + F__n_hi*diag(D_hi)*(R_hi_lo*y_lo + R_hi_hi*y_hi))/sqrt(2*pi);
b_i0 = 2*eta*ctranspose(F__n_lo)*diag(1./z_vec)*s_vec;

disp_flag=1;
if disp_flag;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
subplot(3,3,1);
hold on;
plot(x__s,real(F_al*y__k),'r:',x__s,imag(F_al*y__k),'b:');
plot(x__s,real(F_al_lo*c_ls),'r-',x__s,imag(F_al_lo*c_ls),'b-');
hold off;
xlabel('x'); ylabel('y'); xlim([-pi,pi]); title('y_al and c_ls','Interpreter','none');
subplot(3,3,2);
hold on;
plot(x__s,real(F_al_lo*y_lo),'r:',x__s,imag(F_al_lo*y_lo),'b:');
plot(x__s,real(F_al_lo*c_ls),'r-',x__s,imag(F_al_lo*c_ls),'b-');
hold off;
xlabel('x'); ylabel('y'); xlim([-pi,pi]); title('y_lo and c_ls','Interpreter','none');
subplot(3,3,3);
hold on;
plot(x__s,real(F_al_lo*b_ls),'r:',x__s,imag(F_al_lo*b_ls),'b:');
plot(x__s,real(F_al_lo*c_ls),'r-',x__s,imag(F_al_lo*c_ls),'b-');
hold off;
xlabel('x'); ylabel('y'); xlim([-pi,pi]); title('b_ls and c_ls','Interpreter','none');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
subplot(3,3,4);
hold on;
plot(x__s,real(F_al*y__k),'r:',x__s,imag(F_al*y__k),'b:');
plot(x__s,real(F_al_lo*c_i0),'r-',x__s,imag(F_al_lo*c_i0),'b-');
hold off;
xlabel('x'); ylabel('y'); xlim([-pi,pi]); title('y_al and c_i0','Interpreter','none');
subplot(3,3,5);
hold on;
plot(x__s,real(F_al_lo*y_lo),'r:',x__s,imag(F_al_lo*y_lo),'b:');
plot(x__s,real(F_al_lo*c_i0),'r-',x__s,imag(F_al_lo*c_i0),'b-');
hold off;
xlabel('x'); ylabel('y'); xlim([-pi,pi]); title('y_lo and c_i0','Interpreter','none');
subplot(3,3,6);
hold on;
plot(x__s,real(F_al_lo*b_i0),'r:',x__s,imag(F_al_lo*b_i0),'b:');
plot(x__s,real(F_al_lo*c_i0),'r-',x__s,imag(F_al_lo*c_i0),'b-');
hold off;
xlabel('x'); ylabel('y'); xlim([-pi,pi]); title('b_i0 and c_i0','Interpreter','none');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
subplot(3,3,7);
hold on;
plot(x__s,real(F_al_lo*y_lo),'r:',x__s,imag(F_al_lo*y_lo),'b:');
plot(x__s,real(F_al_lo*c_i1),'r-',x__s,imag(F_al_lo*c_i0),'b-');
hold off;
xlabel('x'); ylabel('y'); xlim([-pi,pi]); title('y_lo and c_i1','Interpreter','none');
subplot(3,3,8);
hold on;
plot(x__s,real(F_al_lo*y_lo),'r:',x__s,imag(F_al_lo*y_lo),'b:');
plot(x__s,real(F_al_lo*c_i2),'r-',x__s,imag(F_al_lo*c_i0),'b-');
hold off;
xlabel('x'); ylabel('y'); xlim([-pi,pi]); title('y_lo and c_i2','Interpreter','none');
subplot(3,3,9);
hold on;
plot(x__s,real(F_al_lo*y_lo),'r:',x__s,imag(F_al_lo*y_lo),'b:');
plot(x__s,real(F_al_lo*c_i3),'r-',x__s,imag(F_al_lo*c_i0),'b-');
hold off;
xlabel('x'); ylabel('y'); xlim([-pi,pi]); title('y_lo and c_i3','Interpreter','none');
end;%if disp_flag;

disp_flag=0;
if disp_flag;
end;%if disp_flag;

disp_flag=0;
if disp_flag;
end;%if disp_flag;


