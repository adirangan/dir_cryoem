%%%%%%%%;
% quick test of half-periodic block-functions in 1d. ;
%%%%%%%%;
nf=0;

k_max = 133;
k_ = transpose(-k_max:2:k_max);
n_k = numel(k_);
nkn1 = efind(k_==-1);
nkp1 = efind(k_==+1);

%%%%%%%%;
% V = 2*cos(2x) ;
%%%%%%%%;
V__ = zeros(n_k,n_k);
V__ = diag(ones(n_k-1,1),+1) + diag(ones(n_k-1,1),-1);
D__ = diag(k_.^2);
n_svd = 6;
[B_kn__,E_nn__] = eigs(V__+D__,n_svd,'smallestabs'); E_n_ = diag(E_nn__);
%%%%;
n_x = 1024*8+1;
x_ = transpose(linspace(0,2*pi,n_x));
V_x_ = cos(2*x_);
F_xk__ = zeros(n_x,n_k);
for nk=0:n_k-1;
k_val = k_(1+nk);
F_xk__(:,1+nk) = exp(i*k_val*x_);
end;%for nk=0:n_k-1;
B_xn__ = F_xk__*B_kn__;
%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,1,1);
plot(k_,B_kn__,'o-','LineWidth',3);
xlim([-k_max-1,+k_max+1]);
xlabel('k_','Interpreter','none');
subplot(2,1,2);
hold on;
plot([pi,pi],[-2,+2],'k-','LineWidth',5);
plot(x_,V_x_,'k-','LineWidth',5);
plot(x_,B_xn__,'-','LineWidth',3);
xlim([0,2*pi]); set(gca,'XTick',[0*pi,1*pi,2*pi],'XTickLabel',{'0','\pi','2\pi'});
xlabel('x_','Interpreter','none');



