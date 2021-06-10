function orthopoly_test_0();

K = input('how many points would you like to use for interpolation (default = 16): ');if isempty(K);K=16;end;
a = input('The interval starts at the point (default t = -1): ');if isempty(a);a=-1;end;
c = input('The interval stops at the point (default t = 4): ');if isempty(c);c=4;end;
M = input('how many points would you like to use for your graphs (default = 1000): ');if isempty(M);M=1000;end;
%K = 16; a = -1; c = +4; M = 1024;
tmid = (c+a)/2; trad = (c-a)/2;
xvals = linspace(-1,+1,M); tvals = linspace(a,c,M);
w = @(x) (1 + x).^2; wv = [1 2 1];
disp(sprintf('Using wv = [1 2 1] --> w(x) = (1+x).^2'));
%w = @(x) 1 + x; wv = [1 1];
%w = @(x) ones(size(x)); wv = [1];
[lx,lw,Lx,Lv] = orthopoly_node_weight_matrix_0(K,wv);
w_at_xv = w(xvals);
w_at_lx = w(lx);

D1 = zeros(K,K); D2 = zeros(K,K);
for nk1=0:K-1;for nk2=0:K-1;
D1(1+nk1,1+nk2) = polyint(conv(wv,conv(Lv(1+nk1,:),Lv(1+nk2,:))),[-1,1]);
D2(1+nk1,1+nk2) = sum(transpose(Lx(1+nk1,:)).*transpose(Lx(1+nk2,:)).*lw);
end;end;%for nk1=0:K-1;for nk2=0:K-1;

lt = lx*trad + tmid;

alpha = 1/10; sigma = 2;
f = @(x) 1./(x.^2 + alpha) + sin(2*pi*x) + cos(4*pi*x).*(1/sqrt(2*pi)/sigma).*exp(-(x).^2/(2*sigma.^2));
f_at_lt = f(lt);
f_at_tv = f(tvals);

orthopoly_coefs = zeros(K,1);
for nk1=0:K-1;
orthopoly_coefs(1+nk1) = sum(f_at_lt.*transpose(Lx(1+nk1,:)).*lw);
end;%for nk1=0:K-1;

orthopoly_p = zeros(1,K+1);
for nk1=0:K-1;
orthopoly_p = orthopoly_p + orthopoly_coefs(1+nk1)*Lv(1+nk1,:);
end;%for nk1=0:K-1;
p_at_xv = polyval(orthopoly_p,xvals);

figure(1);clf;hold on;
plot(tvals,f_at_tv,'r-');
plot(tvals,p_at_xv,'b-');
plot(lt,f_at_lt,'bx');
plot(lt,0,'bx');
xlabel('x');
ylabel('function value');
title('function in red, polynomial in blue, nodes at *');
hold off;
disp(sprintf('the calculated integral of your function is %f',sum(f_at_lt.*lw)*trad));
disp(sprintf('the reiman_sum integral of your function is %f',sum(f_at_tv.*w_at_xv)*2/M*trad));

