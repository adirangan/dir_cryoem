%%%%%%%%;
% testing lsq preconditioner. ;
% using quadrature on simple [-1,+1] interval ;
% Here we assume that sample points (x_rand_) are sampled from the distribution mu(x) = (1+x)/2 ;
% Note that mu has coefficients mu0 = 1/sqrt(2) and mu1 = 1/sqrt(6). ;
% Note that the cdf MU(x) is y=(1/2 + x/2).^2, and the inverse is x=2*sqrt(y)-1;
%%%%%%%%;

setup;
rng(1);
n_quad = 24;
n_rand = 1430;
x_dens_ = linspace(-1,+1,1024);
[x_quad_,w_quad_] = legpts(n_quad);
x_rand_ = 2*sqrt(rand(n_rand,1))-1;
P_quad_ = cell(1+n_quad,1);
for nquad=0:n_quad;
P_quad_{1+nquad} = legpoly(nquad)/sqrt(2/(2*nquad+1));
end;%for nquad=0:n_quad;
%for nquad=0:n_quad;disp(sprintf('nquad %d:',nquad)); disp(rats(sum(P_quad_{1+nquad}*P_quad_{1+nquad}))); end; %<-- verify normality. ;
P_rand_ = cell(1+n_quad,1);
for nquad=0:n_quad;
P_rand_{1+nquad} = P_quad_{1+nquad}(x_rand_);
end;%for nrand=0:n_rand;

%%%%%%%%;
% Now determine integral: ;
% PP0_{m1m2} = \int_{-1}^{+1} P_quad_{m1}(x) P_quad_{m2}(x) P_quad_{1+0}(x) dx ;
% PP1_{m1m2} = \int_{-1}^{+1} P_quad_{m1}(x) P_quad_{m2}(x) P_quad_{1+1}(x) dx ;
%%%%%%%%;
PP0__ = zeros(1+n_quad,1+n_quad);
PP1__ = zeros(1+n_quad,1+n_quad);
PP2__ = zeros(1+n_quad,1+n_quad);
for nquad1=0:n_quad;
for nquad2=nquad1:n_quad;
PP0__(1+nquad1,1+nquad2) = sum(P_quad_{1+nquad1}*P_quad_{1+nquad2}*P_quad_{1+0});
PP0__(1+nquad2,1+nquad1) = PP0__(1+nquad1,1+nquad2);
PP1__(1+nquad1,1+nquad2) = sum(P_quad_{1+nquad1}*P_quad_{1+nquad2}*P_quad_{1+1});
PP1__(1+nquad2,1+nquad1) = PP1__(1+nquad1,1+nquad2);
PP2__(1+nquad1,1+nquad2) = sum(P_quad_{1+nquad1}*P_quad_{1+nquad2}*P_quad_{1+2});
PP2__(1+nquad2,1+nquad1) = PP2__(1+nquad1,1+nquad2);
end;%for nquad2=nquad1:n_quad;
end;%for nquad1=0:n_quad;

%%%%%%%%;
% Now set coefficients of mu. ;
%%%%%%%%;
mu0 = 1/sqrt(2);
mu1 = 1/sqrt(6); 
FF__ = mu0*PP0__(1:n_quad,1:n_quad) + mu1*PP1__(1:n_quad,1:n_quad) ;
[tmp_U,tmp_S,tmp_V] = svd(FF__);
EE__ = tmp_V*inv(sqrt(tmp_S))*transpose(tmp_U);

%%%%%%%%;w
% Now solve lsq for a right-hand-side ;
%%%%%%%%;
y_ = chebfun(@(x) cos(2*pi*x));
y_dens_ = cos(2*pi*x_dens_);
y_quad_ = cos(2*pi*x_quad_);
y_rand_ = cos(2*pi*x_rand_) + randn(size(x_rand_));
%{
f_quad_ = zeros(n_quad,1);
p_quad_ = chebfun(0);
for nquad=0:n_quad-1;
f_quad_(1+nquad) = w_quad_*(P_quad_{1+nquad}(x_quad_).*y_quad_);
p_quad_ = p_quad_ + P_quad_{1+nquad}*f_quad_(1+nquad);
end;%for nquad=0:n_quad-1;
plot(x_rand_,y_rand_,'ko',x_dens_,p_quad_(x_dens_),'r-');
 %}
P_rand__ = zeros(n_rand,n_quad);
P_quad__ = zeros(n_quad,n_quad);
for nquad=0:n_quad-1;
P_rand__(:,1+nquad) = P_quad_{1+nquad}(x_rand_);
P_quad__(:,1+nquad) = P_quad_{1+nquad}(x_quad_);
end;%for nquad=0:n_quad-1;
f_rand_ = P_rand__\y_rand_;
p_rand_ = chebfun(0); for nquad=0:n_quad-1; p_rand_ = p_rand_ + f_rand_(1+nquad)*P_quad_{1+nquad}; end;%for nquad=0:n_quad-1;
plot(x_rand_,y_rand_,'ko',x_rand_,P_rand__*f_rand_,'go',x_dens_,p_rand_(x_dens_),'g-');
%%%%%%%%;
PP__ = transpose(P_rand__)*P_rand__;
PPP__ = inv(PP__)*transpose(P_rand__);
%%%%%%%%;
disp(sprintf(' %% error normal-eq -vs- backslash: %0.16f',norm(f_rand_ - PPP__*y_rand_)));

k_PP = cond(PP__);
k_EEPPEE = cond(EE__*PP__*EE__);
disp(sprintf(' %% cond(PP__) = %0.2f',k_PP));
disp(sprintf(' %% cond(EE__PP__EE__) = %0.2f',k_EEPPEE));
disp(sprintf(' %% ratio: %0.2f = 1/%0.2f',k_PP/k_EEPPEE,k_EEPPEE/k_PP));

%{
tmp_knn_ = knnsearch(x_quad_,x_rand_);
tmp_num_ = zeros(n_quad,1);
for nquad=1:n_quad;tmp_num_(nquad)=length(find(tmp_knn_==nquad)); end;
tmp_div_ = zeros(n_rand,1);
for nquad=1:n_quad;tmp_ij_ = find(tmp_knn_==nquad); tmp_div_(tmp_ij_) = 1/tmp_num_(nquad); end;
I__ = sparse(tmp_knn_,1:n_rand,tmp_div_,n_quad,n_rand);
PWI__ = transpose(P_quad__)*diag(w_quad_)*I__;
 %}


