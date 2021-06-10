%%%%%%%%;
% testing lsq preconditioner. ;
%%%%%%%%;

setup;
rng(1);
n_quad = 16;
n_rand = 43;
x_dens_ = linspace(-1,+1,1024);
[x_quad_,w_quad_] = legpts(n_quad);
x_rand_ = sort(2*rand(n_rand,1)-1);
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

PPP__ = inv(transpose(P_rand__)*P_rand__)*transpose(P_rand__);
Q__ = zeros(n_quad,n_quad,n_quad,n_rand);
for n1=0:n_quad-1;
for n2=0:n_quad-1;
for n3=0:n_quad-1;
for j1=1:n_rand;
Q__(1+n1,1+n2,1+n3,j1) = P_quad__(1+n3,1+n1)*w_quad_(1+n3)*P_rand__(j1,1+n2);
end;%for j1=1:n_rand;
end;%for n3=0:n_quad-1;
end;%for n2=0:n_quad-1;
end;%for n1=0:n_quad-1;
C__ = reshape(reshape(Q__,n_quad^2,n_quad*n_rand)\reshape(eye(n_quad),n_quad^2,1),n_quad,n_rand);
PWCP__ = transpose(P_quad__)*diag(w_quad_)*C__*P_rand__;
PWC__ = transpose(P_quad__)*diag(w_quad_)*C__;
f_cond_ = PWC__*y_rand_;
p_cond_ = chebfun(0); for nquad=0:n_quad-1; p_cond_ = p_cond_ + f_cond_(1+nquad)*P_quad_{1+nquad}; end;%for nquad=0:n_quad-1;
plot(x_rand_,y_rand_,'ko',x_rand_,P_rand__*f_rand_,'go',x_dens_,p_rand_(x_dens_),'g-',x_rand_,P_rand__*f_cond_,'ro',x_dens_,p_cond_(x_dens_),'r-');
disp(sprintf(' %% error lsq %0.4f , cond %0.4f',sum(abs(y_rand_-P_rand__*f_rand_).^2),sum(abs(y_rand_-P_rand__*f_cond_).^2)));

tmp_knn_ = knnsearch(x_quad_,x_rand_);
tmp_num_ = zeros(n_quad,1);
for nquad=1:n_quad;tmp_num_(nquad)=length(find(tmp_knn_==nquad)); end;
tmp_div_ = zeros(n_rand,1);
for nquad=1:n_quad;tmp_ij_ = find(tmp_knn_==nquad); tmp_div_(tmp_ij_) = 1/tmp_num_(nquad); end;
I__ = sparse(tmp_knn_,1:n_rand,tmp_div_,n_quad,n_rand);
PWI__ = transpose(P_quad__)*diag(w_quad_)*I__;


