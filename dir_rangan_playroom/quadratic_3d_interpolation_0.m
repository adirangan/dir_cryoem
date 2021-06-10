function [x_,F,DF,DDF] = quadratic_3d_interpolation_0(D);
% This function creates a quadratic interpolant to the 3^n_d dimensional array D. ;
% D is expected to be a 3-x-3-x-3-x-...-x-3 array of dimension n_d. ;
% D represents the function-evaluations of a n_d-dimensional function evaluated ;
% on the gridpoints [-1,0,+1]^n_d. ;

if nargin<1;
n_d = 3; n_3d = 3.^n_d;
p_test = 27:-1:1; 
for n3d=0:n_3d-1;
x__ = get_x__(n_d,n3d);
D_(1+n3d) = p_eval(n_d,p_test,x__);
end;%for n3d=0:n_3d-1;
D_ = reshape(D_,[3,3,3]);
[x_,F,DF,DDF] = quadratic_3d_interpolation_0(D_);
disp(sprintf(' %% p_(x_) = %0.2f = %0.2f, norm(grad) = %0.16f',p_eval(n_d,p_test,x_),F,norm(DF)));
disp('returning'); return;
end;%if nargin<1;

verbose=0;

n_d = length(size(D)); n_3d = 3^n_d;
if (verbose>0); disp(sprintf(' %% n_d %d n_3d %d',n_d,n_3d)); end;
if (min(size(D))~=3 | max(size(D))~=3) disp(sprintf(' %% Warning! D has dimension of size ~= 3 in quadratic_3d_interpolation')); end;

M___ = zeros(n_3d,n_3d);
nd1_ = zeros(n_d,1); 
nd2_ = zeros(n_d,1); 
for n3d1=0:n_3d-1;
nd1_ = get_nd_(n_d,n3d1); if (verbose>2); disp(sprintf(' %% n3d1 %.3d nd1_ [%d,%d,%d]',n3d1,nd1_)); end;
for n3d2=0:n_3d-1;
nd2_ = get_nd_(n_d,n3d2); if (verbose>3); disp(sprintf(' %% n3d2 %.3d nd2_ [%d,%d,%d]',n3d2,nd2_)); end;
M___(1+n3d1,1+n3d2) = get_M__(n_d,nd1_,nd2_);
end;%for n3d=0:n_3d-1;
end;%for n3d=0:n_3d-1;

if (verbose>1); disp(M___); end;
if (verbose>0); disp(rank(M___)); disp(svds(M___,n_3d)); end;

p___ = M___\D(:);
check_flag=1;
if check_flag;
D0 = zeros(n_3d,1);
for n3d1=0:n_3d-1;
nd1_ = get_nd_(n_d,n3d1);
D0(1+n3d1) = 0;
for n3d2=0:n_3d-1;
nd2_ = get_nd_(n_d,n3d2);
D0(1+n3d1) = D0(1+n3d1) + p___(1+n3d2)*get_M__(n_d,nd1_,nd2_);
end;%for n3d2=0:n_3d-1;
if (verbose>1); disp(sprintf(' %% n3d1 %.3d nd1_ [%d,%d,%d], D(%d) %0.2f vs %0.2f',n3d1,nd1_,n3d1,D(1+n3d1),D0(1+n3d1))); end;
end;%for n3d1=0:n_3d-1;
if (verbose>0); disp(norm(D(:)-D0)); end;
end;%if check_flag;

iteration_max = 16;
x_ = zeros(n_d,1); ni=1; continue_flag=1;
while (continue_flag);
F = p_eval(n_d,p___,x_);
DF = p_grad(n_d,p___,x_);
DDF = p_hess(n_d,p___,x_);
if (rank(DDF)<n_d); dx_ = +DF/32; end;
if (rank(DDF)==n_d); dx_ = -DDF\DF; end;
x_ = max(-1,min(+1,x_ + dx_));
ni = ni+1;
continue_flag = (ni<iteration_max & norm(DF)>1e-9);
end;% while;

if (norm(DF)>1e-9); disp(sprintf(' %% Warning! iteration ni %d, norm(DF) %0.16f in quadratic_3d_interpolation_0',ni,norm(DF))); end;

if (F<p___(1)); 
disp(sprintf(' %% Warning! F %0.2f < on-grid maximum p___(1)=%0.2f',F,p___(1)));
x_ = zeros(n_d,1); F = p___(1); DF = zeros(n_d,1); DDF = zeros(n_d,n_d); 
end;%if (F<p___(1)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

function nd_ = get_nd_(n_d,n3d);
nd_ = zeros(n_d,1);
tmp = n3d; 
for nd=0:n_d-1; 
nd_(1+nd) = mod(tmp,3); 
tmp = tmp-nd_(1+nd); 
tmp = tmp/3; 
end;%for nd=0:n_d-1;

function x__ = get_x__(n_d,n3d);
e___ = [-1,0,1];
nd_ = get_nd_(n_d,n3d);
x__ = zeros(n_d,1);
for nd=0:n_d-1;
x__(1+nd) = e___(1+nd_(1+nd));
end;%for nd=0:n_d-1;

function M__ = get_M__(n_d,nd1_,nd2_);
e___ = [-1,0,1];
tmp = 1;
for nd=0:n_d-1;
tmp = tmp*(e___(1+nd1_(1+nd)).^nd2_(1+nd));
end;%for nd=0:n_d-1;
M__ = tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

function output = p_eval(n_d,p_,x_);
% evaluate quadratic polynomial ;
n_3d = 3.^n_d;
output = 0;
for n3d=0:n_3d-1;
nd_ = get_nd_(n_d,n3d);
output = output + p_(1+n3d) * prod(x_.^nd_);
end;%for n3d=0:n_3d-1;

function output = p_grad(n_d,p_,x_);
% evaluate gradient of quadratic polynomial ;
n_3d = 3.^n_d;
output = zeros(n_d,1);
x__ = zeros(n_d,n_d);
for n3d=0:n_3d-1;
nd_ = get_nd_(n_d,n3d);
ne_ = max(0,repmat(nd_,1,n_d) - eye(n_d));
for nd=0:n_d-1;
x__(:,1+nd) = x_.^(ne_(:,1+nd));
end;%for nd=0:n_d-1;
for nd=0:n_d-1;
output(1+nd) = output(1+nd) + p_(1+n3d) * nd_(1+nd) * prod(x__(:,1+nd));
end;%for nd=0:n_d-1;
end;%for n3d=0:n_3d-1;

function output = p_hess(n_d,p_,x_);
% evaluate hessian of quadratic polynomial ;
n_3d = 3.^n_d;
output = zeros(n_d,n_d);
x___ = zeros(n_d,n_d,n_d);
I2 = repmat(eye(n_d),1,1,n_d);
I3 = permute(repmat(eye(n_d),1,1,n_d),[1,3,2]);
for n3d=0:n_3d-1;
nd_ = get_nd_(n_d,n3d);
ne__ = max(0,repmat(nd_,1,n_d,n_d) - I2 - I3);
for nd1=0:n_d-1; for nd2=0:n_d-1;
x___(:,1+nd1,1+nd2) = x_.^(ne__(:,1+nd1,1+nd2));
end;end;%for nd1=0:n_d-1; for nd2=0:n_d-1;
for nd1=0:n_d-1; for nd2=0:n_d-1;
if (nd1==nd2);
output(1+nd1,1+nd2) = output(1+nd1,1+nd2) + p_(1+n3d) * nd_(1+nd1) * (nd_(1+nd2)-1) * prod(x___(:,1+nd1,1+nd2));
end;%if (nd1==nd2);
if (nd1~=nd2);
output(1+nd1,1+nd2) = output(1+nd1,1+nd2) + p_(1+n3d) * nd_(1+nd1) * nd_(1+nd2) * prod(x___(:,1+nd1,1+nd2));
end;%if (nd1~=nd2);
end;end;%for nd1=0:n_d-1; for nd2=0:n_d-1;
end;%for n3d=0:n_3d-1;


