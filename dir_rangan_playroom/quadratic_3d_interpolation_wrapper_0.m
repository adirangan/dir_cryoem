function [max_d_,x_,F,DF,DDF] = quadratic_3d_interpolation_wrapper_0(X_,n_d_,periodic_flag_);
% This function estimates the maximum of X_ using quadratic_3d_interpolation_0.m ;
% X_ is a (linear) array of dimension n_d = length(n_d_) with prod(n_d_) entries. ;
% n_d_ is an (integer) array of length n_d indicating the size of X_. ;
% periodic_flag_ is an (integer) array of length n_d indicating whether or not the corresponding dimension of X_ is to be treated periodically (1) or not (0). ;

if nargin<1;
n_d = 3;
n_d_ = 2+ceil(4*rand(n_d,1));
X_ = ceil(128*rand(prod(n_d_),1));
periodic_flag_ = 1+zeros(n_d,1);
[max_d_,x_,F,DF,DDF] = quadratic_3d_interpolation_wrapper_0(X_,n_d_,periodic_flag_);
disp(num2str(transpose(n_d_)));
disp(num2str(transpose(max_d_)));
disp(num2str(transpose(x_)));
max_ndd = get_ndd(n_d_,max_d_);
disp(sprintf(' %% max_ndd %d X_(%d) = %0.2f, F = %0.2f',max_ndd,max_ndd,X_(1+max_ndd),F));
disp(reshape(X_,transpose(n_d_)));
disp('returning'); return;
end;%if nargin<1;

verbose=0;

n_d = length(n_d_); n_3d = 3.^n_d;
assert(prod(n_d_)<=numel(X_));
assert(length(periodic_flag_)==n_d);
[max_v,max_ij] = max(X_(:));
max_d_ = get_d__(n_d_,max_ij-1);
D_ = zeros(n_3d,1);
for n3d=0:n_3d-1;
dx_ = get_x__(n_d,n3d);
tmp_d_ = max_d_ + dx_;
for nd=0:n_d-1;
if tmp_d_(1+nd)<0; 
if periodic_flag_(1+nd); tmp_d_(1+nd) = n_d_(1+nd)-1; end;
if ~periodic_flag_(1+nd); tmp_d_(1+nd) = min(1,n_d_(1+nd)-1); end;
end;%if tmp_d_(1+nd)<0; 
if tmp_d_(1+nd)>n_d_(1+nd)-1; 
if periodic_flag_(1+nd); tmp_d_(1+nd) = 0; end;
if ~periodic_flag_(1+nd); tmp_d_(1+nd) = max(0,n_d_(1+nd)-2); end;
end;%if tmp_d_(1+nd)>n_d_(1+nd)-1; 
end;%for nd=0:n_d-1;
tmp_ndd = get_ndd(n_d_,tmp_d_);
D_(1+n3d) = X_(1+tmp_ndd);
end;%for n3d=0:n_3d-1;
D = reshape(D_,3*ones(1,n_d));
[x_,F,DF,DDF] = quadratic_3d_interpolation_0(D);

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

function d_ = get_d__(n_d_,ndd);
n_d = length(n_d_);
d_ = zeros(n_d,1);
tmp = ndd; 
for nd=0:n_d-1; 
d_(1+nd) = mod(tmp,n_d_(1+nd)); 
tmp = tmp-d_(1+nd); 
tmp = tmp/n_d_(1+nd); 
end;%for nd=0:n_d-1;

function ndd = get_ndd(n_d_,d_);
n_d = length(n_d_);
tmp_ = cumprod([1;reshape(n_d_(1:end-1),n_d-1,1)]);
ndd = 0;
for nd=0:n_d-1;
ndd = ndd + d_(1+nd)*tmp_(1+nd);
end;%for nd=0:n_d-1;

