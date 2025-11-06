function v_ = chebval_0(n_t,t_chebcoef_,n_x,x_)
% Evaluates a polynomial given by: ;
% v_(j) = t_(0)*T_{0}(x_{j}) + t_(1)*T_{1}(x_{j}) + t_(2)*T_{2}(x_{j}) + ... + t_(n_t-1)*T_{n_t-1}(x_{j}). ;

if nargin<1;
disp(sprintf(' %% testing chebval_0'));
n_t = 8;
t_chebcoef_ = randn(n_t,1);
t_chebfun_ = chebfun(t_chebcoef_,'coeffs');
n_y = 1024+1;
y_ = linspace(-1,+1,n_y);
n_x = 128;
x_ = sort(2*rand(n_x,1)-1);
plot(y_,t_chebfun_(y_),'k-',x_,chebval_0(n_t,t_chebcoef_,n_x,x_),'ro');
disp('returning'); return;
end;%if nargin<1;

flag_slow = 0;
if flag_slow;
v_ = zeros(size(x_));
for nt=0:n_t-1;
v_ = v_ + t_chebcoef_(1+nt)*cos(nt*acos(max(-1,min(+1,x_))));
end;%for nt=0:n_t-1;
end;%if flag_slow;

flag_fast = 1;
if flag_fast;
v_ = zeros(n_x,1);
nt_ = transpose([0:n_t-1]);
[nt__,x__] = ndgrid(nt_,x_(:));
v_ = reshape(sum(bsxfun(@times,reshape(t_chebcoef_,[n_t,1]),cos(nt__.*acos(max(-1,min(+1,x__))))),1),size(x_));
end;%if flag_fast;
