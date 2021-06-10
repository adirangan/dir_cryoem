function [x,r,n] = conjgrad_0(A, b, x, tol) ;
% taken from wikipedia. ;
if nargin<4; tol=1e-6; end;
if nargin<3; x=b; end;

r = b - A * x;
p = r;
rsold = transpose(r) * r;

n=1;
for i = 1:length(b);
Ap = A * p;
alpha = rsold / (transpose(p) * Ap);
x = x + alpha * p;
r = r - alpha * Ap;
rsnew = transpose(r) * r;
if sqrt(rsnew) < tol;
break;
end;%if sqrt(rsnew) < tol;
p = r + (rsnew / rsold) * p;
rsold = rsnew;
n=n+1;
end;%for i = 1:length(b);

