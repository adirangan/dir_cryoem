function v_ = polyval_r8_reverse_0(n_p,p_,n_x,x_);
% Evaluates a polynomial given by: ;
% v_(j) = p_(0) + p_(1)*x_(j) + p_(2)*x_(j)**2 + .. + p_(n_p-1)*x_(j)**(n_p-1) ;
for nx=0:n_x-1;
x = x_(1+nx);
p = 0.0d0;
nj=n_p-1;
for np=0:n_p-1;
p = p_(1+nj) + p*x;
nj = nj-1;
end%for np=0:n_p-1;
v_(1+nx) = p;
end;%for nx=0:n_x-1;


