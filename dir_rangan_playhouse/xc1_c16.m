function C_ = xc1_c16(n_x,A_,B_) ;
C_=  zeros(size(A_));
na = 0 ;
nb = 0 ;
nc = 0 ;
for nx=0:n_x-1;
C_(1+nc) = A_(1+na)*conj(B_(1+nb)) ;
na = na + 1 ;
nb = nb + 1 ;
nc = nc + 1 ;
end;%for nx=0:n_x-1;

