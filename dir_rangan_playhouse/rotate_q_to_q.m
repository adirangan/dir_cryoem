function M_q_ = rotate_q_to_q(n_r,n_w_,n_A,S_q_,gamma);
% Assumes that M_q_ is the same size and dimensions as S_q_;
% Multiplication performed in place;
M_q_ = zeros(n_A,1);
n_w_max = n_w_(1+n_r-1);
C_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
q = nw - floor(n_w_max/2);
C = exp(-i*q*gamma);
C_(1+nw) = C;
end;%for nw=0:n_w_max-1;
ic=0;
for nr=0:n_r-1;
for nq=0:n_w_(1+nr)-1;
q = nq;
if (q>floor(n_w_(1+nr)/2)-1); q = q-n_w_(1+nr); end;
C = C_(1+floor(n_w_max/2) + q);
M_q_(1+ic) = S_q_(1+ic)*C;
ic = ic+1;
end;%for nq=0:n_w_(1+nr)-1;
end;%for nr=0:n_r-1;



