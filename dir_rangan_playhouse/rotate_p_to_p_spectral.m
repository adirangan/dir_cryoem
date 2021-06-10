function M_p_ = rotate_p_to_p_spectral(n_r,n_w_,n_A,S_p_,gamma);
% Assumes that M_p_ is the same size and dimensions as S_p_;
S_q_ = interp_p_to_q(n_r,n_w_,n_A,S_p_);
M_q_ = zeros(size(S_q_));
ic=0;
icstart=0;
for nr=0:n_r-1         ;
for nq=0:n_w_(1+nr)-1;
q = periodize(nq,-n_w_(1+nr)/2,n_w_(1+nr)/2-1);
M_q_(1+ic) = S_q_(1+ic)*exp(-i*q*gamma);
ic = ic + 1;
end;%for;
icstart = icstart + n_w_(1+nr);
end;%for;
M_p_ = interp_q_to_p(n_r,n_w_,n_A,M_q_);


