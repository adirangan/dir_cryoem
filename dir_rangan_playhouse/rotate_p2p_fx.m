function M_p_ = rotate_p2p_fx(n_r,n_w_,n_A,S_p_,gamma);
% Assumes that M_p_ is the same size and dimensions as S_p_;
M_p_ = zeros(size(S_p_));

ic = 0;
for nr=0:n_r-1;
ic = ic + n_w_(1+nr);
end;%for nr=0:n_r-1;
if (ic~=n_A);
disp(sprintf(' %% Warning! n_A %d, ic %d',n_A,ic));
end;%if (ic~=n_A);

n_w_max = n_w_(1+n_r-1);
C_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
q = nw - floor(n_w_max/2);
C = exp(-i*q*gamma);
C_(1+nw) = C;
end;%for nw=0:n_w_max-1;

ic=0;
for nr=0:n_r-1;
if (n_w_(1+nr)>0);
tmp_S_p_ = S_p_(1+ic+(0:n_w_(1+nr)-1));
tmp_S_q_ = fft(tmp_S_p_);
%disp(sprintf(' %% nr %d; n_w_(1+nr) %d: tmp_S_p_ [%d %d] tmp_S_q_ [%d %d]',nr,n_w_(1+nr),size(tmp_S_p_),size(tmp_S_q_)));
for nq=0:n_w_(1+nr)-1;
q = nq;
if (q>floor(n_w_(1+nr)/2)-1) ;
q = q - n_w_(1+nr);
end;%if (q>floor(n_w_(1+nr)/2)-1) ;
C = C_(1+floor(n_w_max/2) + q);
tmp_S_q_(1+nq) = tmp_S_q_(1+nq) * C;
end;%for nq=0:n_w_(1+nr)-1;
tmp_S_p_ = ifft(tmp_S_q_);
M_p_(1+ic+(0:n_w_(1+nr)-1)) = tmp_S_p_;
ic = ic + n_w_(1+nr);
end;%if (n_w_(1+nr)>0);
end;%for nr=0:n_r-1;
