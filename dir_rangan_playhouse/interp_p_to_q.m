function S_q_ = interp_p_to_q(n_r,n_w_,n_A,S_p_);

if (nargin<4);
test_transforms_dr.m();
disp('returning'); return;
end;%if (nargin<4);

S_q_ = zeros(n_A,1);
n_w_max = n_w_(1+n_r-1);
ic=0;
for nr=0:n_r-1;
if (n_w_(1+nr)>0);
ij_tmp = ic + (0:n_w_(1+nr)-1);
S_q_(1+ij_tmp) = fft(S_p_(1+ij_tmp))/sqrt(n_w_(1+nr));
ic = ic + n_w_(1+nr);
end;%if (n_w_(1+nr)>0);
end;%for nr=0:n_r-1;
