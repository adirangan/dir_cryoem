% test cg_lsq. ;
flag_verbose = 1;
tolerance_cg_lsq = 1e-6;
n_n = 13;
An__ = randn(n_n,n_n);
At__ = transpose(An__);
AtA__ = At__*An__;
lhs_tru_ = randn(n_n,1);
rhs_tru_ = AtA__*lhs_tru_;

% Conjugate Gradient Loop
n_iteration = n_n;
niteration=0;
lhs_tmp_ = zeros(n_n,1);
res_ = rhs_tru_ - AtA__*lhs_tmp_;
pcg_ = res_;
beta_num = sum(res_.^2);
beta_den = 1.0;
beta = beta_num/max(1e-12,beta_den);
flag_continue = 1;
while flag_continue
if flag_verbose; fprintf(' %% niteration = %d beta_num = %f\n', niteration, beta_num); end;
zeta = sum((An__*pcg_).^2,'all');
alph = beta_num / max(1e-12, zeta);
lhs_tmp_ = lhs_tmp_ + alph*pcg_ ;
res_ = res_ - alph * At__*(An__*pcg_);
beta_den = beta_num;
beta_num = sum(res_.^2);
beta = beta_num / max(1e-12,beta_den);
pcg_ = res_ + beta * pcg_ ;
niteration = niteration + 1;
if niteration >= n_iteration; flag_continue=0; end;
if sqrt(beta) < tolerance_cg_lsq; flag_continue=0; end;
end;%while;
disp(sprintf(' %% final error: %0.6f',fnorm(lhs_tru_ - lhs_tmp_)/max(1e-12,fnorm(lhs_tru_))));
