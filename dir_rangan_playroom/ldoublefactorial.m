function output_ = ldoublefactorial(n_);
output_ = zeros(size(n_));
tmp_ij_ = find(mod(n_,2)==0);
k_ = n_(tmp_ij_)/2;
output_(tmp_ij_) = lfactorial(k_) + k_.*log(2);
tmp_ij_ = find(mod(n_,2)==1);
k_ = (n_(tmp_ij_)+1)/2;
output_(tmp_ij_) = lfactorial(2*k_-1) - (k_-1)*log(2) - lfactorial(k_-1);


