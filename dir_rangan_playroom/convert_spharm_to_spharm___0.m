function a_k_Y_mlk___ = convert_spharm_to_spharm___0(n_k_p_r,l_max_,a_k_Y_mlk_);
% converts a linearly-ordered spherical-harmonic-expansion into a stacked array. ;
% More specifically, if a_k_Y_mlk_ stores the components of a such that m_val varies ;
% quickly then l_val varies, and final nk_p_r varies most slowly, ;
% then a_k_Y_mlk___ stores the components so that: ;
% m_val corresponds to the rows, ;
% l_val corresponds to the cols, ;
% nk_p_r corresponds to the layers (dimension 3). ;

l_max_max = max(l_max_);
n_lm_ = (1+l_max_).^2;

ixk=0;
for nk=0:n_k_p_r-1;
ixl=0;
for l_val=0:l_max_(1+nk);
a_k_Y_mlk___(1+l_max_max + [-l_val:+l_val],1+l_val,1+nk) = ...
  a_k_Y_mlk_(1+ixk+ixl+l_val+-[-l_val:+l_val]);
ixl = ixl+(2*l_val+1);
end;%for l_val=0:l_max_(1+nk);
assert(ixl==n_lm_(1+nk));
ixk = ixk+n_lm_(1+nk);
end;%for nk=0:n_k_p_r-1;
assert(ixk==sum(n_lm_));
