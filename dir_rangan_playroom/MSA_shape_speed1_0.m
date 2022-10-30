function B_w_ = MSA_shape_speed1_0(n_w,gamma_w_,A_w_);
%%%%%%%%;
% Reparametrizes a single (complex) image-ring by arclength. ;
% Assumes that gamma_w_ is periodic. ;
%%%%%%%%;
dgamma_w_ = [gamma_w_(2:end);gamma_w_(1+0)] - [gamma_w_(end);gamma_w_(1:end-1)];
dA_w_ = [A_w_(2:end);A_w_(1+0)] - [A_w_(end);A_w_(1:end-1)];
s_w_ = cumsum(abs(dA_w_));
B_w_ = interp1(s_w_,A_w_,transpose(linspace(s_w_(1+0),s_w_(end),n_w)));



