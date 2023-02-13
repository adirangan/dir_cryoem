function A_w_ = MSA_shape_diffuse_0(n_w,gamma_w_,A_w_,t_diffuse);
%%%%%%%%;
% diffuses assuming periodic boundary conditions and mod(n_w,2)==0. ;
%%%%%%%%;
n_q = n_w/2;
q_ = periodize(transpose([0:n_w-1]),-n_w/2,+n_w/2);
A_q_ = fft(A_w_);
A_q_ = A_q_.*exp(-t_diffuse*q_.^2);
A_w_ = ifft(A_q_);
