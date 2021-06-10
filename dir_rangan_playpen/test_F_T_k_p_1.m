function T_k_p_ = test_F_T_k_p_1(r_d_,r_w_,delta_d,delta_w);
verbose=0;
n_r_d = length(r_d_); r_d_ = reshape(r_d_,1,n_r_d);
n_r_w = length(r_w_); r_w_ = reshape(r_w_,1,n_r_w);
phi_k_p_ = repmat(transpose(r_w_),1,n_r_d)-delta_w;
rho_k_p_ = repmat(r_d_,n_r_w,1);
delta_x = delta_d*cos(delta_w); delta_y = delta_d*sin(delta_w);
if (verbose); disp(sprintf(' %% delta_x_c (%f,%f) delta_x_p (%f,%f)',delta_x,delta_y,delta_d,delta_w)); end;
T_k_p_ = exp(-i*delta_d*rho_k_p_.*cos(phi_k_p_));
