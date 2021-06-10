function T_k_c_ = test_F_T_k_c_0(n_point,max_x_c,delta_d,delta_w);
[max_k_c,max_k_p,grid_x_c_,d_x_c,X_x_c_,Y_x_c_,R_x_c_,W_x_c_,grid_k_c_,d_k_c,X_k_c_,Y_k_c_,R_k_c_,W_k_c_,grid_x_r_,grid_x_w_,R_x_p_,W_x_p_,grid_k_r_,d_k_r,grid_k_w_,d_k_w,R_k_p_,W_k_p_,X_k_p_,Y_k_p_] = test_F_grid_0(n_point,max_x_c);
verbose=0;
delta_x = delta_d*cos(delta_w); delta_y = delta_d*sin(delta_w);
if (verbose); disp(sprintf(' %% delta_x_c (%f,%f) delta_x_p (%f,%f)',delta_x,delta_y,delta_d,delta_w)); end;
T_k_c_ = exp(-i*delta_d*2*pi*R_k_c_.*cos(W_k_c_-delta_w));
