function X_k_p_ = test_F_X_k_p_5(n_point,max_x_c,delta_d,delta_w,...
				 n_svd_r,svd_r_,svd_r_m,svd_r_c,svd_r_w_,svd_r_Jv_,n_svd_d,svd_d_,svd_d_m,svd_d_c,svd_d_w_,svd_d_Jv_,n_svd_l,svd_l_,svd_U_d_,svd_s_,svd_V_r_...
				 );
[max_k_c,max_k_p,grid_x_c_,d_x_c,X_x_c_,Y_x_c_,R_x_c_,W_x_c_,grid_k_c_,d_k_c,X_k_c_,Y_k_c_,R_k_c_,W_k_c_,grid_x_r_,grid_x_w_,R_x_p_,W_x_p_,grid_k_r_,d_k_r,grid_k_w_,d_k_w,R_k_p_,W_k_p_,X_k_p_,Y_k_p_] = test_F_grid_0(n_point,max_x_c);

verbose=0;
a_K = length(svd_r_w_); b_K = length(svd_d_w_);
phi_k_p_ = W_k_p_-delta_w;
X_k_p_ = zeros(size(phi_k_p_));
for ns=1:n_svd_l;
l = svd_l_(ns);
S = svd_s_(ns);
U_d = 0;
for nkB=0:b_K-1;
b_tmp = svd_d_Jv_{1+nkB}((delta_d - svd_d_m)/svd_d_c);
U_d = U_d + svd_U_d_(1+nkB,ns)*b_tmp;
end;% for nkB=0:b_K-1;
V_r_ = zeros(1,length(grid_k_r_));
for nkA=0:a_K-1;
a_tmp = svd_r_Jv_{1+nkA}((2*pi*grid_k_r_ - svd_r_m)/svd_r_c);
V_r_ = V_r_ + svd_V_r_(1+nkA,ns)*a_tmp;
end;%for nkA=0:a_K-1;
X_k_p_ = X_k_p_ + exp(-i*l*pi/2).*U_d.*S.*(ones(n_point,1)*V_r_).*exp(+i*l*phi_k_p_);
end;%for ns=1:n_svd_l;

