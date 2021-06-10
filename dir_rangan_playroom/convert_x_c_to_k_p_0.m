function [a_k_p_] = convert_x_c_to_k_p_0(verbosen,n_k_all,k_c_0_all,k_c_1_all,k_c_2_all,X_0_,X_1_,X_2_,a_x_c_);
n_x_c = numel(X_0_); assert(numel(X_1_)<=n_x_c); assert(numel(X_2_)<=n_x_c);
a_k_p_ = nufft3d3(n_x_c,pi*X_0_(:),pi*X_1_(:),pi*X_2_(:),a_x_c_,-1,1e-12,n_k_all,2*pi*k_c_0_all_/pi,2*pi*k_c_1_all_/pi,2*pi*k_c_2_all_/pi) / n_x_c / (2*pi);
