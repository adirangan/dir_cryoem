function a_x_c_ = convert_k_p_to_x_c_0(verbose,n_k_all,k_c_0_all,k_c_1_all,k_c_2_all,weight_k_all_,a_k_all_,k_p_r_max,res);

k_c_0_all_rescaled_ = pi*k_c_0_all_/k_p_r_max ;
k_c_1_all_rescaled_ = pi*k_c_1_all_/k_p_r_max ;
k_c_2_all_rescaled_ = pi*k_c_2_all_/k_p_r_max ;
n_m = res.^3;
m_x_ = linspace(-k_p_r_max/2,+k_p_r_max/2,res+1); m_x_ = m_x_(1:res);
m_y_ = linspace(-k_p_r_max/2,+k_p_r_max/2,res+1); m_y_ = m_y_(1:res);
m_z_ = linspace(-k_p_r_max/2,+k_p_r_max/2,res+1); m_z_ = m_z_(1:res);
[M_X_,M_Y_,M_Z_] = meshgrid(m_x_,m_y_,m_z_);
[a_x_c_,ier] = nufft3d3(n_k_all,k_c_0_all_rescaled_,k_c_1_all_rescaled_,k_c_2_all_rescaled_,a_all_.*weight_all_,+1,1e-12,n_m,M_X_(:),M_Y_(:),M_Z_(:));
a_x_c_ = real(reshape(a_x_c_,res,res,res));
