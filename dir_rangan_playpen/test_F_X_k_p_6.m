function X_k_p_ = test_F_X_k_p_6(r_d_,r_w_,delta_d,delta_w,...
				 n_svd_r,svd_r_,svd_r_m,svd_r_c,svd_r_w_,svd_r_Jv_,n_svd_d,svd_d_,svd_d_m,svd_d_c,svd_d_w_,svd_d_Jv_,n_svd_l,svd_l_,svd_U_d_,svd_s_,svd_V_r_...
				 );
n_r_d = length(r_d_); r_d_ = reshape(r_d_,1,n_r_d);
n_r_w = length(r_w_); r_w_ = reshape(r_w_,1,n_r_w);
verbose=0;
a_K = length(svd_r_w_); b_K = length(svd_d_w_);
phi_k_p_ = repmat(transpose(r_w_),1,n_r_d)-delta_w;
rho_k_p_ = repmat(r_d_,n_r_w,1);
X_k_p_ = zeros(size(phi_k_p_));
for ns=1:n_svd_l;
l = svd_l_(ns);
S = svd_s_(ns);
U_d = 0;
for nkB=0:b_K-1;
b_tmp = svd_d_Jv_{1+nkB}((delta_d - svd_d_m)/svd_d_c);
U_d = U_d + svd_U_d_(1+nkB,ns)*b_tmp;
end;% for nkB=0:b_K-1;
V_r_ = zeros(1,n_r_d);
for nkA=0:a_K-1;
a_tmp = svd_r_Jv_{1+nkA}((r_d_ - svd_r_m)/svd_r_c);
V_r_ = V_r_ + svd_V_r_(1+nkA,ns)*a_tmp;
end;%for nkA=0:a_K-1;
X_k_p_ = X_k_p_ + exp(-i*l*pi/2).*U_d.*S.*(ones(n_r_w,1)*V_r_).*exp(+i*l*phi_k_p_);
end;%for ns=1:n_svd_l;

