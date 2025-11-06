function ...
[ ...
 P_k_p_wk_ ...
,T_k_p_wk_ ...
] = ...
test_F_P_k_p_7( ...
 n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,delta_d ...
,delta_w ...
,FTK ...
);

n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);

index_wk_from_k_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
index_wk_from_k_(1+n_w_csum_(1+nk_p_r)+[0:n_w-1]) = nk_p_r;
end;%for nk_p_r=0:n_k_p_r-1;
k_k_p_wk_ = k_p_r_(1+index_wk_from_k_);

w_k_p_wk_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
w_k_p_ = 2*pi*transpose([0:n_w-1])/max(1,n_w);
w_k_p_wk_(1+n_w_csum_(1+nk_p_r)+[0:n_w-1]) = w_k_p_;
end;%for nk_p_r=0:n_k_p_r-1;

T_k_p_wk_ = exp(-i.*2.*pi.*k_k_p_wk_.*delta_d.*cos(w_k_p_wk_ - delta_w));

%%%%%%%%;
n_r_degree = FTK.n_r_degree;
n_d_degree = FTK.n_d_degree;
n_svd_r = FTK.n_svd_r;
svd_r_ = FTK.svd_r_;
svd_r_m = FTK.svd_r_m;
svd_r_c = FTK.svd_r_c;
svd_r_w_ = FTK.svd_r_w_;
svd_r_Jv_ = FTK.svd_r_Jv_;
n_svd_d = FTK.n_svd_d;
svd_d_ = FTK.svd_d_;
svd_d_m = FTK.svd_d_m;
svd_d_c = FTK.svd_d_c;
svd_d_w_ = FTK.svd_d_w_;
svd_d_Jv_ = FTK.svd_d_Jv_;
n_svd_l = FTK.n_svd_l;
svd_l_ = FTK.svd_l_;
svd_U_d_jacocoef_ = FTK.svd_U_d_jacocoef_;
svd_s_ = FTK.svd_s_;
svd_V_r_jacocoef_ = FTK.svd_V_r_jacocoef_;
svd_U_d_ = FTK.svd_U_d_;
svd_V_r_ = FTK.svd_V_r_;
n_delta_v = FTK.n_delta_v;
delta_x_ = FTK.delta_x_;
delta_y_ = FTK.delta_y_;
svd_d_max = FTK.svd_d_max;
svd_polyval_U_d_ = FTK.svd_polyval_U_d_;
svd_r_max = FTK.svd_r_max;
svd_polyval_V_r_ = FTK.svd_polyval_V_r_;
svd_expiw__ = FTK.svd_expiw__;
svd_U_d_expiw_s__ = FTK.svd_U_d_expiw_s__;
%%%%%%%%;

%%%%%%%%;
n_b = n_svd_d; n_a = n_svd_r;
svd_U_d_jacocoef_bl__ = reshape(svd_U_d_jacocoef_,[n_svd_d,n_svd_l]);
svd_V_r_jacocoef_al__ = reshape(svd_V_r_jacocoef_,[n_svd_r,n_svd_l]);
%%%%%%%%;
b_b_ = zeros(n_b,1);
for nb=0:n_b-1;
b = svd_d_Jv_{1+nb}((delta_d - svd_d_m)/svd_d_c);
b_b_(1+nb) = b;
end;% for nb=0:n_b-1;
%%%%%%%%;
a_ak__ = zeros(n_a,n_k_p_r);
for na=0:n_a-1;
a_k_ = svd_r_Jv_{1+na}((2*pi*k_p_r_ - svd_r_m)/svd_r_c);
a_ak__(1+na,:) = a_k_;
end;%for na=0:n_a-1;
%%%%%%%%;

%%%%%%%%;
phi_k_p_ = w_k_p_wk_-delta_w;
P_k_p_wk_ = zeros(size(phi_k_p_));
for ns=0:n_svd_l-1;
l = svd_l_(1+ns);
S = svd_s_(1+ns);
U_d = dot(svd_U_d_jacocoef_bl__(:,1+ns),b_b_);
V_r_ = transpose(sum(bsxfun(@times,svd_V_r_jacocoef_al__(:,1+ns),a_ak__),1));
V_wk_ = V_r_(1+index_wk_from_k_);
P_k_p_wk_ = P_k_p_wk_ + exp(-i*l*pi/2).*U_d.*S.*V_wk_.*exp(+i*l*phi_k_p_);
end;%for ns=0:n_svd_l-1;
%%%%%%%%;
