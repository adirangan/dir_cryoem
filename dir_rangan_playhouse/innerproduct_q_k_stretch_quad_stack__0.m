function M_q_rwM___ = innerproduct_q_k_stretch_quad_stack__0(n_r,weight_p_r_,n_w_,n_M,M_q__) ;
% Assumes quasi-uniform polar-grid. ;
% Assumes that n_w_ is a list of even integers. ;
n_w_ = n_w_(:);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
M_q_rwM___ = zeros(n_r,n_w_max,n_M);
for nM=0:n_M-1;
for nr=0:n_r-1;
n_w = n_w_(1+nr);
n_w_2 = round(n_w/2);
n_w_csum = n_w_csum_(1+nr);
dAn = weight_p_r_(1+nr)*(2*pi)/max(1,n_w);
M_q_rwM___(1+nr,1+(0:n_w_2-1),1+nM) = M_q__(1+n_w_csum_(1+nr)+(0:n_w_2-1),1+nM)*sqrt(dAn);
M_q_rwM___(1+nr,1+(n_w_2+1:n_w-1)-n_w+n_w_max,1+nM) = M_q__(1+n_w_csum_(1+nr)+(n_w_2+1:n_w-1),1+nM)*sqrt(dAn);
end;%for nr=0:n_r-1;
end;%for nM=0:n_M-1;
