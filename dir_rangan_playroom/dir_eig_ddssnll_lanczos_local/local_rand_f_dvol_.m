function [f_dvol_yk_,f_dvol_yk__] = local_rand_f_dvol_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
scaling_volumetric = (4*pi)^2 * sqrt(pi/2);
f_dvol_yk_ = randn(n_lm_sum,1) + i*randn(n_lm_sum,1);
%%%%%%%%;
na=0;
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
for l_val=0:+l_max;
for m_val=-l_val:+l_val;
if m_val== 0;
if mod(l_val,2)==0; f_dvol_yk_(1+na) = 2*real(f_dvol_yk_(1+na)); end;
if mod(l_val,2)==1; f_dvol_yk_(1+na) = i*2*imag(f_dvol_yk_(1+na)); end;
end;%if m_val== 0;
if m_val> 0;
tmp_index_pos = n_lm_csum_(1+nk_p_r) + l_val^2 + l_val + m_val ;
assert(tmp_index_pos==na);
tmp_index_neg = n_lm_csum_(1+nk_p_r) + l_val^2 + l_val - m_val ;
tmp_factor = +1; if mod(l_val,2)==1; tmp_factor = -1; end;
f_dvol_yk_(1+tmp_index_neg) = tmp_factor*conj(f_dvol_yk_(1+tmp_index_pos));
end;%if m_val> 0;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:+l_max;
end;%for nk_p_r=0:n_k_p_r-1;
assert(na==n_lm_sum);
%%%%%%%%;
tmp_ff_dvol = local_f_dvol_bar_dot_g_dvol_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,f_dvol_yk_,f_dvol_yk_);
f_dvol_yk_ = f_dvol_yk_/max(1e-12,sqrt(tmp_ff_dvol));
f_dvol_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,f_dvol_yk_);
