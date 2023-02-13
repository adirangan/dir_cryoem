
%%%%%%%%;
% test alternative mult when n_M and n_S are full. ;
%%%%%%%%;
nUX_rank = 1;
pm_n_UX_rank = 1+nUX_rank;
if (verbose); disp(sprintf(' %% pm_n_UX_rank %d/%d',pm_n_UX_rank,n_UX_rank)); end;
tmp_UX_weight_kn__ = diag(X_weight_r_)*UX_S_avg_kn__(:,1:pm_n_UX_rank);
tmp_SX_k_ = SX_S_avg_k_(1:pm_n_UX_rank);
UX_weight_CTF_S_k_q_nwS___ = reshape(transpose(tmp_UX_weight_kn__)*reshape(CTF_S_k_q_kwS___,[n_k_p_r,n_w_max*n_S]),[pm_n_UX_rank,n_w_max,n_S]);
UX_weight_conj_CTF_S_k_q_Snw___ = conj(permute(UX_weight_CTF_S_k_q_nwS___,[3,1,2]));
UX_weight_M_k_q_nwM___ = reshape(transpose(tmp_UX_weight_kn__)*reshape(M_k_q_kwM___,[n_k_p_r,n_w_max*n_M]),[pm_n_UX_rank,n_w_max,n_M]);
UX_weight_M_k_q_nMw___ = permute(UX_weight_M_k_q_nwM___,[1,3,2]);


tmp_t = tic();
conj_S_CTF_weight_weight_M_k_q_SMw___ = zeros(n_S,n_M,n_w_max);
for nw=0:n_w_max-1;
conj_S_CTF_weight_weight_M_k_q_SMw___(:,:,1+nw) = UX_weight_conj_CTF_S_k_q_Snw___(:,:,1+nw) * UX_weight_M_k_q_nMw___(:,:,1+nw);
end;%for nw=0:n_w_max-1;
tmp_o = n_w_max * n_S * n_M * pm_n_UX_rank;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% mult: %0.3fs --> %0.2fGHz',tmp_t,tmp_o/tmp_t/1e9)); end;
%%%%;
tmp_t = tic();
conj_S_CTF_weight_weight_M_k_q_SMw___ = zeros(n_S,n_M,n_w_max);
tmp_a_ = reshape(UX_weight_conj_CTF_S_k_q_Snw___,[n_S,pm_n_UX_rank,n_w_max]);
tmp_b_ = reshape(UX_weight_M_k_q_nMw___,[pm_n_UX_rank,n_M,n_w_max]);
tmp_c_ = tmp_a_(:,1,:) .* tmp_b_(1,:,:) + tmp_a_(:,2,:) .* tmp_b_(2,:,:);
tmp_o = n_w_max * n_S * n_M * pm_n_UX_rank;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% mult: %0.3fs --> %0.2fGHz',tmp_t,tmp_o/tmp_t/1e9)); end;
