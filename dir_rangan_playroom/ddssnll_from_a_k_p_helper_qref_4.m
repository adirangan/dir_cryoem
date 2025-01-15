%%%%%%%%;
% Check spherical-grid for consistency. ;
% For efficiency the spherical-grid should be generated using: ;
% flag_uniform_over_n_k_p_r = 1;
% flag_uniform_over_polar_a = 0;
% str_T_vs_L = 'C2';
% For accuracy one might consider increasing the angular-resolution (on each shell) ;
% without necessarily increasing the radial-resolution (i.e., the distance between shells). ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% checking spherical-grid for consistency.')); end;
tmp_t = tic();
if n_k_p_r> 1;
tmp_std = std(diff(n_qk_csum_));
if (flag_verbose>0); disp(sprintf(' %% std(diff(n_qk_csum_)) %0.6f',tmp_std)); end;
if tmp_std> 1e-12; disp(sprintf(' %% Warning, std(diff(n_qk_csum_)) %0.6f in %s',tmp_std,str_thisfunction)); end;
end;%if n_k_p_r> 1;
if n_k_p_r> 1;
tmp_std = std(diff(n_polar_a_k_));
if (flag_verbose>0); disp(sprintf(' %% std(diff(n_polar_a_k_)) %0.6f',tmp_std)); end;
if tmp_std> 1e-12; disp(sprintf(' %% Warning, std(diff(n_polar_a_k_)) %0.6f in %s',tmp_std,str_thisfunction)); end;
end;%if n_k_p_r> 1;
polar_a_single_shell_ = polar_a_ka__{1};
n_azimu_b_single_shell_ = n_azimu_b_ka__{1};
tmp_std = std(diff(polar_a_single_shell_));
if (flag_verbose>0); disp(sprintf(' %% std(diff(polar_a_single_shell_)) %0.6f',tmp_std)); end;
if tmp_std> 1e-12; disp(sprintf(' %% Warning, std(diff(polar_a_single_shell_)) %0.6f in %s',tmp_std,str_thisfunction)); end;
n_polar_a_single_shell = numel(polar_a_single_shell_);
polar_a_single_shell_lim_ = [polar_a_single_shell_(1),polar_a_single_shell_(end)];
tab_error=0; tab_check=0;
for npolar_a_single_shell=0:n_polar_a_single_shell-1;
n_azimu_b_single_shell = n_azimu_b_single_shell_(1+npolar_a_single_shell);
if npolar_a_single_shell==0; if n_azimu_b_single_shell~=1; disp(sprintf(' %% Warning, npolar_a_single_shell %d n_azimu_b_single_shell %d',npolar_a_single_shell,n_azimu_b_single_shell)); tab_error = tab_error+1; else; tab_check = tab_check+1; end; end;
if npolar_a_single_shell==n_polar_a_single_shell-1; if n_azimu_b_single_shell~=1; disp(sprintf(' %% Warning, npolar_a_single_shell %d n_azimu_b_single_shell %d',npolar_a_single_shell,n_azimu_b_single_shell)); tab_error = tab_error+1; else; tab_check = tab_check+1; end; end;
if npolar_a_single_shell> 0 & npolar_a_single_shell< n_polar_a_single_shell-1; if mod(n_azimu_b_single_shell,2)~=0; disp(sprintf(' %% Warning, npolar_a_single_shell %d n_azimu_b_single_shell %d',npolar_a_single_shell,n_azimu_b_single_shell)); tab_error = tab_error+1; else; tab_check = tab_check+1; end; end;
end;%for npolar_a_single_shell=0:n_polar_a_single_shell-1;
assert(tab_check==n_polar_a_single_shell);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% check grid: %0.6fs',tmp_t)); end;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% extracting quadrature-grid on shell.')); end;
tmp_t = tic();
n_q_single_shell = n_qk/n_k_p_r; n_3 = 3; n_9 = 9; n_1 = 1;
nk_p_r = n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
weight_3d_k_p_r = weight_3d_k_p_r_(1+nk_p_r);
tmp_index_ = n_qk_csum_(1+nk_p_r)+[0:n_q_single_shell-1];
k_c_0_single_shell_ = k_c_0_qk_(1+tmp_index_)/max(1e-12,k_p_r);
k_c_1_single_shell_ = k_c_1_qk_(1+tmp_index_)/max(1e-12,k_p_r);
k_c_2_single_shell_ = k_c_2_qk_(1+tmp_index_)/max(1e-12,k_p_r);
k_p_azimu_b_single_shell_ = k_p_azimu_b_qk_(1+tmp_index_);
k_p_polar_a_single_shell_ = k_p_polar_a_qk_(1+tmp_index_);
weight_3d_k_single_shell_ = weight_3d_k_p_qk_(1+tmp_index_); %<-- sum(weight_3d_k_single_shell_) = 4*pi*weight_3d_k_p_r;
if (flag_verbose>0); disp(sprintf(' %% sum(weight_3d_k_single_shell_) %0.16f 4*pi*weight_3d_k_p_r %0.16f',sum(weight_3d_k_single_shell_),4*pi*weight_3d_k_p_r)); end;
weight_shell_qk_single_shell_ = weight_shell_qk_(1+tmp_index_); %<-- sum(weight_shell_qk_single_shell_) = 4*pi*k_p_r^2. ;
if (flag_verbose>0); disp(sprintf(' %% sum(weight_shell_qk_single_shell_) %0.16f 4*pi*k_p_r^2 %0.16f',sum(weight_shell_qk_single_shell_),4*pi*k_p_r^2)); end;
weight_shell_qk_unit_shell_ = weight_shell_qk_single_shell_/max(1e-12,k_p_r)^2; %<-- sum(weight_shell_qk_unit_shell_) = 4*pi. ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% set weight_shell_qk_unit_shell_: %0.6fs',tmp_t)); end;
%%%%%%%%;
str_C2 = 'C2';
flag_tensor_vs_adap = 0;
tmp_t = tic();
[ ...
 qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
] = ...
sample_shell_6( ...
 1.0 ...
,qref_k_eq_d ...
,str_C2 ...
,flag_tensor_vs_adap ...
) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% sample_shell_5: %0.2fs',tmp_t)); end;
%%%%%%%%;
fnorm_disp(flag_verbose,'qref_n_shell',qref_n_shell,'n_q_single_shell',n_q_single_shell);
fnorm_disp(flag_verbose,'qref_azimu_b_shell_',qref_azimu_b_shell_,'k_p_azimu_b_single_shell_',k_p_azimu_b_single_shell_);
fnorm_disp(flag_verbose,'qref_polar_a_shell_',qref_polar_a_shell_,'k_p_polar_a_single_shell_',k_p_polar_a_single_shell_);
fnorm_disp(flag_verbose,'qref_weight_shell_',qref_weight_shell_,'weight_shell_qk_unit_shell_',weight_shell_qk_unit_shell_);
fnorm_disp(flag_verbose,'qref_k_c_0_shell_',qref_k_c_0_shell_,'k_c_0_single_shell_',k_c_0_single_shell_);
fnorm_disp(flag_verbose,'qref_k_c_1_shell_',qref_k_c_1_shell_,'k_c_1_single_shell_',k_c_1_single_shell_);
fnorm_disp(flag_verbose,'qref_k_c_2_shell_',qref_k_c_2_shell_,'k_c_2_single_shell_',k_c_2_single_shell_);
fnorm_disp(flag_verbose,'qref_n_polar_a',qref_n_polar_a,'n_polar_a_single_shell',n_polar_a_single_shell);
fnorm_disp(flag_verbose,'qref_polar_a_',qref_polar_a_,'polar_a_single_shell_',polar_a_single_shell_);
fnorm_disp(flag_verbose,'qref_n_azimu_b_',qref_n_azimu_b_,'n_azimu_b_single_shell_',n_azimu_b_single_shell_);
%%%%%%%%;
qref_n_shell = n_q_single_shell;
qref_azimu_b_shell_ = k_p_azimu_b_single_shell_;
qref_polar_a_shell_ = k_p_polar_a_single_shell_;
qref_weight_shell_ = weight_shell_qk_unit_shell_;
qref_k_c_0_shell_ = k_c_0_single_shell_;
qref_k_c_1_shell_ = k_c_1_single_shell_;
qref_k_c_2_shell_ = k_c_2_single_shell_;
qref_n_polar_a = n_polar_a_single_shell;
qref_polar_a_ = polar_a_single_shell_;
qref_n_azimu_b_ = n_azimu_b_single_shell_;
%%%%%%%%;
