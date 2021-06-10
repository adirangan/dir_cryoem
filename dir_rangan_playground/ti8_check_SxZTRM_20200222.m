% checking ti8_dr_digits for a particular case where svd seems to give very wrong answer (i.e., 30% error). ;
clear;setup;
dir_tmp = '/data/rangan/dir_cryoem/dir_rangan_playground/dir_tmp';
n_delta_v = MDA_read_i4(sprintf('%s/n_delta_v.mda',dir_tmp));
delta_x_ = MDA_read_r8(sprintf('%s/delta_x_.mda',dir_tmp));
delta_y_ = MDA_read_r8(sprintf('%s/delta_y_.mda',dir_tmp));
quad_n_r = MDA_read_i4(sprintf('%s/quad_n_r.mda',dir_tmp));
quad_n_w_ = MDA_read_i4(sprintf('%s/quad_n_w_.mda',dir_tmp));
quad_n_w_sum = MDA_read_i4(sprintf('%s/quad_n_w_sum.mda',dir_tmp));
quad_grid_k_p_r_ = MDA_read_r8(sprintf('%s/quad_grid_k_p_r_.mda',dir_tmp));
quad_weight_k_p_r_ = MDA_read_r8(sprintf('%s/quad_weight_k_p_r_.mda',dir_tmp));
gamma_z = MDA_read_r8(sprintf('%s/gamma_z.mda',dir_tmp));
delta_x = MDA_read_r8(sprintf('%s/delta_x.mda',dir_tmp));
delta_y = MDA_read_r8(sprintf('%s/delta_y.mda',dir_tmp));
quad_S_k_p_ = MDA_read_c16(sprintf('%s/quad_S_k_p_.mda',dir_tmp));
quad_M_k_q_ = MDA_read_c16(sprintf('%s/quad_M_k_q_.mda',dir_tmp));
assert(quad_n_w_sum==numel(quad_S_k_p_));
print_sub(quad_n_w_sum,quad_S_k_p_,' 2 quad_S_k_p_: ');
quad_R_k_p_ = rotate_p_to_p_lagrange(quad_n_r,quad_n_w_,quad_n_w_sum,quad_S_k_p_,gamma_z);
print_sub(quad_n_w_sum,quad_R_k_p_,' quad_R_k_p_: ');
quad_T_k_p_ = rotate_p_to_p_fftw(quad_n_r,quad_n_w_,quad_n_w_sum,quad_S_k_p_,gamma_z);
print_sub(quad_n_w_sum,quad_T_k_p_,' 2 quad_T_k_p_: ');
flag_plot=0;
if flag_plot;
subplot(2,3,1); clim_ = imagesc_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum,real(quad_S_k_p_),[],colormap_beach()); 
axis image; set(gca,'XTick',[],'YTick',[]); title('real S');
subplot(2,3,2); imagesc_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum,real(quad_T_k_p_),clim_,colormap_beach());
axis image; set(gca,'XTick',[],'YTick',[]); title('real T');
subplot(2,3,3); imagesc_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum,real(quad_R_k_p_),clim_,colormap_beach());
axis image; set(gca,'XTick',[],'YTick',[]); title('real R');
subplot(2,3,4); clim_ = imagesc_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum,imag(quad_S_k_p_),[],colormap_beach());
axis image; set(gca,'XTick',[],'YTick',[]); title('imag S');
subplot(2,3,5); imagesc_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum,imag(quad_T_k_p_),clim_,colormap_beach());
axis image; set(gca,'XTick',[],'YTick',[]); title('imag T');
subplot(2,3,6); imagesc_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum,imag(quad_R_k_p_),clim_,colormap_beach());
axis image; set(gca,'XTick',[],'YTick',[]); title('imag R');
figbig;
end;%if flag_plot;
quad_T_k_p_ = transf_p_to_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum,quad_T_k_p_,delta_x,delta_y);
print_sub(quad_n_w_sum,quad_T_k_p_,'  2quad_T_k_p_: ');
quad_T_k_q_ = interp_p_to_q(quad_n_r,quad_n_w_,quad_n_w_sum,quad_T_k_p_);
quad_T_k_q_2_ = quad_T_k_q_;
print_sub(quad_n_w_sum,quad_T_k_q_,' 2 quad_T_k_q_: ');
print_sub(quad_n_w_sum,quad_M_k_q_,' 2 quad_M_k_q_: ');
tmp_ori = 0.5*innerproduct_p_quad(quad_n_r,quad_grid_k_p_r_,quad_weight_k_p_r_,quad_n_w_,quad_n_w_sum,quad_T_k_q_,quad_M_k_q_);
disp(sprintf(' %% S_R_T_x_T_R_CTF_M_q_ori: %0.16f + %0.16fi',real(tmp_ori),imag(tmp_ori)));
%%%%%%%%;
quad_T_k_q_ = interp_p_to_q(quad_n_r,quad_n_w_,quad_n_w_sum,quad_S_k_p_);
print_sub(quad_n_w_sum,quad_T_k_q_,' 1 quad_T_k_q_: ');
quad_T_k_q_ = rotate_q_to_q(quad_n_r,quad_n_w_,quad_n_w_sum,quad_T_k_q_,gamma_z);
print_sub(quad_n_w_sum,quad_T_k_q_,' 1 quad_T_k_q_: ');
n_pixel = 4.5; half_diameter_k_c = 16.0; delta_max = n_pixel/sqrt(2)/half_diameter_k_c;
eps_target = 1e-6; flag_warning = 1; dir_svd = '/data/rangan/dir_cryoem/dir_rangan_playground/dir_gen_Jsvd_6';
FTK = get_svd_FTK_2(eps_target,quad_grid_k_p_r_,quad_n_r,1,delta_max,0,flag_warning,dir_svd);
FTK.n_delta_v = 1;
FTK.delta_x_ = delta_x;
FTK.delta_y_ = delta_y;
FTK.svd_d_max = delta_max;
FTK.svd_polyval_U_d_ = get_svd_polyval_U_d_0(FTK.svd_d_max,FTK.n_svd_d,FTK.svd_d_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_U_d_,FTK.n_delta_v,FTK.delta_x_,FTK.delta_y_);
FTK.svd_r_max = quad_grid_k_p_r_(1+quad_n_r-1);
FTK.svd_polyval_V_r_ = get_svd_polyval_V_r_0(FTK.svd_r_max,FTK.n_svd_r,FTK.svd_r_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_V_r_,quad_n_r,quad_grid_k_p_r_);
quad_T_k_q_ = transf_svd_q_to_q_FTK_5(FTK,quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum,quad_T_k_q_,delta_x,delta_y);
quad_T_k_q_1_ = quad_T_k_q_;
print_sub(quad_n_w_sum,quad_T_k_q_,' 1 quad_T_k_q_: ');
tmp_svd = 0.5*innerproduct_p_quad(quad_n_r,quad_grid_k_p_r_,quad_weight_k_p_r_,quad_n_w_,quad_n_w_sum,quad_T_k_q_,quad_M_k_q_);
disp(sprintf(' %% S_R_T_x_T_R_CTF_M_q_svd: %0.16f + %0.16fi',real(tmp_svd),imag(tmp_svd)));

quad_T_k_p_1_ = interp_q_to_p(quad_n_r,quad_n_w_,quad_n_w_sum,quad_T_k_q_1_);
quad_T_k_p_2_ = interp_q_to_p(quad_n_r,quad_n_w_,quad_n_w_sum,quad_T_k_q_2_);
quad_T_k_p_3_ = quad_T_k_p_1_ - quad_T_k_p_2_;
subplot(2,3,1); clim_ = imagesc_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum,real(quad_T_k_p_1_),[],colormap_beach()); 
axis image; set(gca,'XTick',[],'YTick',[]); title('real T1');
subplot(2,3,2); imagesc_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum,real(quad_T_k_p_2_),clim_,colormap_beach());
axis image; set(gca,'XTick',[],'YTick',[]); title('real T2');
subplot(2,3,3); imagesc_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum,real(quad_T_k_p_3_),clim_,colormap_beach());
axis image; set(gca,'XTick',[],'YTick',[]); title('real T3');
subplot(2,3,4); clim_ = imagesc_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum,imag(quad_T_k_p_1_),[],colormap_beach()); 
axis image; set(gca,'XTick',[],'YTick',[]); title('imag T1');
subplot(2,3,5); imagesc_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum,imag(quad_T_k_p_2_),clim_,colormap_beach());
axis image; set(gca,'XTick',[],'YTick',[]); title('imag T2');
subplot(2,3,6); imagesc_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum,imag(quad_T_k_p_3_),clim_,colormap_beach());
axis image; set(gca,'XTick',[],'YTick',[]); title('imag T3');
figbig;
quad_T_k_q_3_ = interp_p_to_q(quad_n_r,quad_n_w_,quad_n_w_sum,quad_T_k_p_3_);
tmp_dif = 0.5*innerproduct_p_quad(quad_n_r,quad_grid_k_p_r_,quad_weight_k_p_r_,quad_n_w_,quad_n_w_sum,quad_T_k_q_3_,quad_M_k_q_);
quad_M_k_p_ = interp_q_to_p(quad_n_r,quad_n_w_,quad_n_w_sum,quad_M_k_q_);











