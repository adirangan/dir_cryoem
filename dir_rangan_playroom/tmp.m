
for i in $(squeue -u rangan -h -t PD -o %i)
do
scontrol update jobid=$i partition=ccm
done

figure(1);figmed;
x_ = 0:0.01:1;
subplot(1,2,1);
plot(x_,8*x_./(2+6*x_),'-','LineWidth',2);
xlabel('infected fraction where you live');
ylabel('probability you had covid');
title('I had a fever after vaccine!');
grid on; xlim([0,1]); ylim([0,1]); axis square;
subplot(1,2,2);
plot(x_,64*x_./(4+60*x_),'-','LineWidth',2);
xlabel('infected fraction where you live');
ylabel('probability you both had covid');
title('Me and Partner both had fever after vaccine!');
grid on; xlim([0,1]); ylim([0,1]); axis square;


% testing ampmut_wrap_wrap_3.m ;
 
fname_pre = sprintf('%s_mat/ut%.4dn%.2dr%d',dir_pm,floor(1000*delta_r_max_use),dat_n_UX_rank,dat_rseed);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);

if ( exist(fname_mat,'file'));

disp(sprintf(' %% %s found, aligning',fname_mat));
tmp_ = load(fname_mat);
if (~isfield(tmp_.parameter,'fname_pre')); tmp_.parameter.fname_pre = fname_pre; end; %<-- parameter_bookmark. ;
%%%%%%%%;
tmp_fname_pre = sprintf('%s_align_a_CTF_avg_UX_Y_',tmp_.parameter.fname_pre);
tmp_.parameter.fname_align_a_CTF_avg_UX_Y_pre = tmp_fname_pre;
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
disp(sprintf(' %% %s not found, creating',tmp_fname_mat));
ampmut_align_to_a_CTF_avg_UX_Y_0( ...
 tmp_.parameter ...
,l_max_max ...
,dat_n_UX_rank ...
,reshape(a_CTF_avg_UX_Y_quad__(:,1:dat_n_UX_rank),[n_lm_max*dat_n_UX_rank,1]) ...
,n_M ...
,tmp_euler_polar_a_true_ ...
,tmp_euler_azimu_b_true_ ...
,tmp_euler_gamma_z_true_ ...
,tmp_image_delta_x_true_ ...
,tmp_image_delta_y_true_ ...
,[] ...
,tmp_.a_CTF_avg_UX_Y__ ...
,tmp_.euler_polar_a__ ...
,tmp_.euler_azimu_b__ ...
,tmp_.euler_gamma_z__ ...
,tmp_.image_delta_x_acc__ + tmp_.image_delta_x_upd__ ...
,tmp_.image_delta_y_acc__ + tmp_.image_delta_y_upd__ ...
);
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
%%%%%%%%;
tmp_fname_pre = sprintf('%s_align_a_k_Y_',tmp_.parameter.fname_pre);
tmp_.parameter.fname_align_a_k_Y_pre = tmp_fname_pre;
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
disp(sprintf(' %% %s not found, creating',tmp_fname_mat));
ampmut_align_to_a_k_Y_0( ...
 tmp_.parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_M ...
,M_k_p__ ...
,CTF_index_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_quad_ ...
,[] ...
,tmp_.euler_polar_a__ ...
,tmp_.euler_azimu_b__ ...
,tmp_.euler_gamma_z__ ...
,tmp_.image_delta_x_acc__ + tmp_.image_delta_x_upd__ ...
,tmp_.image_delta_y_acc__ + tmp_.image_delta_y_upd__ ...
);
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
%%%%%%%%;

end;%if ( exist(fname_mat,'file'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

load('ampmut_1_debug.mat');

image_delta_x = +0.028;
image_delta_y = -0.012;
nM=0;
%%%%%%%%;
ori_M_k_p_ = ...
transf_p_to_p( ...
 n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,n_w_sum ...
,M_k_p__(:,1+nM) ...
,+image_delta_x ...
,+image_delta_y ...
);
ori_M_k_q_ = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,ori_M_k_p_ ...
);
%%%%%%%%;
ori_svd_VUXM_lwn___ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,1,ori_M_k_q_,pm_n_UX_rank,UX__,X_weight_r_);
%%%%%%%%;
[ori_UX_M_k_q_wn__,ori_UX_M_k_p_wn__] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,1,ori_svd_VUXM_lwn___,zeros(1,1),zeros(1,1));
%%%%%%%%;
pos_M_k_p_ = ...
transf_p_to_p( ...
 n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,n_w_sum ...
,M_k_p__(:,1+nM) ...
,0 ...
,0 ...
);
pos_M_k_q_ = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,pos_M_k_p_ ...
);
%%%%%%%%;
pos_svd_VUXM_lwn___ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,1,pos_M_k_q_,pm_n_UX_rank,UX__,X_weight_r_);
%%%%%%%%;
[pos_UX_M_k_q_wn__,pos_UX_M_k_p_wn__] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,1,pos_svd_VUXM_lwn___,+image_delta_x,+image_delta_y);
%%%%%%%%;
neg_M_k_p_ = ...
transf_p_to_p( ...
 n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,n_w_sum ...
,M_k_p__(:,1+nM) ...
,0 ...
,0 ...
);
neg_M_k_q_ = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,neg_M_k_p_ ...
);
%%%%%%%%;
neg_svd_VUXM_lwn___ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,1,neg_M_k_q_,pm_n_UX_rank,UX__,X_weight_r_);
%%%%%%%%;
[neg_UX_M_k_q_wn__,neg_UX_M_k_p_wn__] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,1,neg_svd_VUXM_lwn___,-image_delta_x,-image_delta_y);
%%%%%%%%;

figure(1);clf;figbig;np=0;
subplot(2,2,1+np);np=np+1;
plot(real(ori_UX_M_k_q_wn__(:)),real(pos_UX_M_k_q_wn__(:)),'kx'); xlabel('ori');ylabel('pos');title('real pos'); axis equal;
subplot(2,2,1+np);np=np+1;
plot(imag(ori_UX_M_k_q_wn__(:)),imag(pos_UX_M_k_q_wn__(:)),'kx'); xlabel('ori');ylabel('pos');title('imag pos'); axis equal;
subplot(2,2,1+np);np=np+1;
plot(real(ori_UX_M_k_q_wn__(:)),real(neg_UX_M_k_q_wn__(:)),'kx'); xlabel('ori');ylabel('neg');title('real neg'); axis equal;
subplot(2,2,1+np);np=np+1;
plot(imag(ori_UX_M_k_q_wn__(:)),imag(neg_UX_M_k_q_wn__(:)),'kx'); xlabel('ori');ylabel('neg');title('imag neg'); axis equal;




%{
sprintf('%s_jpg',dir_trunk);
delta_sigma_use = delta_sigma;
svd_eps = 1e-3;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_UX_rank = 8;
dat_n_iteration = 16;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
%delta_r_max_factor_ = [0,0.125,0.25,0.5,1,sqrt(2),2]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0]; n_delta_r_max_factor = numel(delta_r_max_factor_);
delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0.25]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('ut%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
for ndat_rseed=n_dat_rseed-1;%for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for dat_nUX_rank=0:dat_n_UX_rank-1;%for dat_nUX_rank=0:dat_n_UX_rank-1;
fname_0 = sprintf('tpmutameux_%s_n%.3d',dat_infix,dat_n_M);
fname_2 = sprintf('nUX%.3drng%.3d',dat_nUX_rank,dat_rseed);
MS_fname_mat = sprintf('%s_mat/%s_MS_%s.mat',dir_trunk,fname_0,fname_2);
SM_fname_mat = sprintf('%s_mat/%s_SM_%s.mat',dir_trunk,fname_0,fname_2);
if ( exist(MS_fname_mat,'file'));
disp(sprintf(' %% %s found',MS_fname_mat));
fname_1 = 'MS';
tpmutamf_1( ...
 sprintf('%s_mat',dir_trunk) ...
,sprintf('%s_jpg',dir_trunk) ...
,fname_0 ...
,fname_1 ...
,fname_2 ...
,n_image_sub ...
,M_k_p__ ...
,dat_nUX_rank ...
,dat_n_UX_rank ...
,dat_n_iteration ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,delta_r_max_use ...
,delta_r_max_upb ...
,svd_eps ...
,n_delta_v_requested ...
,n_CTF_rank ...
,CTF_idx_ ...
,CTF_k_p_r__ ...
,CTF_k_p__ ...
,l_max_ ...
,UX__ ...
,X_weight_r_ ...
,a_CTF_UX_Y_quad__ ...
,[] ...
,a_k_Y_quad_ ...
,euler_angle_marina_ ...
,delta_read_x_ ...
,delta_read_y_ ...
,delta_sigma ...
);
end;%if ( exist(MS_fname_mat,'file'));
if ( exist(SM_fname_mat,'file'));
disp(sprintf(' %% %s found',SM_fname_mat));
fname_1 = 'SM';
tpmutamf_1( ...
 sprintf('%s_mat',dir_trunk) ...
,sprintf('%s_jpg',dir_trunk) ...
,fname_0 ...
,fname_1 ...
,fname_2 ...
,n_image_sub ...
,M_k_p__ ...
,dat_nUX_rank ...
,dat_n_UX_rank ...
,dat_n_iteration ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,delta_r_max_use ...
,delta_r_max_upb ...
,svd_eps ...
,n_delta_v_requested ...
,n_CTF_rank ...
,CTF_idx_ ...
,CTF_k_p_r__ ...
,CTF_k_p__ ...
,l_max_ ...
,UX__ ...
,X_weight_r_ ...
,a_CTF_UX_Y_quad__ ...
,[] ...
,a_k_Y_quad_ ...
,euler_angle_marina_ ...
,delta_read_x_ ...
,delta_read_y_ ...
,delta_sigma ...
);
end;%if ( exist(SM_fname_mat,'file'));
end;%for dat_nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
 %}

%{
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_a_k_p_quad_.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_a_k_p_quad_.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_a_k_Y_quad_.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_a_k_Y_quad_.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_c_k_Y_.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_c_k_Y_.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_CTF_k_p__.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_CTF_k_p__.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_M_x_c___.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_M_x_c___.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_S_k_p__.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_S_k_p__.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_M_k_p___.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_M_k_p___.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_S_CTF_k_p___.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_S_CTF_k_p___.mat')
movefile('./dir_principled_marching_mat/test_principled_marching_trpv1_15_X_2d_xcor_d0047.mat','./dir_principled_marching_mat/test_principled_marching_trpv1_16_X_2d_xcor_d0047.mat')
  %}
%{
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_a_k_p_quad_.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_a_k_p_quad_.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_a_k_Y_quad_A.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_a_k_Y_quad_A.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_a_k_Y_quad_.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_a_k_Y_quad_.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_CTF_k_p__.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_CTF_k_p__.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_CTF_k_p_r_xxxx__.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_CTF_k_p_r_xxxx__.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_c_x_u_.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_c_x_u_.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_delta_read_.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_delta_read_.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_euler_angle_marina_.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_euler_angle_marina_.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_M_x_c___center_.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_M_x_c___center_.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_M_x_c___sample.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_M_x_c___sample.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_S_k_p__A.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_S_k_p__A.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_S_k_p__.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_S_k_p__.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_M_k_p___.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_M_k_p___.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_M_k_q___.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_M_k_q___.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_M_k_q___spectrum.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_M_k_q___spectrum.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_S_CTF_k_p___A.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_S_CTF_k_p___A.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_S_CTF_k_p___.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_S_CTF_k_p___.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_S_CTF_k_q___.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_S_CTF_k_q___.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_UX_2d_xcor_d0047_S_CTF_k_q___spectrum.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_UX_2d_xcor_d0047_S_CTF_k_q___spectrum.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_X_2d_xcor_d0047_A.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_X_2d_xcor_d0047_A.jpg')
movefile('dir_principled_marching_jpg/test_principled_marching_trpv1_15_X_2d_xcor_d0047_B.jpg','dir_principled_marching_jpg/test_principled_marching_trpv1_16_X_2d_xcor_d0047_B.jpg')
  %}

%{
  % replicates tmp_X2__ ; 
tmp_X3__ = zeros(tmp_n_delta_v,n_w_max);
%%%%%%%%;
tmp_U_d__ = zeros(tmp_n_delta_v,FTK.n_svd_l);
tmp_delta_r_ = zeros(tmp_n_delta_v,1);
for ndelta_v=0:tmp_n_delta_v-1;
tmp_U_d_ = zeros(FTK.n_svd_l,1);
tmp_delta_r_(1+ndelta_v) = fnorm(tmp_delta__(:,1+ndelta_v));
if (tmp_delta_r>FTK.svd_d_max);
disp(sprintf(' %% Warning, tmp_delta_r %0.6f > svd_d_max %0.6f',tmp_delta_r,FTK.svd_d_max));
end;%if;
tmp_svd_d = (tmp_delta_r_(1+ndelta_v) - FTK.svd_d_m)/FTK.svd_d_c;
for nl=0:FTK.n_svd_l-1;
tmp_U_d__(1+ndelta_v,1+nl) = polyval_r8_reverse_0(FTK.n_svd_d,FTK.svd_U_d_(1+0+nl*FTK.n_svd_d+(0:FTK.n_svd_d-1)),1,tmp_svd_d);
end;%for nl=0:FTK.n_svd_l-1;
end;%for ndelta_v=0:tmp_n_delta_v-1;
%%%%%%%%;
tmp_delta_w_ = atan2(tmp_delta__(1+1,:),tmp_delta__(1+0,:));
tmp_expiw__ = exp(-i*(pi/2 - reshape(tmp_delta_w_,[tmp_n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l]));
%%%%%%%%;
tmp_V_r__ = zeros(FTK.n_svd_l,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
if (2*pi*k_p_r_(1+nk_p_r)>FTK.svd_r_max);
disp(sprintf(' %% Warning, 2*pi*k_p_r_(1+nk_p_r) %0.6f > FTK.svd_r_max %0.6f',2*pi*k_p_r_(1+nk_p_r),FTK.svd_r_max));
end;%if;
tmp_svd_r = (2*pi*k_p_r_(1+nk_p_r) - FTK.svd_r_m)/FTK.svd_r_c;
for nl=0:FTK.n_svd_l-1;
tmp_V_r__(1+nl,1+nk_p_r) = polyval_r8_reverse_0(FTK.n_svd_r,FTK.svd_V_r_(1+0+nl*FTK.n_svd_r+(0:FTK.n_svd_r-1)),1,tmp_svd_r);
end;%for nl=0:FTK.n_svd_l-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_l_max = max(abs(FTK.svd_l_));
tmp_M_k_q__ = zeros(n_w_sum,1+2*tmp_l_max);
ic=0;
for nk_p_r=0:n_k_p_r-1;
for nw=0:n_w_(1+nk_p_r)-1;
tmp_n_w_t = floor(1.0d0*n_w_(1+nk_p_r)/2.0d0);
for tmp_I_l=-tmp_l_max:+tmp_l_max;
nwc = nw;
if (nwc>=tmp_n_w_t); nwc = nwc - n_w_(1+nk_p_r); end;%if;
flag_ict_overflow = 0;
nwd = nwc + tmp_I_l;
if (abs(nwd)<tmp_n_w_t); nwt = periodize(nwd,0,n_w_(1+nk_p_r));  else; nwt = 0; flag_ict_overflow = 1; end;%if;
ict = ic-nw+nwt;
if (flag_ict_overflow==0); tmp_M_k_q__(1+ic,1+tmp_l_max+tmp_I_l) = tmp_M_k_q_(1+ict); end;%if;
end;%for tmp_I_l=-tmp_l_max:+tmp_l_max;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_M_k_q_rwl___ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,weight_2d_k_p_r_/(2*pi),n_w_,1+2*tmp_l_max,tmp_M_k_q__);
tmp_S_k_q_rw__ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,weight_2d_k_p_r_/(2*pi),n_w_,1,tmp_S_k_q_);
%%%%%%%%;
tmp_VSM__ = zeros(FTK.n_svd_l,n_w_max);
for nl=0:FTK.n_svd_l-1;
tmp_VSM__(1+nl,:) = ifft(tmp_V_r__(1+nl,:)*(conj(tmp_S_k_q_rw__(:,:)).*tmp_M_k_q_rwl___(:,:,1+tmp_l_max+FTK.svd_l_(1+nl))))*n_w_max;
end;%for nl=0:FTK.n_svd_l-1;
tmp_USEVSM__ = (tmp_U_d__.*tmp_expiw__)*diag(FTK.svd_s_)*tmp_VSM__;
%%%%%%%%;
tmp_X3__ = tmp_USEVSM__;
disp(sprintf(' %% tmp_X0__ vs tmp_X3__: %0.16f',fnorm(tmp_X0__-tmp_X3__)/fnorm(tmp_X0__)));
 %}

%{
  % replicates tmp_X2__ ; 
tmp_X3__ = zeros(tmp_n_delta_v,n_w_max);
%%%%%%%%;
tmp_U_d__ = zeros(tmp_n_delta_v,FTK.n_svd_l);
tmp_delta_r_ = zeros(tmp_n_delta_v,1);
for ndelta_v=0:tmp_n_delta_v-1;
tmp_U_d_ = zeros(FTK.n_svd_l,1);
tmp_delta_r_(1+ndelta_v) = fnorm(tmp_delta__(:,1+ndelta_v));
if (tmp_delta_r>FTK.svd_d_max);
disp(sprintf(' %% Warning, tmp_delta_r %0.6f > svd_d_max %0.6f',tmp_delta_r,FTK.svd_d_max));
end;%if;
tmp_svd_d = (tmp_delta_r_(1+ndelta_v) - FTK.svd_d_m)/FTK.svd_d_c;
for nl=0:FTK.n_svd_l-1;
tmp_U_d__(1+ndelta_v,1+nl) = polyval_r8_reverse_0(FTK.n_svd_d,FTK.svd_U_d_(1+0+nl*FTK.n_svd_d+(0:FTK.n_svd_d-1)),1,tmp_svd_d);
end;%for nl=0:FTK.n_svd_l-1;
end;%for ndelta_v=0:tmp_n_delta_v-1;
%%%%%%%%;
tmp_delta_w_ = atan2(tmp_delta__(1+1,:),tmp_delta__(1+0,:));
tmp_expiw__ = exp(-i*(pi/2 - reshape(tmp_delta_w_,[tmp_n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l]));
%%%%%%%%;
tmp_V_r__ = zeros(FTK.n_svd_l,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
if (2*pi*k_p_r_(1+nk_p_r)>FTK.svd_r_max);
disp(sprintf(' %% Warning, 2*pi*k_p_r_(1+nk_p_r) %0.6f > FTK.svd_r_max %0.6f',2*pi*k_p_r_(1+nk_p_r),FTK.svd_r_max));
end;%if;
tmp_svd_r = (2*pi*k_p_r_(1+nk_p_r) - FTK.svd_r_m)/FTK.svd_r_c;
for nl=0:FTK.n_svd_l-1;
tmp_V_r__(1+nl,1+nk_p_r) = polyval_r8_reverse_0(FTK.n_svd_r,FTK.svd_V_r_(1+0+nl*FTK.n_svd_r+(0:FTK.n_svd_r-1)),1,tmp_svd_r);
end;%for nl=0:FTK.n_svd_l-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_l_max = max(abs(FTK.svd_l_));
tmp_M_k_q__ = zeros(n_w_sum,1+2*tmp_l_max);
ic=0;
for nk_p_r=0:n_k_p_r-1;
for nw=0:n_w_(1+nk_p_r)-1;
tmp_n_w_t = floor(1.0d0*n_w_(1+nk_p_r)/2.0d0);
for tmp_I_l=-tmp_l_max:+tmp_l_max;
nwc = nw;
if (nwc>=tmp_n_w_t); nwc = nwc - n_w_(1+nk_p_r); end;%if;
flag_ict_overflow = 0;
nwd = nwc + tmp_I_l;
if (abs(nwd)<tmp_n_w_t); nwt = periodize(nwd,0,n_w_(1+nk_p_r));  else; nwt = 0; flag_ict_overflow = 1; end;%if;
ict = ic-nw+nwt;
if (flag_ict_overflow==0); tmp_M_k_q__(1+ic,1+tmp_l_max+tmp_I_l) = tmp_M_k_q_(1+ict); end;%if;
end;%for tmp_I_l=-tmp_l_max:+tmp_l_max;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_M_k_q_rwl___ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,weight_2d_k_p_r_/(2*pi),n_w_,1+2*tmp_l_max,tmp_M_k_q__);
tmp_S_k_q_rw__ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,weight_2d_k_p_r_/(2*pi),n_w_,1,tmp_S_k_q_);
%%%%%%%%;
for ndelta_v=0:tmp_n_delta_v-1;
tmp_Z_q_ = zeros(n_w_sum,1);
tmp_C_q_ = zeros(n_w_(1+n_k_p_r-1),1);
for nw=0:n_w_max-1;
for nk_p_r=0:n_k_p_r-1;
for nl=0:FTK.n_svd_l-1;
tmp_C_q_(1+nw) = tmp_C_q_(1+nw) + ...
  tmp_U_d__(1+ndelta_v,1+nl) * ...
  FTK.svd_s_(1+nl) * ...
  tmp_V_r__(1+nl,1+nk_p_r) * ...
  tmp_expiw__(1+ndelta_v,1+nl) * ...
  tmp_M_k_q_rwl___(1+nk_p_r,1+nw,1+tmp_l_max+FTK.svd_l_(1+nl)) * ...
  conj(tmp_S_k_q_rw__(1+nk_p_r,1+nw));  
end;%for nl=0:FTK.n_svd_l-1;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nw=0:n_w_max-1;
%%%%%%%%;
tmp_X3__(1+ndelta_v,:) = ifft(tmp_C_q_)*n_w_max;
end;%for ndelta_v=0:tmp_n_delta_v-1;
disp(sprintf(' %% tmp_X0__ vs tmp_X3__: %0.16f',fnorm(tmp_X0__-tmp_X3__)/fnorm(tmp_X0__)));
 %}

%{
  % replicates tmp_X2__ ; 
tmp_X3__ = zeros(tmp_n_delta_v,n_w_max);
%%%%%%%%;
tmp_U_d__ = zeros(tmp_n_delta_v,FTK.n_svd_l);
tmp_delta_r_ = zeros(tmp_n_delta_v,1);
for ndelta_v=0:tmp_n_delta_v-1;
tmp_U_d_ = zeros(FTK.n_svd_l,1);
tmp_delta_r_(1+ndelta_v) = fnorm(tmp_delta__(:,1+ndelta_v));
if (tmp_delta_r>FTK.svd_d_max);
disp(sprintf(' %% Warning, tmp_delta_r %0.6f > svd_d_max %0.6f',tmp_delta_r,FTK.svd_d_max));
end;%if;
tmp_svd_d = (tmp_delta_r_(1+ndelta_v) - FTK.svd_d_m)/FTK.svd_d_c;
for nl=0:FTK.n_svd_l-1;
tmp_U_d__(1+ndelta_v,1+nl) = polyval_r8_reverse_0(FTK.n_svd_d,FTK.svd_U_d_(1+0+nl*FTK.n_svd_d+(0:FTK.n_svd_d-1)),1,tmp_svd_d);
end;%for nl=0:FTK.n_svd_l-1;
end;%for ndelta_v=0:tmp_n_delta_v-1;
%%%%%%%%;
tmp_delta_w_ = atan2(tmp_delta__(1+1,:),tmp_delta__(1+0,:));
tmp_expiw__ = exp(-i*(pi/2 - reshape(tmp_delta_w_,[tmp_n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l]));
%%%%%%%%;
tmp_V_r__ = zeros(FTK.n_svd_l,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
if (2*pi*k_p_r_(1+nk_p_r)>FTK.svd_r_max);
disp(sprintf(' %% Warning, 2*pi*k_p_r_(1+nk_p_r) %0.6f > FTK.svd_r_max %0.6f',2*pi*k_p_r_(1+nk_p_r),FTK.svd_r_max));
end;%if;
tmp_svd_r = (2*pi*k_p_r_(1+nk_p_r) - FTK.svd_r_m)/FTK.svd_r_c;
for nl=0:FTK.n_svd_l-1;
tmp_V_r__(1+nl,1+nk_p_r) = polyval_r8_reverse_0(FTK.n_svd_r,FTK.svd_V_r_(1+0+nl*FTK.n_svd_r+(0:FTK.n_svd_r-1)),1,tmp_svd_r);
end;%for nl=0:FTK.n_svd_l-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_l_max = max(abs(FTK.svd_l_));
tmp_M_k_q__ = zeros(n_w_sum,1+2*tmp_l_max);
ic=0;
for nk_p_r=0:n_k_p_r-1;
for nw=0:n_w_(1+nk_p_r)-1;
tmp_n_w_t = floor(1.0d0*n_w_(1+nk_p_r)/2.0d0);
for tmp_I_l=-tmp_l_max:+tmp_l_max;
nwc = nw;
if (nwc>=tmp_n_w_t); nwc = nwc - n_w_(1+nk_p_r); end;%if;
flag_ict_overflow = 0;
nwd = nwc + tmp_I_l;
if (abs(nwd)<tmp_n_w_t); nwt = periodize(nwd,0,n_w_(1+nk_p_r));  else; nwt = 0; flag_ict_overflow = 1; end;%if;
ict = ic-nw+nwt;
if (flag_ict_overflow==0); tmp_M_k_q__(1+ic,1+tmp_l_max+tmp_I_l) = tmp_M_k_q_(1+ict); end;%if;
end;%for tmp_I_l=-tmp_l_max:+tmp_l_max;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
for ndelta_v=0:tmp_n_delta_v-1;
tmp_Z_q_ = zeros(n_w_sum,1);
tmp_C_q_ = zeros(n_w_(1+n_k_p_r-1),1);
ic=0;
for nk_p_r=0:n_k_p_r-1;
tmp_dw = 2*pi/(1.0d0*max(1,n_w_(1+nk_p_r))); tmp_dA = weight_2d_k_p_r_(1+nk_p_r)/(2*pi); tmp_dAn = tmp_dA*tmp_dw;
for nw=0:n_w_(1+nk_p_r)-1;
n_w_t = floor(n_w_(1+nk_p_r)/2);
tmp_Z_q_(1+ic) = 0;
for nl=0:FTK.n_svd_l-1;
tmp_Z_q_(1+ic) = ...
  tmp_Z_q_(1+ic) + ...
  tmp_U_d__(1+ndelta_v,1+nl) * ...
  FTK.svd_s_(1+nl) * ...
  tmp_V_r__(1+nl,1+nk_p_r) * ...
  tmp_expiw__(1+ndelta_v,1+nl) * ...
  tmp_M_k_q__(1+ic,1+tmp_l_max+FTK.svd_l_(1+nl)) ;
end;%for nl=0:FTK.n_svd_l-1;
if (nw>n_w_t); nw_fix = nw - n_w_(1+nk_p_r) + n_w_max; end; 
if (nw<n_w_t); nw_fix = nw; end;
if (nw~=n_w_t); tmp_C_q_(1+nw_fix) = tmp_C_q_(1+nw_fix) + conj(tmp_S_k_q_(1+ic)) * tmp_Z_q_(1+ic) * tmp_dAn ; end;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_X3__(1+ndelta_v,:) = ifft(tmp_C_q_)*n_w_max;
end;%for ndelta_v=0:tmp_n_delta_v-1;
disp(sprintf(' %% tmp_X0__ vs tmp_X3__: %0.16f',fnorm(tmp_X0__-tmp_X3__)/fnorm(tmp_X0__)));
 %}

%{
  % replicates tmp_X2__ ; 
tmp_X3__ = zeros(tmp_n_delta_v,n_w_max);
%%%%%%%%;
tmp_U_d__ = zeros(tmp_n_delta_v,FTK.n_svd_l);
tmp_delta_r_ = zeros(tmp_n_delta_v,1);
for ndelta_v=0:tmp_n_delta_v-1;
tmp_U_d_ = zeros(FTK.n_svd_l,1);
tmp_delta_r_(1+ndelta_v) = fnorm(tmp_delta__(:,1+ndelta_v));
if (tmp_delta_r>FTK.svd_d_max);
disp(sprintf(' %% Warning, tmp_delta_r %0.6f > svd_d_max %0.6f',tmp_delta_r,FTK.svd_d_max));
end;%if;
tmp_svd_d = (tmp_delta_r_(1+ndelta_v) - FTK.svd_d_m)/FTK.svd_d_c;
for nl=0:FTK.n_svd_l-1;
tmp_U_d__(1+ndelta_v,1+nl) = polyval_r8_reverse_0(FTK.n_svd_d,FTK.svd_U_d_(1+0+nl*FTK.n_svd_d+(0:FTK.n_svd_d-1)),1,tmp_svd_d);
end;%for nl=0:FTK.n_svd_l-1;
end;%for ndelta_v=0:tmp_n_delta_v-1;
%%%%%%%%;
tmp_delta_w_ = atan2(tmp_delta__(1+1,:),tmp_delta__(1+0,:));
tmp_expiw__ = exp(-i*(pi/2 - reshape(tmp_delta_w_,[tmp_n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l]));
%%%%%%%%;
tmp_V_r__ = zeros(FTK.n_svd_l,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
if (2*pi*k_p_r_(1+nk_p_r)>FTK.svd_r_max);
disp(sprintf(' %% Warning, 2*pi*k_p_r_(1+nk_p_r) %0.6f > FTK.svd_r_max %0.6f',2*pi*k_p_r_(1+nk_p_r),FTK.svd_r_max));
end;%if;
tmp_svd_r = (2*pi*k_p_r_(1+nk_p_r) - FTK.svd_r_m)/FTK.svd_r_c;
for nl=0:FTK.n_svd_l-1;
tmp_V_r__(1+nl,1+nk_p_r) = polyval_r8_reverse_0(FTK.n_svd_r,FTK.svd_V_r_(1+0+nl*FTK.n_svd_r+(0:FTK.n_svd_r-1)),1,tmp_svd_r);
end;%for nl=0:FTK.n_svd_l-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_l_max = max(abs(FTK.svd_l_));
tmp_M_k_q__ = zeros(n_w_sum,1+2*tmp_l_max);
ic=0;
for nk_p_r=0:n_k_p_r-1;
for nw=0:n_w_(1+nk_p_r)-1;
tmp_n_w_t = floor(1.0d0*n_w_(1+nk_p_r)/2.0d0);
for tmp_I_l=-tmp_l_max:+tmp_l_max;
nwc = nw;
if (nwc>=tmp_n_w_t); nwc = nwc - n_w_(1+nk_p_r); end;%if;
flag_ict_overflow = 0;
nwd = nwc + tmp_I_l;
if (abs(nwd)<tmp_n_w_t); nwt = periodize(nwd,0,n_w_(1+nk_p_r));  else; nwt = 0; flag_ict_overflow = 1; end;%if;
ict = ic-nw+nwt;
if (flag_ict_overflow==0); tmp_M_k_q__(1+ic,1+tmp_l_max+tmp_I_l) = tmp_M_k_q_(1+ict); end;%if;
end;%for tmp_I_l=-tmp_l_max:+tmp_l_max;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
for ndelta_v=0:tmp_n_delta_v-1;
tmp_Z_q_ = zeros(n_w_sum,1);
tmp_C_w_ = zeros(FTK.n_svd_l,1);
ic=0;
for nk_p_r=0:n_k_p_r-1;
for nw=0:n_w_(1+nk_p_r)-1;
tmp_Z_q_(1+ic) = 0;
for nl=0:FTK.n_svd_l-1;
tmp_Z_q_(1+ic) = tmp_Z_q_(1+ic) + ...
                 tmp_U_d__(1+ndelta_v,1+nl) * ...
                 FTK.svd_s_(1+nl) * ...
                 tmp_V_r__(1+nl,1+nk_p_r) * ...
                 tmp_expiw__(1+ndelta_v,1+nl) * ...
                 tmp_M_k_q__(1+ic,1+tmp_l_max+FTK.svd_l_(1+nl)) ;
end;%for nl=0:FTK.n_svd_l-1;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_T_k_q_ = tmp_Z_q_;
%%%%%%%%;
tmp_C_q_ = zeros(n_w_(1+n_k_p_r-1),1);
ic = 0;
for nk_p_r=0:n_k_p_r-1;
tmp_dw = 2*pi/(1.0d0*max(1,n_w_(1+nk_p_r)));
tmp_dA = weight_2d_k_p_r_(1+nk_p_r)/(2*pi);
% We assume that the fourier basis is orthonormal (not merely orthogonal);
tmp_dAn = tmp_dA*tmp_dw;
for nw=0:n_w_(1+nk_p_r)-1;
if (nw>n_w_(1+nk_p_r)/2);
nw_fix = nw - n_w_(1+nk_p_r) + n_w_max;
tmp_C_q = conj(tmp_S_k_q_(1+ic))*tmp_T_k_q_(1+ic);
nw_C = nw_fix;
tmp_C_q_(1+nw_C) = tmp_C_q_(1+nw_C) + tmp_C_q*tmp_dAn;
end;%if (nw>n_w_(1+nk_p_r)/2);
if (nw<n_w_(1+nk_p_r)/2);
nw_fix = nw;
tmp_C_q = conj(tmp_S_k_q_(1+ic))*tmp_T_k_q_(1+ic);
nw_C = nw_fix;
tmp_C_q_(1+nw_C) = tmp_C_q_(1+nw_C) + tmp_C_q*tmp_dAn;
end%if;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_X3__(1+ndelta_v,:) = ifft(tmp_C_q_)*n_w_max;
end;%for ndelta_v=0:tmp_n_delta_v-1;
disp(sprintf(' %% tmp_X0__ vs tmp_X3__: %0.16f',fnorm(tmp_X0__-tmp_X3__)/fnorm(tmp_X0__)));
 %}

%{
  % replicates tmp_X2__ ; 
tmp_X3__ = zeros(tmp_n_delta_v,n_w_max);
for ndelta_v=0:tmp_n_delta_v-1;
%%%%%%%%;
tmp_Z_q_ = zeros(n_w_sum,1);
tmp_U_d_ = zeros(FTK.n_svd_l,1);
tmp_delta_r = fnorm(tmp_delta__(:,1+ndelta_v));
if (tmp_delta_r>FTK.svd_d_max);
disp(sprintf(' %% Warning, tmp_delta_r %0.6f > svd_d_max %0.6f',tmp_delta_r,FTK.svd_d_max));
end;%if;
tmp_delta_w = atan2(tmp_delta__(1+1,1+ndelta_v),tmp_delta__(1+0,1+ndelta_v));
tmp_svd_d = (tmp_delta_r - FTK.svd_d_m)/FTK.svd_d_c;
for nl=0:FTK.n_svd_l-1;
tmp_U_d_(1+nl) = polyval_r8_reverse_0(FTK.n_svd_d,FTK.svd_U_d_(1+0+nl*FTK.n_svd_d+(0:FTK.n_svd_d-1)),1,tmp_svd_d);
end;%for nl=0:FTK.n_svd_l-1;
tmp_V_r_ = zeros(FTK.n_svd_l*n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
if (2*pi*k_p_r_(1+nk_p_r)>FTK.svd_r_max);
disp(sprintf(' %% Warning, 2*pi*k_p_r_(1+nk_p_r) %0.6f > FTK.svd_r_max %0.6f',2*pi*k_p_r_(1+nk_p_r),FTK.svd_r_max));
end;%if;
tmp_svd_r = (2*pi*k_p_r_(1+nk_p_r) - FTK.svd_r_m)/FTK.svd_r_c;
for nl=0:FTK.n_svd_l-1;
tmp_V_r_(1+nl+nk_p_r*FTK.n_svd_l) = polyval_r8_reverse_0(FTK.n_svd_r,FTK.svd_V_r_(1+0+nl*FTK.n_svd_r+(0:FTK.n_svd_r-1)),1,tmp_svd_r);
end;%for nl=0:FTK.n_svd_l-1;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_C_w_ = zeros(FTK.n_svd_l,1);
ic=0;
for nk_p_r=0:n_k_p_r-1;
for nl=0:FTK.n_svd_l-1;
tmp_theta = FTK.svd_l_(1+nl)*(pi/2 - tmp_delta_w);
tmp_C_w = +cos(tmp_theta) - i*sin(tmp_theta);
tmp_D_V_r = tmp_V_r_(1+nl+nk_p_r*FTK.n_svd_l);
tmp_D_U_d = tmp_U_d_(1+nl);
tmp_D_s = FTK.svd_s_(1+nl);
tmp_C_w_(1+nl) = (tmp_D_U_d * tmp_D_s * tmp_D_V_r) * tmp_C_w;
end;%for nl=0:FTK.n_svd_l-1;
for nw=0:n_w_(1+nk_p_r)-1;
tmp_Z_q_(1+ic) = 0;
tmp_n_w_t = floor(1.0d0*n_w_(1+nk_p_r)/2.0d0);
for nl=0:FTK.n_svd_l-1;
tmp_I_l = FTK.svd_l_(1+nl);
tmp_C_q = tmp_C_w_(1+nl);
nwc = nw;
if (nwc>=tmp_n_w_t);
nwc = nwc - n_w_(1+nk_p_r);
end;%if;
flag_ict_overflow = 0;
nwd = nwc + tmp_I_l;
if (abs(nwd)<tmp_n_w_t);
nwt = periodize(nwd,0,n_w_(1+nk_p_r));
 else;
nwt = 0;
flag_ict_overflow = 1;
end;%if;
ict = ic-nw+nwt;
if (flag_ict_overflow==0);
tmp_Z_q_(1+ic) = tmp_Z_q_(1+ic) + tmp_C_q*tmp_M_k_q_(1+ict);
end;%if;
end;%for nl=0:FTK.n_svd_l-1;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_T_k_q_ = tmp_Z_q_;
%%%%%%%%;
tmp_C_q_ = zeros(n_w_(1+n_k_p_r-1),1);
ic = 0;
for nk_p_r=0:n_k_p_r-1;
tmp_dw = 2*pi/(1.0d0*max(1,n_w_(1+nk_p_r)));
tmp_dA = weight_2d_k_p_r_(1+nk_p_r)/(2*pi);
% We assume that the fourier basis is orthonormal (not merely orthogonal);
tmp_dAn = tmp_dA*tmp_dw;
for nw=0:n_w_(1+nk_p_r)-1;
if (nw>n_w_(1+nk_p_r)/2);
nw_fix = nw - n_w_(1+nk_p_r) + n_w_max;
tmp_C_q = conj(tmp_S_k_q_(1+ic))*tmp_T_k_q_(1+ic);
nw_C = nw_fix;
tmp_C_q_(1+nw_C) = tmp_C_q_(1+nw_C) + tmp_C_q*tmp_dAn;
end;%if (nw>n_w_(1+nk_p_r)/2);
if (nw<n_w_(1+nk_p_r)/2);
nw_fix = nw;
tmp_C_q = conj(tmp_S_k_q_(1+ic))*tmp_T_k_q_(1+ic);
nw_C = nw_fix;
tmp_C_q_(1+nw_C) = tmp_C_q_(1+nw_C) + tmp_C_q*tmp_dAn;
end%if;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_X3__(1+ndelta_v,:) = ifft(tmp_C_q_)*n_w_max;
end;%for ndelta_v=0:tmp_n_delta_v-1;
 %}

%{
tmp_U__ = 2*pi * reshape(FTK.svd_polyval_U_d_,[FTK.n_delta_v,FTK.n_svd_l]) .* exp(-i*atan2(FTK.delta_y_(:),FTK.delta_x_(:))*transpose(FTK.svd_l_(:))) * diag((-i).^FTK.svd_l_);
tmp_S__ = diag(FTK.svd_s_);
tmp_V__ = reshape(FTK.svd_polyval_V_r_,[n_k_p_r,FTK.n_svd_l]);
tmp_s__ = zeros(n_w_max,n_k_p_r); %<-- new zero-mode is at n_w_max/2 instead of 0;
tmp_m__ = zeros(n_w_max,n_k_p_r);
tmp_n_w_max_2 = round(n_w_max/2);
for nk_p_r=0:n_k_p_r-1;
tmp_n_w_2 = round(n_w_(1+nk_p_r)/2)-1; 
tmp_ij_sub_ = 0:tmp_n_w_2-1;
tmp_ij_set_ = tmp_ij_sub_ + tmp_n_w_max_2;
tmp_ij_ = n_w_csum_(1+nk_p_r) + tmp_ij_sub_;
tmp_s__(1+tmp_ij_set_,1+nk_p_r) = tmp_S_k_q_(1+tmp_ij_);
tmp_m__(1+tmp_ij_set_,1+nk_p_r) = tmp_M_k_q_(1+tmp_ij_);
tmp_ij_sub_ = n_w_(1+nk_p_r)-tmp_n_w_2:n_w_(1+nk_p_r)-1;
tmp_ij_set_ = tmp_ij_sub_ + tmp_n_w_max_2 - n_w_(1+nk_p_r);
tmp_ij_ = n_w_csum_(1+nk_p_r) + tmp_ij_sub_;
tmp_s__(1+tmp_ij_set_,1+nk_p_r) = tmp_S_k_q_(1+tmp_ij_);
tmp_m__(1+tmp_ij_set_,1+nk_p_r) = tmp_M_k_q_(1+tmp_ij_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_svwm__ = zeros(FTK.n_svd_l,n_w_max);
for nsvd_l=0:FTK.n_svd_l-1;
tmp_l = FTK.svd_l_(1+nsvd_l);
if (tmp_l==0);
tmp_ij_ = 0:n_w_max-1; 
tmp_ = sum(conj(tmp_s__).*repmat(transpose(weight_2d_k_p_r_(:).*tmp_V__(:,1+nsvd_l)),[n_w_max,1]).*tmp_m__(1+tmp_ij_,:),2);
tmp_ij_ = [[tmp_n_w_max_2:n_w_max-1],[0:tmp_n_w_max_2-1]];
tmp_svwm__(1+nsvd_l,:) = ifft(tmp_(1+tmp_ij_))*n_w_max;
end;%if (tmp_l==0);
if (tmp_l> 0); 
tmp_ij_ = 0:n_w_max-1-tmp_l; 
tmp_ = sum(conj(tmp_s__).*repmat(transpose(weight_2d_k_p_r_(:).*tmp_V__(:,1+nsvd_l)),[n_w_max,1]).*[zeros(tmp_l,n_k_p_r);tmp_m__(1+tmp_ij_,:)],2);
tmp_ij_ = [[tmp_n_w_max_2:n_w_max-1],[0:tmp_n_w_max_2-1]];
tmp_svwm__(1+nsvd_l,:) = ifft(tmp_(1+tmp_ij_))*n_w_max;
end;%if (tmp_l> 0); 
if (tmp_l< 0); 
tmp_ij_ = -tmp_l:n_w_max-1;
tmp_ = sum(conj(tmp_s__).*repmat(transpose(weight_2d_k_p_r_(:).*tmp_V__(:,1+nsvd_l)),[n_w_max,1]).*[tmp_m__(1+tmp_ij_,:);zeros(-tmp_l,n_k_p_r)],2);
tmp_ij_ = [[tmp_n_w_max_2:n_w_max-1],[0:tmp_n_w_max_2-1]];
tmp_svwm__(1+nsvd_l,:) = ifft(tmp_(1+tmp_ij_))*n_w_max;
end;%if (tmp_l< 0); 
end;%for nsvd_l=0:FTK.n_svd_l-1;
%%%%%%%%;
tmp_X3__ = tmp_U__*tmp_S__*tmp_svwm__;
%%%%%%%%;
flag_plot=0;
if flag_plot;
clf;
colormap(colormap_beach());
subplot(2,2,1+0); imagesc(real(tmp_X0__)); axisnotick; title('X0');
subplot(2,2,1+1); imagesc(real(tmp_X1__)); axisnotick; title('X1');
subplot(2,2,1+2); imagesc(real(tmp_X2__)); axisnotick; title('X2');
subplot(2,2,1+3); imagesc(real(tmp_X3__)); axisnotick; title('X3');
figbig;
end;%if flag_plot;
 %}


%{
function M_q_ = transf_svd_q_to_q_FTK_5(FTK,n_r,grid_p_,n_w_,n_A,S_q_,delta_x,delta_y);

warning_flag = 1;

svd_r_max = FTK.svd_r_max;
n_svd_r = FTK.n_svd_r;
svd_r_ = FTK.svd_r_;
svd_d_max = FTK.svd_d_max;
n_svd_d = FTK.n_svd_d;
svd_d_ = FTK.svd_d_;
n_svd_l = FTK.n_svd_l;
svd_l_ = FTK.svd_l_;
svd_U_d_ = FTK.svd_U_d_;
svd_s_ = FTK.svd_s_;
svd_V_r_ = FTK.svd_V_r_;

Z_q_ = zeros(n_A,1);
U_d_ = zeros(n_svd_l,1);
svd_d_m = svd_d_max / 2.0;
svd_d_c = svd_d_m;
delta = sqrt(delta_x^2 + delta_y^2);
omega = atan2(delta_y,delta_x);
if (delta>svd_d_max & warning_flag);
disp(sprintf(' %% Warning, delta %0.6f > svd_d_max %0.6f',delta,svd_d_max));
end;%if;
svd_d = (delta - svd_d_m)/svd_d_c;
for nl=0:n_svd_l-1;
U_d_(1+nl) = polyval_r8_reverse_0(n_svd_d,svd_U_d_(1+0+nl*n_svd_d+(0:n_svd_d-1)),1,svd_d);
end;%for nl=0:n_svd_l-1;
V_r_ = zeros(n_svd_l*n_r,1);
svd_r_m = svd_r_max / 2.0;
svd_r_c = svd_r_m;
for nr=0:n_r-1;
if (grid_p_(1+nr)>svd_r_max & warning_flag);
disp(sprintf(' %% Warning, grid_p_(1+nr) %0.6f > svd_r_max %0.6f',grid_p_(1+nr),svd_r_max));
end;%if;
svd_r = (grid_p_(1+nr) - svd_r_m)/svd_r_c;
for nl=0:n_svd_l-1;
V_r_(1+nl+nr*n_svd_l) = polyval_r8_reverse_0(n_svd_r,svd_V_r_(1+0+nl*n_svd_r+(0:n_svd_r-1)),1,svd_r);
end;%for nl=0:n_svd_l-1;
end;%for nr=0:n_r-1;
C_w_ = zeros(n_svd_l,1);
ic=0;
for nr=0:n_r-1;
for nl=0:n_svd_l-1;
theta = svd_l_(1+nl)*(pi/2 - omega);
C_w = +cos(theta) - i*sin(theta);
D_V_r = V_r_(1+nl+nr*n_svd_l);
D_U_d = U_d_(1+nl);
D_s = svd_s_(1+nl);
C_w_(1+nl) = (D_U_d * D_s * D_V_r) * C_w;
end;%for nl=0:n_svd_l-1;
for nw=0:n_w_(1+nr)-1;
Z_q_(1+ic) = 0;
n_w_t = floor(1.0d0*n_w_(1+nr)/2.0d0);
for nl=0:n_svd_l-1;
I_l = svd_l_(1+nl);
C_q = C_w_(1+nl);
nwc = nw;
if (nwc>=n_w_t);
nwc = nwc - n_w_(1+nr);
end;%if;
flag_ict_overflow = 0;
nwd = nwc + I_l;
if (abs(nwd)<n_w_t);
nwt = periodize(nwd,0,n_w_(1+nr));
 else;
nwt = 0;
flag_ict_overflow = 1;
end;%if;
ict = ic-nw+nwt;
if (flag_ict_overflow==0);
Z_q_(1+ic) = Z_q_(1+ic) + C_q*S_q_(1+ict);
end;%if;
end;%for nl=0:n_svd_l-1;
ic = ic + 1;
end;%for nw=0:n_w_(1+nr)-1;
end;%for nr=0:n_r-1;
M_q_ = Z_q_;
 %}

%{
function C_q_ = innerproduct_q_k_stretch_quad_0(n_r,grid_p_,weight_p_,n_w_,n_A,T_q_,M_q_) ;
% Assumes that M_q_ is the same size and dimensions as T_q_. ;
% Assumes quasi-uniform polar-grid. ;
% Assumes that C_q_ is large enough to hold all n_w_(1+n_r-1) modes ;
% (assuming of course that n_w_(1+n_r-1) is the largest value within n_w_). ;
% Stores C_q_ in fft ordering (mode 0 , mode 1 , ... , mode -1) ;
% Upgraded to ignore frequencies of magnitude n_w_(1+nr)/2 or larger. ;
% Upgraded to include radial weight. ;
verbose=0;
C_q_ = zeros(n_w_(1+n_r-1),1);
if (verbose>0); disp(sprintf(' %% [entering innerproduct_q_k_stretch_quad_0] n_r %d',n_r)); end%if;
n_w_max = n_w_(1+n_r-1);
if (verbose>0); disp(sprintf(' %% n_w_max %d',n_w_max)); end%if;
for nw=0:n_w_max-1; C_q_(1+nw) = 0.0; end;%for nw=0:n_w_max-1; C_q = 0.0;
ic = 0;
for nr=0:n_r-1;
dw = 2*pi/(1.0d0*max(1,n_w_(1+nr)));
dA = weight_p_(1+nr);
% We assume that the fourier basis is orthonormal (not merely orthogonal);
dAn = dA*dw;
for nw=0:n_w_(1+nr)-1;
if (nw>n_w_(1+nr)/2);
nw_fix = nw - n_w_(1+nr) + n_w_max;
if (verbose>2 & nr<5); disp(sprintf(' %% nw %d nw_fix %d (full loop)',nw,nw_fix)); end%if;
C_q = conj(T_q_(1+ic))*M_q_(1+ic);
nw_C = nw_fix;
C_q_(1+nw_C) = C_q_(1+nw_C) + C_q*dAn;
 elseif (nw==n_w_(1+nr)/2);
nw_fix = nw;
if (verbose>2 & nr<5); disp(sprintf(' %% nw %d nw_fix %d (first orig)',nw,nw_fix)); end%if;
C_q = 0.0d0*0.5d0*conj(T_q_(1+ic))*M_q_(1+ic);
nw_C = nw_fix;
C_q_(1+nw_C) = C_q_(1+nw_C) + C_q*dAn;
nw_fix = nw - n_w_(1+nr) + n_w_max;
if (verbose>2 & nr<5); disp(sprintf(' %% nw %d nw_fix %d (then loop)',nw,nw_fix)); end%if;
C_q = 0.0d0*0.5d0*conj(T_q_(1+ic))*M_q_(1+ic);
nw_C = nw_fix;
C_q_(1+nw_C) = C_q_(1+nw_C) + C_q*dAn;
 else;
nw_fix = nw;
if (verbose>2 & nr<5); disp(sprintf(' %% nw %d nw_fix %d (full orig)',nw,nw_fix)); end;%if;
C_q = conj(T_q_(1+ic))*M_q_(1+ic);
nw_C = nw_fix;
C_q_(1+nw_C) = C_q_(1+nw_C) + C_q*dAn;
end%if;
ic = ic + 1;
end;%for nw=0:n_w_(1+nr)-1;
end;%for nr=0:n_r-1;
if (verbose>0); disp(sprintf(' %% [finished innerproduct_q_k_stretch_quad_0]')); end;%if;
 %}

%{


 %tmp_ms_ = load('dir_principled_marching_mat/tpmhtameux_2d_xcor_d0047_le3v128_n1024_MS_nUX002rng000.mat');
 %tmp_sm_ = load('dir_principled_marching_mat/tpmhtameux_2d_xcor_d0047_le3v128_n1024_SM_nUX002rng000.mat');
 tmp_ms_ = load('dir_principled_marching_mat/tpmhtamsux_2d_xcor_d0047_le7v128_n1024s0125_MS_nUX003rng000.mat');
 tmp_sm_ = load('dir_principled_marching_mat/tpmhtamsux_2d_xcor_d0047_le7v128_n1024s0125_SM_nUX003rng000.mat');
  % plots histograms of image parameters. ;
   figure(1);clf;
  subplot(2,3,1); hist(tmp_ms_.euler_polar_a_MS__(:,end),linspace(0,2*pi,64)); title('ms polar');
  subplot(2,3,2); hist(tmp_ms_.euler_azimu_b_MS__(:,end),linspace(0,2*pi,64)); title('ms azimu');
  subplot(2,3,3); hist(tmp_ms_.euler_gamma_z_MS__(:,end),linspace(0,2*pi,64)); title('ms gamma');
  subplot(2,3,4); hist(tmp_sm_.euler_polar_a_SM__(:,end),linspace(0,2*pi,64)); title('sm polar');
  subplot(2,3,5); hist(tmp_sm_.euler_azimu_b_SM__(:,end),linspace(0,2*pi,64)); title('sm azimu');
  subplot(2,3,6); hist(tmp_sm_.euler_gamma_z_SM__(:,end),linspace(0,2*pi,64)); title('sm gamma');
  figbig;
  % plots histograms of image parameters. ;
   figure(2);clf;
  subplot(2,2,1); hist(tmp_ms_.image_delta_x_MS__(:,end),linspace(-0.11,+0.11,16)); title('ms polar');
  subplot(2,2,2); hist(tmp_ms_.image_delta_y_MS__(:,end),linspace(-0.11,+0.11,16)); title('ms azimu');
  subplot(2,2,3); hist(tmp_sm_.image_delta_x_SM__(:,end),linspace(-0.11,+0.11,16)); title('sm polar');
  subplot(2,2,4); hist(tmp_sm_.image_delta_y_SM__(:,end),linspace(-0.11,+0.11,16)); title('sm azimu');
  figbig;
  % plots scatterplots of image parameters. ;
   figure(3);clf;
   ms_ij = 32;
   sm_ij = 32;
   subplot(2,3,1); plot(tmp_ms_.euler_polar_a_MS__(:,ms_ij),tmp_sm_.euler_polar_a_SM__(:,sm_ij),'x'); axis equal;
   subplot(2,3,2); plot(tmp_ms_.euler_azimu_b_MS__(:,ms_ij),tmp_sm_.euler_azimu_b_SM__(:,sm_ij),'x'); axis equal;
   subplot(2,3,3); plot(tmp_ms_.euler_gamma_z_MS__(:,ms_ij),tmp_sm_.euler_gamma_z_SM__(:,sm_ij),'x'); axis equal;
   subplot(2,3,4); plot(tmp_ms_.image_delta_x_MS__(:,ms_ij),tmp_sm_.image_delta_x_SM__(:,sm_ij),'x'); axis equal;
   subplot(2,3,5); plot(tmp_ms_.image_delta_y_MS__(:,ms_ij),tmp_sm_.image_delta_y_SM__(:,sm_ij),'x'); axis equal;
  figbig;


  %}
