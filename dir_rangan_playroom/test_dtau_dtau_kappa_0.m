%%%%%%%%;
% test derivatives of kappa. ;
%%%%%%%%;
rng(0);

%%%%%%%%;
% First assume only one spatial dimension. ;
%%%%%%%%;
tmp_c_ = randn(1+16,1);
tmp_chebfun_kernel_ = chebfun(leg2cheb(tmp_c_,'norm'),'coeffs');
n_n = 3;
n_w = 5;
n_M = 7;
qref_wMn___ = 0.1*randn(n_w,n_M,n_n);
tmp_HwM33____ = 0.1*randn(n_w,n_M,3,3);
tmp_HwM33____ = 0.5*(tmp_HwM33____ + permute(tmp_HwM33____,[1,2,4,3]));
tmp_JwM3___ = 0.1*randn(n_w,n_M,3);
tmp_KwM__ = 0.1*randn(n_w,n_M);
%%%%;
tmp_point_wM__ = @(da_M_,db_M_,dc_M_) ...
  + tmp_KwM__(:,:) ...
  + bsxfun(@times,tmp_JwM3___(:,:,1+0),reshape(da_M_,[1,n_M,1])) ...
  + bsxfun(@times,tmp_JwM3___(:,:,1+1),reshape(db_M_,[1,n_M,1])) ...
  + bsxfun(@times,tmp_JwM3___(:,:,1+2),reshape(dc_M_,[1,n_M,1])) ...
  + 0.5*bsxfun(@times,tmp_HwM33____(:,:,1+0,1+0),reshape(da_M_,[1,n_M,1,1]).*reshape(da_M_,[1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HwM33____(:,:,1+0,1+1),reshape(da_M_,[1,n_M,1,1]).*reshape(db_M_,[1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HwM33____(:,:,1+0,1+2),reshape(da_M_,[1,n_M,1,1]).*reshape(dc_M_,[1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HwM33____(:,:,1+1,1+0),reshape(db_M_,[1,n_M,1,1]).*reshape(da_M_,[1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HwM33____(:,:,1+1,1+1),reshape(db_M_,[1,n_M,1,1]).*reshape(db_M_,[1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HwM33____(:,:,1+1,1+2),reshape(db_M_,[1,n_M,1,1]).*reshape(dc_M_,[1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HwM33____(:,:,1+2,1+0),reshape(dc_M_,[1,n_M,1,1]).*reshape(da_M_,[1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HwM33____(:,:,1+2,1+1),reshape(dc_M_,[1,n_M,1,1]).*reshape(db_M_,[1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HwM33____(:,:,1+2,1+2),reshape(dc_M_,[1,n_M,1,1]).*reshape(dc_M_,[1,n_M,1,1])) ...
  ;
%%%%;
dtau = 1e-4;
da_M_ = randn(n_M,1);
db_M_ = randn(n_M,1);
dc_M_ = randn(n_M,1);
dtau_fnorm = fnorm([da_M_,db_M_,dc_M_]);
da_M_ = da_M_ / max(1e-12,dtau_fnorm);
db_M_ = db_M_ / max(1e-12,dtau_fnorm);
dc_M_ = dc_M_ / max(1e-12,dtau_fnorm);
%%%%;
z0_M_ = zeros(n_M,1);
tmp_point_mid_wM__ = tmp_point_wM__(z0_M_ + 0*da_M_,z0_M_ + 0*db_M_,z0_M_ + 0*dc_M_);
tmp_point_pos_wM__ = tmp_point_wM__(z0_M_ + 1*dtau*da_M_,z0_M_ + 1*dtau*db_M_,z0_M_ + 1*dtau*dc_M_);
tmp_point_neg_wM__ = tmp_point_wM__(z0_M_ - 1*dtau*da_M_,z0_M_ - 1*dtau*db_M_,z0_M_ - 1*dtau*dc_M_);
tmp_dtau_point_dif_wM__ = (tmp_point_pos_wM__ - tmp_point_neg_wM__)/max(1e-12,2*dtau);
tmp_dtau_point_mid_wM__ = ...
 + bsxfun(@times,tmp_JwM3___(:,:,1+0),reshape(da_M_,[1,n_M,1])) ...
 + bsxfun(@times,tmp_JwM3___(:,:,1+1),reshape(db_M_,[1,n_M,1])) ...
 + bsxfun(@times,tmp_JwM3___(:,:,1+2),reshape(dc_M_,[1,n_M,1])) ...
;
disp(sprintf(' %% tmp_dtau_point_dif_wM__ vs tmp_dtau_point_mid_wM__: %0.16f',fnorm(tmp_dtau_point_dif_wM__-tmp_dtau_point_mid_wM__)/fnorm(tmp_dtau_point_dif_wM__)));
tmp_dtau_dtau_point_dif_wM__ = (tmp_point_pos_wM__ - 2*tmp_point_mid_wM__ + tmp_point_neg_wM__)/max(1e-12,dtau*dtau);
tmp_dtau_dtau_point_mid_wM__ = ...
  + 1.0*bsxfun(@times,tmp_HwM33____(:,:,1+0,1+0),reshape(da_M_,[1,n_M,1,1]).*reshape(da_M_,[1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HwM33____(:,:,1+0,1+1),reshape(da_M_,[1,n_M,1,1]).*reshape(db_M_,[1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HwM33____(:,:,1+0,1+2),reshape(da_M_,[1,n_M,1,1]).*reshape(dc_M_,[1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HwM33____(:,:,1+1,1+0),reshape(db_M_,[1,n_M,1,1]).*reshape(da_M_,[1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HwM33____(:,:,1+1,1+1),reshape(db_M_,[1,n_M,1,1]).*reshape(db_M_,[1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HwM33____(:,:,1+1,1+2),reshape(db_M_,[1,n_M,1,1]).*reshape(dc_M_,[1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HwM33____(:,:,1+2,1+0),reshape(dc_M_,[1,n_M,1,1]).*reshape(da_M_,[1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HwM33____(:,:,1+2,1+1),reshape(dc_M_,[1,n_M,1,1]).*reshape(db_M_,[1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HwM33____(:,:,1+2,1+2),reshape(dc_M_,[1,n_M,1,1]).*reshape(dc_M_,[1,n_M,1,1])) ...
;
disp(sprintf(' %% tmp_dtau_dtau_point_dif_wM__ vs tmp_dtau_dtau_point_mid_wM__: %0.16f',fnorm(tmp_dtau_dtau_point_dif_wM__-tmp_dtau_dtau_point_mid_wM__)/fnorm(tmp_dtau_dtau_point_dif_wM__)));
%%%%;
tmp_distsquared_wMn___ = @(da_M_,db_M_,dc_M_) ...
  bsxfun(@minus,qref_wMn___,tmp_point_wM__(da_M_,db_M_,dc_M_)).^2 ...
  ;
%%%%;
tmp_distsquared_mid_wMn___ = tmp_distsquared_wMn___(z0_M_ + 0*da_M_,z0_M_ + 0*db_M_,z0_M_ + 0*dc_M_);
tmp_distsquared_pos_wMn___ = tmp_distsquared_wMn___(z0_M_ + 1*dtau*da_M_,z0_M_ + 1*dtau*db_M_,z0_M_ + 1*dtau*dc_M_);
tmp_distsquared_neg_wMn___ = tmp_distsquared_wMn___(z0_M_ - 1*dtau*da_M_,z0_M_ - 1*dtau*db_M_,z0_M_ - 1*dtau*dc_M_);
tmp_dtau_distsquared_dif_wMn___ = (tmp_distsquared_pos_wMn___ - tmp_distsquared_neg_wMn___)/max(1e-12,2*dtau);
tmp_d_distsquared_wMn___ = -2*bsxfun(@minus,qref_wMn___,tmp_point_mid_wM__);
tmp_dd_distsquared_wMn___ = 2*ones(n_w,n_M,n_n);
tmp_dtau_distsquared_mid_wMn___ = bsxfun(@times,tmp_d_distsquared_wMn___,tmp_dtau_point_mid_wM__);
disp(sprintf(' %% tmp_dtau_distsquared_dif_wMn___ vs tmp_dtau_distsquared_mid_wMn___: %0.16f',fnorm(tmp_dtau_distsquared_dif_wMn___-tmp_dtau_distsquared_mid_wMn___)/fnorm(tmp_dtau_distsquared_dif_wMn___)));
tmp_dtau_dtau_distsquared_dif_wMn___ = (tmp_distsquared_pos_wMn___ - 2*tmp_distsquared_mid_wMn___ + tmp_distsquared_neg_wMn___)/max(1e-12,dtau*dtau);
tmp_dtau_dtau_distsquared_mid_wMn___ = ...
 + bsxfun(@times,tmp_dd_distsquared_wMn___,bsxfun(@times,tmp_dtau_point_mid_wM__,tmp_dtau_point_mid_wM__)) ...
 + bsxfun(@times,tmp_d_distsquared_wMn___,tmp_dtau_dtau_point_mid_wM__) ...
;
disp(sprintf(' %% tmp_dtau_dtau_distsquared_dif_wMn___ vs tmp_dtau_dtau_distsquared_mid_wMn___: %0.16f',fnorm(tmp_dtau_dtau_distsquared_dif_wMn___-tmp_dtau_dtau_distsquared_mid_wMn___)/fnorm(tmp_dtau_dtau_distsquared_dif_wMn___)));
%%%%;
tmp_kappa_wMn___ = @(da_M_,db_M_,dc_M_) ...
  tmp_chebfun_kernel_(1-tmp_distsquared_wMn___(da_M_,db_M_,dc_M_)/2) ...
  ;
%%%%;
tmp_dchebfun_kernel_ = diff(tmp_chebfun_kernel_);
tmp_ddchebfun_kernel_ = diff(tmp_dchebfun_kernel_);
tmp_kappa_mid_wMn___ = tmp_kappa_wMn___(z0_M_ + 0*da_M_,z0_M_ + 0*db_M_,z0_M_ + 0*dc_M_);
tmp_kappa_pos_wMn___ = tmp_kappa_wMn___(z0_M_ + 1*dtau*da_M_,z0_M_ + 1*dtau*db_M_,z0_M_ + 1*dtau*dc_M_);
tmp_kappa_neg_wMn___ = tmp_kappa_wMn___(z0_M_ - 1*dtau*da_M_,z0_M_ - 1*dtau*db_M_,z0_M_ - 1*dtau*dc_M_);
tmp_dtau_kappa_dif_wMn___ = (tmp_kappa_pos_wMn___ - tmp_kappa_neg_wMn___)/max(1e-12,2*dtau);
tmp_dtau_kappa_mid_wMn___ = tmp_dchebfun_kernel_(1-tmp_distsquared_mid_wMn___/2).*(-0.5).*tmp_dtau_distsquared_mid_wMn___;
disp(sprintf(' %% tmp_dtau_kappa_dif_wMn___ vs tmp_dtau_kappa_mid_wMn___: %0.16f',fnorm(tmp_dtau_kappa_dif_wMn___-tmp_dtau_kappa_mid_wMn___)/fnorm(tmp_dtau_kappa_dif_wMn___)));
tmp_dtau_dtau_kappa_dif_wMn___ = (tmp_kappa_pos_wMn___ - 2*tmp_kappa_mid_wMn___ + tmp_kappa_neg_wMn___)/max(1e-12,dtau*dtau);
tmp_dtau_dtau_kappa_mid_wMn___ = ...
 + tmp_dchebfun_kernel_(1-tmp_distsquared_mid_wMn___/2).*(-0.5).*tmp_dtau_dtau_distsquared_mid_wMn___ ...
 + tmp_ddchebfun_kernel_(1-tmp_distsquared_mid_wMn___/2).*(-0.5).*tmp_dtau_distsquared_mid_wMn___.*(-0.5).*tmp_dtau_distsquared_mid_wMn___ ...
;
disp(sprintf(' %% tmp_dtau_dtau_kappa_dif_wMn___ vs tmp_dtau_dtau_kappa_mid_wMn___: %0.16f',fnorm(tmp_dtau_dtau_kappa_dif_wMn___-tmp_dtau_dtau_kappa_mid_wMn___)/fnorm(tmp_dtau_dtau_kappa_dif_wMn___)));

%%%%%%%%;
% Now try with 3 spatial dimensions. ;
%%%%%%%%;
tmp_c_ = randn(1+16,1);
tmp_chebfun_kernel_ = chebfun(leg2cheb(tmp_c_,'norm'),'coeffs');
n_c = 3;
n_n = 11;
n_w = 5;
n_M = 7;
qref_cwMn____ = 0.1*randn(n_c,n_w,n_M,n_n);
tmp_HcwM33_____ = 0.1*randn(n_c,n_w,n_M,3,3);
tmp_HcwM33_____ = 0.5*(tmp_HcwM33_____ + permute(tmp_HcwM33_____,[1,2,3,5,4]));
tmp_JcwM3____ = 0.1*randn(n_c,n_w,n_M,3);
tmp_KcwM___ = 0.1*randn(n_c,n_w,n_M);
%%%%;
tmp_point_cwM___ = @(da_M_,db_M_,dc_M_) ...
  + tmp_KcwM___(:,:,:) ...
  + bsxfun(@times,tmp_JcwM3____(:,:,:,1+0),reshape(da_M_,[1,1,n_M,1])) ...
  + bsxfun(@times,tmp_JcwM3____(:,:,:,1+1),reshape(db_M_,[1,1,n_M,1])) ...
  + bsxfun(@times,tmp_JcwM3____(:,:,:,1+2),reshape(dc_M_,[1,1,n_M,1])) ...
  + 0.5*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+0,1+0),reshape(da_M_,[1,1,n_M,1,1]).*reshape(da_M_,[1,1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+0,1+1),reshape(da_M_,[1,1,n_M,1,1]).*reshape(db_M_,[1,1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+0,1+2),reshape(da_M_,[1,1,n_M,1,1]).*reshape(dc_M_,[1,1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+1,1+0),reshape(db_M_,[1,1,n_M,1,1]).*reshape(da_M_,[1,1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+1,1+1),reshape(db_M_,[1,1,n_M,1,1]).*reshape(db_M_,[1,1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+1,1+2),reshape(db_M_,[1,1,n_M,1,1]).*reshape(dc_M_,[1,1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+2,1+0),reshape(dc_M_,[1,1,n_M,1,1]).*reshape(da_M_,[1,1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+2,1+1),reshape(dc_M_,[1,1,n_M,1,1]).*reshape(db_M_,[1,1,n_M,1,1])) ...
  + 0.5*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+2,1+2),reshape(dc_M_,[1,1,n_M,1,1]).*reshape(dc_M_,[1,1,n_M,1,1])) ...
  ;
%%%%;
dtau = 1e-4;
da_M_ = randn(n_M,1);
db_M_ = randn(n_M,1);
dc_M_ = randn(n_M,1);
dtau_fnorm = fnorm([da_M_,db_M_,dc_M_]);
da_M_ = da_M_ / max(1e-12,dtau_fnorm);
db_M_ = db_M_ / max(1e-12,dtau_fnorm);
dc_M_ = dc_M_ / max(1e-12,dtau_fnorm);
%%%%;
z0_M_ = zeros(n_M,1);
tmp_point_mid_cwM___ = tmp_point_cwM___(z0_M_ + 0*da_M_,z0_M_ + 0*db_M_,z0_M_ + 0*dc_M_);
tmp_point_pos_cwM___ = tmp_point_cwM___(z0_M_ + 1*dtau*da_M_,z0_M_ + 1*dtau*db_M_,z0_M_ + 1*dtau*dc_M_);
tmp_point_neg_cwM___ = tmp_point_cwM___(z0_M_ - 1*dtau*da_M_,z0_M_ - 1*dtau*db_M_,z0_M_ - 1*dtau*dc_M_);
tmp_dtau_point_dif_cwM___ = (tmp_point_pos_cwM___ - tmp_point_neg_cwM___)/max(1e-12,2*dtau);
tmp_dtau_point_mid_cwM___ = ...
 + bsxfun(@times,tmp_JcwM3____(:,:,:,1+0),reshape(da_M_,[1,1,n_M,1])) ...
 + bsxfun(@times,tmp_JcwM3____(:,:,:,1+1),reshape(db_M_,[1,1,n_M,1])) ...
 + bsxfun(@times,tmp_JcwM3____(:,:,:,1+2),reshape(dc_M_,[1,1,n_M,1])) ...
;
disp(sprintf(' %% tmp_dtau_point_dif_cwM___ vs tmp_dtau_point_mid_cwM___: %0.16f',fnorm(tmp_dtau_point_dif_cwM___-tmp_dtau_point_mid_cwM___)/fnorm(tmp_dtau_point_dif_cwM___)));
tmp_dtau_dtau_point_dif_cwM___ = (tmp_point_pos_cwM___ - 2*tmp_point_mid_cwM___ + tmp_point_neg_cwM___)/max(1e-12,dtau*dtau);
tmp_dtau_dtau_point_mid_cwM___ = ...
  + 1.0*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+0,1+0),reshape(da_M_,[1,1,n_M,1,1]).*reshape(da_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+0,1+1),reshape(da_M_,[1,1,n_M,1,1]).*reshape(db_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+0,1+2),reshape(da_M_,[1,1,n_M,1,1]).*reshape(dc_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+1,1+0),reshape(db_M_,[1,1,n_M,1,1]).*reshape(da_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+1,1+1),reshape(db_M_,[1,1,n_M,1,1]).*reshape(db_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+1,1+2),reshape(db_M_,[1,1,n_M,1,1]).*reshape(dc_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+2,1+0),reshape(dc_M_,[1,1,n_M,1,1]).*reshape(da_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+2,1+1),reshape(dc_M_,[1,1,n_M,1,1]).*reshape(db_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,tmp_HcwM33_____(:,:,:,1+2,1+2),reshape(dc_M_,[1,1,n_M,1,1]).*reshape(dc_M_,[1,1,n_M,1,1])) ...
;
disp(sprintf(' %% tmp_dtau_dtau_point_dif_cwM___ vs tmp_dtau_dtau_point_mid_cwM___: %0.16f',fnorm(tmp_dtau_dtau_point_dif_cwM___-tmp_dtau_dtau_point_mid_cwM___)/fnorm(tmp_dtau_dtau_point_dif_cwM___)));
%%%%;
tmp_distsquared_1wMn___ = @(da_M_,db_M_,dc_M_) ...
  reshape(sum(bsxfun(@minus,qref_cwMn____,tmp_point_cwM___(da_M_,db_M_,dc_M_)).^2,1),[n_w,n_M,n_n]) ...
  ;
%%%%;
tmp_distsquared_mid_1wMn___ = tmp_distsquared_1wMn___(z0_M_ + 0*da_M_,z0_M_ + 0*db_M_,z0_M_ + 0*dc_M_);
tmp_distsquared_pos_1wMn___ = tmp_distsquared_1wMn___(z0_M_ + 1*dtau*da_M_,z0_M_ + 1*dtau*db_M_,z0_M_ + 1*dtau*dc_M_);
tmp_distsquared_neg_1wMn___ = tmp_distsquared_1wMn___(z0_M_ - 1*dtau*da_M_,z0_M_ - 1*dtau*db_M_,z0_M_ - 1*dtau*dc_M_);
tmp_dtau_distsquared_dif_1wMn___ = (tmp_distsquared_pos_1wMn___ - tmp_distsquared_neg_1wMn___)/max(1e-12,2*dtau);
tmp_d_distsquared_cwMn____ = -2*bsxfun(@minus,qref_cwMn____,tmp_point_mid_cwM___);
tmp_dtau_distsquared_mid_1wMn___ = reshape(sum(bsxfun(@times,tmp_d_distsquared_cwMn____,tmp_dtau_point_mid_cwM___),1),[n_w,n_M,n_n]);
disp(sprintf(' %% tmp_dtau_distsquared_dif_1wMn___ vs tmp_dtau_distsquared_mid_1wMn___: %0.16f',fnorm(tmp_dtau_distsquared_dif_1wMn___-tmp_dtau_distsquared_mid_1wMn___)/fnorm(tmp_dtau_distsquared_dif_1wMn___)));
tmp_dtau_dtau_distsquared_dif_1wMn___ = (tmp_distsquared_pos_1wMn___ - 2*tmp_distsquared_mid_1wMn___ + tmp_distsquared_neg_1wMn___)/max(1e-12,dtau*dtau);
tmp_dtau_dtau_distsquared_mid_1wMn___ = ...
 + reshape(sum(bsxfun(@times,2*ones(1,1,1,n_n),bsxfun(@times,tmp_dtau_point_mid_cwM___,tmp_dtau_point_mid_cwM___)),1),[n_w,n_M,n_n]) ...
 + reshape(sum(bsxfun(@times,tmp_d_distsquared_cwMn____,tmp_dtau_dtau_point_mid_cwM___),1),[n_w,n_M,n_n]) ...
;
disp(sprintf(' %% tmp_dtau_dtau_distsquared_dif_1wMn___ vs tmp_dtau_dtau_distsquared_mid_1wMn___: %0.16f',fnorm(tmp_dtau_dtau_distsquared_dif_1wMn___-tmp_dtau_dtau_distsquared_mid_1wMn___)/fnorm(tmp_dtau_dtau_distsquared_dif_1wMn___)));
%%%%;
tmp_kappa_1wMn___ = @(da_M_,db_M_,dc_M_) ...
  tmp_chebfun_kernel_(1-tmp_distsquared_1wMn___(da_M_,db_M_,dc_M_)/2) ...
  ;
%%%%;
tmp_dchebfun_kernel_ = diff(tmp_chebfun_kernel_);
tmp_ddchebfun_kernel_ = diff(tmp_dchebfun_kernel_);
tmp_kappa_mid_1wMn___ = tmp_kappa_1wMn___(z0_M_ + 0*da_M_,z0_M_ + 0*db_M_,z0_M_ + 0*dc_M_);
tmp_kappa_pos_1wMn___ = tmp_kappa_1wMn___(z0_M_ + 1*dtau*da_M_,z0_M_ + 1*dtau*db_M_,z0_M_ + 1*dtau*dc_M_);
tmp_kappa_neg_1wMn___ = tmp_kappa_1wMn___(z0_M_ - 1*dtau*da_M_,z0_M_ - 1*dtau*db_M_,z0_M_ - 1*dtau*dc_M_);
tmp_dtau_kappa_dif_1wMn___ = (tmp_kappa_pos_1wMn___ - tmp_kappa_neg_1wMn___)/max(1e-12,2*dtau);
tmp_dtau_kappa_mid_1wMn___ = tmp_dchebfun_kernel_(1-tmp_distsquared_mid_1wMn___/2).*(-0.5).*tmp_dtau_distsquared_mid_1wMn___;
disp(sprintf(' %% tmp_dtau_kappa_dif_1wMn___ vs tmp_dtau_kappa_mid_1wMn___: %0.16f',fnorm(tmp_dtau_kappa_dif_1wMn___-tmp_dtau_kappa_mid_1wMn___)/fnorm(tmp_dtau_kappa_dif_1wMn___)));
tmp_dtau_dtau_kappa_dif_1wMn___ = (tmp_kappa_pos_1wMn___ - 2*tmp_kappa_mid_1wMn___ + tmp_kappa_neg_1wMn___)/max(1e-12,dtau*dtau);
tmp_dtau_dtau_kappa_mid_1wMn___ = ...
 + tmp_dchebfun_kernel_(1-tmp_distsquared_mid_1wMn___/2).*(-0.5).*tmp_dtau_dtau_distsquared_mid_1wMn___ ...
 + tmp_ddchebfun_kernel_(1-tmp_distsquared_mid_1wMn___/2).*(-0.5).*tmp_dtau_distsquared_mid_1wMn___.*(-0.5).*tmp_dtau_distsquared_mid_1wMn___ ...
;
disp(sprintf(' %% tmp_dtau_dtau_kappa_dif_1wMn___ vs tmp_dtau_dtau_kappa_mid_1wMn___: %0.16f',fnorm(tmp_dtau_dtau_kappa_dif_1wMn___-tmp_dtau_dtau_kappa_mid_1wMn___)/fnorm(tmp_dtau_dtau_kappa_dif_1wMn___)));