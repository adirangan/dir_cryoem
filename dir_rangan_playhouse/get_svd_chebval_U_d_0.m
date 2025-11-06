function svd_chebval_U_d_ = get_svd_chebval_U_d_0(svd_d_max,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_chebcoef_,n_delta_v,delta_x_,delta_y_);
flag_verbose=0; flag_warning=1;
if (flag_verbose>0); disp(sprintf(' %% [entering get_svd_chebval_U_d_0]')); end;
svd_d_m = svd_d_max / 2.0d0;
svd_d_c = svd_d_m;

%%%%%%%%;
% loop version. ;
%%%%%%%%;
%{
svd_chebval_U_d_ = zeros(n_svd_l*n_delta_v,1);
for ndv=0:n_delta_v-1;
delta_x = delta_x_(1+ndv);
delta_y = delta_y_(1+ndv);
delta_r = sqrt(delta_x^2 + delta_y^2);
if (delta_r>svd_d_max & flag_warning);
disp(sprintf(' %% Warning, delta_r %0.2f > svd_d_max %0.2f',delta_r,svd_d_max));
end;%if (delta_r>svd_d_max & flag_warning);
delta_w = atan2(delta_y,delta_x);
svd_d = (delta_r - svd_d_m)/svd_d_c; if (~isfinite(svd_d)); svd_d=0; end;
for nl=0:n_svd_l-1;
svd_chebval_U_d_(1+nl+ndv*n_svd_l) = chebval_0(n_svd_d,svd_U_d_chebcoef_(1+0+nl*n_svd_d+(0:n_svd_d-1)),1,svd_d);
end;%for nl=0:n_svd_l-1;
end;%for ndv=0:n_delta_v-1;
%}

%%%%%%%%;
% vect version. ;
%%%%%%%%;
tolerance_margin = 1e-6;
delta_r_ = sqrt(delta_x_(:).^2 + delta_y_(:).^2);
%delta_w_ = atan2(delta_y_(:),delta_x_(:));
if sum(delta_r_>svd_d_max+tolerance_margin)>0 & flag_warning;
disp(sprintf(' %% Warning, delta_r %0.2f > svd_d_max %0.2f',delta_r,svd_d_max));
end;%if sum(delta_r_>svd_d_max+tolerance_margin)>0 & flag_warning;
tmp_svd_d_ = (delta_r_ - svd_d_m)/max(1e-12,svd_d_c);
svd_chebval_U_d_ld__ = zeros(n_svd_l,n_delta_v);
for nl=0:n_svd_l-1;
svd_chebval_U_d_ld__(1+nl,:) = chebval_0(n_svd_d,svd_U_d_chebcoef_(1+0+nl*n_svd_d+(0:n_svd_d-1)),n_delta_v,tmp_svd_d_);
end;%for nl=0:n_svd_l-1;
svd_chebval_U_d_ = reshape(svd_chebval_U_d_ld__,[n_svd_l*n_delta_v,1]);

if (flag_verbose>0); disp(sprintf(' %% [finished get_svd_chebval_U_d_0]')); end;
