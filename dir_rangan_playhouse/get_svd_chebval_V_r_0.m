function svd_chebval_V_r_ = get_svd_chebval_V_r_0(svd_r_max,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_V_r_chebcoef_,n_r,grid_p_);
flag_verbose=0; flag_warning=1;
if (flag_verbose>0); disp(sprintf(' %% [entering get_svd_chebval_V_r_0]')); end;
svd_r_m = svd_r_max / 2.0d0;
svd_r_c = svd_r_m;

%%%%%%%%;
% loop version. ;
%%%%%%%%;
%{
for nr=0:n_r-1;
if (grid_p_(1+nr)>svd_r_max & flag_warning);
disp(sprintf(' %% Warning, grid_p_(1+nr) %0.2f > svd_r_max %0.2f',grid_p_(1+nr),svd_r_max));
end;%if (grid_p_(1+nr)>svd_r_max & flag_warning);
svd_r = (grid_p_(1+nr) - svd_r_m)/svd_r_c; if (~isfinite(svd_r)); svd_r=0; end;
for nl=0:n_svd_l-1;
svd_chebval_V_r_(1+nl+nr*n_svd_l) = chebval_0(n_svd_r,svd_V_r_chebcoef_(1+0+nl*n_svd_r+(0:n_svd_r-1)),1,svd_r);
end;%for nl=0:n_svd_l-1;
end;%for nr=0:n_r-1;
%}

%%%%%%%%;
% vect version. ;
%%%%%%%%;
tolerance_margin = 1e-6;
if sum(grid_p_>svd_r_max+tolerance_margin & flag_warning);
disp(sprintf(' %% Warning, grid_p_ > svd_r_max %0.2f',svd_r_max));
end;%if sum(grid_p_>svd_r_max+tolerance_margin & flag_warning);
tmp_svd_r_ = (grid_p_ - svd_r_m)/max(1e-12,svd_r_c);
svd_chebval_V_r_lr__ = zeros(n_svd_l,n_r);
for nl=0:n_svd_l-1;
svd_chebval_V_r_lr__(1+nl,:) = chebval_0(n_svd_r,svd_V_r_chebcoef_(1+0+nl*n_svd_r+(0:n_svd_r-1)),n_r,tmp_svd_r_);
end;%for nl=0:n_svd_l-1;
svd_chebval_V_r_ = reshape(svd_chebval_V_r_lr__,[n_svd_l*n_r,1]);

if (flag_verbose>0); disp(sprintf(' %% [finished get_svd_chebval_V_r_0]')); end;
