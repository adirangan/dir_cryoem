%%%%%%%%;
% We construct the riesz integration-weights on the sphere. ;
% These are associated with the riesz-potential 1/k^2.5, ;
% or a weighting-function (for the squared-L2-norm) of 1/k. ;
%%%%%%%%;
weight_3d_riesz_k_p_r_ = weight_3d_k_p_r_;
weight_3d_riesz_k_all_ = weight_3d_k_all_;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
weight_3d_k_p_r = weight_3d_k_p_r_(1+nk_p_r);
weight_2d_k_p_r = weight_2d_k_p_r_(1+nk_p_r);
weight_3d_riesz_k_p_r_(1+nk_p_r) = weight_3d_k_p_r_(1+nk_p_r) * weight_2d_k_p_r / max(1e-16,weight_3d_k_p_r);
tmp_index_ = n_k_all_csum_(1+nk_p_r):n_k_all_csum_(1+nk_p_r+1)-1;
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(weight_3d_k_all_(1+tmp_index_))/(weight_3d_k_p_r * 4*pi): %0.16f',k_p_r,sum(weight_3d_k_all_(1+tmp_index_))/(4*pi*weight_3d_k_p_r))); end;
weight_3d_riesz_k_all_(1+tmp_index_) = weight_3d_k_all_(1+tmp_index_) * weight_2d_k_p_r / max(1e-16,weight_3d_k_p_r);
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(weight_3d_riesz_k_all_(1+tmp_index_))/(weight_2d_k_p_r * 4*pi): %0.16f',k_p_r,sum(weight_3d_riesz_k_all_(1+tmp_index_))/(4*pi*weight_2d_k_p_r))); end;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;

%%%%%%%%;
% Calibrate scaling factor. ;
%%%%%%%%;
term_deltafunc = sqrt(2*pi);
term_2 = (pi*k_p_r_max^2)/(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% sum(weight_2d_wk_) vs (pi*k_p_r_max^2)/(4*pi^2): %0.16f',fnorm(sum(weight_2d_wk_) - term_2))); end;
term_3 = (4/3)*pi*k_p_r_max^3;
if (flag_verbose>0); disp(sprintf(' %% sum(weight_3d_k_all_) vs (4/3)*pi*k_p_r_max^3: %0.16f',fnorm(sum(weight_3d_k_all_) - term_3))); end;
term_3r = (4*pi^2*k_p_r_max^2);
if (flag_verbose>0); disp(sprintf(' %% sum(weight_3d_riesz__all_) vs 4*pi^2*k_p_r_max^2: %0.16f',fnorm(sum(weight_3d_riesz_k_all_) - term_3r))); end;
scaling_volumetric = term_3r / term_2 / term_deltafunc ;
if (flag_verbose>0); disp(sprintf(' %% scaling_volumetric: %+0.6f',scaling_volumetric)); end;
if (flag_verbose>0); disp(sprintf(' %% (4*pi)^2 * sqrt(pi/2): %+0.6f',(4*pi)^2 * sqrt(pi/2))); end;
%%%%%%%%;
