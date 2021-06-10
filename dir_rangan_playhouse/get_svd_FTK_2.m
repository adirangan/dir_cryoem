function FTK = get_svd_FTK_2(eps_target,grid_p_,n_r,n_delta_v,delta_x_,delta_y_,flag_warning,dir_svd);
verbose=0;
if (verbose>0); disp(sprintf(' %% [entering get_svd_FTK_2]')); end;
pi = 4.0d0*atan(1.0d0);
R_max = 2.0d0*pi*grid_p_(1+n_r-1);
K_max = grid_p_(1+n_r-1);
delta_max = 0.0d0;
for ndv=0:n_delta_v-1;
delta = sqrt(delta_x_(1+ndv)^2 + delta_y_(1+ndv)^2);
if (delta>delta_max); delta_max = delta; end;%if
end;%for ndv=0:n_delta_v-1;
n_pixel = delta_max/sqrt(2.0d0)*2.0d0*K_max;
if (verbose>0); disp(sprintf(' %% R_max: %0.2f, delta_max: %0.2f, K_max: %0.2f, n_pixel: %0.2f',R_max,delta_max,K_max,n_pixel)); end;%if
if (n_pixel>5 & flag_warning); disp(sprintf(' %% Warning, n_pixel %0.2f too large in get_svd_FTK_2',n_pixel)); end;%if
if (eps_target<1.0d-6 & flag_warning); disp(sprintf(' %% Warning, eps_target %0.8f too small in get_svd_FTK_2',eps_target)); end;%if
svd_fname = sprintf('%s/gen_Jsvd_N%.2dl10e6.txt',dir_svd,round(10*n_pixel));
if (~exist(svd_fname,'file')); disp(sprintf(' %% Warning, %s not found',svd_fname)); end;
if (verbose>0); disp(sprintf(' %% svd_fname: %s',svd_fname)); end;% if;
FTK = gen_Jsvd_svdload_FTK(svd_fname);
if (verbose>0); disp(sprintf(' %% n_svd_r %d, n_svd_d %d, n_svd_l %d',FTK.n_svd_r,FTK.n_svd_d,FTK.n_svd_l)); end;%if
flag_s_ = zeros(FTK.n_svd_l,1);
n_svd_l_out_max = length(find(FTK.svd_s_>=eps_target));
svd_l_out_ = zeros(n_svd_l_out_max,1);
svd_U_d_out_ = zeros(n_svd_l_out_max*FTK.n_svd_d,1);
svd_s_out_ = zeros(n_svd_l_out_max,1);
svd_V_r_out_ = zeros(n_svd_l_out_max*FTK.n_svd_r,1);
n_svd_l_out = 0;
for nl=0:FTK.n_svd_l-1;
flag_s_(1+nl) = 0;
if (FTK.svd_s_(1+nl)<eps_target);
if (verbose>1); disp(sprintf(' %% nl %d, svd_s %0.2f, skipping',nl,FTK.svd_s_(1+nl))); end;%if (verbose>1);
end;%if (FTK.svd_s_(1+nl)<eps_target);
if (FTK.svd_s_(1+nl)>=eps_target);
flag_s_(1+nl) = 1;
if (verbose>1); disp(sprintf(' %% nl %d, svd_s %0.2f, retaining: n_svd_l_out %d',nl,FTK.svd_s_(1+nl),n_svd_l_out)); end;%if (verbose>1);
svd_l_out_(1+n_svd_l_out) = FTK.svd_l_(1+nl);
svd_U_d_out_(1+n_svd_l_out*FTK.n_svd_d + (0:FTK.n_svd_d-1)) = FTK.svd_U_d_(1+nl*FTK.n_svd_d + (0:FTK.n_svd_d-1));
svd_s_out_(1+n_svd_l_out) = FTK.svd_s_(1+nl);
svd_V_r_out_(1+n_svd_l_out*FTK.n_svd_r + (0:FTK.n_svd_r-1)) = FTK.svd_V_r_(1+nl*FTK.n_svd_r + (0:FTK.n_svd_r-1));
n_svd_l_out = n_svd_l_out + 1;
end;%if (FTK.svd_s_(1+nl)>=eps_target);
end;%for nl=0:FTK.n_svd_l-1;
if (verbose>0); disp(sprintf(' %% sum(flag_s_) %d',sum(flag_s_))); end;%if (verbose>0);
FTK.eps_target = eps_target;
FTK.R_max = R_max;
FTK.K_max = K_max;
FTK.n_pixel = n_pixel;
FTK.n_svd_l = n_svd_l_out;
FTK.svd_l_ = svd_l_out_;
FTK.svd_U_d_ = svd_U_d_out_;
FTK.svd_s_ = svd_s_out_;
FTK.svd_V_r_ = svd_V_r_out_;
if (verbose>0); disp(sprintf(' %% [finished get_svd_FTK_2]')); end;%if (verbose>0);
