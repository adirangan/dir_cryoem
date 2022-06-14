%%%%%%%%;
% makes vfigs (i.e., volume-figs) for Principled_Marching_2?.tex. ;
%%%%%%%%;

%%%%%%%%;
% first recapitulates test_pm_trpv1_2. ;
%%%%%%%%;
test_pm_trpv1c_9b;
flag_invert=0;

%%%%%%%%;
% Now generate figures for paper. ;
%%%%%%%%;
% rsync -avum rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_trpv1c/dir_pm_mat/X_2d_xcor_d0_a1t????n18r?_a_k_Y_.mat /home/rangan/dir_cryoem/dir_trpv1c/dir_pm_mat/ ;

fname_mat = sprintf('%s_mat/pm_vfig_b_x_u_xi__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
fname_b_k_Y_ = { ...
 'X_2d_xcor_d0_a1t0014n18r0_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0029n18r0_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0059n18r0_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0089n18r0_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0014n18r1_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0029n18r1_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0059n18r1_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0089n18r1_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0014n18r2_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0029n18r2_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0059n18r2_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0089n18r2_a_k_Y_.mat' ...
};
n_fname_b_k_Y = numel(fname_b_k_Y_);
b_k_Y_yi__ = zeros(n_lm_sum,n_fname_b_k_Y);
for nfname_b_k_Y=0:n_fname_b_k_Y-1;
fname_b_k_Y = fname_b_k_Y_{1+nfname_b_k_Y};
tmp_ = load(sprintf('%s_mat/%s',dir_pm,fname_b_k_Y));
b_k_Y_yi__(:,1+nfname_b_k_Y) = tmp_.a_k_Y_reco_;
clear tmp_;
end;%for nfname_b_k_Y=0:n_fname_b_k_Y-1;
%%%%%%%%;
b_k_p_ki__ = zeros(n_k_all,n_fname_b_k_Y);
b_x_u_xi__ = zeros(n_x_u_pack^3,n_fname_b_k_Y);
for nfname_b_k_Y=0:n_fname_b_k_Y-1;
tmp_b_k_Y_ = b_k_Y_yi__(:,1+nfname_b_k_Y);
tmp_t = tic;
[tmp_b_k_p_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,tmp_b_k_Y_);
tmp_t = toc(tmp_t); disp(sprintf(' %% tmp_b_k_Y_ --> tmp_b_k_p_ time %0.2fs',tmp_t));
eta = pi/k_p_r_max; tmp_t = tic;
tmp_b_x_u_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,tmp_b_k_p_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: tmp_b_x_u_ time %0.2fs',tmp_t));
tmp_b_x_u_ = reshape(tmp_b_x_u_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]);
b_k_p_ki__(:,1+nfname_b_k_Y) = tmp_b_k_p_(:);
b_x_u_xi__(:,1+nfname_b_k_Y) = tmp_b_x_u_(:);
end;%for nfname_b_k_Y=0:n_fname_b_k_Y-1;
save(fname_mat ...
     ,'fname_b_k_Y_','n_fname_b_k_Y' ...
     ,'b_k_Y_yi__' ...
     ,'b_k_p_ki__' ...
     ,'b_x_u_xi__' ...
     );
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

fname_fig = sprintf('%s_jpg/pm_vfig_volume_FIGA_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;
p_row = 1; p_col = 4; ns=0;
tmp_n_x = round(0.5352*n_x_u_pack/2);
tmp_x_ij = n_x_u_pack/2 + [-tmp_n_x:+tmp_n_x];
tmp_prct = 90;
alim_ = mean(real(a_x_u_pack_(tmp_x_ij,tmp_x_ij,tmp_x_ij)),'all') + std(real(a_x_u_pack_(tmp_x_ij,tmp_x_ij,tmp_x_ij)),1,'all')*2.5*[-1,+1];
for pcol=0:p_col-1;
tmp_a_x_u_ = reshape(a_x_u_pack_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]);
tmp_str = 'ground truth';
if (pcol> 0);
tmp_a_x_u_ = reshape(b_x_u_xi__(:,1+pcol-1),[n_x_u_pack,n_x_u_pack,n_x_u_pack]);
tmp_str = sprintf('sample %d',1+pcol-1);
end;%if (pcol> 0);
subplot(p_row,p_col,1+0*p_col+pcol);cla;
isosurface_f_x_u_0(tmp_a_x_u_(tmp_x_ij,tmp_x_ij,tmp_x_ij),tmp_prct);
axisnotick; xlabel([]); ylabel([]); zlabel([]); title(tmp_str);
%subplot(p_row,p_col,1+1*p_col+pcol);cla;
%imagesc(squeeze(real(tmp_a_x_u_(tmp_x_ij,tmp_x_ij,round(n_x_u_pack/2)))),alim_);
%subplot(p_row,p_col,1+2*p_col+pcol);cla;
%imagesc(squeeze(real(tmp_a_x_u_(tmp_x_ij,round(n_x_u_pack/2),tmp_x_ij))),alim_);
end;%for pcol=0:p_col-1;
set(gcf,'Position',1+[0,0,1024+512,512]);
%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% Now find principal-modes for ground-truth volume. ;
%%%%%%%%;
a_k_Y_mlk_ = a_k_Y_quad_;
a_k_Y_mlk___ = convert_spharm_to_spharm___0(n_k_p_r,l_max_,a_k_Y_mlk_);
n_l_max = 1+l_max_max;
%%%%%%%%;
% calculate radial principal-modes. ;
%%%%%%%%;
tmp_t = tic();
X_kk__ = zeros(n_k_p_r,n_k_p_r);
X_weight_r_ = sqrt(weight_3d_k_p_r_);
for nk_p_r_0=0:n_k_p_r-1;
a_k_Y_ml0__ = a_k_Y_mlk___(:,:,1+nk_p_r_0);
X_weight_0 = X_weight_r_(1+nk_p_r_0);
for nk_p_r_1=0:n_k_p_r-1;
a_k_Y_ml1__ = a_k_Y_mlk___(:,:,1+nk_p_r_1);
X_weight_1 = X_weight_r_(1+nk_p_r_1);
X_weight_01 = X_weight_0*X_weight_1;
tmp_X = sum(conj(a_k_Y_ml0__(:,2:end)).*a_k_Y_ml1__(:,2:end),'all');
X_kk__(1+nk_p_r_0,1+nk_p_r_1) = real(tmp_X)*X_weight_01*(4*pi)^2;
end;%for nk_p_r_1=0:n_k_p_r-1;
end;%for nk_p_r_0=0:n_k_p_r-1;
tmp_X_kk__ = principled_marching_cost_matrix_2(n_k_p_r,weight_3d_k_p_r_,l_max_max,a_k_Y_quad__);
disp(sprintf(' %% X_kk__ vs tmp_X_kk__: %0.16f',fnorm(X_kk__-tmp_X_kk__)/fnorm(X_kk__)));
clear tmp_X_kk__;
tmp_t = toc(tmp_t);
if (verbose); disp(sprintf(' %% X_kk__: %0.3fs',tmp_t)); end;
tot_t_X_kk = tmp_t;
%%%%%%%%;
% Now compress k. ;
%%%%%%%%;
tmp_t = tic();
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
[UX__,SX__,VX__] = svds(X_kk__,n_UX_rank); SX_ = diag(SX__);
a_UX_Y_mln___ = zeros(n_m_max,n_l_max,n_UX_rank);
a_UX_Y_mln___ = reshape(reshape(a_k_Y_mlk___,[n_m_max*n_l_max,n_k_p_r])*diag(X_weight_r_)*UX__,[n_m_max,n_l_max,n_UX_rank]);
tmp_t = toc(tmp_t);
if (verbose); disp(sprintf(' %% UX__: %0.3fs',tmp_t)); end;
tot_t_UX = tmp_t;
%%%%%%%%;
% Now calculate order principal-modes. ;
%%%%%%%%;
tmp_t = tic();
pm_n_UX_rank = n_UX_rank;
Z_ll__ = zeros(n_l_max,n_l_max);
for nl_max_0=1:n_l_max-1;
a_UX_Y_m0n__ = a_UX_Y_mln___(:,1+nl_max_0,1:pm_n_UX_rank);
for nl_max_1=1:n_l_max-1;
a_UX_Y_m1n__ = a_UX_Y_mln___(:,1+nl_max_1,1:pm_n_UX_rank);
tmp_Z = sum(conj(a_UX_Y_m0n__).*a_UX_Y_m1n__,'all');
Z_ll__(1+nl_max_0,1+nl_max_1) = real(tmp_Z)*(4*pi)^2;
end;%for nl_max_1=1:n_l_max-1;
end;%for nl_max_0=1:n_l_max-1;
tmp_t = toc(tmp_t);
if (verbose); disp(sprintf(' %% Z_ll__: %0.3fs',tmp_t)); end;
tot_t_Z_ll = tmp_t;
%%%%%%%%;
tmp_t = tic();
n_UZ_rank = n_l_max-0; %<-- to match norm. ;
[UZ__,SZ__,VZ__] = svds(Z_ll__,n_UZ_rank); SZ_ = diag(SZ__);
tmp_t = toc(tmp_t);
if (verbose); disp(sprintf(' %% UZ__: %0.3fs',tmp_t)); end;
tot_t_UZ = tmp_t;
%%%%%%%%;
flag_plot=0;
if flag_plot;
llim_ = [-20,-5];
%%%%;
figure(1);clf;figbig;figbeach;
for nk_p_r=0:n_k_p_r-1;
subplot(7,7,1+nk_p_r);
imagesc(log(abs(a_k_Y_mlk___(:,:,1+nk_p_r))),llim_)
end;%for nk_p_r=0:16-1;
title('a_k_Y_','Interpreter','none');
%%%%;
llim_ = [-20,-5];
figure(2);clf;figbig;figbeach;
for nUX_rank=0:n_UX_rank-1;
subplot(7,7,1+nUX_rank);
imagesc(log(abs(a_UX_Y_mln___(:,:,1+nUX_rank))),llim_);
end;%for nUX_rank=0:16-1;
title('a_UX_Y_','Interpreter','none');
%%%%;
end;%if flag_plot;

%%%%%%%%;
% Now perform stripped down calculation of innerproducts. ;
%%%%%%%%;
n_UZ_rank = n_l_max;
%%%%%%%%;
% calculate d_. ;
%%%%%%%%;
n_beta = n_w_max; beta_ = transpose(linspace(-pi,+pi,n_beta+1)); beta_ = beta_(1:end-1);
t_0in = tic;
d_mmlb____ = zeros(n_m_max,n_m_max,n_l_max,n_beta);
d_mmzb____ = zeros(n_m_max,n_m_max,n_UZ_rank,n_beta);
for nbeta=0:n_beta-1;
if (verbose); if (mod(nbeta,16)==0); disp(sprintf(' %% nbeta %d/%d',nbeta,n_beta)); end;
beta = beta_(1+nbeta);
W_ = wignerd_b(l_max_max,-beta);
n_W_ = zeros(1,1+l_max_max); for (l_val=0:l_max_max); n_W_(1+l_val) = numel(W_{1+l_val}); end;
d_mml___ = zeros(n_m_max,n_m_max,1+l_max_max);
for l_val=0:l_max_max;
d_mml___(1+l_max_max + [-l_val:+l_val],1+l_max_max + [-l_val:+l_val],1+l_val) = W_{1+l_val};
end;%for l_val=0:l_max_max;
d_mmz___ = reshape(reshape(d_mml___,[n_m_max*n_m_max,n_l_max])*UZ__,[n_m_max,n_m_max,n_UZ_rank]);
d_mmlb____(:,:,:,1+nbeta) = d_mml___;
d_mmzb____(:,:,:,1+nbeta) = d_mmz___;
end;%for nbeta=0:n_beta-1;
t_out = toc(t_0in);
if (verbose>1); disp(sprintf(' %% calculate d_: t %0.6f',t_out)); end;
%%%%%%%%;
% load b_. ;
%%%%%%%%;
tmp_ = load('/home/rangan/dir_cryoem/dir_trpv1_x0/dir_pm_mat/a_k_Y_0qbp_reco_.mat');
b_k_Y_ = tmp_.a_k_Y_0qbp_reco_;
b_x_u_ = tmp_.a_x_u_0qbp_reco_;
clear tmp_;
b_k_Y_mlk_ = b_k_Y_;
b_k_Y_mlk___ = convert_spharm_to_spharm___0(n_k_p_r,l_max_,b_k_Y_mlk_);
%%%%%%%%;

%%%%%%%%;
% compress b_. ;
%%%%%%%%;
t_0in = tic;
b_UX_Y_mln___ = reshape(reshape(b_k_Y_mlk___,[n_m_max*n_l_max,n_k_p_r])*diag(X_weight_r_)*UX__,[n_m_max,n_l_max,n_UX_rank]);
t_out = toc(t_0in);
if (verbose>1); disp(sprintf(' %% compress b: t %0.6f',t_out)); end;

%%%%%%%%;
% Now test. ;
%%%%%%%%;
flag_test=0;
if flag_test;
n_test = 8;
for ntest=0:n_test-1;
nbeta=max(0,min(n_beta-1,floor(n_beta*ntest/n_test)));
beta = beta_(1+nbeta);
d_mml___ = d_mmlb____(:,:,:,1+nbeta);
d_mmz___ = d_mmzb____(:,:,:,1+nbeta);
tmp_X0_mm__ = register_spharm_to_spharm_single_beta_3(0,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_mlk_,b_k_Y_mlk_,beta);
%%%%;
a_k_Y_mkl___ = permute(a_k_Y_mlk___,[1,3,2]); b_k_Y_mkl___ = permute(b_k_Y_mlk___,[1,3,2]);
tmp_X1_mm__ = register_spharm_to_spharm_single_beta_3_stripped_0(0,n_m_max,n_k_p_r,n_l_max,weight_3d_k_p_r_,a_k_Y_mkl___,b_k_Y_mkl___,[],[],1,d_mml___);
disp(sprintf(' %% tmp_X0_mm__ vs tmp_X1_mm__: %0.16f',fnorm(tmp_X0_mm__-tmp_X1_mm__)/fnorm(tmp_X0_mm__)));
disp(sprintf(' %% corr(tmp_X1_mm__,tmp_X0_mm__): %0.16f',corr(tmp_X1_mm__(:),tmp_X0_mm__(:))));
clear tmp_X1_mm__;
%%%%;
a_UX_Y_mnl___ = permute(a_UX_Y_mln___,[1,3,2]); b_UX_Y_mnl___ = permute(b_UX_Y_mln___,[1,3,2]);
tmp_X1_mm__ = register_spharm_to_spharm_single_beta_3_stripped_0(0,n_m_max,n_UX_rank,n_l_max,[],a_UX_Y_mnl___,b_UX_Y_mnl___,[],[],1,d_mml___);
disp(sprintf(' %% tmp_X0_mm__ vs tmp_X1_mm__: %0.16f',fnorm(tmp_X0_mm__-tmp_X1_mm__)/fnorm(tmp_X0_mm__)));
disp(sprintf(' %% corr(tmp_X1_mm__,tmp_X0_mm__): %0.16f',corr(tmp_X1_mm__(:),tmp_X0_mm__(:))));
clear tmp_X1_mm__;
%%%%;
tmp_X1_mm__ = register_spharm_to_spharm_single_beta_3_stripped_0(0,n_m_max,n_UX_rank,n_l_max,[],a_UX_Y_mnl___,b_UX_Y_mnl___,n_UZ_rank,UZ__,1,d_mmz___);
disp(sprintf(' %% tmp_X0_mm__ vs tmp_X1_mm__: %0.16f',fnorm(tmp_X0_mm__-tmp_X1_mm__)/fnorm(tmp_X0_mm__)));
disp(sprintf(' %% corr(tmp_X1_mm__,tmp_X0_mm__): %0.16f',corr(tmp_X1_mm__(:),tmp_X0_mm__(:))));
clear tmp_X1_mm__;
%%%%;
end;%for ntest=0:n_test-1;
end;%if flag_test;

fname_mat = sprintf('%s_mat/pm_vfig_timing_FIGA__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
% Now time. ;
%%%%%%%%;
drank = 1;
[ ...
 X0_mmb___ ...
 t_sumn_0 ...
 t_abab_0 ...
 t_sumz_0 ...
 t_fft2_0 ...
 n_mult_0 ...
] = ...
register_spharm_to_spharm_single_beta_3_stripped_0( ...
 verbose ...
,n_m_max ...
,n_k_p_r ...
,n_l_max ...
,weight_3d_k_p_r_ ...
,a_k_Y_mkl___ ...
,b_k_Y_mkl___ ...
,[] ...
,[] ...
,n_beta ...
,d_mmlb____ ...
);
X0_mmb___ = real(X0_mmb___);
%%%%%%%%;
t_sumn_nz__ = zeros(n_UX_rank,n_UZ_rank);
t_abab_nz__ = zeros(n_UX_rank,n_UZ_rank);
t_sumz_nz__ = zeros(n_UX_rank,n_UZ_rank);
t_fft2_nz__ = zeros(n_UX_rank,n_UZ_rank);
n_mult_nz__ = zeros(n_UX_rank,n_UZ_rank);
l2er_X_X_nz__ = zeros(n_UX_rank,n_UZ_rank);
corr_X_X_nz__ = zeros(n_UX_rank,n_UZ_rank);
prct_X_X_nz__ = zeros(n_UX_rank,n_UZ_rank);
%%%%%%%%;
for nUZ_rank=n_UZ_rank-1:-drank:0;
if (verbose); disp(sprintf(' %% nUZ_rank %d/%d',nUZ_rank,n_UZ_rank)); end;
for nUX_rank=n_UX_rank-1:-drank:0;
%%%%%%%%;
[ ...
 X1_mmb___ ...
 t_sumn_1 ...
 t_abab_1 ...
 t_sumz_1 ...
 t_fft2_1 ...
 n_mult_1 ...
] = ...
register_spharm_to_spharm_single_beta_3_stripped_0( ...
 0*verbose ...
,n_m_max ...
,1+nUX_rank ...
,n_l_max ...
,[] ...
,a_UX_Y_mnl___ ...
,b_UX_Y_mnl___ ...
,1+nUZ_rank ...
,UZ__ ...
,n_beta ...
,d_mmzb____ ...
);
X1_mmb___ = real(X1_mmb___);
t_sumn_nz__(1+nUX_rank,1+nUZ_rank) = t_sumn_1;
t_abab_nz__(1+nUX_rank,1+nUZ_rank) = t_abab_1;
t_sumz_nz__(1+nUX_rank,1+nUZ_rank) = t_sumz_1;
t_fft2_nz__(1+nUX_rank,1+nUZ_rank) = t_fft2_1;
n_mult_nz__(1+nUX_rank,1+nUZ_rank) = n_mult_1;
%%%%;
l2er_X_X = fnorm(X0_mmb___ - X1_mmb___)/fnorm(X0_mmb___);
corr_X_X = corr(X0_mmb___(:),X1_mmb___(:));
[~,tmp_index_nmmb] = max(X1_mmb___,[],'all','linear'); tmp_index_nmmb = tmp_index_nmmb-1;
tmp_X0 = X0_mmb___(1+tmp_index_nmmb);
prct_X_X = sum(X0_mmb___(:)>tmp_X0)/numel(X0_mmb___);
%%%%;
l2er_X_X_nz__(1+nUX_rank,1+nUZ_rank) = l2er_X_X;
corr_X_X_nz__(1+nUX_rank,1+nUZ_rank) = corr_X_X;
prct_X_X_nz__(1+nUX_rank,1+nUZ_rank) = prct_X_X;
%%%%%%%%;
end;%for nUX_rank=n_UX_rank-1:-1:0;
end;%for nUZ_rank=n_UZ_rank-1:-1:0;
%%%%%%%%;
save(fname_mat ...
     ,'n_UX_rank','n_UZ_rank' ...
     ,'t_sumn_0','t_abab_0','t_sumz_0','t_fft2_0' ...
     ,'t_sumn_nz__' ...
     ,'t_abab_nz__' ...
     ,'t_sumz_nz__' ...
     ,'t_fft2_nz__' ...
     ,'n_mult_nz__' ...
     ,'l2er_X_X_nz__' ...
     ,'corr_X_X_nz__' ...
     ,'prct_X_X_nz__' ...
     );
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

flag_plot=0;
if flag_plot;
figure(1);clf;figbig;figbeach();
tmp_ij_n_ = 1+flip(n_UX_rank-1:-drank:0);
tmp_ij_z_ = 1+flip(n_UZ_rank-1:-drank:0);
p_row = 2; p_col = 4; ns=0;
t_comp_0 = t_sumn_0 + t_abab_0 + t_sumz_0 + t_fft2_0;
t_comp_nz__ = t_sumn_nz__ + t_abab_nz__ + t_sumz_nz__ + t_fft2_nz__;
subplot(p_row,p_col,1+ns);ns=ns+1; imagesc(t_sumn_nz__(tmp_ij_n_,tmp_ij_z_)); colorbar;
subplot(p_row,p_col,1+ns);ns=ns+1; imagesc(t_abab_nz__(tmp_ij_n_,tmp_ij_z_)); colorbar;
subplot(p_row,p_col,1+ns);ns=ns+1; imagesc(t_sumz_nz__(tmp_ij_n_,tmp_ij_z_)); colorbar;
subplot(p_row,p_col,1+ns);ns=ns+1; imagesc(n_mult_nz__(tmp_ij_n_,tmp_ij_z_)); colorbar;
subplot(p_row,p_col,1+ns);ns=ns+1; imagesc(t_comp_nz__(tmp_ij_n_,tmp_ij_z_)); colorbar;
subplot(p_row,p_col,1+ns);ns=ns+1; imagesc(0+l2er_X_X_nz__(tmp_ij_n_,tmp_ij_z_)); colorbar;
subplot(p_row,p_col,1+ns);ns=ns+1; imagesc(0+corr_X_X_nz__(tmp_ij_n_,tmp_ij_z_)); colorbar;
subplot(p_row,p_col,1+ns);ns=ns+1; imagesc(1-prct_X_X_nz__(tmp_ij_n_,tmp_ij_z_)); colorbar;
end;%if flag_plot;












