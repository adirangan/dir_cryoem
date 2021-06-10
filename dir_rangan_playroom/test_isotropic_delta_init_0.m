% test isotropic delta estimation for the images. ;
% Verdict: does not seem to be better than center of mass for rib80s.; 
% Obviously this would be better if the average mass is 0. ;
% But then we can try calculating the center of the absolute-value of the mass. ;

fname_euler_angle = sprintf('%s/euler_angles',dir_data);
fp = fopen(fname_euler_angle,'r');
euler_angle_load_ = textscan(fp,'%f%f%f\n%f%f\n');
fclose(fp);
euler_angle_marina_ = zeros(3,n_image); %<-- [polar_a,azimu_b,gamma_z]. ;
for nimage=0:n_image-1;
euler_angle_marina_(:,1+nimage) = convert_euler_relion_to_marina([euler_angle_load_{1}(1+nimage),euler_angle_load_{2}(1+nimage),euler_angle_load_{3}(1+nimage)]);
end;%for nimage=0:n_image-1;
delta_read_x_ = euler_angle_load_{4}*(2/n_x_u);
delta_read_y_ = euler_angle_load_{5}*(2/n_x_u);
% Note that delta_read_r_ is not too large. ;
% prctile(sqrt(delta_read_x_(:).^2 + delta_read_y_(:).^2),95) ; %<-- 0.1182. ;
grid_x_u_ = linspace(-1,+1,n_x_u+1); grid_x_u_ = x_p_r_max*grid_x_u_(1:end-1);
fname_image_mda = sprintf('%s/images_mda',dir_data);
M_x_c___ = MDA_read_r8(fname_image_mda);
n_image_sub = min(n_image_sub,size(M_x_c___,3));
%%%%%%%%;
flag_plot=0;
if flag_plot;
figure(1);clf;
imagesc(log10(1+hist2d_0(delta_read_x_,delta_read_y_,n_x_u,n_x_u,[-1,+1],[-1,+1])),[0,3]); 
axis image; xlabel('x'); ylabel('y'); title('delta read');
colormap(colormap_beach());
figbig;
end;%if flag_plot;
%%%%%%%%;
% Simple center of mass centering. ;
%%%%%%%%;
x_c_0_ = -half_diameter_x_c + transpose([0:n_x_u-1]/n_x_u)*diameter_x_c;
x_c_1_ = -half_diameter_x_c + transpose([0:n_x_u-1]/n_x_u)*diameter_x_c;
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_);
M_abs_x_c_0_avg_ = zeros(n_image_sub,1);
M_abs_x_c_1_avg_ = zeros(n_image_sub,1);
for nimage_sub=0:n_image_sub-1;
M_abs_x_c_ = abs(squeeze(M_x_c___(:,:,1+nimage_sub)));
M_abs_avg = mean(M_abs_x_c_,'all');
M_abs_x_c_0_avg = mean(M_abs_x_c_/M_abs_avg.*x_c_0__,'all');
M_abs_x_c_1_avg = mean(M_abs_x_c_/M_abs_avg.*x_c_1__,'all');
M_abs_x_c_0_avg_(1+nimage_sub) = M_abs_x_c_0_avg;
M_abs_x_c_1_avg_(1+nimage_sub) = M_abs_x_c_1_avg;
clear M_abs_x_c_;
end;%for nimage_sub=0:n_image_sub-1;
%%%%%%%%;
% Now visualize. ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
nimage_sub=4;
M_x_c__ = squeeze(M_x_c___(:,:,1+nimage_sub));
M_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,M_x_c__,n_k_p_r,k_p_r_,n_w_) ;
N_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,-1*M_abs_x_c_0_avg_(1+nimage_sub),-1*M_abs_x_c_1_avg_(1+nimage_sub));
N_x_c__ = interp_k_p_to_x_c_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,N_k_p_.*weight_2d_k_all_);
O_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,+1*delta_read_x_(1+nimage_sub),+1*delta_read_y_(1+nimage_sub));
O_x_c__ = interp_k_p_to_x_c_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,O_k_p_.*weight_2d_k_all_);
figure(3);clf;figbeach();
subplot(1,3,1); imagesc(real(M_x_c__)); axis image; title('M');
subplot(1,3,2); imagesc(real(N_x_c__)); axis image; title('N');
subplot(1,3,3); imagesc(real(O_x_c__)); axis image; title('O');
end;%if flag_plot;
%%%%%%%%;
N_x_c___ = M_x_c___;
O_x_c___ = M_x_c___;
for nimage_sub=0:n_image_sub-1;
if (mod(nimage_sub,128)==0); disp(sprintf(' %% nimage_sub %d/%d',nimage_sub,n_image_sub)); end;
M_x_c__ = squeeze(M_x_c___(:,:,1+nimage_sub));
M_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,M_x_c__,n_k_p_r,k_p_r_,n_w_) ;
N_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,-1*M_abs_x_c_0_avg_(1+nimage_sub),-1*M_abs_x_c_1_avg_(1+nimage_sub));
N_x_c___(:,:,1+nimage_sub) = interp_k_p_to_x_c_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,N_k_p_.*weight_2d_k_all_);
O_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,+1*delta_read_x_(1+nimage_sub),+1*delta_read_y_(1+nimage_sub));
O_x_c___(:,:,1+nimage_sub) = interp_k_p_to_x_c_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,O_k_p_.*weight_2d_k_all_);
end;%for nimage_sub=0:n_image_sub-1;
M_x_c_avg__ = real(mean(M_x_c___,3));
N_x_c_avg__ = real(mean(N_x_c___,3));
O_x_c_avg__ = real(mean(O_x_c___,3));
%%%%%%%%;
flag_plot=1;
if flag_plot;
figure(2);clf;
subplot(2,2,1); imagesc(M_x_c_avg__); axis image; axisnotick; colormap(colormap_beach()); title('M (orig)');
subplot(2,2,2); imagesc(N_x_c_avg__); axis image; axisnotick; colormap(colormap_beach()); title('N (-M_abs_x_c)','Interpreter','none');
subplot(2,2,3); imagesc(O_x_c_avg__); axis image; axisnotick; colormap(colormap_beach()); title('O (+delta)');
subplot(2,2,4); imagesc(P_x_c_avg__); axis image; axisnotick; colormap(colormap_beach()); title('P (none)');
figbig;
end;%if flag_plot;
%%%%%%%%;
delta_read_plus_M_abs_x_c_0_avg_ = delta_read_x_(1:n_image_sub) + M_abs_x_c_0_avg_;
delta_read_plus_M_abs_x_c_1_avg_ = delta_read_y_(1:n_image_sub) + M_abs_x_c_1_avg_;
%%%%%%%%;
% Now set up annular kernels ;
%%%%%%%%;
n_annulus = ceil(0.5*k_p_r_max);
annulus_x_c___ = zeros(n_x_u,n_x_u,n_annulus);
x_p_r__ = sqrt(x_c_0__.^2 + x_c_1__.^2);
x_p_w__ = atan2(x_c_1__,x_c_0__);
for nannulus=0:n_annulus-1;
tmp_x_p_r_lo = 0.50*half_diameter_x_c*(nannulus+0)/n_annulus;
tmp_x_p_r_hi = 0.50*half_diameter_x_c*(nannulus+1)/n_annulus;
annulus_x_c__ = (x_p_r__>=tmp_x_p_r_lo) & (x_p_r__< tmp_x_p_r_hi);
annulus_x_c__ = annulus_x_c__/sum(annulus_x_c__,'all') ;
annulus_x_c___(:,:,1+nannulus) = annulus_x_c__;
end;%for nannulus=0:n_annulus-1;
%%%%%%%%;
% Now convolve. ;
% We expect the convolution evaluated at (y0,y1) to return: ;
% \integral_{dx} M_x_c__(x) * annulus_x_c__(x-y) ;
%%%%%%%%;
fft2_recenter2_annulus_x_c___ = zeros(n_x_u,n_x_u,n_annulus);
for nannulus=0:n_annulus-1;
annulus_x_c__ = annulus_x_c___(:,:,1+nannulus);
fft2_recenter2_annulus_x_c___(:,:,1+nannulus) = fft2(recenter2(annulus_x_c__));
end;%for nannulus=0:n_annulus-1;
L_x_c_xyaM____ = zeros(n_x_u,n_x_u,n_annulus,n_image_sub);
for nimage_sub=0:n_image_sub-1;
if (mod(nimage_sub,128)==0); disp(sprintf(' %% nimage_sub %d/%d',nimage_sub,n_image_sub)); end;
M_x_c__ = squeeze(M_x_c___(:,:,1+nimage_sub));
fft2_recenter2_M_x_c__ = fft2(recenter2(M_x_c__));
for nannulus=0:n_annulus-1;
L_x_c_xyaM____(:,:,1+nannulus,1+nimage_sub) = decenter2(ifft2(fft2_recenter2_M_x_c__.*fft2_recenter2_annulus_x_c___(:,:,1+nannulus)));
end;%for nannulus=0:n_annulus-1;
end;%for nimage_sub=0:n_image_sub-1;
%%%%%%%%;
% Now use message passing to determine center for each image. ;
%%%%%%%%;
%%%%%%%%;
% initialize. ;
%%%%%%%%;
index_center_ini_ = zeros(n_image_sub,1);
for nimage_sub=0:n_image_sub-1;
L_x_c_xy__ = sum(L_x_c_xyaM____(:,:,:,1+nimage_sub).^2,3);
[~,index_center] = max(L_x_c_xy__(:)); index_center = index_center-1;
index_center_ini_(1+nimage_sub) = index_center;
end;%for nimage_sub=0:n_image_sub-1;
annular_ini_x_c_0_ = zeros(n_image_sub,1);
annular_ini_x_c_1_ = zeros(n_image_sub,1);
for nimage_sub=0:n_image_sub-1;
index_center = index_center_ini_(1+nimage_sub);
nx0 = mod(index_center,n_x_u);
nx1 = (index_center-nx0)/n_x_u;
annular_ini_x_c_0_(1+nimage_sub) = x_c_0_(1+nx0);
annular_ini_x_c_1_(1+nimage_sub) = x_c_1_(1+nx1);
end;%for nimage_sub=0:n_image_sub-1;
%%%%%%%%;
% iterate. ;
%%%%%%%%;
index_center_pre_ = index_center_ini_;
n_iteration = 16;
niteration=0; flag_continue=1;
while (flag_continue);
if (verbose); disp(sprintf(' %% niteration %d/%d',niteration,n_iteration)); end;
%%%%;
L_x_c_a_ = zeros(n_annulus,1);
for nimage_sub=0:n_image_sub-1;
index_center = index_center_pre_(1+nimage_sub);
nx0 = mod(index_center,n_x_u);
nx1 = (index_center-nx0)/n_x_u;
L_x_c_a_ = L_x_c_a_ + reshape(L_x_c_xyaM____(1+nx0,1+nx1,:,1+nimage_sub),[n_annulus,1])/n_image_sub;
end;%for nimage_sub=0:n_image_sub-1;
L_x_c_a_rep___ = repmat(reshape(L_x_c_a_,[1,1,n_annulus]),[n_x_u,n_x_u,1]);
%%%%;
for nimage_sub=0:n_image_sub-1;
L_x_c_xy__ = sum((L_x_c_xyaM____(:,:,:,1+nimage_sub) - L_x_c_a_rep___).^2,3);
[~,index_center] = min(L_x_c_xy__(:)); index_center = index_center-1;
index_center_pos_(1+nimage_sub) = index_center;
end;%for nimage_sub=0:n_image_sub-1;
%%%%;
niteration=niteration+1;
flag_continue = (niteration<n_iteration) & (fnorm(index_center_pos_-index_center_pre_)>0);
index_center_pre_ = index_center_pos_;
%%%%;
end;%while (flag_continue);
%%%%;
annular_align_x_c_0_ = zeros(n_image_sub,1);
annular_align_x_c_1_ = zeros(n_image_sub,1);
for nimage_sub=0:n_image_sub-1;
index_center = index_center_pos_(1+nimage_sub);
nx0 = mod(index_center,n_x_u);
nx1 = (index_center-nx0)/n_x_u;
annular_align_x_c_0_(1+nimage_sub) = x_c_0_(1+nx0);
annular_align_x_c_1_(1+nimage_sub) = x_c_1_(1+nx1);
end;%for nimage_sub=0:n_image_sub-1;
%%%%%%%%;
% Now compare. ;
%%%%%%%%;
figure(4);clf;
subplot(2,2,1); plot(delta_read_x_(1:n_image_sub),M_abs_x_c_0_avg_,'b.',delta_read_x_(1:n_image_sub),annular_ini_x_c_0_,'r.');
xlabel('delta_read_x','Interpreter','none'); ylabel('approx'); legend({'M abs','annular ini'});
subplot(2,2,2); plot(delta_read_y_(1:n_image_sub),M_abs_x_c_1_avg_,'b.',delta_read_y_(1:n_image_sub),annular_ini_x_c_1_,'r.');
xlabel('delta_read_y','Interpreter','none'); ylabel('approx'); legend({'M abs','annular ini'});
subplot(2,2,3); plot(delta_read_x_(1:n_image_sub),M_abs_x_c_0_avg_,'b.',delta_read_x_(1:n_image_sub),annular_align_x_c_0_,'r.');
xlabel('delta_read_x','Interpreter','none'); ylabel('approx'); legend({'M abs','annular align'});
subplot(2,2,4); plot(delta_read_y_(1:n_image_sub),M_abs_x_c_1_avg_,'b.',delta_read_y_(1:n_image_sub),annular_align_x_c_1_,'r.');
xlabel('delta_read_y','Interpreter','none'); ylabel('approx'); legend({'M abs','annular align'});
figbig;

%%%%%%%%;
% Now try the same thing with an isotropic molecule. ;
% Note how well it works! ;
%%%%%%%%;
snr = 0.25;
iso_S_x_c__ = exp(-x_p_r__.^2/(2*0.25^2)).*sin(2*pi*8*x_p_r__);
iso_S_l2 = fnorm(iso_S_x_c__);
iso_nx0_ = max(0,min(n_x_u-1,floor(n_x_u*rand(n_image_sub,1))));
iso_nx1_ = max(0,min(n_x_u-1,floor(n_x_u*rand(n_image_sub,1))));
iso_delta_true_x_c_0_ = x_c_0_(1+iso_nx0_);
iso_delta_true_x_c_1_ = x_c_1_(1+iso_nx1_);
iso_M_x_c___ = zeros(n_x_u,n_x_u,n_image_sub);
for nimage_sub=0:n_image_sub-1;
iso_M_x_c__ = circshift(circshift(iso_S_x_c__,1+iso_nx0_(1+nimage_sub)-n_x_u/2,1+0),1+iso_nx1_(1+nimage_sub)-n_x_u/2,1+1);
iso_M_x_c___(:,:,1+nimage_sub) = iso_M_x_c__ + (1.0/snr)*randn(n_x_u,n_x_u)*iso_S_l2/n_x_u;
clear iso_M_x_c__;
end;%for nimage_sub=0:n_image_sub-1;
%%%%%%%%;
iso_M_abs_x_c_0_avg_ = zeros(n_image_sub,1);
iso_M_abs_x_c_1_avg_ = zeros(n_image_sub,1);
for nimage_sub=0:n_image_sub-1;
iso_M_abs_x_c_ = abs(squeeze(iso_M_x_c___(:,:,1+nimage_sub)));
iso_M_abs_avg = mean(iso_M_abs_x_c_,'all');
iso_M_abs_x_c_0_avg = mean(iso_M_abs_x_c_/iso_M_abs_avg.*x_c_0__,'all');
iso_M_abs_x_c_1_avg = mean(iso_M_abs_x_c_/iso_M_abs_avg.*x_c_1__,'all');
iso_M_abs_x_c_0_avg_(1+nimage_sub) = iso_M_abs_x_c_0_avg;
iso_M_abs_x_c_1_avg_(1+nimage_sub) = iso_M_abs_x_c_1_avg;
clear iso_M_abs_x_c_;
end;%for nimage_sub=0:n_image_sub-1;
%%%%%%%%;
n_annulus = ceil(1.0*k_p_r_max);
annulus_x_c___ = zeros(n_x_u,n_x_u,n_annulus);
x_p_r__ = sqrt(x_c_0__.^2 + x_c_1__.^2);
x_p_w__ = atan2(x_c_1__,x_c_0__);
for nannulus=0:n_annulus-1;
tmp_x_p_r_lo = 0.25*half_diameter_x_c*(nannulus+0)/n_annulus;
tmp_x_p_r_hi = 0.25*half_diameter_x_c*(nannulus+1)/n_annulus;
annulus_x_c__ = (x_p_r__>=tmp_x_p_r_lo) & (x_p_r__< tmp_x_p_r_hi);
annulus_x_c__ = annulus_x_c__/sum(annulus_x_c__,'all') ;
annulus_x_c___(:,:,1+nannulus) = annulus_x_c__;
end;%for nannulus=0:n_annulus-1;
%%%%%%%%;
% Now convolve. ;
% We expect the convolution evaluated at (y0,y1) to return: ;
% \integral_{dx} iso_M_x_c__(x) * annulus_x_c__(x-y) ;
%%%%%%%%;
fft2_recenter2_annulus_x_c___ = zeros(n_x_u,n_x_u,n_annulus);
for nannulus=0:n_annulus-1;
annulus_x_c__ = annulus_x_c___(:,:,1+nannulus);
fft2_recenter2_annulus_x_c___(:,:,1+nannulus) = fft2(recenter2(annulus_x_c__));
end;%for nannulus=0:n_annulus-1;
L_x_c_xyaM____ = zeros(n_x_u,n_x_u,n_annulus,n_image_sub);
for nimage_sub=0:n_image_sub-1;
if (mod(nimage_sub,128)==0); disp(sprintf(' %% nimage_sub %d/%d',nimage_sub,n_image_sub)); end;
iso_M_x_c__ = squeeze(iso_M_x_c___(:,:,1+nimage_sub));
fft2_recenter2_M_x_c__ = fft2(recenter2(iso_M_x_c__));
for nannulus=0:n_annulus-1;
L_x_c_xyaM____(:,:,1+nannulus,1+nimage_sub) = decenter2(ifft2(fft2_recenter2_M_x_c__.*fft2_recenter2_annulus_x_c___(:,:,1+nannulus)));
end;%for nannulus=0:n_annulus-1;
end;%for nimage_sub=0:n_image_sub-1;
%%%%%%%%;
% Now use message passing to determine center for each image. ;
%%%%%%%%;
%%%%%%%%;
% initialize. ;
%%%%%%%%;
index_center_ini_ = zeros(n_image_sub,1);
for nimage_sub=0:n_image_sub-1;
L_x_c_xy__ = sum(L_x_c_xyaM____(:,:,:,1+nimage_sub).^2,3);
[~,index_center] = max(L_x_c_xy__(:)); index_center = index_center-1;
index_center_ini_(1+nimage_sub) = index_center;
end;%for nimage_sub=0:n_image_sub-1;
annular_ini_x_c_0_ = zeros(n_image_sub,1);
annular_ini_x_c_1_ = zeros(n_image_sub,1);
for nimage_sub=0:n_image_sub-1;
index_center = index_center_ini_(1+nimage_sub);
nx0 = mod(index_center,n_x_u);
nx1 = (index_center-nx0)/n_x_u;
annular_ini_x_c_0_(1+nimage_sub) = x_c_0_(1+nx0);
annular_ini_x_c_1_(1+nimage_sub) = x_c_1_(1+nx1);
end;%for nimage_sub=0:n_image_sub-1;
%%%%%%%%;
% iterate. ;
%%%%%%%%;
index_center_pre_ = index_center_ini_;
n_iteration = 16;
niteration=0; flag_continue=1;
while (flag_continue);
if (verbose); disp(sprintf(' %% niteration %d/%d',niteration,n_iteration)); end;
%%%%;
L_x_c_a_ = zeros(n_annulus,1);
for nimage_sub=0:n_image_sub-1;
index_center = index_center_pre_(1+nimage_sub);
nx0 = mod(index_center,n_x_u);
nx1 = (index_center-nx0)/n_x_u;
L_x_c_a_ = L_x_c_a_ + reshape(L_x_c_xyaM____(1+nx0,1+nx1,:,1+nimage_sub),[n_annulus,1])/n_image_sub;
end;%for nimage_sub=0:n_image_sub-1;
L_x_c_a_rep___ = repmat(reshape(L_x_c_a_,[1,1,n_annulus]),[n_x_u,n_x_u,1]);
%%%%;
for nimage_sub=0:n_image_sub-1;
L_x_c_xy__ = sum((L_x_c_xyaM____(:,:,:,1+nimage_sub) - L_x_c_a_rep___).^2,3);
[~,index_center] = min(L_x_c_xy__(:)); index_center = index_center-1;
index_center_pos_(1+nimage_sub) = index_center;
end;%for nimage_sub=0:n_image_sub-1;
%%%%;
niteration=niteration+1;
flag_continue = (niteration<n_iteration) & (fnorm(index_center_pos_-index_center_pre_)>0);
index_center_pre_ = index_center_pos_;
%%%%;
end;%while (flag_continue);
%%%%;
annular_align_x_c_0_ = zeros(n_image_sub,1);
annular_align_x_c_1_ = zeros(n_image_sub,1);
for nimage_sub=0:n_image_sub-1;
index_center = index_center_pos_(1+nimage_sub);
nx0 = mod(index_center,n_x_u);
nx1 = (index_center-nx0)/n_x_u;
annular_align_x_c_0_(1+nimage_sub) = x_c_0_(1+nx0);
annular_align_x_c_1_(1+nimage_sub) = x_c_1_(1+nx1);
end;%for nimage_sub=0:n_image_sub-1;
%%%%%%%%;
% Now compare. ;
%%%%%%%%;
figure(4);clf;
subplot(2,2,1); plot(iso_delta_true_x_c_0_(1:n_image_sub),iso_M_abs_x_c_0_avg_,'b.',iso_delta_true_x_c_0_(1:n_image_sub),annular_ini_x_c_0_,'r.');
xlabel('delta_read_x','Interpreter','none'); ylabel('approx'); legend({'M abs','annular ini'});
subplot(2,2,2); plot(iso_delta_true_x_c_1_(1:n_image_sub),iso_M_abs_x_c_1_avg_,'b.',iso_delta_true_x_c_1_(1:n_image_sub),annular_ini_x_c_1_,'r.');
xlabel('delta_read_y','Interpreter','none'); ylabel('approx'); legend({'M abs','annular ini'});
subplot(2,2,3); plot(iso_delta_true_x_c_0_(1:n_image_sub),iso_M_abs_x_c_0_avg_,'b.',iso_delta_true_x_c_0_(1:n_image_sub),annular_align_x_c_0_,'r.');
xlabel('delta_read_x','Interpreter','none'); ylabel('approx'); legend({'M abs','annular align'});
subplot(2,2,4); plot(iso_delta_true_x_c_1_(1:n_image_sub),iso_M_abs_x_c_1_avg_,'b.',iso_delta_true_x_c_1_(1:n_image_sub),annular_align_x_c_1_,'r.');
xlabel('delta_read_y','Interpreter','none'); ylabel('approx'); legend({'M abs','annular align'});
figbig;

