ctfw2___ = MDA_read_c16('./dir_test33b/ctfw2___.mda'); 
all_refim___ = MDA_read_c16('./dir_test33b/all_refim___.mda'); 
[ngridc,ngridr,nrefim]=size(all_refim___); 
xnodesr_ = MDA_read_r8('./dir_test33b/xnodesr_.mda');
figure(1);clf;
for ni=1:nrefim;
subplot(4,3,ni+0); 
imagesc_p(ngridr,xnodesr_,ones(ngridr,1)*ngridc,ngridr*ngridc,reshape(real(all_refim___(:,:,ni)),ngridr*ngridc,1),[-1,+1],colormap_beach()); set(gca,'XTick',[],'YTick',[]);
axis image;
subplot(4,3,ni+3); 
imagesc_p(ngridr,xnodesr_,ones(ngridr,1)*ngridc,ngridr*ngridc,reshape(imag(all_refim___(:,:,ni)),ngridr*ngridc,1),[-1,+1],colormap_beach()); set(gca,'XTick',[],'YTick',[]);
axis image;
end;%for ni=1:nrefim;
figbig;
for ni=1:nrefim;
subplot(4,3,ni+6); 
imagesc_p(ngridr,xnodesr_,ones(ngridr,1)*ngridc,ngridr*ngridc,reshape(real(ctfw2___(:,:,ni)),ngridr*ngridc,1),[-1,+1],colormap_beach()); set(gca,'XTick',[],'YTick',[]);
axis image;
subplot(4,3,ni+9); 
imagesc_p(ngridr,xnodesr_,ones(ngridr,1)*ngridc,ngridr*ngridc,reshape(imag(ctfw2___(:,:,ni)),ngridr*ngridc,1),[-1,+1],colormap_beach()); set(gca,'XTick',[],'YTick',[]);
axis image;
end;%for ni=1:nrefim;
figbig;

dirname = './dir_Mode2_k20_dRTTRx7p100m10_s1e100_nF_cN_tF_hF_rF_M3_S0_s0_MS_r1_m3G';
CTF_k_p__ = MDA_read_c16(sprintf('%s/CTF_k_p__.mda',dirname)); 
Y_slice__ = MDA_read_c16(sprintf('%s/Y_slice__.mda',dirname)); 
[n_A,n_M]=size(Y_slice__); 
grid_k_p_r_ = MDA_read_r8(sprintf('%s/grid_k_p_r_.mda',dirname));
n_k_p_r_max = length(grid_k_p_r_);
n_w_ = MDA_read_i4(sprintf('%s/n_w_.mda',dirname));
figure(2);clf;
for ni=1:n_M;
subplot(4,3,ni+0); 
imagesc_p(n_k_p_r_max,grid_k_p_r_,n_w_,n_A,real(Y_slice__(:,ni)),[-1,+1],colormap_beach()); set(gca,'XTick',[],'YTick',[]);
axis image;
subplot(4,3,ni+3); 
imagesc_p(n_k_p_r_max,grid_k_p_r_,n_w_,n_A,imag(Y_slice__(:,ni)),[-1,+1],colormap_beach()); set(gca,'XTick',[],'YTick',[]);
axis image;
end;%for ni=1:n_M;
figbig;
for ni=1:n_M;
subplot(4,3,ni+6); 
imagesc_p(n_k_p_r_max,grid_k_p_r_,n_w_,n_A,real(CTF_k_p__(:,ni)),[-1,+1],colormap_beach()); set(gca,'XTick',[],'YTick',[]);
axis image;
subplot(4,3,ni+9); 
imagesc_p(n_k_p_r_max,grid_k_p_r_,n_w_,n_A,imag(CTF_k_p__(:,ni)),[-1,+1],colormap_beach()); set(gca,'XTick',[],'YTick',[]);
axis image;
end;%for ni=1:n_M;
figbig;

dirname_relion = '/data/rangan/dir_cryoem/dir_cryoEM-SCDA/marina_work/trpv1_nosym/try40_relion'; n_ = -1; infix = '';
isosurface_ver8(n_,dirname_relion,infix);

dirname_adi = './dir_Mode2_k20_dRTTRx7p100m10_s1e100_nF_cN_tF_hF_rF_M4096_S0_s0_MS_r1_m3G'; n_ = 20; infix = 'tru';
isosurface_ver8(n_,dirname_adi,infix);

dirname_adi = './dir_Mode2_k20_dRTTRx7p100m10_s1e100_nF_cN_tF_hF_rF_M512_S0_s0_MS_r1_m3G'; n_ = 20; infix = 'tru';
isosurface_ver8(n_,dirname_adi,infix);

dirname_adi = './dir_Mode2_k20_dRTTRx7p100m10_s1e100_nF_cN_tF_hTa100b0c0d0_rTl2i7_M512_S0_s0_MS_r1_m3G'; n_ = 20; infix = 'tru';
isosurface_ver8(n_,dirname_adi,infix);

dirname_adi = './dir_Mode2_k20_dRTTRx7p100m10_s1e100_nF_cN_tF_hTa100b0c0d0_rF_M512_S0_s0_MS_r1_m3G'; n_ = 20; infix = 'tru';
isosurface_ver8(n_,dirname_adi,infix);

dirname_adi = './dir_Mode2_k20_dRTTRx7p100m10_s1e100_nF_cN_tF_hF_rTl2i7_M512_S0_s0_MS_r1_m3G'; n_ = 20; infix = 'tru';
isosurface_ver8(n_,dirname_adi,infix);
residual_scatter_ver0(n_,dirname_adi,infix);

dirname_adi = './dir_Mode2_k20_dRTTRx7p100m10_s1e100_nF_cN_tF_hF_rTl2i7_M4096_S0_s0_MS_r1_m3G'; n_ = 20; infix = 'tru';
isosurface_ver8(n_,dirname_adi,infix);
residual_scatter_ver0(n_,dirname_adi,infix);

dirname_adi = './dir_Mode2_k20_dRTTRx3p50m25_s1e100_nF_cN_tF_hF_rF_M512_S0_s0_MS_r1_m3G'; n_ = 2:5;
isosurface_ver8(n_,dirname_adi,'tru');
isosurface_ver8(n_,dirname_adi,'est');
alphadist_ver3(n_,dirname_adi,'est',0.25);

dirname_adi = './dir_Mode2_k20_dRTTRx3p50m25_s1e100_nF_cN_tF_hF_rF_M512_S0_s0_MS_r1_m3G'; n_ = 5:5:20;
isosurface_ver8(n_,dirname_adi,'tru');
isosurface_ver8(n_,dirname_adi,'est');
alphadist_ver3(n_,dirname_adi,'est',0.25);

dirname_adi = './dir_Mode2_k25_dRTTRx3p50m25_s1e100_nF_cN_tF_hF_rF_M1024_S0_s0_MS_r1_m4G'; n_ = 5:3:20;
isosurface_ver8(n_,dirname_adi,'tru');
isosurface_ver8(n_,dirname_adi,'est');
alphadist_ver3(n_,dirname_adi,'est',0.25);

dirname_adi = './dir_Mode2_k20_dRTTRx3p50m25_s1e100_nF_cN_tF_hF_rF_M1536_S0_s0_MS_r1_m4G'; n_ = 5:3:20;
isosurface_ver8(n_,dirname_adi,'tru');
isosurface_ver8(n_,dirname_adi,'est');
alphadist_ver3(n_,dirname_adi,'est',0.25);

dirname_adi = './dir_Mode2_k24_dRTTRx3p50m25_s1e100_nF_cN_tF_hF_rF_M1536_S0_s0_MS_r1_m4G'; n_ = 4:4:24;
%isosurface_ver8(20,dirname_adi,'est');
isosurface_ver8(n_,dirname_adi,'est');
isosurface_ver8(n_,dirname_adi,'tru');
alphadist_ver3(n_,dirname_adi,'est',0.25);
%n_k = n_(end)
%tmp_a_ = MDA_read_c16(sprintf('%s/Y_tru_%d_.mda',dirname_adi,n_k)); 
%tmp_b_ = MDA_read_c16(sprintf('%s/Y_est_%d_.mda',dirname_adi,n_k)); 
%k_ = MDA_read_r8(sprintf('%s/grid_k_p_r_.mda',dirname_adi)); 

dirname_adi = './dir_Mode2_k20_dF_s2_nF_cN_tF_hF_rF_M1536_S0_s0_MS_r1_m4G/'; n_ = 5:3:20;
%isosurface_ver8(20,dirname_adi,'est');
isosurface_ver8(n_,dirname_adi,'est');
isosurface_ver8(n_,dirname_adi,'tru');
alphadist_ver3(n_,dirname_adi,'est',0.25);
n_k = n_(end);
n_l_ = MDA_read_i4(sprintf('%s/n_Y_l__.mda',dirname_adi)); 
a_ = MDA_read_c16(sprintf('%s/Y_tru_%d_.mda',dirname_adi,n_k)); 
b_ = MDA_read_c16(sprintf('%s/Y_est_%d_.mda',dirname_adi,n_k)); 
k_ = MDA_read_r8(sprintf('%s/grid_k_p_r_.mda',dirname_adi)); 
X_best = test_spharm_13(n_k,k_,n_l_,a_,b_);

n_lm_ = (n_l_+1).^2;
X_best = test_spharm_13(n_k-1,k_(2:end),n_l_(2:end),a_(n_lm_(1)+1:end),b_(n_lm_(1)+1:end));

dirname_adi = './dir_Mode2_k24_dRTTRx3p50m25_s1e100_nF_cN_tF_hF_rF_M1536_S0_s0_SMa0_r1_m4G'; n_ = 4:4:24;
isosurface_ver8(n_,dirname_adi,'est');
isosurface_ver8(n_,dirname_adi,'tru');
alphadist_ver3(n_,dirname_adi,'est',0.25);

dirname_adi = './dir_Mode2_k24_dRTTRx5p50m25_s1e100_nF_cN_tF_hF_rF_M1536_S0_s0_MS_r1_m4G'; n_ = 4:4:24;
alphadist_ver3(n_,dirname_adi,'est',0.25);
isosurface_ver8(n_,dirname_adi,'est');
%isosurface_ver8(n_,dirname_adi,'tru');

dirname_adi = './dir_Mode2_k20_dF_s2_nF_cN_tF_hF_rF_M128_S0_s0_SMa0_r1_m4G'; n_ = 5:3:20;
alphadist_ver3(n_,dirname_adi,'est',0.25);
isosurface_ver8(n_,dirname_adi,'est');
isosurface_ver8(n_,dirname_adi,'tru');

dirname_adi = './dir_Mode2_k20_dF_s2_nF_cN_tF_hF_rF_M16_S0_s0_SMa0_r1_m4G'; n_ = 3:3:12;
alphadist_ver3(n_,dirname_adi,'est',0.25);
isosurface_ver8(n_,dirname_adi,'est');
isosurface_ver8(n_,dirname_adi,'tru');

dirname_adi = './dir_Mode2_k64_dRTTRx7p50m25_s1e100_nF_cN_tF_hF_rF_M1536_S0_s0_MS_r1_m4G/'; n_ = 4:12:64;
alphadist_ver3(n_,dirname_adi,'est',0.25);
isosurface_ver8(n_,dirname_adi,'est');
isosurface_ver8(n_,dirname_adi,'tru');


./test_ver18.out 2 48 1 0 0 0 0 1 7 1

dirname_adi = './dir_Mode2_k48_dRTTRx7p50m25_s1e100_nF_cN_tF_hF_rF_M4096_S0_s0_MS_r1_m4G'; n_ = 48;
isosurface_ver8(n_,dirname_adi,'tru');
isosurface_ver8(n_,dirname_adi,'est'); 
dirname_adi = './dir_Mode2_k48_dRTTRx15p50m25_s1e100_nF_cN_tF_hF_rF_M4096_S0_s0_MS_r1_m4G'; n_ = 48;
isosurface_ver8(n_,dirname_adi,'est'); 
dirname_adi = './dir_Mode2_k48_dRTTRx15p50m25_s1e100_nF_cN_tF_hF_rF_M4096_S0_s0_SMa0_r1_m4G'; n_ = 48; n_ = 39;
isosurface_ver8(n_,dirname_adi,'est'); 
dirname_adi = './dir_Mode2_k48_dF_s2_nF_cN_tF_hF_rF_M4096_S0_s0_SMa0_r1_m4G'; n_ = 48;
isosurface_ver8(n_,dirname_adi,'est'); 
dirname_adi = './dir_Mode2_k48_dF_s2_nF_cN_tF_hF_rF_M4096_S0_s0_MS_r1_m4G'; n_ = 48;
isosurface_ver8(n_,dirname_adi,'est'); 
