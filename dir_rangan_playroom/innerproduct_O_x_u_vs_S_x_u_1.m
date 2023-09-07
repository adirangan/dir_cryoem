function ...
[ ...
 parameter ...
,OCS_x_u_xxw___ ...
] = ...
innerproduct_O_x_u_vs_S_x_u_1( ...
 parameter ...
,n_O_x_0 ...
,n_O_x_1 ...
,O_x_u_xx__ ...
,n_S_x_0 ...
,n_S_x_1 ...
,S_x_u_xx__ ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,S_k_p_wk_ ...
,weight_2d_wk_ ...
,CTF_k_p_r_k_ ...
);
%%%%%%%%;
% calculates the convolution (i.e., innerproduct) between ;
% micrograph O_x_u_xx__ and template S_x_u_xx__. ;
% The array OCS_x_u_xxw___(:,:,1+nw) ;
% corresponds to the conv2 between O_x_u_xx__ ;
% and the template S_k_p_wk_ ;
% rotated by tmp_gamma_z_(1+nw) = (2*pi*nw)/max(1,n_w_max). ;
%%%%%%%%;
str_thisfunction = 'innerproduct_O_x_u_vs_S_x_u_1';

if nargin<1;
flag_verbose=1;
if (flag_verbose); disp(sprintf(' %% testing %s',str_thisfunction)); end;
%flag_micrograph_scale = 15; %<-- typical (960-x-682). ;
flag_micrograph_scale = 15;
nf=0;
%%%%;
k_p_r_max = 2*48/(2*pi); k_eq_d = 1.0/(2*pi);
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
);
%%%%;
l_max_upb = round(2*pi*k_p_r_max); l_max_max = l_max_upb;
n_w_0in_ = 2*(l_max_max+1)*ones(n_k_p_r,1);
template_k_eq_d = -1;
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum = cumsum([0;n_w_(:)]);
%%%%;
CTF_k_p_r_k_ = 1+sin(2*pi*3*k_p_r_/k_p_r_max);
%%%%;
sigma_O = 1/256.0;
n_source_O = 1024;
rng(0);
source_O_2s__ = 2*rand(2,n_source_O)-1;
n_O_x_0 = round(960/15*flag_micrograph_scale);
n_O_x_1 = round(682/15*flag_micrograph_scale);
O_x_0_ = linspace(-1,+1,n_O_x_0);
O_x_1_ = linspace(-1,+1,n_O_x_1);
[O_x_0__,O_x_1__] = ndgrid(O_x_0_,O_x_1_);
O_x_u_xx__ = zeros(n_O_x_0,n_O_x_1);
for nsource_O=0:n_source_O-1;
O_x_r2__ = (O_x_0__-source_O_2s__(1+0,1+nsource_O)).^2 + (O_x_1__-source_O_2s__(1+1,1+nsource_O)).^2;
O_x_u_xx__ = O_x_u_xx__ + exp(-O_x_r2__/2/sigma_O^2);
end;%for nsource_O=0:n_source_O-1;
%%%%;
sigma_S = 0.0625;
dS00 = -0.162; dS01 = +0.170;
dS10 = +0.185; dS11 = -0.170;
%f_S = @(x0,x1) ...
%  1/(2*pi)/sigma_S^2 * exp(-((x0-dS00).^2+(x1-dS01).^2)/(2*sigma_S^2)) ...
%+ 1/(2*pi)/sigma_S^2 * exp(-((x0-dS10).^2+(x1-dS11).^2)/(2*sigma_S^2)) ...
%;
f_S = @(x0,x1) ...
  1/(2*pi)/sigma_S^2 * exp(-max(abs(x0-dS00),abs(x1-dS01)).^2/(2*sigma_S^2)) ...
+ 1/(2*pi)/sigma_S^2 * exp(-(abs(x0-dS10)+abs(x1-dS11)).^2/(2*sigma_S^2)) ...
;
n_S_x_0 = 64;
n_S_x_1 = 64;
S_x_0_ = linspace(-1,+1,n_S_x_0);
S_x_1_ = linspace(-1,+1,n_S_x_1);
[S_x_0__,S_x_1__] = ndgrid(S_x_0_,S_x_1_);
S_x_u_xx__ = f_S(S_x_0__,S_x_1__);
%%%%;
diameter_x_c = 2.0; dx_0 = diameter_x_c/max(1,n_S_x_0); dx_1 = diameter_x_c/max(1,n_S_x_1);
S_k_p_wk_ = interp_x_c_to_k_p_xxnufft(n_S_x_0,diameter_x_c,n_S_x_1,diameter_x_c,S_x_u_xx__,n_k_p_r,k_p_r_,n_w_)*sqrt(n_S_x_0*n_S_x_1)*dx_0*dx_1;
CTF_S_k_p_wk_ = zeros(n_w_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
CTF_k_p_r = CTF_k_p_r_k_(1+nk_p_r);
for nw=0:n_w-1;
CTF_S_k_p_wk_(1+na) = S_k_p_wk_(1+na)*CTF_k_p_r;
na=na+1;
end;%for nw=0:n_w-1;
end;%for nk_p_r=0:n_k_p_r-1;
CTF_S_x_c_xx_ = interp_k_p_to_x_c_xxnufft(n_S_x_0,diameter_x_c,n_S_x_1,diameter_x_c,n_k_p_r,k_p_r_,n_w_,CTF_S_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_S_x_0*n_S_x_1) * n_w_sum;
CTF_S_x_c_xx__ = real(reshape(CTF_S_x_c_xx_,[n_S_x_0,n_S_x_1]));
if (flag_verbose); disp(sprintf(' %% S_x_u_xx__ vs conv2(PSF,S_x_c_xx__): %0.16f',fnorm(S_x_u_xx__-CTF_S_x_c_xx__)/fnorm(S_x_u_xx__))); end;
%%%%;
flag_check=0;
if flag_check;
if (flag_verbose); disp(sprintf(' %% checking disk_k_p_scatter_from_tensor_interpolate_n_4')); end;
tmp_n_k = 128;
tmp_k_0_ = linspace(-k_p_r_max,+k_p_r_max,tmp_n_k);
tmp_k_1_ = linspace(-k_p_r_max,+k_p_r_max,tmp_n_k);
[tmp_k_0__,tmp_k_1__] = ndgrid(tmp_k_0_,tmp_k_1_);
gamma_z_scatter__ = atan2(tmp_k_1__,tmp_k_0__);
k_p_rad_scatter__ = sqrt(tmp_k_0__.^2+tmp_k_1__.^2);
n_scatter = tmp_n_k^2; n_order = 5;
tmp_t = tic();
[ ...
 scatter_from_tensor_swk__ ...
] = ...
disk_k_p_scatter_from_tensor_interpolate_n_4( ...
 n_order ...
,n_w_max ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_scatter ...
,gamma_z_scatter__(:) ...
,k_p_rad_scatter__(:) ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% disk_k_p_scatter_from_tensor_interpolate_n_4: %0.2fs',tmp_t)); end;
if (flag_verbose); disp(sprintf(' %% scatter_from_tensor_swk__: numel %d, nnz %d --> %0.4f GB',numel(scatter_from_tensor_swk__),nnz(scatter_from_tensor_swk__),nnz(scatter_from_tensor_swk__)*3*8/1e9)); end;
S_k_p_kk__ = reshape(scatter_from_tensor_swk__*S_k_p_wk_,[tmp_n_k,tmp_n_k]);
tmp_nw = floor(n_w_max/6); tmp_gamma_z = (2*pi*tmp_nw)/n_w_max;
T_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,+tmp_gamma_z);
tmp_t = tic();
[tmp_ij_row_,tmp_ij_col_,tmp_val_] = find(scatter_from_tensor_swk__);
tmp_index_row_ = tmp_ij_row_-1;
tmp_index_col_ = tmp_ij_col_-1;
tmp_index_nw_rot_ = mod(tmp_index_col_,n_w_max);
tmp_index_nk_rot_ = floor(tmp_index_col_/n_w_max);
tmp_index_nw_rot_ = periodize(tmp_index_nw_rot_-tmp_nw,0,n_w_max);
tmp_index_col_rot_ = tmp_index_nw_rot_ + n_w_max*tmp_index_nk_rot_;
scatter_from_tensor_rot_swk__ = sparse(1+tmp_index_row_,1+tmp_index_col_rot_,tmp_val_,n_scatter,n_w_sum);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% scatter_from_tensor_rot_swk__: %0.2fs',tmp_t)); end;
T_k_p_kk__ = reshape(scatter_from_tensor_rot_swk__*S_k_p_wk_,[tmp_n_k,tmp_n_k]);
figure(1+nf);nf=nf+1;clf;figbig; 
p_row=2; p_col=4; np=0;
tmp_real_clim__ = prctile(real(S_k_p_wk_), [5,95]);
tmp_imag_clim__ = prctile(imag(S_k_p_wk_), [5,95]);
tmp_c_use__ = colormap_80s;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wk_),tmp_real_clim__,tmp_c_use__); axis image; axisnotick;
title('real(S_k_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(tmp_n_k,tmp_k_0_,tmp_n_k,tmp_k_1_,real(S_k_p_kk__),tmp_real_clim__,tmp_c_use__); axis image; axisnotick;
title('real(S_k_p_kk__)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(S_k_p_wk_),tmp_imag_clim__,tmp_c_use__); axis image; axisnotick;
title('imag(S_k_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(tmp_n_k,tmp_k_0_,tmp_n_k,tmp_k_1_,imag(S_k_p_kk__),tmp_imag_clim__,tmp_c_use__); axis image; axisnotick;
title('imag(S_k_p_kk__)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_wk_),tmp_real_clim__,tmp_c_use__); axis image; axisnotick;
title('real(T_k_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(tmp_n_k,tmp_k_0_,tmp_n_k,tmp_k_1_,real(T_k_p_kk__),tmp_real_clim__,tmp_c_use__); axis image; axisnotick;
title('real(T_k_p_kk__)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(T_k_p_wk_),tmp_imag_clim__,tmp_c_use__); axis image; axisnotick;
title('imag(T_k_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(tmp_n_k,tmp_k_0_,tmp_n_k,tmp_k_1_,imag(T_k_p_kk__),tmp_imag_clim__,tmp_c_use__); axis image; axisnotick;
title('imag(T_k_p_kk__)','Interpreter','none');
if (flag_verbose); disp(sprintf(' %% finished check')); end;
error('stopping early');
end;%if flag_check;
%%%%;
tmp_t = tic();
[~,OCS_x_u_xxw___] = ...
innerproduct_O_x_u_vs_S_x_u_1( ...
 [] ...
,n_O_x_0 ...
,n_O_x_1 ...
,O_x_u_xx__ ...
,n_S_x_0 ...
,n_S_x_1 ...
,[] ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,S_k_p_wk_ ...
,weight_2d_wk_ ...
,CTF_k_p_r_k_ ...
);
tmp_t = toc(tmp_t);
if (flag_verbose); disp(sprintf(' %% innerproduct_O_x_u_vs_S_x_u_1: %0.2fs',tmp_t)); end;
%%%%;
tmp_t = tic();
for nw=0:n_w_max-1;
tmp_gamma_z = (2*pi*nw)/max(1,n_w_max);
tmp_CTF_S_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,CTF_S_k_p_wk_,+tmp_gamma_z);
tmp_CTF_S_x_c_xx_ = interp_k_p_to_x_c_xxnufft(n_S_x_0,diameter_x_c,n_S_x_1,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_CTF_S_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_S_x_0*n_S_x_1) * n_w_sum;
tmp_CTF_S_x_c_xx__ = real(reshape(tmp_CTF_S_x_c_xx_,[n_S_x_0,n_S_x_1]));
OCS_x_c_xxw___(:,:,1+nw) = conv2(O_x_u_xx__,tmp_CTF_S_x_c_xx__,'same');
end;%for nw=0:n_w_max-1;
tmp_t = toc(tmp_t);
if (flag_verbose); disp(sprintf(' %% conv2: %0.2fs',tmp_t)); end;
figure(1+nf);nf=nf+1;clf;figbig; fig80s; p_row=2;p_col=4;np=0;
subplot(p_row,p_col,1+np);np=np+1; imagesc(O_x_u_xx__); title('O'); axisnotick; colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(S_x_u_xx__); title('S'); axisnotick; colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(CTF_S_x_c_xx__); title('conv2(PSF,S)'); axisnotick; colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(OCS_x_u_xxw___(:,:,1+0)); title('ifft2 nw 0'); axisnotick; colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(OCS_x_c_xxw___(:,:,1+0)); title('conv2 nw 0'); axisnotick; colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(max(OCS_x_u_xxw___,[],3)); title('max(OCS_x_u_xxw___,[],3)','Interpreter','none'); axisnotick; colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(max(OCS_x_c_xxw___,[],3)); title('max(OCS_x_c_xxw___,[],3)','Interpreter','none'); axisnotick; colorbar;
subplot(p_row,p_col,1+np);np=np+1; 
hlim_ = [ ...
 min(min(OCS_x_u_xxw___,[],'all'),min(OCS_x_c_xxw___,[],'all')) ...
,max(max(OCS_x_u_xxw___,[],'all'),max(OCS_x_c_xxw___,[],'all')) ...
];
n_h = 128;
h2__ = hist2d_0(OCS_x_u_xxw___,OCS_x_c_xxw___,n_h,n_h,hlim_,hlim_);
imagesc(log2(1+h2__)); axis image; axisnotick; colorbar; set(gca,'ydir','normal');
title('scatterplot');
disp(sprintf(' %% OCS_x_u_xxw___ vs OCS_x_c_xxw___: %0.16f',fnorm(OCS_x_u_xxw___ - OCS_x_c_xxw___)/max(1e-12,fnorm(OCS_x_u_xxw___))));
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_O_x_0=[]; end; na=na+1;
if (nargin<1+na); n_O_x_1=[]; end; na=na+1;
if (nargin<1+na); O_x_u_xx__=[]; end; na=na+1;
if (nargin<1+na); n_S_x_0=[]; end; na=na+1;
if (nargin<1+na); n_S_x_1=[]; end; na=na+1;
if (nargin<1+na); S_x_u_xx__=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); S_k_p_wk_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_wk_=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_r_k_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%;
O_x_u_xx__ = reshape(O_x_u_xx__,[n_O_x_0,n_O_x_1]);
if ~isempty(S_x_u_xx__);
S_x_u_xx__ = reshape(S_x_u_xx__,[n_S_x_0,n_S_x_1]);
end;%if ~isempty(S_x_u_xx__);
%%%%;
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum = cumsum([0;n_w_(:)]);
%%%%;
diameter_x_c = 2.0; dx_0 = diameter_x_c/max(1,n_S_x_0); dx_1 = diameter_x_c/max(1,n_S_x_1);
if  isempty(S_k_p_wk_);
S_k_p_wk_ = interp_x_c_to_k_p_xxnufft(n_S_x_0,diameter_x_c,n_S_x_1,diameter_x_c,S_x_u_xx__,n_k_p_r,k_p_r_,n_w_)*sqrt(n_S_x_0*n_S_x_1)*dx_0*dx_1;
end;%if  isempty(S_k_p_wk_);
CTF_S_k_p_wk_ = zeros(n_w_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
CTF_k_p_r = CTF_k_p_r_k_(1+nk_p_r);
for nw=0:n_w-1;
CTF_S_k_p_wk_(1+na) = S_k_p_wk_(1+na)*CTF_k_p_r;
na=na+1;
end;%for nw=0:n_w-1;
end;%for nk_p_r=0:n_k_p_r-1;
CTF_S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,CTF_S_k_p_wk_);
CTF_S_x_c_xx_ = interp_k_p_to_x_c_xxnufft(n_S_x_0,diameter_x_c,n_S_x_1,diameter_x_c,n_k_p_r,k_p_r_,n_w_,CTF_S_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_S_x_0*n_S_x_1) * n_w_sum;
CTF_S_x_c_xx__ = real(reshape(CTF_S_x_c_xx_,[n_S_x_0,n_S_x_1]));
%%%%;

%%%%%%%%;
n_N_x_0 = max(n_O_x_0,n_S_x_0); n_N_x_0 = n_N_x_0 + mod(n_N_x_0,2);
n_N_x_1 = max(n_O_x_1,n_S_x_1); n_N_x_1 = n_N_x_1 + mod(n_N_x_1,2);
N_x_u_xx__ = zeros(n_N_x_0,n_N_x_1);
ij_O_0_ =  n_N_x_0/2 - floor(n_O_x_0/2) + [1:n_O_x_0] ;
ij_O_1_ =  n_N_x_1/2 - floor(n_O_x_1/2) + [1:n_O_x_1] ;
N_x_u_xx__(ij_O_0_,ij_O_1_) = O_x_u_xx__;
N_k_u_kk__ = fft2(fftshift(N_x_u_xx__));
%N_k_u_kk__ = fft2(N_x_u_xx__);
%%%%%%%%;
n_R_x_0 = n_N_x_0;
n_R_x_1 = n_N_x_1;
ij_S_0_ =  n_N_x_0/2 - floor(n_S_x_0/2) + [1:n_S_x_0] ;
ij_S_1_ =  n_N_x_1/2 - floor(n_S_x_1/2) + [1:n_S_x_1] ;
R_x_u_xxw___ = zeros(n_R_x_0,n_R_x_1,n_w_max);
tmp_CTF_S_k_p_wk_ = CTF_S_k_p_wk_; tmp_CTF_S_x_c_xx__ = CTF_S_x_c_xx__;
for nw=0:n_w_max-1;
tmp_gamma_z = (2*pi*nw)/max(1,n_w_max);
%tmp_CTF_S_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,CTF_S_k_p_wk_,+tmp_gamma_z);
tmp_CTF_S_k_q_wk_ = rotate_q_to_q(n_k_p_r,n_w_,n_w_sum,CTF_S_k_q_wk_,+tmp_gamma_z);
tmp_CTF_S_k_p_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,tmp_CTF_S_k_q_wk_);
tmp_CTF_S_x_c_xx_ = interp_k_p_to_x_c_xxnufft(n_S_x_0,diameter_x_c,n_S_x_1,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_CTF_S_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_S_x_0*n_S_x_1) * n_w_sum;
tmp_CTF_S_x_c_xx__ = real(reshape(tmp_CTF_S_x_c_xx_,[n_S_x_0,n_S_x_1]));
R_x_u_xxw___(ij_S_0_,ij_S_1_,1+nw) = tmp_CTF_S_x_c_xx__;
end;%for nw=0:n_w_max-1;
R_k_u_kkw___ = fft2(fftshift(fftshift(R_x_u_xxw___,1),2));
%R_k_u_kkw___ = fft2(R_x_u_xxw___);
NR_k_u_kkw___ = bsxfun(@times,N_k_u_kk__,R_k_u_kkw___);
NR_x_u_xxw___ = real(fftshift(fftshift(ifft2(NR_k_u_kkw___),1),2));
%NR_x_u_xxw___ = real(ifft2(NR_k_u_kkw___));
OCS_x_u_xxw___ = NR_x_u_xxw___(ij_O_0_,ij_O_1_,:);
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

%{
%%%%%%%%;
% example of 'upsampling' via padding fft. ;
%%%%%%%%;
tmp_n_x = 32;
tmp_x0_ = linspace(-1,+1,tmp_n_x);
tmp_x1_ = linspace(-1,+1,tmp_n_x); 
[tmp_x0__,tmp_x1__] = ndgrid(tmp_x0_,tmp_x1_);
tmp_xr__ = sqrt(tmp_x0__.^2 + tmp_x1__.^2);
tmp_A_xx__ = sin(2*pi*tmp_x0__).*cos(2*pi*tmp_x1__).*exp(-tmp_xr__.^2/(2*0.125.^2));
tmp_B_kk__ = fftshift(fft2(fftshift(tmp_A_xx__)));
tmp_C_kk__ = zeros(3*tmp_n_x,3*tmp_n_x);
tmp_C_kk__(tmp_n_x+[1:tmp_n_x],tmp_n_x+[1:tmp_n_x]) = tmp_B_kk__;
tmp_D_xx__ = fftshift(ifft2(fftshift(tmp_C_kk__)));
%}
