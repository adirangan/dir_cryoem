%function test_transforms_dr();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% Generating initial S_x_c_ ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
n_r = 32;
n_x = 2*n_r; n_y = 2*n_r;
max_x_c = 1.0; max_y_c = 1.0;
grid_x_c_ = linspace(-max_x_c/2,+max_x_c/2,n_x);
x_c_ = linspace(-max_x_c/2,+max_x_c/2,n_x);
y_c_ = linspace(-max_y_c/2,+max_y_c/2,n_y);
X_c_ = ones(n_y,1)*x_c_;
Y_c_ = transpose(y_c_)*ones(1,n_x);
R_ = sqrt(X_c_.^2 + Y_c_.^2); W_ = atan2(Y_c_,X_c_);
P_ = 0.125*(1 + 0.5*cos(5*W_) + 0.25*sin(3*W_) + 0.125*cos(1*W_));
G_x_c_ = 0.25*exp(-((X_c_-0.125).^2 + (Y_c_-0.125).^2)/2/0.125^2); 
S_x_c_ = exp(-R_./P_) + G_x_c_;
S_x_c_ = recenter2(S_x_c_);
tmp_kx_ij = n_x/2 + (-n_r/2:n_r/2-1); tmp_ky_ij = n_x/2 + (-n_r/2:n_r/2-1);
S_k_c_ = fft2(S_x_c_); S_k_c_ = recenter2(S_k_c_); S_k_c_ = S_k_c_(1+tmp_kx_ij,1+tmp_ky_ij); S_k_c_ = decenter2(S_k_c_);
max_k_c = 1.0*n_r/max_x_c;
k_c_ = linspace(-max_k_c/2,+max_k_c/2,n_r);
grid_k_p_ = linspace(1,max_k_c/2,n_r);
n_w_ = max(6,ceil(2*pi*grid_k_p_)); n_A = sum(n_w_); n_w_max = max(n_w_);
S_k_p_ = interp_c_to_p(n_r,max_k_c,n_r,max_k_c,S_k_c_,n_r,grid_k_p_,n_w_,n_A);
S_k_q_ = interp_p_to_q(n_r,n_w_,n_A,S_k_p_);
C_x_c = innerproduct_c(n_x,x_c_,n_y,y_c_,S_x_c_,S_x_c_);
disp(sprintf(' %% <S_x_c_,S_x_c_> = (%0.4f,%0.4f)',real(C_x_c),imag(C_x_c)));
C_k_c = innerproduct_c(n_r,k_c_,n_r,k_c_,S_k_c_,S_k_c_);
disp(sprintf(' %% <S_k_c_,S_k_c_> = (%0.4f,%0.4f)',real(C_k_c),imag(C_k_c)));
C_k_p = innerproduct_p(n_r,grid_k_p_,n_w_,n_A,S_k_p_,S_k_p_);
disp(sprintf(' %% <S_k_p_,S_k_p_> = (%0.4f,%0.4f)',real(C_k_p),imag(C_k_p)));
C_k_q = innerproduct_p(n_r,grid_k_p_,n_w_,n_A,S_k_q_,S_k_q_);
disp(sprintf(' %% <S_k_q_,S_k_q_> = (%0.4f,%0.4f)',real(C_k_q),imag(C_k_q)));
flag_disp=0;
if flag_disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
figure(1); pcols = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
subplot(2,pcols,1+0*pcols); 
imagesc_c(n_x,x_c_,n_y,y_c_,real(decenter2(S_x_c_)),[],colormap('jet'));
xlim(max_x_c/2*[-1,1]);ylim(max_y_c/2*[-1,1]);
title('real(S_x_c_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[],'Ydir','normal'); axis square;
subplot(2,pcols,1+1*pcols); 
imagesc_c(n_x,x_c_,n_y,y_c_,imag(decenter2(S_x_c_)),[],colormap('jet'));
xlim(max_x_c/2*[-1,1]);ylim(max_y_c/2*[-1,1]);
title('imag(S_x_c_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[],'Ydir','normal'); axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
subplot(2,pcols,2+0*pcols); 
imagesc_c(n_r,k_c_,n_r,k_c_,real(decenter2(S_k_c_)),[],colormap('jet'));
xlim(max_k_c/2*[-1,1]);ylim(max_k_c/2*[-1,1]);
title('real(S_k_c_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[],'Ydir','normal'); axis square;
subplot(2,pcols,2+1*pcols); 
imagesc_c(n_r,k_c_,n_r,k_c_,imag(decenter2(S_k_c_)),[],colormap('jet'));
xlim(max_k_c/2*[-1,1]);ylim(max_k_c/2*[-1,1]);
title('imag(S_k_c_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[],'Ydir','normal'); axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
subplot(2,pcols,3+0*pcols);
imagesc_p(n_r,grid_k_p_,n_w_,n_A,real(S_k_p_),[],colormap('jet'));
xlim(max_k_c/2*[-1,1]);ylim(max_k_c/2*[-1,1]);
title('real(S_k_p_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[]); axis square;
subplot(2,pcols,3+1*pcols);
imagesc_p(n_r,grid_k_p_,n_w_,n_A,imag(S_k_p_),[],colormap('jet'));
xlim(max_k_c/2*[-1,1]);ylim(max_k_c/2*[-1,1]);
title('imag(S_k_p_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[]); axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
subplot(2,pcols,4+0*pcols);
imagesc_q(n_r,grid_k_p_,n_w_,n_A,real(S_k_q_),[],colormap('jet'));
xlim(max_k_c/2*[0,1]);ylim(n_w_max/2*[-1,1]);
title('real(S_k_q_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[]); axis square;
subplot(2,pcols,4+1*pcols);
imagesc_q(n_r,grid_k_p_,n_w_,n_A,imag(S_k_q_),[],colormap('jet'));
xlim(max_k_c/2*[0,1]);ylim(n_w_max/2*[-1,1]);
title('imag(S_k_q_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[]); axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
end;%if flag_disp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% Generating initial M_x_c_ by rotating S_x_c_ by pi/3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
gamma_z_tru = 1.0*pi/3.0;
M_x_c_ = rotate_c_to_c(n_x,max_x_c,n_y,max_y_c,S_x_c_,gamma_z_tru);
M_k_c_ = fft2(M_x_c_); M_k_c_ = recenter2(M_k_c_); M_k_c_ = M_k_c_(1+tmp_kx_ij,1+tmp_ky_ij); M_k_c_ = decenter2(M_k_c_);
M_k_p_ = interp_c_to_p(n_r,max_k_c,n_r,max_k_c,M_k_c_,n_r,grid_k_p_,n_w_,n_A);
N_k_p_ = rotate_p_to_p(n_r,n_w_,n_A,S_k_p_,gamma_z_tru);
flag_disp=0;
if flag_disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
figure(2); pcols=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
subplot(2,pcols,1+0*pcols); 
imagesc_c(n_x,x_c_,n_y,y_c_,real(decenter2(S_x_c_)),[],colormap('jet'));
xlim(max_x_c/2*[-1,1]);ylim(max_y_c/2*[-1,1]);
title('real(S_x_c_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[],'Ydir','normal'); axis square;
subplot(2,pcols,1+1*pcols); 
imagesc_c(n_x,x_c_,n_y,y_c_,real(decenter2(M_x_c_)),[],colormap('jet'));
xlim(max_x_c/2*[-1,1]);ylim(max_y_c/2*[-1,1]);
title('real(M_x_c_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[],'Ydir','normal'); axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
subplot(2,pcols,2+0*pcols);
imagesc_p(n_r,grid_k_p_,n_w_,n_A,real(M_k_p_),[],colormap('jet'));
xlim(max_k_c/2*[-1,1]);ylim(max_k_c/2*[-1,1]);
title('real(M_k_p_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[]); axis square;
subplot(2,pcols,2+1*pcols);
imagesc_p(n_r,grid_k_p_,n_w_,n_A,imag(M_k_p_),[],colormap('jet'));
xlim(max_k_c/2*[-1,1]);ylim(max_k_c/2*[-1,1]);
title('imag(M_k_p_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[]); axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
subplot(2,pcols,3+0*pcols);
imagesc_p(n_r,grid_k_p_,n_w_,n_A,real(N_k_p_),[],colormap('jet'));
xlim(max_k_c/2*[-1,1]);ylim(max_k_c/2*[-1,1]);
title('real(N_k_p_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[]); axis square;
subplot(2,pcols,3+1*pcols);
imagesc_p(n_r,grid_k_p_,n_w_,n_A,imag(N_k_p_),[],colormap('jet'));
xlim(max_k_c/2*[-1,1]);ylim(max_k_c/2*[-1,1]);
title('imag(N_k_p_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[]); axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
end;%if flag_disp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% Calculating Z_(delta==0,gamma_z) = < R_{gamma_z}(S) , M >
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
n_gamma_z = 2*n_w_max;
gamma_z_ = get_gamma_0(n_gamma_z);
%%%%%%%%%%%%%%%% ;
D_ = zeros(n_gamma_z,1);
for ngz=0:n_gamma_z-1;
gamma_z = gamma_z_(1+ngz);
D_(1+ngz) = innerproduct_p(n_r,grid_k_p_,n_w_,n_A,rotate_p_to_p_lagrange(n_r,n_w_,n_A,S_k_p_,gamma_z),M_k_p_);
end;%for ngz=0:n_gamma_z-1;
%%%%%%%%%%%%%%%% ;
M_k_q_ = interp_p_to_q(n_r,n_w_,n_A,M_k_p_);
C_q_ = ifft(innerproduct_q_k_only(n_r,grid_k_p_,n_w_,n_A,S_k_q_,M_k_q_))*n_r*pi;
C_ = zeros(n_gamma_z,1);
for ngz=0:n_gamma_z-1;
gamma_z = gamma_z_(1+ngz);
C_(1+ngz) = interp1_0(n_w_max,0,2*pi,C_q_,gamma_z);
end;%for ngz=0:n_gamma_z-1;
%%%%%%%%%%%%%%%% ;
flag_disp=0;
if flag_disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
figure(3); hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
plot(gamma_z_,real(D_),'r.-',gamma_z_,imag(D_),'b.-');
plot(gamma_z_,real(C_),'ro',gamma_z_,imag(C_),'bo');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
end;%if flag_disp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% Generating initial M_x_c_ by translating S_x_c_ by (0.1,0.2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
delta_x_tru = 0.1; delta_y_tru = 0.2;
M_x_c_ = transl_c_to_c(n_x,max_x_c,n_y,max_y_c,S_x_c_,delta_x_tru,delta_y_tru);
M_k_c_ = fft2(M_x_c_); M_k_c_ = recenter2(M_k_c_); M_k_c_ = M_k_c_(1+tmp_kx_ij,1+tmp_ky_ij); M_k_c_ = decenter2(M_k_c_);
N_k_c_ = transf_c_to_c(n_r,max_k_c,n_r,max_k_c,S_k_c_,delta_x_tru,delta_y_tru);
C_k_c = innerproduct_c(n_r,k_c_,n_r,k_c_,M_k_c_,M_k_c_);
disp(sprintf(' %% <M_k_c_,M_k_c_> = (%0.4f,%0.4f)',real(C_k_c),imag(C_k_c)));
C_k_c = innerproduct_c(n_r,k_c_,n_r,k_c_,N_k_c_,N_k_c_);
disp(sprintf(' %% <N_k_c_,N_k_c_> = (%0.4f,%0.4f)',real(C_k_c),imag(C_k_c)));
flag_disp=0;
if flag_disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
figure(4); pcols=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
subplot(2,pcols,1+0*pcols); 
imagesc_c(n_x,x_c_,n_y,y_c_,real(decenter2(S_x_c_)),[],colormap('jet'));
xlim(max_x_c/2*[-1,1]);ylim(max_y_c/2*[-1,1]);
title('real(S_x_c_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[],'Ydir','normal'); axis square;
subplot(2,pcols,1+1*pcols); 
imagesc_c(n_x,x_c_,n_y,y_c_,real(decenter2(M_x_c_)),[],colormap('jet'));
xlim(max_x_c/2*[-1,1]);ylim(max_y_c/2*[-1,1]);
title('real(M_x_c_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[],'Ydir','normal'); axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
subplot(2,pcols,2+0*pcols); 
imagesc_c(n_r,k_c_,n_r,k_c_,real(decenter2(S_k_c_)),[],colormap('jet'));
xlim(max_k_c/2*[-1,1]);ylim(max_k_c/2*[-1,1]);
title('real(S_k_c_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[],'Ydir','normal'); axis square;
subplot(2,pcols,2+1*pcols); 
imagesc_c(n_r,k_c_,n_r,k_c_,imag(decenter2(S_k_c_)),[],colormap('jet'));
xlim(max_k_c/2*[-1,1]);ylim(max_k_c/2*[-1,1]);
title('imag(S_k_c_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[],'Ydir','normal'); axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
subplot(2,pcols,3+0*pcols); 
imagesc_c(n_r,k_c_,n_r,k_c_,real(decenter2(M_k_c_)),[],colormap('jet'));
xlim(max_k_c/2*[-1,1]);ylim(max_k_c/2*[-1,1]);
title('real(M_k_c_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[],'Ydir','normal'); axis square;
subplot(2,pcols,3+1*pcols); 
imagesc_c(n_r,k_c_,n_r,k_c_,imag(decenter2(M_k_c_)),[],colormap('jet'));
xlim(max_k_c/2*[-1,1]);ylim(max_k_c/2*[-1,1]);
title('imag(M_k_c_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[],'Ydir','normal'); axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
subplot(2,pcols,4+0*pcols); 
imagesc_c(n_r,k_c_,n_r,k_c_,real(decenter2(N_k_c_)),[],colormap('jet'));
xlim(max_k_c/2*[-1,1]);ylim(max_k_c/2*[-1,1]);
title('real(N_k_c_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[],'Ydir','normal'); axis square;
subplot(2,pcols,4+1*pcols); 
imagesc_c(n_r,k_c_,n_r,k_c_,imag(decenter2(N_k_c_)),[],colormap('jet'));
xlim(max_k_c/2*[-1,1]);ylim(max_k_c/2*[-1,1]);
title('imag(N_k_c_)','Interpreter','none');
set(gca,'XTick',[],'YTick',[],'Ydir','normal'); axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
end;%if flag_disp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% Calculating Z_(delta,gamma_z) = < R_{gamma_z}(T_{gamma_z}(S)) , M >
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
N_pixels = 3.0; N_pixels_in = N_pixels/2.0; K_max = 40; eps_target = 0.001; l_max = 35;
[n_svd_r,svd_r_,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,svd_s_,svd_V_r_] = gen_Jsvd_2(K_max,N_pixels,eps_target,l_max);
half_diameter_x_c = max_x_c/2; 
n_delta_x = 7; n_delta_y = 7;
[delta_x_,delta_y_] = get_delta_0(N_pixels_in,n_r,half_diameter_x_c,n_delta_x,n_delta_y);
Z_S_svdd_ = innerproduct_q_k_svdd_1(1,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,n_delta_x,delta_x_,n_delta_y,delta_y_,l_max);
Z_M_svdd_ = innerproduct_q_k_svdd_1(0,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,n_delta_x,delta_x_,n_delta_y,delta_y_,l_max);
M_k_q_ = interp_p_to_q(n_r,n_w_,n_A,M_k_p_);
Z_S_svds_ = innerproduct_q_k_svds_1(1,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_,svd_V_r_,n_r,grid_k_p_,n_w_,n_A,S_k_q_,M_k_q_);
Z_M_svds_ = innerproduct_q_k_svds_1(0,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_,svd_V_r_,n_r,grid_k_p_,n_w_,n_A,S_k_q_,M_k_q_);
for nl=0:n_svd_l-1;
Z_S_svds_(1+nl,:) = ifft(Z_S_svds_(1+nl,:))*n_r*pi;
Z_M_svds_(1+nl,:) = ifft(Z_M_svds_(1+nl,:))*n_r*pi;
end;%for nl=0:n_svd_l-1;
n_gamma_z = 1;
gamma_z_ = get_gamma_0(n_gamma_z) + 1*3*pi/4;
%%%%%%%%%%%%%%%% ;
E_S__ = zeros(n_gamma_z,n_delta_x,n_delta_y);
for ndx=0:n_delta_x-1; delta_x = delta_x_(1+ndx);
for ndy=0:n_delta_y-1; delta_y = delta_y_(1+ndy);
for ngz=0:n_gamma_z-1; gamma_z = gamma_z_(1+ngz);
T_k_q_ = S_k_q_;
T_k_q_ = transf_svd_q_to_q(n_svd_r,svd_r_,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,svd_s_,svd_V_r_,n_r,grid_k_p_,n_w_,n_A,T_k_q_,delta_x,delta_y);
T_k_q_ = rotate_q_to_q(n_r,n_w_,n_A,T_k_q_,gamma_z);
E_S__(1+ngz,1+ndx,1+ndy) = innerproduct_p(n_r,grid_k_p_,n_w_,n_A,T_k_q_,M_k_q_);
end;%for ngz=0:n_gamma_z-1; gamma_z = gamma_z_(1+ngz);
end;%for ndy=0:n_delta_y-1; delta_y = delta_y_(1+ndy);
end;%for ndx=0:n_delta_x-1; delta_x = delta_x_(1+ndx);
%%%%%%%%%%%%%%%% ;
D_S__ = zeros(n_gamma_z,n_delta_x,n_delta_y);
for ndx=0:n_delta_x-1; delta_x = delta_x_(1+ndx);
for ndy=0:n_delta_y-1; delta_y = delta_y_(1+ndy);
for ngz=0:n_gamma_z-1; gamma_z = gamma_z_(1+ngz);
T_k_p_ = S_k_p_;
T_k_p_ = transf_p_to_p(n_r,grid_k_p_,n_w_,n_A,T_k_p_,delta_x,delta_y);
T_k_p_ = rotate_p_to_p_lagrange(n_r,n_w_,n_A,T_k_p_,gamma_z);
D_S__(1+ngz,1+ndx,1+ndy) = innerproduct_p(n_r,grid_k_p_,n_w_,n_A,T_k_p_,M_k_p_);
end;%for ngz=0:n_gamma_z-1; gamma_z = gamma_z_(1+ngz);
end;%for ndy=0:n_delta_y-1; delta_y = delta_y_(1+ndy);
end;%for ndx=0:n_delta_x-1; delta_x = delta_x_(1+ndx);
%%%%%%%%%%%%%%%% ;
C_q__ = transpose(Z_S_svdd_)*Z_S_svds_;
C_S__ = zeros(n_gamma_z,n_delta_x,n_delta_y);
for ndx=0:n_delta_x-1; delta_x = delta_x_(1+ndx);
for ndy=0:n_delta_y-1; delta_y = delta_y_(1+ndy);
T_q_ = C_q__(1 + ndx + ndy*n_delta_x,:);
for ngz=0:n_gamma_z-1; gamma_z = gamma_z_(1+ngz);
C_S__(1+ngz,1+ndx,1+ndy) = interp1_0(n_w_max,0,2*pi,T_q_,gamma_z);
end;%for ngz=0:n_gamma_z-1; gamma_z = gamma_z_(1+ngz);
end;%for ndy=0:n_delta_y-1; delta_y = delta_y_(1+ndy);
end;%for ndx=0:n_delta_x-1; delta_x = delta_x_(1+ndx);
%%%%%%%%%%%%%%%% ;
E_M__ = zeros(n_gamma_z,n_delta_x,n_delta_y);
for ndx=0:n_delta_x-1; delta_x = delta_x_(1+ndx);
for ndy=0:n_delta_y-1; delta_y = delta_y_(1+ndy);
for ngz=0:n_gamma_z-1; gamma_z = gamma_z_(1+ngz);
T_k_q_ = S_k_q_;
T_k_q_ = rotate_q_to_q(n_r,n_w_,n_A,T_k_q_,gamma_z);
T_k_q_ = transf_svd_q_to_q(n_svd_r,svd_r_,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,svd_s_,svd_V_r_,n_r,grid_k_p_,n_w_,n_A,T_k_q_,delta_x,delta_y);
E_M__(1+ngz,1+ndx,1+ndy) = innerproduct_p(n_r,grid_k_p_,n_w_,n_A,T_k_q_,M_k_q_);
end;%for ngz=0:n_gamma_z-1; gamma_z = gamma_z_(1+ngz);
end;%for ndy=0:n_delta_y-1; delta_y = delta_y_(1+ndy);
end;%for ndx=0:n_delta_x-1; delta_x = delta_x_(1+ndx);
%%%%%%%%%%%%%%%% ;
D_M__ = zeros(n_gamma_z,n_delta_x,n_delta_y);
for ndx=0:n_delta_x-1; delta_x = delta_x_(1+ndx);
for ndy=0:n_delta_y-1; delta_y = delta_y_(1+ndy);
for ngz=0:n_gamma_z-1; gamma_z = gamma_z_(1+ngz);
T_k_p_ = S_k_p_;
T_k_p_ = rotate_p_to_p_lagrange(n_r,n_w_,n_A,T_k_p_,gamma_z);
T_k_p_ = transf_p_to_p(n_r,grid_k_p_,n_w_,n_A,T_k_p_,delta_x,delta_y);
D_M__(1+ngz,1+ndx,1+ndy) = innerproduct_p(n_r,grid_k_p_,n_w_,n_A,T_k_p_,M_k_p_);
end;%for ngz=0:n_gamma_z-1; gamma_z = gamma_z_(1+ngz);
end;%for ndy=0:n_delta_y-1; delta_y = delta_y_(1+ndy);
end;%for ndx=0:n_delta_x-1; delta_x = delta_x_(1+ndx);
%%%%%%%%%%%%%%%% ;
C_q__ = transpose(Z_M_svdd_)*Z_M_svds_;
C_M__ = zeros(n_gamma_z,n_delta_x,n_delta_y);
for ndx=0:n_delta_x-1; delta_x = delta_x_(1+ndx);
for ndy=0:n_delta_y-1; delta_y = delta_y_(1+ndy);
T_q_ = C_q__(1 + ndx + ndy*n_delta_x,:);
for ngz=0:n_gamma_z-1; gamma_z = gamma_z_(1+ngz);
C_M__(1+ngz,1+ndx,1+ndy) = interp1_0(n_w_max,0,2*pi,T_q_,gamma_z);
end;%for ngz=0:n_gamma_z-1; gamma_z = gamma_z_(1+ngz);
end;%for ndy=0:n_delta_y-1; delta_y = delta_y_(1+ndy);
end;%for ndx=0:n_delta_x-1; delta_x = delta_x_(1+ndx);
%%%%%%%%%%%%%%%%;
disp_flag=1;
figure(5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
subplot(1,2,1);
hold on;
plot((E_S__(:)-C_S__(:))./abs(E_S__(:)),'gx');title('error for C_S__','Interpreter','none');
plot((D_S__(:)-C_S__(:))./abs(D_S__(:)),'ms');title('error for C_S__','Interpreter','none');
subplot(1,2,2);
hold on;
plot((E_M__(:)-C_M__(:))./abs(E_M__(:)),'gx');title('error for C_M__','Interpreter','none');
plot((D_M__(:)-C_M__(:))./abs(D_M__(:)),'ms');title('error for C_M__','Interpreter','none');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

