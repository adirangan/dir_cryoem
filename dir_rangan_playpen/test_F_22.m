function test_F_22(N_pixel);
% comparing svd-expansion with functional lsq. ;
% analytically calculate the overlap integral between two plane-waves on a disc, ;
% as well as the derivative of that overlap integral. : ;
% Can explicitly write down function that can be used to determine the derivative ;
% of terms like f_{j}^{\dagger} f_{j^{\prime}}, ;
% as well as the error E. ;
% We can calculate the total (average) error across disc using either svd-expansion or lsq-interpolation. ;
% Using quadrature grid. ;
% Switching to pseudoinverse. ;
% Pulling out delta_t independent term from derivative. ;
% Using roughly equispaced initial nodes. ;
% Attempting to vectorize derivative calculation. ;
% Upgrading to jacpts instead of using orthopoly_node_weight_matrix_0. ;

if (nargin<1); N_pixel = 1.5; end;

verbose=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% generate grids. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
n_point = 128; max_x_c = 1; 
[max_k_c,max_k_p,grid_x_c_,d_x_c,X_x_c_,Y_x_c_,R_x_c_,W_x_c_,grid_k_c_,d_k_c,X_k_c_,Y_k_c_,R_k_c_,W_k_c_,grid_x_r_,grid_x_w_,R_x_p_,W_x_p_,grid_k_r_,d_k_r,grid_k_w_,d_k_w,R_k_p_,W_k_p_,X_k_p_,Y_k_p_] = test_F_grid_0(n_point,max_x_c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% Determine svd expansion ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
%N_pixel = 0.5; % number of wavelengths allowed in displacement grid. ;
%eps_target = 0.1; % tolerance used for svd-expansion. ;
l_max = 32; % maximum order of bessel-functions to retain. ;
n_r_degree = 63; % degree of orthonormal-polynomial to use for r = |k|. ;
n_d_degree = 65; % degree of orthonormal-polynomial to use for d = |delta|. ;
R_target = max_k_p-0.5; z_target = N_pixel*pi*sqrt(2); D_target = z_target/(2*pi*R_target);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% use to generate svd-expansion. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
%[n_svd_r,svd_r_,svd_r_m,svd_r_c,svd_r_w_,svd_r_Jv_,n_svd_d,svd_d_,svd_d_m,svd_d_c,svd_d_w_,svd_d_Jv_,n_svd_l,svd_l_,svd_U_d_,svd_s_,svd_V_r_] = gen_Jsvd_5(max_k_p,N_pixel,eps_target,l_max,n_r_degree,n_d_degree);

flag_skip=1;
if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% compare T with svd-expansion on regular polar grid. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
delta_d = D_target; delta_w = 1*pi/6; 
[T_k_p_] = test_F_T_k_p_0(n_point,max_x_c,delta_d,delta_w);
Tlim_k = [-1,+1]; 
for eps_target=[0.01];%for eps_target=[0.1,0.01,0.001];
clear n_svd_r svd_r_ svd_r_m svd_r_c svd_r_w_ svd_r_Jv_ n_svd_d svd_d_ svd_d_m svd_d_c svd_d_w_ svd_d_Jv_ n_svd_l svd_l_ svd_U_d_ svd_s_ svd_V_r_;
[n_svd_r,svd_r_,svd_r_m,svd_r_c,svd_r_w_,svd_r_Jv_,n_svd_d,svd_d_,svd_d_m,svd_d_c,svd_d_w_,svd_d_Jv_,n_svd_l,svd_l_,svd_U_d_,svd_s_,svd_V_r_] = gen_Jsvd_5(max_k_p,N_pixel,eps_target,l_max,n_r_degree,n_d_degree);
[X_k_p_] = test_F_X_k_p_5(n_point,max_x_c,delta_d,delta_w,n_svd_r,svd_r_,svd_r_m,svd_r_c,svd_r_w_,svd_r_Jv_,n_svd_d,svd_d_,svd_d_m,svd_d_c,svd_d_w_,svd_d_Jv_,n_svd_l,svd_l_,svd_U_d_,svd_s_,svd_V_r_);
Xlim_k = [-1,+1];
tmp_ = abs(T_k_p_-X_k_p_).^2; tmp_ = tmp_*diag(2*pi*grid_k_r_); tmp = sum(tmp_(:)).*(2*pi*d_k_r)*d_k_w; tmp = sqrt(tmp/(pi*max_k_p.^2)); Dlim_k = 0.25*tmp*[-1,+1];
cra = colormap_pm(64);
figure(1); clf;
subplot(2,3,1); clim = polarpatch(R_k_p_,W_k_p_,real(T_k_p_),Tlim_k,0,0,1,cra);  title(sprintf('real(T) [%0.1f,%0.1f]',Tlim_k)); axis equal;
%colorbar_handle = colorbar('v'); tmp = get(colorbar_handle,'Limits'); set(colorbar_handle,'XTick',[min(tmp),max(tmp)],'XTickLabel',Tlim_k); 
subplot(2,3,4); clim = polarpatch(R_k_p_,W_k_p_,imag(T_k_p_),Tlim_k,0,0,1,cra);  title(sprintf('imag(T) [%0.1f,%0.1f]',Tlim_k)); axis equal;
%colorbar_handle = colorbar('v'); tmp = get(colorbar_handle,'Limits'); set(colorbar_handle,'XTick',[min(tmp),max(tmp)],'XTickLabel',Tlim_k); 
subplot(2,3,2); clim = polarpatch(R_k_p_,W_k_p_,real(X_k_p_),Xlim_k,0,0,1,cra);  title(sprintf('real(T^{svd}) [%0.1f,%0.1f]',Xlim_k)); axis equal;
%colorbar_handle = colorbar('v'); tmp = get(colorbar_handle,'Limits'); set(colorbar_handle,'XTick',[min(tmp),max(tmp)],'XTickLabel',Xlim_k); 
subplot(2,3,5); clim = polarpatch(R_k_p_,W_k_p_,imag(X_k_p_),Xlim_k,0,0,1,cra);  title(sprintf('imag(T^{svd}) [%0.1f,%0.1f]',Xlim_k)); axis equal;
%colorbar_handle = colorbar('v'); tmp = get(colorbar_handle,'Limits'); set(colorbar_handle,'XTick',[min(tmp),max(tmp)],'XTickLabel',Xlim_k); 
subplot(2,3,3); clim = polarpatch(R_k_p_,W_k_p_,real(T_k_p_-X_k_p_),Dlim_k,0,0,1,cra);  title(sprintf('real(T-T^{svd}) [%f,%f]',Dlim_k)); axis equal;
%colorbar_handle = colorbar('v'); tmp = get(colorbar_handle,'Limits'); set(colorbar_handle,'XTick',[min(tmp),max(tmp)],'XTickLabel',Dlim_k); 
subplot(2,3,6); clim = polarpatch(R_k_p_,W_k_p_,imag(T_k_p_-X_k_p_),Dlim_k,0,0,1,cra);  title(sprintf('imag(T-T^{svd}) [%f,%f]',Dlim_k)); axis equal;
%colorbar_handle = colorbar('v'); tmp = get(colorbar_handle,'Limits'); set(colorbar_handle,'XTick',[min(tmp),max(tmp)],'XTickLabel',Dlim_k); 
%suptitle(sprintf('n_{point} %d  K_{max} 2*pi*%d  N^{pixel} %0.1f  eps %0.3f(%f)',n_point,max_k_p,N_pixel,eps_target,eps_target/svd_s_(1)));
set(gcf,'Position',1+[0,0,1024*1.5,1024]); 
fname_base = sprintf('X_n%d_K%d_N%.2d_e%.3d',n_point,max_k_p,round(10*N_pixel),round(1000*eps_target));
print('-djpeg',sprintf('./dir_svd/%s.jpg',fname_base));
print('-depsc',sprintf('./dir_svd/%s.eps',fname_base));
%%%%%%%%%%%%%%%% ;
end;%for eps_target=[0.1,0.01,0.001];
disp('returning'); return;
end;%if ~flag_skip;

flag_skip=0;
if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% compare T with svd-expansion on ray from a jacobi polar grid. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
eps_target_ = 10.^[-1:-0.5:-4]; l_max=32; n_r_degree = 63; n_d_degree = 65;
E_abs__ = zeros(n_d_degree,length(eps_target_));
E_rel__ = zeros(n_d_degree,length(eps_target_));
n_svd_l_ = zeros(length(eps_target_),1);
for neps = 1:length(eps_target_);
eps_target = eps_target_(neps);
disp(sprintf(' %% nesp %d/%d; eps_target %f',neps,length(eps_target_),eps_target));
clear n_svd_r svd_r_ svd_r_m svd_r_c svd_r_w_ svd_r_Jv_ n_svd_d svd_d_ svd_d_m svd_d_c svd_d_w_ svd_d_Jv_ n_svd_l svd_l_ svd_U_d_ svd_s_ svd_V_r_;
[n_svd_r,svd_r_,svd_r_m,svd_r_c,svd_r_w_,svd_r_Jv_,n_svd_d,svd_d_,svd_d_m,svd_d_c,svd_d_w_,svd_d_Jv_,n_svd_l,svd_l_,svd_U_d_,svd_s_,svd_V_r_] = gen_Jsvd_5(max_k_p,N_pixel,eps_target,l_max,n_r_degree,n_d_degree);
n_svd_l_(neps) = n_svd_l;
n_w = max(3,ceil(2*pi*n_r_degree)); n_w_ = n_w*ones(n_svd_r,1);
r_quad_weight_ = ( (2*pi/n_w)*ones(n_w,1) * transpose(svd_r_w_) ) .* ((2*pi*R_target).^2/4);
r_d_ = svd_r_; r_w_ = linspace(0,2*pi,n_w+1); r_w_ = r_w_(1:end-1); 
delta_w = 1*pi/6; 
for ndelta = 1:n_svd_d;
delta_d = svd_d_(ndelta);
disp(sprintf(' %% %% ndelta %d/%d; delta_d %f',ndelta,n_svd_d,delta_d));
[T_k_p_] = test_F_T_k_p_1(r_d_,r_w_,delta_d,delta_w);
Tlim_k = [-1,+1]; 
[X_k_p_] = test_F_X_k_p_6(r_d_,r_w_,delta_d,delta_w,n_svd_r,svd_r_,svd_r_m,svd_r_c,svd_r_w_,svd_r_Jv_,n_svd_d,svd_d_,svd_d_m,svd_d_c,svd_d_w_,svd_d_Jv_,n_svd_l,svd_l_,svd_U_d_,svd_s_,svd_V_r_);
Xlim_k = [-1,+1];
Dlim_k = sqrt(10)*eps_target*[-1,+1];
flag_disp=0;
if flag_disp;
subplot(2,3,1); polarpatch_adaptive(r_d_,n_w_,real(T_k_p_(:)),colormap(colormap_beach(64)),Tlim_k,-1.25*2*pi*R_target,0,1); title('real(T)'); set(gca,'Xtick',[],'Ytick',[]); axis equal;
subplot(2,3,4); polarpatch_adaptive(r_d_,n_w_,imag(T_k_p_(:)),colormap(colormap_beach(64)),Tlim_k,-1.25*2*pi*R_target,0,1); title('imag(T)'); set(gca,'Xtick',[],'Ytick',[]); axis equal;
subplot(2,3,2); polarpatch_adaptive(r_d_,n_w_,real(X_k_p_(:)),colormap(colormap_beach(64)),Xlim_k,-1.25*2*pi*R_target,0,1); title('real(X)'); set(gca,'Xtick',[],'Ytick',[]); axis equal;
subplot(2,3,5); polarpatch_adaptive(r_d_,n_w_,imag(X_k_p_(:)),colormap(colormap_beach(64)),Xlim_k,-1.25*2*pi*R_target,0,1); title('imag(X)'); set(gca,'Xtick',[],'Ytick',[]); axis equal;
subplot(2,3,3); polarpatch_adaptive(r_d_,n_w_,real(T_k_p_(:)-X_k_p_(:)),colormap(colormap_beach(64)),Dlim_k,-1.25*2*pi*R_target,0,1); title(sprintf('real(T-X) %0.2f',Dlim_k(2))); set(gca,'Xtick',[],'Ytick',[]); axis equal;
subplot(2,3,6); polarpatch_adaptive(r_d_,n_w_,imag(T_k_p_(:)-X_k_p_(:)),colormap(colormap_beach(64)),Dlim_k,-1.25*2*pi*R_target,0,1); title(sprintf('imag(T-X) %0.2f',Dlim_k(2))); set(gca,'Xtick',[],'Ytick',[]); axis equal;
drawnow();
end;%if flag_disp;
tmp_1 = sum(sum(abs(T_k_p_ - X_k_p_).*r_quad_weight_));
tmp_2 = sum(sum(abs(T_k_p_).*r_quad_weight_));
tmp_E_abs = tmp_1; tmp_E_rel = tmp_1./tmp_2;
E_abs__(ndelta,neps) = tmp_E_abs;
E_rel__(ndelta,neps) = tmp_E_rel;
end;%for ndelta = 1:n_svd_d;
end;%for neps = 1:length(eps_target_);
E_abs_ = transpose(E_abs__)*svd_d_w_;
E_rel_ = transpose(E_rel__)*svd_d_w_;
%%%%%%%%;
figure(3);clf;
subplot(2,2,1); plot(svd_d_,log10(E_abs__)); xlabel('d'); ylabel('log10(E)'); title('E abs');
subplot(2,2,2); plot(svd_d_,log10(E_rel__)); xlabel('d'); ylabel('log10(E)'); title('E rel');
subplot(2,2,3); plot(-log10(eps_target_),log10(E_rel_),'ro-'); xlabel('-log10(eps)'); ylabel('log10(E)'); title('E rel');
subplot(2,2,4); plot(n_svd_l_,log10(E_rel_),'ro-'); xlabel('n terms'); ylabel('log10(E)'); title('E rel');
fname_base = sprintf('E_N%.2d',round(10*N_pixel));
print('-djpeg',sprintf('./dir_svd/%s.jpg',fname_base));
print('-depsc',sprintf('./dir_svd/%s.eps',fname_base));
%disp('returning'); return;
end;%if ~flag_skip;

flag_skip=0;
if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% compare T with brute-force + linear-interpolation on ray from a jacobi polar grid. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% First svd-expansion; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
verbose=1;
eps_target = 1.0; l_max=0; n_r_degree = 63; n_d_degree = 128;
clear n_svd_r svd_r_ svd_r_m svd_r_c svd_r_w_ svd_r_Jv_ n_svd_d svd_d_ svd_d_m svd_d_c svd_d_w_ svd_d_Jv_ n_svd_l svd_l_ svd_U_d_ svd_s_ svd_V_r_;
[n_svd_r,svd_r_,svd_r_m,svd_r_c,svd_r_w_,svd_r_Jv_,n_svd_d,svd_d_,svd_d_m,svd_d_c,svd_d_w_,svd_d_Jv_,n_svd_l,svd_l_,svd_U_d_,svd_s_,svd_V_r_] = gen_Jsvd_5(max_k_p,N_pixel,eps_target,l_max,n_r_degree,n_d_degree);
%%%%%%%%;
n_w = max(3,ceil(2*pi*n_r_degree)); n_w_ = n_w*ones(n_svd_r,1);
r_quad_weight_ = ( (2*pi/n_w)*ones(n_w,1) * transpose(svd_r_w_) ) .* ((2*pi*R_target).^2/4);
r_d_ = svd_r_; r_w_ = linspace(0,2*pi,n_w+1); r_w_ = r_w_(1:end-1); 
delta_w = 1*pi/6; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% Now brute-force + linear-interpolation;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
verbose=1;
n_node_ = [8,16,24,32];
F_abs__ = zeros(n_svd_d,length(n_node_));
F_rel__ = zeros(n_svd_d,length(n_node_));
for nnode_=1:length(n_node_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
n_node = n_node_(nnode_);
if (verbose); disp(sprintf(' %% n_node %f',n_node)); end;
T_k_p__ = cell(n_node,1);
delta_d_ = zeros(n_node,1);
for nnode=1:n_node;
delta_d = D_target * ((nnode-1)/(n_node-1)); delta_w = 1*pi/6;
delta_d_(nnode) = delta_d;
T_k_p__{nnode} = test_F_T_k_p_1(r_d_,r_w_,delta_d,delta_w);
end;%for nnode=1:n_node;
%%%%%%%%%%%%%%%%;
for nd=1:n_svd_d;
%%%%%%%%%%%%%%%% ;
delta_d = svd_d_(nd); delta_w = 1*pi/6; 
[T_k_p_] = test_F_T_k_p_1(r_d_,r_w_,delta_d,delta_w);
%%%%%%%%%%%%%%%% ;
ij_par = find(abs(delta_d_ - delta_d)<1e-7);
if (~isempty(ij_par)); 
if (verbose>1); disp(sprintf(' %% nd %d delta_d %f ij_par %d',nd,delta_d,ij_par)); end;
B_k_p_ = T_k_p__{ij_par};
end;%if (~isempty(ij_par)); 
if (isempty(ij_par));
ij_pre = max(find(delta_d_<delta_d));
ij_pos = min(find(delta_d_>delta_d));
if (verbose>2); disp(sprintf(' %% nd %d delta_d %f ij_pre %d ij_pos %d',nd,delta_d,ij_pre,ij_pos)); end;
d_pre = delta_d - delta_d_(ij_pre); d_pos = delta_d_(ij_pos) - delta_d;
w_pre = d_pos/(d_pos+d_pre); w_pos = d_pre/(d_pos+d_pre);
B_k_p_ = w_pre*T_k_p__{ij_pre} + w_pos*T_k_p__{ij_pos};
end;%if (isempty(ij_par));
%%%%%%%%%%%%%%%% ;
tmp_1 = sum(sum(abs(T_k_p_ - B_k_p_).*r_quad_weight_));
tmp_2 = sum(sum(abs(T_k_p_).*r_quad_weight_));
tmp_F_abs = tmp_1; tmp_F_rel = tmp_1./tmp_2;
F_abs__(nd,nnode_) = tmp_F_abs;
F_rel__(nd,nnode_) = tmp_F_rel;
%%%%%%%%%%%%%%%% ;
end;%for nd=1:n_d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
end;%for nnode_=1:length(n_node_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
F_abs_ = transpose(F_abs__)*svd_d_w_;
F_rel_ = transpose(F_rel__)*svd_d_w_;
%%%%%%%%;
figure(4);clf;
subplot(2,2,1); plot(svd_d_,log10(F_abs__)); xlabel('d'); ylabel('log10(F)'); title('F abs');
subplot(2,2,2); plot(svd_d_,log10(F_rel__)); xlabel('d'); ylabel('log10(F)'); title('F rel');
subplot(2,2,3); plot(n_node_,log10(F_rel_),'bs-'); xlabel('n nodes'); ylabel('log10(F)'); title('F rel');
subplot(2,2,4); plot(ceil(pi*n_node_.^2),log10(F_rel_),'bs-'); xlabel('n terms'); ylabel('log10(F)'); title('F rel');
fname_base = sprintf('F_N%.2d',round(10*N_pixel));
print('-djpeg',sprintf('./dir_svd/%s.jpg',fname_base));
print('-depsc',sprintf('./dir_svd/%s.eps',fname_base));
%disp('returning'); return;
end;%if ~flag_skip;

flag_skip=0;
if ~flag_skip;
figure(5);clf;
subplot(1,1,1);
hold on;
plot(log10(n_svd_l_),log10(E_rel_),'ro-'); 
plot(log10(ceil(pi*n_node_.^2)),log10(F_rel_),'bs-'); 
xlabel('log10(nterms)'); ylabel('log10(error)'); title('E rel');
legend('svd-expansion','linear-interpolation');
fname_base = sprintf('EvsF_N%.2d',round(10*N_pixel));
print('-djpeg',sprintf('./dir_svd/%s.jpg',fname_base));
print('-depsc',sprintf('./dir_svd/%s.eps',fname_base));
end;%if ~flag_skip;

fname_base = sprintf('EvsF_N%.2d',round(10*N_pixel));
fname_mat = sprintf('./dir_svd/%s.mat',fname_base);
save(fname_mat,'n_svd_l_','E_rel_','n_node_','F_rel_','N_pixel');
