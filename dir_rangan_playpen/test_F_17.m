% test functional svd-expansion. ;
% testing functional lsq. ;
% determining the overlap integral between two plane-waves on a disc, ;
% as well as the derivative of that overlap integral. : ;
% Now defining functions that can be used to determine the derivative ;
% of terms like f_{j}^{\dagger} f_{j^{\prime}}, ;
% as well as the error E. ;

verbose=1;
max_k_p = 64; N_pixel = 3.0;
R_target = max_k_p-0.5; z_target = N_pixel*pi*sqrt(2); D_target = z_target/(2*pi*R_target);
C = 2*pi*max_k_p;
h_ = @(kd) 4*pi*besselj(1,2*pi*kd) ./ kd ; 
dh_ = @(kd) -4*pi*besselj(1,2*pi*kd)./kd.^2 + 4*pi^2*(besselj(0,2*pi*kd)-besselj(2,2*pi*kd))./kd ; 

%{
% check derivative of ff;
delta_1 = rand(2,1)*D_target/sqrt(2); delta_2 = rand(2,1)*D_target/sqrt(2);
for abseps = 10.^[-1:-1:-6];
epsilon = rand(2,1)*norm(delta_1)*abseps;
ff0 = h_(norm(delta_2-delta_1)*max_k_p);
ff1 = h_(norm(delta_2-(delta_1 + epsilon))*max_k_p);
df1 = (ff1 - ff0)/norm(epsilon);
dd = delta_2-delta_1;
df2 = ( dh_(norm(dd)*max_k_p)*max_k_p * dot(-dd/norm(dd),epsilon) )/norm(epsilon) ;
disp(sprintf(' %% abseps %0.6f; df1 %f df2 %f relative difference %f',abseps,df1,df2,abs(df1-df2)/abs(df2)));
end;%for abseps = 10.^[-1:-1:-6];
 %}

%{
% check derivative of inv(A+eps) ;
A = randn(4); B = inv(A);
for abseps = 10.^[-1:-1:-6];
epsilon = randn(4)*norm(A)*abseps;
B1 = B;
B2 = inv(A+epsilon);
dB1 = (B2 - B1)/norm(epsilon);
dB2 = -(B*epsilon*B)/norm(epsilon);
disp(sprintf(' %% abseps %0.6f; relative difference %f',abseps,norm(dB1-dB2)/norm(dB2)));
end;%for abseps = 10.^[-1:-1:-6];
 %}

% check derivative of error term (residual l2-norm). ;
n_point = 128; max_x_c = 1; 
[max_k_c,max_k_p,grid_x_c_,d_x_c,X_x_c_,Y_x_c_,R_x_c_,W_x_c_,grid_k_c_,d_k_c,X_k_c_,Y_k_c_,R_k_c_,W_k_c_,grid_x_r_,grid_x_w_,R_x_p_,W_x_p_,grid_k_r_,d_k_r,grid_k_w_,d_k_w,R_k_p_,W_k_p_,X_k_p_,Y_k_p_] = test_F_grid_0(n_point,max_x_c);
N_pixel = 3.0; % number of wavelengths allowed in displacement grid. ;
R_target = max_k_p-0.5; z_target = N_pixel*pi*sqrt(2); D_target = z_target/(2*pi*R_target);
%%%%%%%%%%%%%%%%;
n_node = 8;
T_k_p__ = cell(n_node,1);
delta1_ = zeros(2,n_node);
for nnode=1:n_node; 
delta_d = D_target*rand(); delta_w = 2*pi*rand();
delta1_(:,nnode) = delta_d*[cos(delta_w);sin(delta_w)];
T_k_p__{nnode} = test_F_T_k_p_0(n_point,max_x_c,delta_d,delta_w);
end;%for nnode=1:n_node;
%%%%%%%%%%%%%%%%;
F_ = zeros(n_node,n_node);
H1_ = zeros(n_node,n_node);
for nnodeA=1:n_node; for nnodeB=1:n_node;
tmp_ = conj(T_k_p__{nnodeA}) .* T_k_p__{nnodeB}; tmp_ = tmp_*diag(2*pi*grid_k_r_); tmp = sum(tmp_(:)).*(2*pi*d_k_r)*d_k_w; tmp = tmp/(pi*max_k_p.^2); 
F_(nnodeA,nnodeB) = tmp;
d = (delta1_(:,nnodeB) - delta1_(:,nnodeA));
if (norm(d)>1e-6); H1_(nnodeA,nnodeB) = h_(norm(d)*max_k_p); end;%if (dj~=0); 
if (norm(d)<=1e-6); H1_(nnodeA,nnodeB) = (2*pi).^2; end;%if (dj==0); 
end;end;%for nnodeA=1:n_node; for nnodeB=1:n_node;
G1_ = inv(H1_);
%%%%%%%%%%%%%%%%;
delta_d = D_target*rand(); delta_w = 2*pi*rand(); delta_t = delta_d*[cos(delta_w);sin(delta_w)];
ft1_ = zeros(n_node,1);
for nnode=1:n_node; d = (delta_t - delta1_(:,nnode)); ft1_(nnode) = h_(norm(d)*max_k_p); end;
alpha1_ = G1_*ft1_;
E1 = (2*pi).^2 - transpose(ft1_)*G1_*ft1_;
%%%%%%%%%%%%%%%%;
k = 3;
for abseps = 10.^(-1:-1:-6);
epsilon = randn(2,1)*norm(delta1_(:,k))*abseps;
delta2_ = delta1_; delta2_(:,k) = delta1_(:,k) + epsilon;
H2_ = zeros(n_node,n_node);
for nnodeA=1:n_node; for nnodeB=1:n_node;
d = (delta2_(:,nnodeB) - delta2_(:,nnodeA));
if (norm(d)>1e-6); H2_(nnodeA,nnodeB) = h_(norm(d)*max_k_p); end;%if (dj~=0); 
if (norm(d)<=1e-6); H2_(nnodeA,nnodeB) = (2*pi).^2; end;%if (dj==0); 
end;end;%for nnodeA=1:n_node; for nnodeB=1:n_node;
G2_ = inv(H2_);
ft2_ = zeros(n_node,1);
for nnode=1:n_node; d = (delta_t - delta2_(:,nnode)); ft2_(nnode) = h_(norm(d)*max_k_p); end;
alpha2_ = G2_*ft2_;
E2 = (2*pi).^2 - transpose(ft2_)*G2_*ft2_;
%%%%%%%%%%%%%%%%;
%{
dH1_ = (H2_-H1_)/norm(epsilon);
dh2_ = zeros(n_node,1); 
for nnode=1:n_node; 
if (nnode~=k);
dd = delta1_(:,nnode)-delta1_(:,k);
dh2_(nnode) = ( dh_(norm(dd)*max_k_p)*max_k_p * dot(-dd/norm(dd),epsilon) )/norm(epsilon) ;
end;%if (nnode~=k);
end;%for nnode=1:n_node; 
dH2_ = zeros(n_node); dH2_(k,:) = transpose(dh2_); dH2_(:,k) = dh2_;
disp(sprintf(' %% abseps %0.6f; relative difference in dH %f',abseps,norm(dH2_-dH1_)/norm(dH2_)));
dG1_ = (G2_-G1_)/norm(epsilon);
dG2_ = -G1_*dH2_*G1_;
disp(sprintf(' %% abseps %0.6f; relative difference in dG %f',abseps,norm(dG2_-dG1_)/norm(dG2_)));
dft1_ = (ft2_ - ft1_)/norm(epsilon);
dft2_ = zeros(n_node,1); 
dd = delta_t-delta1_(:,k);
dft2_(k) = ( dh_(norm(dd)*max_k_p)*max_k_p * dot(-dd/norm(dd),epsilon) )/norm(epsilon) ;
disp(sprintf(' %% abseps %0.6f; relative difference in dft %f',abseps,norm(dft2_-dft1_)/norm(dft2_)));
dalpha1_ = (alpha2_-alpha1_)/norm(epsilon);
dalpha2_ = dG2_*ft1_ + G1_*dft2_;
disp(sprintf(' %% abseps %0.6f; relative difference in dalpha %f',abseps,norm(dalpha2_-dalpha1_)/norm(dalpha2_)));
dE1 = (E2-E1)/norm(epsilon);
dE2 = 2*alpha1_(k)*dot(dh2_,alpha1_) - 2*alpha1_(k)*dft2_(k);
disp(sprintf(' %% abseps %0.6f; relative difference in dE %f',abseps,norm(dE2-dE1)/norm(dE2)));
 %}
%%%%%%%%%%%%%%%%;
%{
dH1_ = H2_-H1_ ;
dhk_ = zeros(n_node,1); 
for nnode=1:n_node; 
if (nnode~=k);
dd = delta1_(:,nnode)-delta1_(:,k);
dhk_(nnode) = dh_(norm(dd)*max_k_p)*max_k_p * dot(-dd/norm(dd),epsilon);
end;%if (nnode~=k);
end;%for nnode=1:n_node; 
dH2_ = zeros(n_node); dH2_(k,:) = transpose(dhk_); dH2_(:,k) = dhk_;
disp(sprintf(' %% abseps %0.6f; relative difference in dH %f',abseps,norm(dH2_-dH1_)/norm(dH2_)));
dG1_ = G2_-G1_ ;
dG2_ = -G1_*dH2_*G1_;
disp(sprintf(' %% abseps %0.6f; relative difference in dG %f',abseps,norm(dG2_-dG1_)/norm(dG2_)));
dd = delta_t-delta1_(:,k);
dftk = dh_(norm(dd)*max_k_p)*max_k_p * dot(-dd/norm(dd),epsilon);
dft1_ = ft2_-ft1_;
dft2_ = zeros(n_node,1); dft2_(k) = dftk;
disp(sprintf(' %% abseps %0.6f; relative difference in dft %f',abseps,norm(dft2_-dft1_)/norm(dft2_)));
dalpha1_ = alpha2_-alpha1_;
dalpha2_ = dG2_*ft1_ + G1_*dft2_;
disp(sprintf(' %% abseps %0.6f; relative difference in dalpha %f',abseps,norm(dalpha2_-dalpha1_)/norm(dalpha2_)));
dE1 = E2-E1;
dE2 = 2*alpha1_(k)*dot(dhk_,alpha1_) - 2*alpha1_(k)*dftk;
disp(sprintf(' %% abseps %0.6f; relative difference in dE %f',abseps,norm(dE2-dE1)/norm(dE2)));
 %}
%%%%%%%%%%%%%%%%;
dH1_ = H2_-H1_ ;
dd_ = delta1_ - repmat(delta1_(:,k),1,n_node); hypot_dd_ = hypot(dd_(1,:),dd_(2,:));
dhk__ = - repmat(dh_(hypot_dd_*max_k_p)*max_k_p,2,1) .* (dd_*diag(1./hypot_dd_));
dhk__(:,find(hypot_dd_==0)) = [0;0];
dhk_ = transpose(dhk__)*epsilon;
dH2_ = zeros(n_node); dH2_(k,:) = transpose(dhk_); dH2_(:,k) = dhk_;
disp(sprintf(' %% abseps %0.6f; relative difference in dH %f',abseps,norm(dH2_-dH1_)/norm(dH2_)));
dG1_ = G2_-G1_ ;
dG2_ = -G1_*dH2_*G1_;
disp(sprintf(' %% abseps %0.6f; relative difference in dG %f',abseps,norm(dG2_-dG1_)/norm(dG2_)));
dd = delta_t-delta1_(:,k); hypot_dd = hypot(dd(1),dd(2));
dftk_ = - dh_(hypot_dd*max_k_p)*max_k_p * (dd/hypot_dd); if (hypot_dd==0); dftk_ = [0;0]; end;
dftk = transpose(dftk_)*epsilon;
dft1_ = ft2_-ft1_;
dft2_ = zeros(n_node,1); dft2_(k) = dftk;
disp(sprintf(' %% abseps %0.6f; relative difference in dft %f',abseps,norm(dft2_-dft1_)/norm(dft2_)));
dalpha1_ = alpha2_-alpha1_;
dalpha2_ = dG2_*ft1_ + G1_*dft2_;
disp(sprintf(' %% abseps %0.6f; relative difference in dalpha %f',abseps,norm(dalpha2_-dalpha1_)/norm(dalpha2_)));
dE1 = E2-E1;
dE2 = 2*alpha1_(k)*dot(dhk_,alpha1_) - 2*alpha1_(k)*dftk;
disp(sprintf(' %% abseps %0.6f; relative difference in dE %f',abseps,norm(dE2-dE1)/norm(dE2)));
%%%%%%%%%%%%%%%%;
end;%for abseps = 10.^(-1:-1:-6);

% compact representation of various derivatives. ;
%dd_ = delta_node_ - repmat(delta_node_(:,k),1,n_node); hypot_dd_ = hypot(dd_(1,:),dd_(2,:));
%dhk__ = - repmat(dh_(hypot_dd_*max_k_p)*max_k_p,2,1) .* (dd_*diag(1./hypot_dd_));
%dhk__(:,find(hypot_dd_==0)) = [0;0];
%dd = delta_t-delta_node_(:,k); hypot_dd = hypot(dd(1),dd(2));
%dftk_ = - dh_(hypot_dd*max_k_p)*max_k_p * (dd/hypot_dd); if (hypot_dd==0); dftk_ = [0;0]; end;
%dEk_ = 2*alpha_(k)*(dhk__*alpha_) - 2*alpha_(k)*dftk_;
