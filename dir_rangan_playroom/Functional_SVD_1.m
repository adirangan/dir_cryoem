function Functional_SVD_1();
% Tests functional SVD. ;
% This time the test involves a non-separable function. ;

disp(sprintf(''));
disp(sprintf(' %% First we set up a non-separable 2-dimensional function:'));
disp(sprintf(' %% F(a,b) = besselj(1,a.*b); '));
F = @(a,b) besselj(2,a.*b);

a_N = 128; b_N = 129;
a_m = 30; a_r = 30; % [0,60] ;
b_m = 6/60; b_r = 6/60; % [0,12/60] ;
a_x_ = linspace(-1,+1,a_N);
a_t_ = linspace(a_m-a_r,a_m+a_r,a_N); 
b_x_ = linspace(-1,+1,b_N);
b_t_ = linspace(b_m-b_r,b_m+b_r,b_N); 
[A_x_,B_x_] = meshgrid(a_x_,b_x_);
[A_t_,B_t_] = meshgrid(a_t_,b_t_);
F_t_ = F(A_t_,B_t_);

a_w = @(a) ones(size(a)); a_w_ = [1 1]; % k ;
b_w = @(b) (1 + b).^2; b_w_ = [1 1]; % delta ;
A_w_x_ = a_w(A_x_); B_w_x_ = b_w(B_x_);

a_K = 16; b_K = 17;
[a_lx,a_lw,a_Lx,a_Lv] = orthopoly_node_weight_matrix_0(a_K,a_w_);
a_lt = a_lx*a_r + a_m;
[b_lx,b_lw,b_Lx,b_Lv] = orthopoly_node_weight_matrix_0(b_K,b_w_);
b_lt = b_lx*b_r + b_m;

[A_lt_,B_lt_] = meshgrid(a_lt,b_lt);
F_lt_ = F(A_lt_,B_lt_);

F_ = zeros(b_K,a_K);
lw = b_lw*transpose(a_lw);
for nkA=0:a_K-1;for nkB=0:b_K-1;
L_tmp = transpose(b_Lx(1+nkB,:))*a_Lx(1+nkA,:);
S_tmp = F_lt_.*L_tmp.*lw;
F_(1+nkB,1+nkA) = sum(S_tmp(:));
end;end;%for nkA=0:a_K-1;for nkB=0:b_K-1;

G_t_ = zeros(b_N,a_N);
for nkA=0:a_K-1;
a_tmp = polyval(a_Lv(1+nkA,:),A_x_);
for nkB=0:b_K-1;
b_tmp = polyval(b_Lv(1+nkB,:),B_x_);
S_tmp = a_tmp.*b_tmp;
G_t_ = G_t_ + S_tmp.*F_(1+nkB,1+nkA);
end;end;%for nkA=0:a_K-1;for nkB=0:b_K-1;

[U_,S_,V_] = svds(F_,5); S_ = diag(S_);
disp(sprintf(' %% Note that S_ = %0.4f,%0.4f,%0.4f,%0.4f,%0.4f',S_));
U1_ = U_(:,1); S1_ = S_(1); V1_ = V_(:,1); 
U2_ = U_(:,2); S2_ = S_(2); V2_ = V_(:,2); 
U3_ = U_(:,3); S3_ = S_(3); V3_ = V_(:,3); 
U4_ = U_(:,4); S4_ = S_(4); V4_ = V_(:,4); 

U1_t_ = zeros(1,b_N);
U2_t_ = zeros(1,b_N);
U3_t_ = zeros(1,b_N);
U4_t_ = zeros(1,b_N);
for nkB=0:b_K-1;
b_tmp = polyval(b_Lv(1+nkB,:),b_x_);
U1_t_ = U1_t_ + U1_(1+nkB)*b_tmp;
U2_t_ = U2_t_ + U2_(1+nkB)*b_tmp;
U3_t_ = U3_t_ + U3_(1+nkB)*b_tmp;
U4_t_ = U4_t_ + U4_(1+nkB)*b_tmp;
end;%for nkB=0:b_K-1;

V1_t_ = zeros(1,a_N);
V2_t_ = zeros(1,a_N);
V3_t_ = zeros(1,a_N);
V4_t_ = zeros(1,a_N);
for nkA=0:a_K-1;
a_tmp = polyval(a_Lv(1+nkA,:),a_x_);
V1_t_ = V1_t_ + V1_(1+nkA)*a_tmp;
V2_t_ = V2_t_ + V2_(1+nkA)*a_tmp;
V3_t_ = V3_t_ + V3_(1+nkA)*a_tmp;
V4_t_ = V4_t_ + V4_(1+nkA)*a_tmp;
end;%for nkA=0:a_K-1;

H1_t_ = (transpose(U1_t_)*V1_t_)*S1_;
H2_t_ = H1_t_ + (transpose(U2_t_)*V2_t_)*S2_;
H3_t_ = H2_t_ + (transpose(U3_t_)*V3_t_)*S3_;
H4_t_ = H3_t_ + (transpose(U4_t_)*V4_t_)*S4_;

disp_flag=1;
if disp_flag;
%subplot(2,3,1); surfl(A_t_,B_t_,F_t_); xlabel('a'); ylabel('b'); zlabel('f'); title('F original');
%subplot(2,3,2); surfl(A_t_,B_t_,H1_t_); xlabel('a'); ylabel('b'); zlabel('f'); title('H1 factorized');
%subplot(2,3,3); surfl(A_t_,B_t_,H2_t_); xlabel('a'); ylabel('b'); zlabel('f'); title('H2 factorized');
%subplot(2,3,4); surfl(A_t_,B_t_,H3_t_); xlabel('a'); ylabel('b'); zlabel('f'); title('H3 factorized');
%subplot(2,3,5); surfl(A_t_,B_t_,H3_t_); xlabel('a'); ylabel('b'); zlabel('f'); title('H4 factorized');
subplot(2,3,1); imagesc(F_t_); xlabel('a'); ylabel('b'); zlabel('f'); title('F original');
subplot(2,3,2); imagesc(H1_t_); xlabel('a'); ylabel('b'); zlabel('f'); title('H1 factorized');
subplot(2,3,3); imagesc(H2_t_); xlabel('a'); ylabel('b'); zlabel('f'); title('H2 factorized');
subplot(2,3,4); imagesc(H3_t_); xlabel('a'); ylabel('b'); zlabel('f'); title('H3 factorized');
subplot(2,3,5); imagesc(H4_t_); xlabel('a'); ylabel('b'); zlabel('f'); title('H4 factorized');
subplot(2,3,6); plot(log10(S_),'.-','MarkerSize',25); title('log10(S)');
end;%if disp_flag;
