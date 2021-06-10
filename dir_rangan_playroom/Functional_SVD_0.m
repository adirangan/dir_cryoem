function Functional_SVD_0();
% Tests functional SVD. ;

disp(sprintf(''));
disp(sprintf(' %% First we set up a separable 2-dimensional function:'));
disp(sprintf(' %% F(a,b) = sin(a).*(b.^2); '));
F = @(a,b) sin(a).*(b.^2);

a_N = 128; b_N = 129;
a_m = pi; a_r = pi;
b_m = 0.5; b_r = 1.5;
a_x_ = linspace(-1,+1,a_N);
a_t_ = linspace(a_m-a_r,a_m+a_r,a_N); 
b_x_ = linspace(-1,+1,b_N);
b_t_ = linspace(b_m-b_r,b_m+b_r,b_N); 
[A_x_,B_x_] = meshgrid(a_x_,b_x_);
[A_t_,B_t_] = meshgrid(a_t_,b_t_);
F_t_ = F(A_t_,B_t_);

a_w = @(a) ones(size(a)); a_w_ = [1]; 
b_w = @(b) (1 + b).^2; b_w_ = [1 2 1];
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

[U_,S_,V_] = svds(F_,3); S_ = diag(S_);
disp(sprintf(' %% Note that S_ = %0.4f,%0.4f,%0.4f',S_));
U_ = U_(:,1); S_ = S_(1,1); V_ = V_(:,1); 

U_t_ = zeros(1,b_N);
for nkB=0:b_K-1;
b_tmp = polyval(b_Lv(1+nkB,:),b_x_);
U_t_ = U_t_ + U_(1+nkB)*b_tmp;
end;%for nkB=0:b_K-1;

V_t_ = zeros(1,a_N);
for nkA=0:a_K-1;
a_tmp = polyval(a_Lv(1+nkA,:),a_x_);
V_t_ = V_t_ + V_(1+nkA)*a_tmp;
end;%for nkA=0:a_K-1;

H_t_ = (transpose(U_t_)*V_t_)*S_;

disp_flag=1;
if disp_flag;
subplot(1,3,1); surfl(A_t_,B_t_,F_t_); xlabel('a'); ylabel('b'); zlabel('f'); title('F original');
subplot(1,3,2); surfl(A_t_,B_t_,G_t_); xlabel('a'); ylabel('b'); zlabel('f'); title('G interpolant');
subplot(1,3,3); surfl(A_t_,B_t_,H_t_); xlabel('a'); ylabel('b'); zlabel('f'); title('H factorized');
end;%if disp_flag;
