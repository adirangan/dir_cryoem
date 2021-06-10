function Functional_SVD_2(K_max,N_pixel);
% Tests functional SVD. ;
% This time the test involves the non-separable bessel function. ;
% We go out to N_pixel. ;

ni=1;
if (nargin<ni); K_max = 128; end; ni = ni+1; % Largest value of K (not including 2*pi). ;
if (nargin<ni); N_pixel = 3.0; end; ni = ni+1; % Maximum pixels for displacement. ;

disp(sprintf(''));
disp(sprintf(' %% First we set up a non-separable 2-dimensional function:'));
disp(sprintf(' %% F(a,b) = besselj(1,a.*b); '));
F = @(a,b) besselj(1,a.*b);

R_target = K_max-0.5; % Target maximum R (i.e., K). ;
z_target = N_pixel*pi*sqrt(2); % Target maximum product of 2*pi*R*D. ;
D_target = z_target/(2*pi*R_target); % Target maximum D. ;
r_max = 2*pi*R_target;
d_max = D_target;

a_N = 128; b_N = 129;
a_m = r_max/2; a_r = a_m;
b_m = d_max/2; b_r = b_m;
a_x_ = linspace(-1,+1,a_N);
a_t_ = linspace(a_m-a_r,a_m+a_r,a_N); 
b_x_ = linspace(-1,+1,b_N);
b_t_ = linspace(b_m-b_r,b_m+b_r,b_N); 
[A_x_,B_x_] = meshgrid(a_x_,b_x_);
[A_t_,B_t_] = meshgrid(a_t_,b_t_);
F_t_ = F(A_t_,B_t_);

a_w = @(a) 1+a; a_w_ = [1 1]; % k ;
b_w = @(b) 1+b; b_w_ = [1 1]; % delta ;
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

n_svd = min([b_K,a_K,10]);
[U_,S_,V_] = svds(F_,n_svd); S_ = diag(S_);
disp(sprintf(' %% Note that S_ = %s',num2str(transpose(S_))));
U__  = cell(n_svd,1);
S__  = cell(n_svd,1);
V__  = cell(n_svd,1);
for ns=1:n_svd;
U__{ns} = U_(:,ns); 
S__{ns} = S_(ns);
V__{ns} = V_(:,ns);
end;%for ns=1:n_svd;

U_t__ = cell(n_svd,1);
for ns=1:n_svd;
U_t__{ns} = zeros(1,b_N);
for nkB=0:b_K-1;
b_tmp = polyval(b_Lv(1+nkB,:),b_x_);
U_t__{ns} = U_t__{ns} + U__{ns}(1+nkB)*b_tmp;
end;%for nkB=0:b_K-1;
end;%for ns=1:n_svd;

V_t__ = cell(n_svd,1);
for ns=1:n_svd;
V_t__{ns} = zeros(1,a_N);
for nkA=0:a_K-1;
a_tmp = polyval(a_Lv(1+nkA,:),a_x_);
V_t__{ns} = V_t__{ns} + V__{ns}(1+nkA)*a_tmp;
end;%for nkA=0:a_K-1;
end;%for ns=1:n_svd;

H_t__ = cell(n_svd,1);
H_t__{1} = transpose(U_t__{1})*V_t__{1}*S__{1};
for ns=2:n_svd;
H_t__{ns} = H_t__{ns-1} + transpose(U_t__{ns})*V_t__{ns}*S__{ns};
end;%for ns=2:n_svd;

disp_flag=1;
if disp_flag;
prows = 3; pcols = ceil((n_svd + 2)/prows);
%colormap(colormap_beach(64));
colormap(colormap_pm(64));
subplot(prows,pcols,1); imagesc(F_t_); xlabel('a'); ylabel('b'); zlabel('f'); title('F original'); 
subplot(prows,pcols,2); plot(1:n_svd,log10(S_),'.-','MarkerSize',25); title('log10(S)');
for ns=1:n_svd;
subplot(prows,pcols,2+ns); imagesc(H_t__{ns}); xlabel('a'); ylabel('b'); zlabel('f'); title(sprintf('H%d factorized',ns),'Interpreter','none');
end;%for ns=1:n_svd;
end;%if disp_flag;

