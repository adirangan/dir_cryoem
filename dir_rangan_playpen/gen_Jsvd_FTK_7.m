function [FTK] = gen_Jsvd_FTK_7(K_max,N_pixel,eps_target,l_max,a_K,b_K);
% Generating functional svd-expansion of various bessel-functions. ;
%
% Inputs: ;
% K_max: maximum K-value; K_target will be set to K_max. ;
% N_pixel: maximum pixel width (in wavelengths) allowed for square grid of displacements. ;
% eps_target: tolerance of svd-expansion. ;
% l_max: maximum order of bessel-functions in svd-expansion. ;
% a_K: degree of jacobi-polynomials used to represent r = 2*pi*|k|. ;
% b_K: degree of jacobi-polynomials used to represent d = |delta|. ;
% ;
% Outputs: a structure FTK with the following fields: ;
% n_svd_r: number of coefficients used to represent polynomials of r = 2*pi*k. ;
% svd_r_: vector of quadrature nodes associated with polynomials of r = 2*pi*k. ;
% svd_r_m: midpoint of interval for r = 2*pi*k. ;
% svd_r_c: half-width of interval for r = 2*pi*k. ;
% svd_r_w_: vector of quadrature weights for r = 2*pi*k. ;
% svd_r_Jv_: cell array of chebfuns for r = 2*pi*k. ; 
% n_svd_d: number of coefficients used to represent polynomials of d = delta. ;
% svd_d_: vector of quadrature nodes associated with polynomials of d = delta. ;
% svd_d_m: midpoint of interval for d = delta. ;
% svd_d_c: half-width of interval for d = delta. ;
% svd_d_w_: vector of quadrature weights for d = delta. ;
% svd_d_Jv_: cell array of chebfuns for d = delta. ; 
% n_svd_l: number of terms in svd-expansion. ;
% svd_l_: vector of bessel-orders for the terms in the svd-expansion. ;
% svd_U_d_chebcoef_: array of size b_K-x-n_svd_l storing the delta-side of the svd-expansion. ;
% svd_s_: vector of length n_svd_l storing the singular-values of the svd-expansion. ;
% svd_V_r_chebcoef_: array of size a_K-x-n_svd_l storing the r-side of the svd-expansion. ;
%     Note: if we define ij_tmp_ to be the list of indices for which svd_l_ = l_tmp,: ;
%           then we can use the svd-expansion to reconstruct the bessel-function of order l_tmp. ;
%           For example, try: ;
%           gen_Jsvd_FTK_7();
% svd_U_d_: array of size b_K-x-n_svd_l storing the delta-side of the svd-expansion evaluated at svd_d_. ;
% svd_V_r_: array of size a_K-x-n_svd_l storing the r-side of the svd-expansion evaluated at svd_r_. ;

if (nargin<3);
K_max = 48; N_pixel = 3.0; l_max = 24; n_r_degree = 31; n_d_degree = 33; verbose=1;
eps_target_ = 10.^[-1:-0.5:-4];
n_svd_l_ = zeros(length(eps_target_),1);
E_abs_ = zeros(1+l_max,length(eps_target_));
E_rel_ = zeros(1+l_max,length(eps_target_));
%%%%%%%%;
for neps=0:length(eps_target_)-1;
eps_target = eps_target_(1+neps);
disp(sprintf(' %% testing gen_Jsvd_FTK_7 with K_max %d N_pixel %0.2f eps_target %0.2f and radial weight function',K_max,N_pixel,eps_target));
[FTK] = gen_Jsvd_FTK_7(K_max,N_pixel,eps_target,l_max,n_r_degree,n_d_degree);
n_svd_l_(1+neps) = FTK.n_svd_l;
a_K = n_r_degree; b_K = n_d_degree;
if (verbose>1); figure(1);clf; end;
for l_tmp = 0:l_max;
ij_tmp_ = find(FTK.svd_l_==l_tmp)-1;
n_delta = FTK.n_svd_d; delta_value_ = reshape(FTK.svd_d_,1,n_delta);
n_r = FTK.n_svd_r; r_value_ = reshape(FTK.svd_r_,1,n_r);
J_ = besselj(l_tmp,transpose(delta_value_)*r_value_);
H_ = zeros(size(J_));
for ij=ij_tmp_;
U_d_ = zeros(1,length(delta_value_));
for nkB=0:b_K-1;
b_tmp = FTK.svd_d_Jv_{1+nkB}((delta_value_ - FTK.svd_d_m)/FTK.svd_d_c);
U_d_ = U_d_ + FTK.svd_U_d_chebcoef_(1+nkB,1+ij)*b_tmp;
end;% for nkB=0:b_K-1;
V_r_ = zeros(1,length(r_value_));
a_K = size(FTK.svd_V_r_chebcoef_,1);
for nkA=0:a_K-1;
a_tmp = FTK.svd_r_Jv_{1+nkA}((r_value_ - FTK.svd_r_m)/FTK.svd_r_c);
V_r_ = V_r_ + FTK.svd_V_r_chebcoef_(1+nkA,1+ij)*a_tmp;
end;%for nkA=0:a_K-1;
H_ = H_ + transpose(U_d_)*FTK.svd_s_(1+ij)*V_r_;
end;%for ij=ij_tmp_;
if (verbose>1 & l_tmp<=5);
subplot(6,3,1 + 3*l_tmp); colormap(colormap_beach(64)); imagesc(J_); set(gca,'Xtick',[],'Ytick',[]); axis square; title(sprintf('J%d',l_tmp));
subplot(6,3,2 + 3*l_tmp); colormap(colormap_beach(64)); imagesc(H_); set(gca,'Xtick',[],'Ytick',[]); axis square; title(sprintf('H%d',l_tmp));
subplot(6,3,3 + 3*l_tmp); colormap(colormap_beach(64)); imagesc(J_-H_,eps_target*[-1,1]); set(gca,'Xtick',[],'Ytick',[]); axis square;  title('J-H'); colorbar();  
suptitle(sprintf('max delta %0.3f, eps %0.3f',max(delta_value_),eps_target));
set(gcf,'Position',1+[0,0,768,1024]);
drawnow(); figbig;
end;%if (verbose>1 & l_tmp<=5);
tmp_1 = sum(sum(abs(J_-H_).*(FTK.svd_d_w_*transpose(FTK.svd_r_w_))));
tmp_2 = sum(sum(abs(H_).*(FTK.svd_d_w_*transpose(FTK.svd_r_w_))));
E_abs_(1+l_tmp,1+neps) = tmp_1; E_rel_(1+l_tmp,1+neps) = tmp_1./tmp_2;
end;%for l_tmp = 0:l_max;
end;%for neps=0:length(eps_target_)-1;
%%%%%%%%;
E_abs_(find(~isfinite(E_abs_)))=0;
E_rel_(find(~isfinite(E_rel_)))=0;
if (verbose);
figure(2);clf; 
subplot(1,3,1); plot(0:l_max,log10(E_rel_),'o-','LineWidth',4); xlim([0,l_max]);
xlabel('l');ylabel('log10((J-H)/H)');title('log10(relative error)');
subplot(1,3,2); plot(0:l_max,log10(E_abs_),'o-','LineWidth',4); xlim([0,l_max]);
xlabel('l');ylabel('log10((J-H)/1)');title('log10(absolute error)');
subplot(1,3,3); 
hold on;
plot(n_svd_l_,log10(mean(E_abs_,1)),'o-','LineWidth',4);
plot(n_svd_l_,log10(mean(E_rel_,1)),'x-','LineWidth',4);
hold off;
xlabel('n terms'); ylabel('log10(error)'); legend('E abs','E rel');
end;%if (verbose);
disp('returning'); return;
end;%if (nargin<3);

verbose=1;
ni = 4;
if (nargin<ni); l_max = 24; end; ni = ni+1;
if (nargin<ni); a_K = 32; end; ni = ni+1;
if (nargin<ni); b_K = 32; end; ni = ni+1;

if isempty(l_max); l_max = 24; end;
if isempty(a_K); a_K = 32; end;
if isempty(b_K); b_K = 32; end;
if (verbose>1); 
disp(sprintf(' %% [entering gen_Jsvd_FTK_7] N_pixel %0.2f l_max %d, a_K %d, b_K %d, eps_target %0.6f',N_pixel,l_max,a_K,b_K,eps_target)); 
end;%if (verbose>1);

%%%%%%%%;
K_target = K_max;
z_target = N_pixel*pi*sqrt(2);
D_target = z_target/(2*pi*K_target);
%%%%%%%%;
r_max = 2*pi*K_target;
d_max = D_target;
a_m = r_max/2; a_r = a_m;
b_m = d_max/2; b_r = b_m;
%%%%%%%%;
[a_jx,a_jw] = jacpts(a_K,0,1); a_jw=transpose(a_jw);
a_Jv_ = cell(1+a_K,1+a_K); a_Jx = zeros(a_K,a_K);
for nkA=0:a_K;
aj = jacpoly(nkA,0,1)*sqrt(nkA+1)/sqrt(2);
if (nkA<a_K); a_Jx(1+nkA,:) = aj(a_jx); end;
a_Jv_{1+nkA} = aj;
end;%for nkA=1:a_K;
a_jt = a_jx*a_r + a_m;
%%%%%%%%;
[b_jx,b_jw] = jacpts(b_K,0,1); b_jw=transpose(b_jw);
b_Jv_ = cell(1+b_K,1+b_K); b_Jx = zeros(b_K,b_K);
for nkB=0:b_K;
bj = jacpoly(nkB,0,1)*sqrt(nkB+1)/sqrt(2);
if (nkB<b_K); b_Jx(1+nkB,:) = bj(b_jx); end;
b_Jv_{1+nkB} = bj;
end;%for nkB=1:b_K;
b_jt = b_jx*b_r + b_m;
%%%%%%%%;
[A_jt_,B_jt_] = ndgrid(a_jt,b_jt);
%%%%%%%%;

clear S_l_ S_u_ S_s_ S_v_ ;
l=0; n_S=0; continue_flag=1;
while (continue_flag);
if (l==0); l_ = [0]; else l_ = [-l,+l]; end;
for l_tmp = l_;
F = @(a,b) besselj(l_tmp,a.*b);
F_jt_ = F(A_jt_,B_jt_);
F_ = zeros(a_K,b_K);
jw_ = a_jw*transpose(b_jw);
for nkA=0:a_K-1;for nkB=0:b_K-1;
J_tmp_ = transpose(a_Jx(1+nkA,:))*b_Jx(1+nkB,:);
S_tmp = F_jt_.*J_tmp_.*jw_;
F_(1+nkA,1+nkB) = sum(S_tmp(:));
end;end;%for nkA=0:a_K-1;for nkB=0:b_K-1;
n_svds = min([a_K,b_K]);
[U_,S_,V_] = svds(transpose(F_),n_svds); S_ = diag(S_); [ij_ret_] = find(S_>eps_target)-1 ;
if ~isempty(ij_ret_);
if (verbose>1); disp(sprintf(' %% l %+.2d, found %d terms [%0.2f,..,%0.2f];',l_tmp,length(1+ij_ret_),S_(1+ij_ret_(1)),S_(1+ij_ret_(end)))); end;%if
for ij = 0:length(ij_ret_)-1;
S_l_(1+n_S) = l_tmp;
S_u_(:,1+n_S) = U_(:,1+ij_ret_(1+ij));
S_s_(1,1+n_S) = S_(1+ij_ret_(1+ij),1);
S_v_(:,1+n_S) = V_(:,1+ij_ret_(1+ij));
n_S = n_S + 1;
end;%for ij = 0:length(ij_ret_)-1;
end;%if ~isempty(ij_ret_);
end;%for l_tmp = l_;
l=l+1;
if (l>l_max); continue_flag=0; else continue_flag=1; end;
end;%while (continue_flag);
if (verbose>1); disp(sprintf(' %% total of n_S %d terms found;',n_S)); end%if;

FTK = struct('type','FTK');
FTK.n_r_degree = a_K; FTK.n_d_degree = b_K;
FTK.n_svd_r = length(a_jt); FTK.svd_r_ = a_jt;
FTK.svd_r_m = a_m; FTK.svd_r_c = a_r; FTK.svd_r_w_ = a_jw; FTK.svd_r_Jv_ = a_Jv_;
FTK.n_svd_d = length(b_jt); FTK.svd_d_ = b_jt;
FTK.svd_d_m = b_m; FTK.svd_d_c = b_r; FTK.svd_d_w_ = b_jw; FTK.svd_d_Jv_ = b_Jv_;
FTK.n_svd_l = n_S;
FTK.svd_l_ = []; if (exist('S_l_','var')); FTK.svd_l_ = S_l_; end;
FTK.svd_U_d_chebcoef_ = []; if (exist('S_u_','var')); FTK.svd_U_d_chebcoef_ = S_u_; end;
FTK.svd_s_ = []; if (exist('S_s_','var')); FTK.svd_s_ = S_s_; end;
FTK.svd_V_r_chebcoef_ = []; if (exist('S_v_','var')); FTK.svd_V_r_chebcoef_ = S_v_; end;

svd_j_d_ = zeros(FTK.n_d_degree,FTK.n_svd_l);
svd_j_r_ = zeros(FTK.n_r_degree,FTK.n_svd_l);
for nl=0:FTK.n_svd_l-1;
%%%%%%%%;
tmp_p_j_ = zeros(FTK.n_d_degree,1);
for nd_degree=0:FTK.n_d_degree-1;
tmp_j_ = flip(poly(FTK.svd_d_Jv_{1+nd_degree}*FTK.svd_U_d_chebcoef_(1+nd_degree,1+nl)));
if (length(tmp_j_)>FTK.n_d_degree); disp(sprintf(' %% Warning! chebfun d-degree exceeds %d',FTK.n_d_degree)); end;
tmp_ij_ = 0 : min(FTK.n_d_degree,length(tmp_j_)) - 1;
tmp_p_j_(1+tmp_ij_) = tmp_p_j_(1+tmp_ij_) + transpose(tmp_j_(1+tmp_ij_));
end;% for nd_degree=0:FTK.n_d_degree-1;
svd_j_d_(:,1+nl) = tmp_p_j_;
%%%%%%%%;
tmp_p_j_ = zeros(FTK.n_r_degree,1);
for nr_degree=0:FTK.n_r_degree-1;
tmp_j_ = flip(poly(FTK.svd_r_Jv_{1+nr_degree}*FTK.svd_V_r_chebcoef_(1+nr_degree,1+nl)));
if (length(tmp_j_)>FTK.n_r_degree); disp(sprintf(' %% Warning! chebfun r-degree exceeds %d',FTK.n_r_degree)); end;
tmp_ij_ = 0 : min(FTK.n_r_degree,length(tmp_j_)) - 1;
tmp_p_j_(1+tmp_ij_) = tmp_p_j_(1+tmp_ij_) + transpose(tmp_j_(1+tmp_ij_));
end;%for nr_degree=0:FTK.n_r_degree-1;
svd_j_r_(:,1+nl) = tmp_p_j_;
%%%%%%%%;
end;%for nl=0:FTK.n_svd_l-1;

FTK.svd_U_d_ = svd_j_d_(:);
FTK.svd_V_r_ = svd_j_r_(:);
