function [ FA_x_c , X_1_ , X_2_ , X_3_ ]  = spharm_to_fig_0(n_k,k_,n_l_,a_,res,alpha_,delta_);
% plots the spherical harmonic representation passed in by modsph ;
% ;
% n_k = integer maximum k ;
% k_ = integer array of length n_k; k_(nk) = k_value for shell nk ;
% n_l_ = integer array of length n_k; n_l_(nk) = spherical harmonic order on shell nk; n_l_(nk) corresponds to n_lm_(nk) = (n_l_(nk)+1)^2 coefficients ;
% a_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% res = resolution to plot ;
% alpha_ = real vector of euler-angles to apply for rotation [ alpha , beta , gamma ];
% delta_ = real vector of displacements to apply [delta_x , delta_y , delta_z ];
% test with: ;
%{
  spharm_to_fig_0();
  %}

verbose=1;

if nargin<7; delta_ = 0*0.125*[1,2,3]; end;
if nargin<6; alpha_ = [0*pi/4,0*pi/4,0*pi/6]; end;
if nargin<5; res = 20; end;
if nargin<4;
isph_start_ = MDA_read_i4('./dir_mda6/isph_start_.mda');
nterms_sph_ = MDA_read_i4('./dir_mda6/nterms_sph_.mda');
modsph_A_ori_ = MDA_read_c16('./dir_mda6/modsph_A_ori_.mda');
n_k = length(isph_start_);
k_ = 1:n_k;
n_l_ = nterms_sph_;
a_ = modsph_A_ori_;
end;%if nargin<4;

n_lm_ = (n_l_+1).^2;

if (norm(alpha_)>0);
if (verbose); disp(sprintf(' %% rotating molecule by: [%0.2f %0.2f %0.2f]',+alpha_)); end;
if (verbose); disp(sprintf(' %% rotating coordinate_frame by: [%0.2f %0.2f %0.2f]',-alpha_(3),-alpha_(2),-alpha_(1))); end;
b_ = zeros(size(a_));
for nk=1:n_k;
n_l = n_l_(nk); n_lm = n_lm_(nk); ix_base = sum(n_lm_(1:nk-1));
a_k_ = a_(ix_base + (1:n_lm)); b_k_ = zeros(size(a_k_));
W_beta_ = wignerd_b(n_l,-alpha_(2));
for nl=0:n_l;
W_alpha = diag(exp(+i*[-nl:+nl]*-alpha_(3)));
W_gamma = diag(exp(+i*[-nl:+nl]*-alpha_(1)));
a_k_tmp = a_k_(1+nl*(nl+1) + (-nl:+nl));
a_k_tmp = reshape(a_k_tmp,2*nl+1,1);
b_k_tmp = W_gamma*W_beta_{1+nl}*W_alpha*a_k_tmp;
b_k_(1+nl*(nl+1) + (-nl:+nl)) = b_k_tmp;
end;%for nl=0:n_l;
b_(ix_base + (1:n_lm)) = b_k_;
end;%for nk=1:n_k;
a_ = b_;
end;%if (norm(alpha_)>0);

k_max = k_(n_k);
k_1_ = linspace(-k_max,+k_max,res); k_2_ = linspace(-k_max,+k_max,res); k_3_ = linspace(-k_max,+k_max,res);
[K_2_,K_1_,K_3_] = meshgrid(k_2_,k_1_,k_3_);

K_RAD_min = 1; K_RAD_max = n_k;
K_THETA_ = atan2(K_2_,K_1_);
K_RXY_ = max(1e-6,sqrt(K_2_.^2 + K_1_.^2));
%K_PHI_ = atan2(K_RXY_,K_3_);
K_RADIUS_ = min(K_RAD_max+1,max(1e-6,sqrt(K_3_.^2 + K_2_.^2 + K_1_.^2)));
K_PHI_ = acos(K_3_./K_RADIUS_);
K_RAD_PRE_ = min(K_RAD_max+1,max(K_RAD_min,floor(K_RADIUS_)));
K_RAD_POS_ = min(K_RAD_max+1,max(K_RAD_min,ceil(K_RADIUS_)));
DK_PRE_ = abs(K_RAD_PRE_ - K_RADIUS_);
DK_POS_ = abs(K_RAD_POS_ - K_RADIUS_);
DK_A_ = DK_PRE_;
DK_B_ = DK_POS_;
DK_C_ = DK_A_ + DK_B_;
DK_ALPHA_ = zeros(size(DK_A_)); DK_BETA_ = ones(size(DK_B_));
ij_tmp = find(DK_C_>0);
DK_ALPHA_(ij_tmp) = DK_A_(ij_tmp)./DK_C_(ij_tmp);
DK_BETA_(ij_tmp) = DK_B_(ij_tmp)./DK_C_(ij_tmp);
MOL_A_K_ = zeros(res,res,res);

%{
% uncomment to test gaussian blob in fourier space ;
c_ = 5*[1,2,3]; % center;
for nk1=1:length(k_1_); k_1 = k_1_(nk1);
for nk2=1:length(k_2_); k_2 = k_2_(nk2);
for nk3=1:length(k_3_); k_3 = k_3_(nk3);
r2_tmp = (k_1 - c_(1)).^2 + (k_2 - c_(2)).^2 + (k_3 - c_(3)).^2;
MOL_A_K_(nk1,nk2,nk3) = 1/sqrt(2*pi)^3 * exp(-r2_tmp/2);
end;end;end;%for nk1,nk2,nk3;
 %}

%{
% uncomment to test gaussian blob in real space ;
c_ = 0.2*[1,2,3]; % center;
for nk1=1:length(k_1_); k_1 = k_1_(nk1);
for nk2=1:length(k_2_); k_2 = k_2_(nk2);
for nk3=1:length(k_3_); k_3 = k_3_(nk3);
k2_tmp = (k_1).^2 + (k_2).^2 + (k_3).^2;
kc = k_1*c_(1) + k_2*c_(2) + k_3*c_(3);
sigma_tmp = 10.0;
MOL_A_K_(nk1,nk2,nk3) = exp(i*kc) * 1/sqrt(2*pi*sigma_tmp.^2)^3 * exp(-k2_tmp/2/sigma_tmp.^2);
end;end;end;%for nk1,nk2,nk3;
 %}

for nk=1:n_k;
k = k_(nk); 
K_RAD_PRE_IJ = find(K_RAD_PRE_==k); K_RAD_POS_IJ = find(K_RAD_POS_==k);
length_pre = length(K_RAD_PRE_IJ); length_pos = length(K_RAD_POS_IJ);
K_PHI_PRE_ = K_PHI_(K_RAD_PRE_IJ); K_PHI_POS_ = K_PHI_(K_RAD_POS_IJ);
K_THETA_PRE_ = transpose(K_THETA_(K_RAD_PRE_IJ)); K_THETA_POS_ = transpose(K_THETA_(K_RAD_POS_IJ));
DK_BETA_PRE_ = transpose(DK_BETA_(K_RAD_PRE_IJ)); DK_BETA_POS_ = transpose(DK_BETA_(K_RAD_POS_IJ));
DK_ALPHA_PRE_ = transpose(DK_ALPHA_(K_RAD_PRE_IJ)); DK_ALPHA_POS_ = transpose(DK_ALPHA_(K_RAD_POS_IJ));
n_l = n_l_(nk); n_lm = n_lm_(nk); ix_base = sum(n_lm_(1:nk-1));
a_k_ = a_(ix_base + (1:n_lm)); 
l_ = []; m_ = []; for nl=0:n_l; l_ = [l_ , nl*ones(1,2*nl+1) ]; m_ = [m_ , [-nl:+nl] ]; end;%for nl=0:n_l;
if (verbose>0); disp(sprintf(' %% nk %d k %d/%d: length_pre %d length_pos %d; n_l %d n_lm %d ix_base %d',nk,k,k_max,length_pre,length_pos,n_l,n_lm,ix_base)); end;
for nl=0:n_l;
l_val = nl;
if (verbose>2); disp(sprintf(' %% nk %d k %d/%d: nl %d l_val %d',nk,k,k_max,nl,l_val)); end;
if (length_pre>0); Llm_pre__=legendre(l_val,cos(K_PHI_PRE_),'unnorm'); end;
if (length_pos>0); Llm_pos__=legendre(l_val,cos(K_PHI_POS_),'unnorm'); end;
A_pre_ = zeros(1,length_pre); A_pos_ = zeros(1,length_pos);
for m_val = -l_val:+l_val;
ix = 1+l_val*(l_val+1)+m_val;
m_abs = abs(m_val);
if (length_pre>0); if (l_val>0); Llm_pre_ = squeeze(Llm_pre__(1+m_abs,:,:)); end; if (l_val==0); Llm_pre_ = Llm_pre__(:,:); end; end;
if (length_pos>0); if (l_val>0); Llm_pos_ = squeeze(Llm_pos__(1+m_abs,:,:)); end; if (l_val==0); Llm_pos_ = Llm_pos__(:,:); end; end;
a1=((2*l_val+1)/(4*pi));
a2=exp(lfactorial(l_val-m_abs) - lfactorial(l_val+m_abs));
c=sqrt(a1*a2);
s=(-1)^((m_val<0)*m_val); % needed to preserve condon-shortley phase. ;
%s=1; % original phase ;
if (length_pre>0); 
Ylm_pre_ = s*c*Llm_pre_.*exp(+i*m_val*K_THETA_PRE_); 
A_pre_ = A_pre_ + s*a_k_(ix)*Ylm_pre_.*DK_BETA_PRE_;
end;%if (length_pre>0); 
if (length_pos>0); 
Ylm_pos_ = s*c*Llm_pos_.*exp(+i*m_val*K_THETA_POS_); 
A_pos_ = A_pos_ + s*a_k_(ix)*Ylm_pos_.*DK_ALPHA_POS_;
end;%if (length_pos>0); 
end;%for m_val = -l_val:+l_val;
if (length_pre>0); 
MOL_A_K_(K_RAD_PRE_IJ) = MOL_A_K_(K_RAD_PRE_IJ) + transpose(A_pre_);
end;%if (length_pre>0); 
if (length_pos>0);
MOL_A_K_(K_RAD_POS_IJ) = MOL_A_K_(K_RAD_POS_IJ) + transpose(A_pos_);
end;%if (length_pos>0);
end;%for nl=0:n_l;
end;%for nk=1:n_k;

if (norm(delta_)>0);
if (verbose); disp(sprintf(' %% translating by [%0.2f %0.2f %0.2f]',delta_)); end;
for nk1=1:length(k_1_);
k_1 = k_1_(nk1);
for nk2=1:length(k_2_);
k_2 = k_2_(nk2);
for nk3=1:length(k_3_);
k_3 = k_3_(nk3);
kd = k_1*delta_(1) + k_2*delta_(2) + k_3*delta_(3);
MOL_A_K_(nk1,nk2,nk3) = MOL_A_K_(nk1,nk2,nk3) * exp(i*kd);
end;end;end;% for nk3,nk2,nk1;
end;%if (norm(delta_)>0);

FA_k_c = MOL_A_K_; FA_l_k_c = abs(FA_k_c); FA_w_k_c = angle(FA_k_c);

plot_flag=1;
if plot_flag;
% plot molecule_A in fourier space;
subplot(1,2,1);
cra = colormap('spring'); cra = cra(end:-1:1,:);  ncra = size(cra,1);
F_tmp = permute(abs(FA_k_c),[2,1,3]);
%F_tmp = abs(FA_k_c);
v_avg = mean(F_tmp(:)); v_std = std(F_tmp(:)); v_min = min(F_tmp(:)); v_max = max(F_tmp(:));
vlim = v_avg + 2.5*v_std*[-1,1];
v_ = linspace(vlim(1),vlim(2),7);
for nv=length(v_):-1:1;
v = v_(nv);
hpatch = patch(isosurface(K_2_,K_1_,K_3_,F_tmp,v)); isonormals(K_2_,K_1_,K_3_,F_tmp,hpatch);
nc = max(1,min(ncra,floor(ncra*nv/length(v_))));
hpatch.FaceColor = cra(nc,:);
hpatch.EdgeColor = 'none';
hpatch.FaceAlpha = (nv/length(v_)).^4 * 1.0;
xlim(k_max*[-1,1]);ylim(k_max*[-1,1]);zlim(k_max*[-1,1]);
xlabel('x');ylabel('y');zlabel('z');
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for nv=1:length(v_);
title(sprintf('FA_k, r_max %d res %d',k_max,res));
% end plot;
end;%if plot_flag;

sample_rate = res / k_max;
k_max = k_(n_k);
k_1_p_ = linspace(-sample_rate*k_max,+sample_rate*k_max,sample_rate*res); k_2_p_ = linspace(-sample_rate*k_max,+sample_rate*k_max,sample_rate*res); k_3_p_ = linspace(-sample_rate*k_max,+sample_rate*k_max,sample_rate*res);
[K_2_p_,K_1_p_,K_3_p_] = meshgrid(k_2_p_,k_1_p_,k_3_p_);
ij_tmp = (1:res) - floor(res/2) + floor(sample_rate*res/2);
FA_k_c_p = zeros(sample_rate*res,sample_rate*res,sample_rate*res);
FA_k_c_p(ij_tmp,ij_tmp,ij_tmp) = FA_k_c;

x_max = +1.0*sample_rate;
x_1_ = linspace(-x_max,+x_max,sample_rate*res); x_2_ = linspace(-x_max,+x_max,sample_rate*res); x_3_ = linspace(-x_max,+x_max,sample_rate*res);
[X_2_,X_1_,X_3_] = meshgrid(x_2_,x_1_,x_3_);

FA_x_c = real(recenter3(fftn(recenter3(FA_k_c_p))));
FA_x_c = permute(FA_x_c,[2,1,3]);
%FA_x_c = FA_x_c;

plot_flag=1;
if plot_flag;
% plot molecule_A in real space ;
subplot(1,2,2);
cra = colormap('spring'); cra = cra(end:-1:1,:);  ncra = size(cra,1);
F_tmp = FA_x_c(ij_tmp,ij_tmp,ij_tmp);
v_avg = mean(F_tmp(:)); v_std = std(F_tmp(:)); v_min = min(F_tmp(:)); v_max = max(F_tmp(:));
vlim = v_avg + 2.5*v_std*[-1,1];
v_ = linspace(vlim(1),vlim(2),7);
for nv=length(v_):-1:1;
v = v_(nv);
hpatch = patch(isosurface(X_2_,X_1_,X_3_,FA_x_c,v)); isonormals(X_2_,X_1_,X_3_,FA_x_c,hpatch);
nc = max(1,min(ncra,floor(ncra*nv/length(v_))));
hpatch.FaceColor = cra(nc,:);
hpatch.EdgeColor = 'none';
hpatch.FaceAlpha = (nv/length(v_)).^4 * 1.0;
xlim(1*[-1,1]);ylim(1*[-1,1]);zlim(1*[-1,1]);
xlabel('x');ylabel('y');zlabel('z');
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for nv=1:length(v_);
title(sprintf('FA_x, r_max %d',k_max));
% end plot;
end;%if plot_flag;


