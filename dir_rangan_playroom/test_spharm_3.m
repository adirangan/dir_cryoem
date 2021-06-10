function test_spharm_3(n_k,k_,n_l_,a_,b_,res);
% n_k = integer maximum k ;
% k_ = integer array of length n_k; k_(nk) = k_value for shell nk ;
% n_l_ = integer array of length n_k; n_l_(nk) = spherical harmonic order on shell nk; n_l_(nk) corresponds to n_lm_(nk) = (n_l_(nk)+1)^2 coefficients ;
% a_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% b_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% a_ corresponds to molecule_a, b_ to molecule_b ;
% res = resolution ;
% test with: ;
%{
  test_spharm_3();
  %}

verbose=1;

if nargin<6; res_ = 60; end;
if nargin<5;
res=20; 
isph_start_ = MDA_read_i4('./dir_mda6/isph_start_.mda');
nterms_sph_ = MDA_read_i4('./dir_mda6/nterms_sph_.mda');
modsph_A_ori_ = MDA_read_c16('./dir_mda6/modsph_A_ori_.mda');
modsph_B_ori_ = MDA_read_c16('./dir_mda6/modsph_B_ori_.mda');
n_k = length(isph_start_);
k_ = 1:n_k;
n_l_ = nterms_sph_;
n_lm_ = (n_l_+1).^2;
a_ = modsph_A_ori_;
b_ = modsph_B_ori_;
end;%if nargin<4;

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
MOL_B_K_ = zeros(res,res,res);

for nk=1:n_k;
k = k_(nk); 
K_RAD_PRE_IJ = find(K_RAD_PRE_==k); K_RAD_POS_IJ = find(K_RAD_POS_==k);
length_pre = length(K_RAD_PRE_IJ); length_pos = length(K_RAD_POS_IJ);
K_PHI_PRE_ = K_PHI_(K_RAD_PRE_IJ); K_PHI_POS_ = K_PHI_(K_RAD_POS_IJ);
K_THETA_PRE_ = transpose(K_THETA_(K_RAD_PRE_IJ)); K_THETA_POS_ = transpose(K_THETA_(K_RAD_POS_IJ));
DK_BETA_PRE_ = transpose(DK_BETA_(K_RAD_PRE_IJ)); DK_BETA_POS_ = transpose(DK_BETA_(K_RAD_POS_IJ));
DK_ALPHA_PRE_ = transpose(DK_ALPHA_(K_RAD_PRE_IJ)); DK_ALPHA_POS_ = transpose(DK_ALPHA_(K_RAD_POS_IJ));
n_l = n_l_(nk); n_lm = n_lm_(nk); ix_base = sum(n_lm_(1:nk-1));
a_k_ = a_(ix_base + (1:n_lm)); b_k_ = b_(ix_base + (1:n_lm));
l_ = []; m_ = []; for nl=0:n_l; l_ = [l_ , nl*ones(1,2*nl+1) ]; m_ = [m_ , [-nl:+nl] ]; end;%for nl=0:n_l;
if (verbose>0); disp(sprintf(' %% nk %d k %d/%d: length_pre %d length_pos %d; n_l %d n_lm %d ix_base %d',nk,k,k_max,length_pre,length_pos,n_l,n_lm,ix_base)); end;
for nl=0:n_l;
l_val = nl;
if (verbose>2); disp(sprintf(' %% nk %d k %d/%d: nl %d l_val %d',nk,k,k_max,nl,l_val)); end;
if (length_pre>0); Llm_pre__=legendre(l_val,cos(K_PHI_PRE_),'unnorm'); end;
if (length_pos>0); Llm_pos__=legendre(l_val,cos(K_PHI_POS_),'unnorm'); end;
A_pre_ = zeros(1,length_pre); A_pos_ = zeros(1,length_pos);
B_pre_ = zeros(1,length_pre); B_pos_ = zeros(1,length_pos);
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
B_pre_ = B_pre_ + s*b_k_(ix)*Ylm_pre_.*DK_BETA_PRE_;
end;%if (length_pre>0); 
if (length_pos>0); 
Ylm_pos_ = s*c*Llm_pos_.*exp(+i*m_val*K_THETA_POS_); 
A_pos_ = A_pos_ + s*a_k_(ix)*Ylm_pos_.*DK_ALPHA_POS_;
B_pos_ = B_pos_ + s*b_k_(ix)*Ylm_pos_.*DK_ALPHA_POS_;
end;%if (length_pos>0); 
end;%for m_val = -l_val:+l_val;
if (length_pre>0); 
MOL_A_K_(K_RAD_PRE_IJ) = MOL_A_K_(K_RAD_PRE_IJ) + transpose(A_pre_);
MOL_B_K_(K_RAD_PRE_IJ) = MOL_B_K_(K_RAD_PRE_IJ) + transpose(B_pre_);
end;%if (length_pre>0); 
if (length_pos>0);
MOL_A_K_(K_RAD_POS_IJ) = MOL_A_K_(K_RAD_POS_IJ) + transpose(A_pos_);
MOL_B_K_(K_RAD_POS_IJ) = MOL_B_K_(K_RAD_POS_IJ) + transpose(B_pos_);
end;%if (length_pos>0);
end;%for nl=0:n_l;
end;%for nk=1:n_k;

x_max = +1.0;
x_1_ = linspace(-x_max,+x_max,res); x_2_ = linspace(-x_max,+x_max,res); x_3_ = linspace(-x_max,+x_max,res);
[X_2_,X_1_,X_3_] = meshgrid(x_2_,x_1_,x_3_);

X_1_p_ = permute(X_1_,[2,1,3]);
X_2_p_ = permute(X_2_,[2,1,3]);
X_3_p_ = permute(X_3_,[2,1,3]);
FA_k_c = MOL_A_K_; FA_l_k_c = abs(FA_k_c); FA_w_k_c = angle(FA_k_c);
FA_x_c = real(recenter3(fftn(recenter3(FA_k_c))));
FA_x_c = permute(FA_x_c,[2,1,3]);
FB_k_c = MOL_B_K_; FB_l_k_c = abs(FB_k_c); FB_w_k_c = angle(FB_k_c);
FB_x_c = real(recenter3(fftn(recenter3(FB_k_c))));
FB_x_c = permute(FB_x_c,[2,1,3]);

save('tmp.mat');

cra = colormap('spring'); cra = cra(end:-1:1,:);  ncra = size(cra,1);
v_avg = mean(FB_x_c(:)); v_std = std(FB_x_c(:)); v_min = min(FB_x_c(:)); v_max = max(FB_x_c(:));
vlim = v_avg + 2.5*v_std*[-1,1];
v_ = linspace(vlim(1),vlim(2),7);
for nv=length(v_):-1:1;
v = v_(nv);
hpatch = patch(isosurface(X_1_p_,X_2_p_,X_3_p_,FB_x_c,v)); isonormals(X_1_p_,X_2_p_,X_3_p_,FB_x_c,hpatch);
nc = max(1,min(ncra,floor(ncra*nv/length(v_))));
hpatch.FaceColor = cra(nc,:);
hpatch.EdgeColor = 'none';
hpatch.FaceAlpha = (nv/length(v_)).^4 * 1.0;
%xlim(1*[-1,1]);ylim(1*[-1,1]);zlim(1*[-1,1]);
xlabel('x');ylabel('y');zlabel('z');
view([-65,20]); 
axis vis3d; %camlight left; lighting gouraud;
end;%for nv=1:length(v_);
title(sprintf('FB_x, r_max %d',k_max));

%{

  dirname = './dir_mda6';
  F2_x_c = MDA_read_c16(sprintf('%s/S_B_ori_.mda',dirname)); F2_x_c = real(permute(F2_x_c,[2,1,3]));
  X = MDA_read_r8(sprintf('%s/x_pre_.mda',dirname)); X = permute(X,[2,1,3]);
  Y = MDA_read_r8(sprintf('%s/y_pre_.mda',dirname)); Y = permute(Y,[2,1,3]);
  Z = MDA_read_r8(sprintf('%s/z_pre_.mda',dirname)); Z = permute(Z,[2,1,3]);
  F2_k_c = recenter3(fftn(recenter3(F2_x_c))); F2_l_k_c = abs(F2_k_c); F2w_k_c = angle(F2_k_c);

  %{
    cra = colormap('spring'); cra = cra(end:-1:1,:);  ncra = size(cra,1);
    v_avg = mean(F2_x_c(:)); v_std = std(F2_x_c(:)); v_min = min(F2_x_c(:)); v_max = max(F2_x_c(:));
    vlim = v_avg + 2.5*v_std*[-1,1];
    v_ = linspace(vlim(1),vlim(2),7);
    for nv=length(v_):-1:1;
    v = v_(nv);
    hpatch = patch(isosurface(X,Y,Z,F2_x_c,v)); isonormals(X,Y,Z,F2_x_c,hpatch);
    %hpatch = patch(isosurface(X_1_p_,X_2_p_,X_3_p_,F2_x_c,v)); isonormals(X_1_p_,X_2_p_,X_3_p_,F2_x_c,hpatch);
    nc = max(1,min(ncra,floor(ncra*nv/length(v_))));
    hpatch.FaceColor = cra(nc,:);
    hpatch.EdgeColor = 'none';
    hpatch.FaceAlpha = (nv/length(v_)).^4 * 1.0;
    xlim(1*[-1,1]);ylim(1*[-1,1]);zlim(1*[-1,1]);
    xlabel('x');ylabel('y');zlabel('z');
    view([-65,20]); 
    axis vis3d; %camlight left; lighting gouraud;
    end;%for nv=1:length(v_);
    title(sprintf('F2_x, r_max %d',k_max));
    %}

  pcols = 9;
  pz_ = round(linspace(1,60,pcols));
  for np=1:pcols;
  pz = pz_(np);
  subplot(2,pcols,np + 0*pcols); imagesc(squeeze(F2_l_k_c(:,:,pz))); title(sprintf('F2 %d',pz));
  subplot(2,pcols,np + 1*pcols); imagesc(squeeze(FA_l_k_c(:,:,pz))); title(sprintf('FA %d',pz));
  end;%for np=1:pcols;

  %}




return


