function [W_,pp_,m2_,l2_] = wignert_lsq(n_l,delta,n_s) ;
% generates wigner-t matrices up to n_l;
% uses least-squares solve ;
% test with: ;
%{

  n_l = 20; 
  K = 100; 
  d_ = linspace(0,0.5/K,32); n_dk = length(d_); k = K;
  %delta = [0,0,1/K]; k_ = linspace(0,K,32); n_k = length(k_);
  W_ = zeros((n_l+1).^2,(n_l+1).^2,n_dk);
  pp_ = zeros((n_l+1).^2,n_dk);
  m2_ = zeros((n_l+1).^2,n_dk);
  l2_ = zeros((n_l+1).^2,n_dk);
  for ndk=1:n_dk;
  [W_tmp,pp_tmp,m2_tmp,l2_tmp] = wignert_lsq(n_l,[0,0,d_(ndk)*k]);
  W_(:,:,ndk) = W_tmp;
  pp_(:,ndk) = pp_tmp;
  m2_(:,ndk) = m2_tmp;
  l2_(:,ndk) = l2_tmp;
  end;% for ndk=1:n_dk;

  %{
  % examine t for m==0 as function of k ;
  n_l_cut = 10;
  m2_tmp = m2_(:,1); l2_tmp = l2_(:,1); m2_ij = find(m2_tmp==0 & l2_tmp<=n_l_cut); m2_length = length(m2_ij);
  W00 = zeros(m2_length,m2_length, ndk);
  for ndk=1:n_dk;
  W_tmp = W_(pp_(:,ndk)+1,pp_(:,ndk)+1,ndk); 
  W00(:,:,ndk) = W_tmp(m2_ij,m2_ij);
  end;%for ndk=1:n_dk;
  %}
  %{
  % examine t for various values of m as function of k ;
  n_s = 5;
  S_ = zeros(n_s,n_dk);
  n_l_cut = 10;
  m_val_ = -n_l_cut:n_l_cut;
  for nm=1:length(m_val_);
  m_val = m_val_(nm);
  m2_tmp = m2_(:,1); l2_tmp = l2_(:,1); m2_ij = find(m2_tmp==m_val & l2_tmp<=n_l_cut); m2_length = length(m2_ij);
  W_m = zeros(m2_length,m2_length, ndk);
  for ndk=1:n_dk;
  W_tmp = W_(pp_(:,ndk)+1,pp_(:,ndk)+1,ndk); 
  W_m(:,:,ndk) = W_tmp(m2_ij,m2_ij);
  end;%for ndk=1:n_dk;
  S_tmp = svds(reshape(W_m,m2_length*m2_length,n_dk),n_s);
  S_(1:length(S_tmp),nm) = S_tmp;
  end;% for nm=1:length(m_val_);
  %}

  % examine t for all values of m as function of delta*k ;
  clear WXX_;
  WXX_ = zeros(0,n_dk);
  n_l_cut = 10;
  m_val_ = -n_l_cut:n_l_cut;
  for nm=1:length(m_val_);
  m_val = m_val_(nm);
  m2_tmp = m2_(:,1); l2_tmp = l2_(:,1); m2_ij = find(m2_tmp==m_val & l2_tmp<=n_l_cut); m2_length = length(m2_ij);
  W_m = zeros(m2_length,m2_length, ndk);
  for ndk=1:n_dk;
  W_tmp = W_(pp_(:,ndk)+1,pp_(:,ndk)+1,ndk); 
  W_m(:,:,ndk) = W_tmp(m2_ij,m2_ij);
  end;%for ndk=1:n_dk;
  WXX_ = [ WXX_ ; reshape(W_m,m2_length*m2_length,n_dk) ];
  end;% for nm=1:length(m_val_);
  n_s = 6;
  [UXX_,SXX_,VXX_] = svds(WXX_,n_s);

  % estimate errors associated with WXX_ factorization ;
  n_l_cut = 10;
  m_val_ = -n_l_cut:n_l_cut;
  S_error_ = zeros(length(m_val_),n_dk,n_s);
  n_sum = 0;
  for nm=1:length(m_val_);
  m_val = m_val_(nm);
  m2_tmp = m2_(:,1); l2_tmp = l2_(:,1); m2_ij = find(m2_tmp==m_val & l2_tmp<=n_l_cut); m2_length = length(m2_ij);
  for ndk=1:n_dk;
  W_tmp = W_(pp_(:,ndk)+1,pp_(:,ndk)+1,ndk); 
  W_tmp = W_tmp(m2_ij,m2_ij);
  W_XXX = WXX_(n_sum + (1:m2_length.^2),ndk);
  W_XXX = reshape(W_XXX,m2_length,m2_length);
  U_XXX = UXX_(n_sum + (1:m2_length.^2),1:n_s);
  USV_X = zeros(m2_length,m2_length); 
  for ns=1:n_s; USV_X = USV_X + reshape(U_XXX(:,ns),m2_length,m2_length)*SXX_(ns,ns)*VXX_(ndk,ns); S_error_(nm,ndk,ns) = norm(W_XXX-USV_X); end;
  %disp(sprintf(' %% nm %d m_val %d ndk %d e0 %f e1 %f',nm,m_val,ndk,norm(W_XXX-W_tmp),norm(W_XXX-USV_X)));
  end;%for ndk=1:n_dk;
  n_sum = n_sum + m2_length.^2;
  end;% for nm=1:length(m_val_);
  figure; for ns=1:n_s; subplot(2,3,ns); imagesc(log10(squeeze(S_error_(:,:,ns))),[-5,0]); colorbar; end; title(sprintf('kd max %d',max(d_)*K));


  %}

verbose=0;

warning on verbose;
warning('error','MATLAB:rankDeficientMatrix');
warning('error','MATLAB:nearlySingularMatrix');

if nargin<3; n_s = 3; end; % n_s = oversampling parameter ;

n_lm = (n_l+1).^2;
l_ = []; m_ = [];
for nl=0:n_l;
l_ = [l_ , nl*ones(1,2*nl+1) ];
m_ = [m_ , [-nl:+nl] ];
end;%for nl=0:n_l;

W_ = zeros(n_lm,n_lm);
%theta_ = linspace( 0 , 2*pi , 2+2*nl );  % Azimuthal/Longitude/Circumferential ;
%phi_   = linspace( 0 ,   pi , 2+2*nl );  % Altitude /Latitude /Elevation ;

% try and obtain W_ matrix ;

continue_flag = 1; 
n_s_current = n_s;
while (continue_flag & n_s_current<4*n_s);

theta_ = sort(2*pi*rand(1,ceil(sqrt(n_s_current*n_lm))));
phi_   = sort(1*pi*rand(1,ceil(sqrt(n_s_current*n_lm))));
[THETA_,PHI_] = meshgrid(theta_,phi_);
n_x = length(theta_)*length(phi_);
Y_orig_ = zeros(n_lm,n_x);
Y_tran_ = zeros(n_lm,n_x);

for nl=0:n_l;
n_m = 1+2*nl;
Ylm_orig_ = ylm(nl,THETA_,PHI_);
Xn_ = sin(PHI_).*cos(THETA_);
Yn_ = sin(PHI_).*sin(THETA_);
Zn_ = cos(PHI_);
Ylm_tran_ = Ylm_orig_;
for nm=1:n_m;
Ylm_tran_{nm} = Ylm_tran_{nm}.*exp(+i*(delta(1)*Xn_+delta(2)*Yn_+delta(3)*Zn_));
end;%for nm=1:n_m;
for nm=1:n_m;
ix = nl.^2 + nm;
if (verbose); disp(sprintf(' %% nl %d nm %d ix %d',nl,nm,ix)); end;
Y_orig_(ix,:) = reshape(Ylm_orig_{nm},1,n_x);
Y_tran_(ix,:) = reshape(Ylm_tran_{nm},1,n_x);
end;%for nm=1:n_m;
end;% for nl=0:n_l;

try; 
W_ = transpose(Y_tran_ / Y_orig_) ;
%disp(sprintf(' %% success with: n_s_current %d',n_s_current));
continue_flag = 0;
catch;
%disp(sprintf(' %% n_s_current %d --> %d',n_s_current,n_s_current+1));
continue_flag = 1;
n_s_current = n_s_current + 1;
end;%try;

end;%while (continue_flag);

if (n_s_current>=4*n_s); 
disp(sprintf(' %% Warning! n_s_current %d / %d; could not find W_ in wignert_lsq.m',n_s_current,n_s)); 
end;%if (n_s_current>=2*n_s); 

Wn_ = W_; Wt_ = ctranspose(Wn_);
for nlm1=1:n_lm; for nlm2=1:n_lm;
if (abs(Wn_(nlm1,nlm2))< abs(Wt_(nlm1,nlm2))); W_(nlm1,nlm2) = Wn_(nlm1,nlm2); end;
if (abs(Wn_(nlm1,nlm2))>=abs(Wt_(nlm1,nlm2))); W_(nlm1,nlm2) = conj(Wt_(nlm1,nlm2)); end;
end;end;%for nlm1=1:n_lm; for nlm2=1:n_lm;

pp_ = [];
m2_ = [];
l2_ = [];
for nm=0:n_l;
l1_ = nm:n_l;
ll_ = l1_.*(l1_+1);
if (nm==0); pp_ = [pp_,ll_]; m2_ = [m2_,zeros(1,1+n_l)]; l2_ = [l2_,l1_]; end;
if (nm>0); 
%pp_tmp = [ll_-nm ; ll_+nm]; pp_tmp = reshape(pp_tmp,1,2*(n_l-nm+1));
%m2_tmp = [-nm*ones(1,n_l-nm+1) ; +nm*ones(1,n_l-nm+1)]; m2_tmp = reshape(m2_tmp,1,2*(n_l-nm+1));
%l2_tmp = [l1_ ; l1_]; l2_tmp = reshape(l2_tmp,1,2*(n_l-nm+1));
pp_tmp = [ll_-nm ; ll_+nm]; pp_tmp = reshape(transpose(pp_tmp),1,2*(n_l-nm+1));
m2_tmp = [-nm*ones(1,n_l-nm+1) ; +nm*ones(1,n_l-nm+1)]; m2_tmp = reshape(transpose(m2_tmp),1,2*(n_l-nm+1));
l2_tmp = [l1_ ; l1_]; l2_tmp = reshape(transpose(l2_tmp),1,2*(n_l-nm+1));
pp_ = [pp_,pp_tmp]; m2_ = [m2_,m2_tmp]; l2_ = [l2_,l2_tmp];
end;%if (nm>0); 
end;%for nm=0:n_l;

plot_flag=0;
if plot_flag;
figure;
subplot(1,2,1);
imagesc(abs(W_),[0,1]); colorbar;
set(gca,'XTick',1:n_lm,'XTickLabel',m_);
set(gca,'YTick',1:n_lm,'YTickLabel',m_);
title(sprintf('m,l'));
subplot(1,2,2);
imagesc(abs(W_(1+pp_,1+pp_)),[0,1]); colorbar;
set(gca,'XTick',1:n_lm,'XTickLabel',m2_);
set(gca,'YTick',1:n_lm,'YTickLabel',m2_);
title(sprintf('l,m'));
figure;
subplot(1,2,1);
imagesc(log10(abs(W_)),[-2,0]); colorbar;
set(gca,'XTick',1:n_lm,'XTickLabel',m_);
set(gca,'YTick',1:n_lm,'YTickLabel',m_);
title(sprintf('m,l'));
subplot(1,2,2);
imagesc(log10(abs(W_(1+pp_,1+pp_))),[-2,0]); colorbar;
set(gca,'XTick',1:n_lm,'XTickLabel',m2_);
set(gca,'YTick',1:n_lm,'YTickLabel',m2_);
title(sprintf('l,m'));
end;%if plot_flag;
