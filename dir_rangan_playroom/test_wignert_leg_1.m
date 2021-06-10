% generate wignert_leg;
% testing to see if omega is separable from delta;
N_pixels = 1.5;
n_l = 20; 
K = 100; d = N_pixels/sqrt(2)/(K-1); k = K; omega_ = linspace(0,pi,32); n_omega = length(omega_); 
W_ = zeros((n_l+1).^2,(n_l+1).^2,n_omega);
pp_ = zeros((n_l+1).^2,n_omega);
m2_ = zeros((n_l+1).^2,n_omega);
l2_ = zeros((n_l+1).^2,n_omega);
for nomega=1:n_omega;
disp(sprintf(' %% nomega %d/%d',nomega,n_omega));
omega = omega_(nomega);
[W_tmp,pp_tmp,m2_tmp,l2_tmp] = wignert_leg(n_l,d*k*[sin(omega),0,cos(omega)]);
W_(:,:,nomega) = W_tmp;
pp_(:,nomega) = pp_tmp;
m2_(:,nomega) = m2_tmp;
l2_(:,nomega) = l2_tmp;
end;% for nomega=1:n_omega;

% examine s for all values of n,m as function of omega ;
clear WXX_;
pp_tmp = pp_(:,1); m2_tmp = m2_(:,1); l2_tmp = l2_(:,1); 
WXX_ = zeros(length(pp_tmp).^2,n_omega);
for nomega=1:n_omega;
W_tmp = W_(1+pp_tmp,1+pp_tmp,nomega);
WXX_(:,nomega) = reshape(W_tmp,length(pp_tmp).^2,1);
end;%for nomega=1:n_omega;
n_s = max(n_l,nomega);
[UXX_,SXX_,VXX_] = svds(WXX_,n_s);
disp(transpose(diag(SXX_)));

% examine s for all values of n,m as function of omega ;
clear WXX_;
WXX_ = zeros(0,n_omega);
n_l_cut = 10;
m_val_ = -n_l_cut:n_l_cut;
for nm=1:length(m_val_);
m_val = m_val_(nm);
m2_tmp = m2_(:,1); l2_tmp = l2_(:,1); m2_ij = find(m2_tmp==m_val & l2_tmp<=n_l_cut); m2_length = length(m2_ij);
W_m = zeros(m2_length,m2_length, nomega);
for nomega=1:n_omega;
W_tmp = W_(pp_(:,nomega)+1,pp_(:,nomega)+1,nomega); 
W_m(:,:,nomega) = W_tmp(m2_ij,m2_ij);
end;%for nomega=1:n_omega;
WXX_ = [ WXX_ ; reshape(W_m,m2_length*m2_length,n_omega) ];
end;% for nm=1:length(m_val_);
n_s = 6;
[UXX_,SXX_,VXX_] = svds(WXX_,n_s);

% estimate errors associated with WXX_ factorization ;
n_l_cut = 10;
m_val_ = -n_l_cut:n_l_cut;
S_error_ = zeros(length(m_val_),n_omega,n_s);
n_sum = 0;
for nm=1:length(m_val_);
m_val = m_val_(nm);
m2_tmp = m2_(:,1); l2_tmp = l2_(:,1); m2_ij = find(m2_tmp==m_val & l2_tmp<=n_l_cut); m2_length = length(m2_ij);
for nomega=1:n_omega;
W_tmp = W_(pp_(:,nomega)+1,pp_(:,nomega)+1,nomega); 
W_tmp = W_tmp(m2_ij,m2_ij);
W_XXX = WXX_(n_sum + (1:m2_length.^2),nomega);
W_XXX = reshape(W_XXX,m2_length,m2_length);
U_XXX = UXX_(n_sum + (1:m2_length.^2),1:n_s);
USV_X = zeros(m2_length,m2_length); 
for ns=1:n_s; USV_X = USV_X + reshape(U_XXX(:,ns),m2_length,m2_length)*SXX_(ns,ns)*VXX_(nomega,ns); S_error_(nm,nomega,ns) = norm(W_XXX-USV_X); end;
%disp(sprintf(' %% nm %d m_val %d nomega %d e0 %f e1 %f',nm,m_val,nomega,norm(W_XXX-W_tmp),norm(W_XXX-USV_X)));
end;%for nomega=1:n_omega;
n_sum = n_sum + m2_length.^2;
end;% for nm=1:length(m_val_);
figure; 
for ns=1:n_s; 
subplot(2,3,ns); 
S_tmp = squeeze(S_error_(:,:,ns));
imagesc(log10(S_tmp),[-5,0]); colorbar; 
title(sprintf('emax %f',max(S_tmp(:))));
end; %for ns=1:n_s; 

