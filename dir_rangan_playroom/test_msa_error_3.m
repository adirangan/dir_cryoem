%%%%%%%%;
% setting up simple multi-slice-alignment for testing errors. ;
% assume a 2d-volume x, comprising a small number of rings. ;
% Assume that x is restricted to just 2*Q+1 bessel-coefficients. ;
%%%%%%%%;

%clear;
%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); string_root = 'home'; end;
if (strcmp(platform,'eval1')); string_root = 'home'; end;
if (strcmp(platform,'rusty')); string_root = 'mnt/home'; end;
%%%%%%%%;

flag_verbose = 1;
flag_disp = flag_verbose;
flag_recalc = 0;
flag_replot = 1;
nf = 0;

n_k_p_r = 2; str_k_p_r = sprintf('k%d',n_k_p_r);
q_max = 16; str_q_max = sprintf('q%d',q_max);
n_M = 2048; str_n_M = sprintf('M%.4d',round(n_M));
nls_true = 5; str_nls_true = sprintf('nls%.3d',round(10*nls_true));
n_iteration = 128*1; str_n_iteration = sprintf('i%.3d',n_iteration);
n_rseed = 64; str_n_rseed = sprintf('r%.2d',n_rseed);
str_xfix = sprintf('%s%s%s%s%s%s',str_k_p_r,str_q_max,str_n_M,str_nls_true,str_n_iteration,str_n_rseed);
%%%%%%%%;
n_newt = 2;
k_p_r_ = transpose(1:n_k_p_r);
q_ = transpose(-q_max:1:+q_max);
n_q = numel(q_);
n_gamma = max(n_q*32,1024);
gamma_ = linspace(0,2*pi,1+n_gamma);
gamma_ = reshape(gamma_(1:n_gamma),[n_gamma,1]);
dgamma = mean(diff(gamma_));
n_w = n_gamma;
F_wq__ = exp(+i*gamma_*transpose(q_));
F_inv_qw__ = ctranspose(F_wq__)/n_gamma;
%%%%%%%%;

rseed_ = 0:n_rseed-1;
nrseed=0;
rseed = rseed_(1+nrseed);
rng(rseed);rseed=rseed+1;
gamma_ini_M_ = linspace(0,2*pi,1+n_M); gamma_ini_M_ = reshape(gamma_ini_M_(1:n_M),[n_M,1]);
eps_Mk__ = crandn(n_M,n_k_p_r);

L2_X00_qk__ = @(X00_qk__) sqrt(2*pi*sum(abs(X00_qk__).^2,'all'));
L2_X00_wk__ = @(X00_wk__) sqrt(dgamma*sum(abs(X00_wk__).^2,'all'));
L2_A00_qk__vs_X00_qki___ = @(A00_qk__,X00_qki___) reshape( sqrt(2*pi*sum(abs(bsxfun(@minus,A00_qk__,X00_qki___)).^2,[1,2])) , [size(X00_qki___,3),1] );

rng(rseed);rseed=rseed+1;
A00_qk__ = randn(n_q,n_k_p_r) + i*randn(n_q,n_k_p_r);
A00_qk__ = A00_qk__/L2_X00_qk__(A00_qk__);
Ap0_qk__ = bsxfun(@times,A00_qk__,i*reshape(q_,n_q,1));
App_qk__ = bsxfun(@times,Ap0_qk__,i*reshape(q_,n_q,1));

A00_wk__ = F_wq__*A00_qk__;
Ap0_wk__ = F_wq__*Ap0_qk__;
App_wk__ = F_wq__*App_qk__;

A00_qkr___(:,:,1+nrseed) = A00_qk__;
Ap0_qkr___(:,:,1+nrseed) = Ap0_qk__;
App_qkr___(:,:,1+nrseed) = App_qk__;

f_A00_w = @(gamma) sum(bsxfun(@times,A00_qk__,exp(+i*gamma*q_)),1);
f_Ap0_w = @(gamma) sum(bsxfun(@times,Ap0_qk__,exp(+i*gamma*q_)),1);
f_App_w = @(gamma) sum(bsxfun(@times,App_qk__,exp(+i*gamma*q_)),1);

f_A00_w_ = @(gamma_) exp(+i*gamma_(:)*transpose(q_))*A00_qk__;
f_Ap0_w_ = @(gamma_) exp(+i*gamma_(:)*transpose(q_))*Ap0_qk__;
f_App_w_ = @(gamma_) exp(+i*gamma_(:)*transpose(q_))*App_qk__;

f_R00_w_ = @(gamma_,M_k_) sum(abs(bsxfun(@minus,f_A00_w_(gamma_),M_k_)).^2,2);
f_Rp0_w_ = @(gamma_,M_k_) 2*real(sum(conj(bsxfun(@minus,f_A00_w_(gamma_),M_k_)).*f_Ap0_w_(gamma_),2));
f_Rpp_w_ = @(gamma_,M_k_) 2*real(sum(conj(f_Ap0_w_(gamma_)).*f_Ap0_w_(gamma_),2)) + 2*real(sum(conj(bsxfun(@minus,f_A00_w_(gamma_),M_k_)).*f_App_w_(gamma_),2));

f_X00_w_ = @(gamma_,X00_qk__) exp(+i*gamma_(:)*transpose(q_))*X00_qk__;
f_Xp0_w_ = @(gamma_,Xp0_qk__) exp(+i*gamma_(:)*transpose(q_))*Xp0_qk__;
f_Xpp_w_ = @(gamma_,Xpp_qk__) exp(+i*gamma_(:)*transpose(q_))*Xpp_qk__;

f_RX00_w_ = @(gamma_,M_k_,X00_qk__) sum(abs(bsxfun(@minus,f_X00_w_(gamma_),M_k_)).^2,2);
f_RXp0_w_ = @(gamma_,M_k_,X00_qk__,Xp0_qk__) 2*real(sum(conj(bsxfun(@minus,f_X00_w_(gamma_),M_k_)).*f_Xp0_w_(gamma_),2));
f_RXpp_w_ = @(gamma_,M_k_,X00_qk__,Xp0_qk__,Xpp_qk__) 2*real(sum(conj(f_Xp0_w_(gamma_)).*f_Xp0_w_(gamma_),2)) + 2*real(sum(conj(bsxfun(@minus,f_X00_w_(gamma_),M_k_)).*f_Xpp_w_(gamma_),2));

flag_check=0;
if flag_check;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
plot(gamma_,real(A00_wk__(:,1)),'k-',gamma_,real(Ap0_wk__(:,1)),'rx',gamma_,real(circshift(A00_wk__(:,1),-1) - circshift(A00_wk__(:,1),+1))/(2*dgamma),'go'); xlabel('gamma'); ylabel('A00');
subplot(1,2,2);
plot(gamma_,real(Ap0_wk__(:,1)),'k-',gamma_,real(App_wk__(:,1)),'rx',gamma_,real(circshift(Ap0_wk__(:,1),-1) - circshift(Ap0_wk__(:,1),+1))/(2*dgamma),'go'); xlabel('gamma'); ylabel('Ap0');
end;%if flag_check;

flag_check=0;
if flag_check;
ngamma = n_gamma/4;
gamma_pre = gamma_(1+ngamma);
sigma = 1.0;
M_k_ = f_A00_w(gamma_pre) + sigma*crandn(1,n_k_p_r);
tmp_R_w_ = sum(abs(bsxfun(@minus,A00_wk__,M_k_)).^2,2);
[~,ij_gamma_upd] = min(tmp_R_w_); index_gamma_upd = ij_gamma_upd - 1;
gamma_pos = gamma_(1+index_gamma_upd);
gamma_upd = gamma_pos + dgamma*quadratic_1d_interpolation_0(-tmp_R_w_(1+periodize(index_gamma_upd+[-1:+1],0,n_gamma)));
gamma_new = gamma_upd;
for nnewt=0:n_newt-1;
tmp_Rpp = f_Rpp_w_(gamma_new,M_k_); if abs(tmp_Rpp)<1e-16; tmp_Rpp = 1e-16*sign(tmp_Rpp); end;
gamma_new = gamma_new - f_Rp0_w_(gamma_new,M_k_)/tmp_Rpp;
end;%for nnewt=0:n_newt-1;
%%%%;
if (flag_verbose);
disp(sprintf(' %% f_R00_w_(gamma_pre): %+0.16f, f_Rp0_w_(gamma_pre): %+0.16f',f_R00_w_(gamma_pre,M_k_),f_Rp0_w_(gamma_pre,M_k_)));
disp(sprintf(' %% f_R00_w_(gamma_pos): %+0.16f, f_Rp0_w_(gamma_pos): %+0.16f',f_R00_w_(gamma_pos,M_k_),f_Rp0_w_(gamma_pos,M_k_)));
disp(sprintf(' %% f_R00_w_(gamma_upd): %+0.16f, f_Rp0_w_(gamma_upd): %+0.16f',f_R00_w_(gamma_upd,M_k_),f_Rp0_w_(gamma_upd,M_k_)));
disp(sprintf(' %% f_R00_w_(gamma_new): %+0.16f, f_Rp0_w_(gamma_new): %+0.16f',f_R00_w_(gamma_new,M_k_),f_Rp0_w_(gamma_new,M_k_)));
end;%if (flag_verbose);
figure(1+nf);nf=nf+1;clf;figbig;p_row=1;p_col=3;np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot(A00_wk__,'-');
plot(M_k_,'o');
plot(f_A00_w(gamma_pre),'x');
plot(f_A00_w(gamma_pos),'^');
plot(f_A00_w(gamma_upd),'s');
plot(f_A00_w(gamma_new),'h');
hold off;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot(gamma_,f_R00_w_(gamma_,M_k_),'k-');
plot(gamma_,f_Rp0_w_(gamma_,M_k_),'r-');
plot(gamma_,(circshift(f_R00_w_(gamma_,M_k_),-1) - circshift(f_R00_w_(gamma_,M_k_),+1))/(2*dgamma),'g-');
plot(gamma_pre,f_R00_w_(gamma_pre,M_k_),'ko');
plot(gamma_pos,f_R00_w_(gamma_pos,M_k_),'k^');
plot(gamma_upd,f_R00_w_(gamma_upd,M_k_),'ks');
plot(gamma_new,f_R00_w_(gamma_new,M_k_),'kh');
plot(gamma_pre,f_Rp0_w_(gamma_pre,M_k_),'ro');
plot(gamma_pos,f_Rp0_w_(gamma_pos,M_k_),'r^');
plot(gamma_upd,f_Rp0_w_(gamma_upd,M_k_),'rs');
plot(gamma_new,f_Rp0_w_(gamma_new,M_k_),'rh');
hold off;
xlim([0,2*pi]); xlabel('gamma'); ylabel('R00j');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot(gamma_,f_Rp0_w_(gamma_,M_k_),'k-');
plot(gamma_,f_Rpp_w_(gamma_,M_k_),'r-');
plot(gamma_,(circshift(f_Rp0_w_(gamma_,M_k_),-1) - circshift(f_Rp0_w_(gamma_,M_k_),+1))/(2*dgamma),'g-');
hold off;
xlim([0,2*pi]); xlabel('gamma'); ylabel('Rp0');
%%%%;
end;%if flag_check;

flag_check=1;
if flag_check;

%%%%%%%%;
% Check the terms involved in a single step. ;
%%%%%%%%;
sigma=1e-3;
F_ini_Mq__ = exp(+i*gamma_ini_M_*transpose(q_));
F_ini_inv_qM__ = pinv(F_ini_Mq__);
N_ini_kM__ = transpose(f_A00_w_(gamma_ini_M_));
n_iteration = 1;
L2_A00_vs_N00_i_ = zeros(n_iteration,1);
%L2_A00_vs_N00_predict = sqrt(2*pi)*norm(F_ini_inv_qM__)*sqrt(n_q)*sqrt(n_k_p_r)*(sigma);
L2_A00_vs_N00_predict = sqrt(2*pi)*sqrt(n_q/n_M)*sqrt(n_k_p_r)*(sigma);
L2_A00_vs_N00_new_i_ = zeros(n_iteration,1);
F_new_inv_norm_i_ = zeros(n_iteration,1);
L2_A00_vs_N00_new_predict_i_ = zeros(n_iteration,1);
gamma_new_Mi__ = zeros(n_M,n_iteration);
L2_A00_vs_N00_two_i_ = zeros(n_iteration,1);
F_two_inv_norm_i_ = zeros(n_iteration,1);
L2_A00_vs_N00_two_predict_i_ = zeros(n_iteration,1);
gamma_two_Mi__ = zeros(n_M,n_iteration);
%%%%%%%%;
for niteration=0:n_iteration-1;
rng(64*niteration)
N_kM__ = N_ini_kM__ + sigma*crandn(n_k_p_r,n_M);
N00_qk__ = F_ini_inv_qM__*transpose(N_kM__);
N00_wk__ = F_wq__*N00_qk__;
flag_disp=0;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
hold on;
plot(A00_wk__,'-','LineWidth',3);
plot(N00_wk__,'-','LineWidth',1);
hold off;
end;%if flag_disp;
L2_A00_vs_N00 = L2_A00_qk__vs_X00_qki___(A00_qk__,N00_qk__);
if (flag_verbose>1); disp(sprintf(' %% niteration %d/%d: A00_qk__ vs N00_qk__: %0.16f',niteration,n_iteration,L2_A00_vs_N00)); end;
L2_A00_vs_N00_i_(1+niteration) = L2_A00_vs_N00;
%%%%;
[ ...
 ~ ...
,N00_new_qk__ ...
,N00_new_wk__ ...
,gamma_new_M_ ...
,Rp0_new_M_ ...
,F_new_inv_qM__ ...
,F_new_Mq__ ...
] = ...
test_msa_update_3( ...
 [] ...	   
,A00_qk__ ...
,N_kM__ ...
);
gamma_new_Mi__(:,1+niteration) = gamma_new_M_;
F_new_inv_norm = norm(F_new_inv_qM__,2);
F_new_inv_norm_i_(1+niteration) = F_new_inv_norm;
L2_A00_vs_N00_new_predict = sqrt(2*pi)*F_new_inv_norm*sqrt(n_q)*sqrt(n_k_p_r)*(sigma);
L2_A00_vs_N00_new_predict_i_(1+niteration) = L2_A00_vs_N00_new_predict;
L2_A00_vs_N00_new = L2_A00_qk__vs_X00_qki___(A00_qk__,N00_new_qk__);
if (flag_verbose>1); disp(sprintf(' %% niteration %d/%d: A00_qk__ vs N00_new_qk__: %0.16f',niteration,n_iteration,L2_A00_vs_N00_new)); end;
L2_A00_vs_N00_new_i_(1+niteration) = L2_A00_vs_N00_new;
%%%%;
[ ...
 ~ ...
,N00_two_qk__ ...
,N00_two_wk__ ...
,gamma_two_M_ ...
,Rp0_two_M_ ...
,F_two_inv_qM__ ...
,F_two_Mq__ ...
] = ...
test_msa_update_3( ...
 [] ...	   
,N00_new_qk__ ...
,N_kM__ ...
);
gamma_two_Mi__(:,1+niteration) = gamma_two_M_;
F_two_inv_norm = norm(F_two_inv_qM__,2);
F_two_inv_norm_i_(1+niteration) = F_two_inv_norm;
L2_A00_vs_N00_two_predict = sqrt(2*pi)*F_two_inv_norm*sqrt(n_q)*sqrt(n_k_p_r)*(sigma);
L2_A00_vs_N00_two_predict_i_(1+niteration) = L2_A00_vs_N00_two_predict;
L2_A00_vs_N00_two = L2_A00_qk__vs_X00_qki___(A00_qk__,N00_two_qk__);
if (flag_verbose>1); disp(sprintf(' %% niteration %d/%d: A00_qk__ vs N00_two_qk__: %0.16f',niteration,n_iteration,L2_A00_vs_N00_two)); end;
L2_A00_vs_N00_two_i_(1+niteration) = L2_A00_vs_N00_two;
%%%%;
end;%for niteration=0:n_iteration-1;
if (flag_verbose); disp(sprintf(' %% mean(L2_A00_vs_N00_i_): %0.9f, L2_A00_vs_N00_predict: %0.9f',mean(L2_A00_vs_N00_i_),L2_A00_vs_N00_predict)); end;
if (flag_verbose); disp(sprintf(' %% mean(L2_A00_vs_N00_new_i_): %0.9f, mean(L2_A00_vs_N00_new_predict_i_): %0.9f',mean(L2_A00_vs_N00_new_i_),mean(L2_A00_vs_N00_new_predict_i_))); end;
if (flag_verbose); disp(sprintf(' %% mean(L2_A00_vs_N00_two_i_): %0.9f, mean(L2_A00_vs_N00_two_predict_i_): %0.9f',mean(L2_A00_vs_N00_two_i_),mean(L2_A00_vs_N00_two_predict_i_))); end;

%%%%%%%%;
% Measure typical angle between Ap0_wk__ and Np0_new_wk__. ;
%%%%%%%%;
Np0_new_qk__ = bsxfun(@times,N00_new_qk__,i*reshape(q_,n_q,1));
Np0_new_wk__ = F_wq__*Np0_new_qk__;
Np0_two_qk__ = bsxfun(@times,N00_two_qk__,i*reshape(q_,n_q,1));
Np0_two_wk__ = F_wq__*Np0_two_qk__;

tmp_A00_vs_N00_new = sqrt(2*pi*mean(sum(abs(A00_wk__ - N00_new_wk__).^2,2)));
tmp_A00_vs_N00_new = L2_A00_qk__vs_X00_qki___(A00_qk__,N00_new_qk__);
tmp_A00_vs_N00_new = L2_A00_vs_N00_new;
disp(sprintf(' %% tmp_A00_vs_N00_new/sqrt(2*pi): %0.9f',tmp_A00_vs_N00_new/sqrt(2*pi))); %<-- typical distance between points on A00_wk__ and N00_new_wk__. ;
tmp_A00_vs_N00_new_predict = L2_A00_vs_N00_new_predict;
disp(sprintf(' %% tmp_A00_vs_N00_new_predict/sqrt(2*pi): %0.9f',tmp_A00_vs_N00_new_predict/sqrt(2*pi))); %<-- predicted distance between A00_wk__ and N00_new_wk__. ;
tmp_Ap0_vs_Np0_new = sqrt(2*pi*mean(sum(abs(Ap0_wk__ - Np0_new_wk__).^2,2)));
tmp_Ap0_vs_Np0_new = L2_A00_qk__vs_X00_qki___(Ap0_qk__,Np0_new_qk__);
disp(sprintf(' %% tmp_Ap0_vs_Np0_new/sqrt(2*pi): %0.9f',tmp_Ap0_vs_Np0_new/sqrt(2*pi))); %<-- typical distance between points on Ap0_wk__ and Np0_new_wk__. ;


Ap0_unit_wk__ = bsxfun(@rdivide,Ap0_wk__,sqrt(sum(abs(Ap0_wk__).^2,2)));
Np0_new_unit_wk__ = bsxfun(@rdivide,Np0_new_wk__,sqrt(sum(abs(Np0_new_wk__).^2,2)));
Np0_two_unit_wk__ = bsxfun(@rdivide,Np0_two_wk__,sqrt(sum(abs(Np0_two_wk__).^2,2)));
theta_Ap0_vs_Np0_new_w_ = acos(real(sum(conj(Ap0_unit_wk__).*Np0_new_unit_wk__,2)));
theta_Ap0_vs_Np0_two_w_ = acos(real(sum(conj(Ap0_unit_wk__).*Np0_two_unit_wk__,2)));
disp(sprintf(' %% typical distance between A00_wk__ and N00_new_wk__: %0.9f <-- %0.9f',L2_A00_vs_N00_new/sqrt(2*pi),sqrt(mean(sum(abs(A00_wk__ - N00_new_wk__).^2,2)))));
disp(sprintf(' %% typical theta between Ap0_wk__ and Np0_new_wk__: %0.6f',mean(theta_Ap0_vs_Np0_new_w_)));
A00_at_new_Mk__ = f_X00_w_(gamma_new_M_,A00_qk__);
N00_new_at_new_Mk__ = f_X00_w_(gamma_new_M_,N00_new_qk__);
N00_new_at_two_Mk__ = f_X00_w_(gamma_two_M_,N00_new_qk__);
A00_at_new_vs_N00_new_at_new_Mk__ = A00_at_new_Mk__ - N00_new_at_new_Mk__;
N00_new_at_new_vs_N00_new_at_two_Mk__ = N00_new_at_two_Mk__ - N00_new_at_new_Mk__;
Ap0_at_new_Mk__ = f_Xp0_w_(gamma_new_M_,Ap0_qk__);
Np0_new_at_new_Mk__ = f_Xp0_w_(gamma_new_M_,Np0_new_qk__);
Np0_new_at_two_Mk__ = f_Xp0_w_(gamma_two_M_,Np0_new_qk__);
Ap0_unit_at_new_Mk__ = bsxfun(@rdivide,Ap0_at_new_Mk__,sqrt(sum(abs(Ap0_at_new_Mk__).^2,2)));
Np0_new_unit_at_new_Mk__ = bsxfun(@rdivide,Np0_new_at_new_Mk__,sqrt(sum(abs(Np0_new_at_new_Mk__).^2,2)));
Np0_new_unit_at_two_Mk__ = bsxfun(@rdivide,Np0_new_at_two_Mk__,sqrt(sum(abs(Np0_new_at_two_Mk__).^2,2)));
theta_Ap0_vs_Np0_new_at_new_Mk__ = acos(real(sum(conj(Ap0_unit_at_new_Mk__).*Np0_new_unit_at_new_Mk__,2)));
theta_Ap0_vs_Np0_new_at_two_Mk__ = acos(real(sum(conj(Ap0_unit_at_new_Mk__).*Np0_new_unit_at_two_Mk__,2)));
N_diff_Mk__ = transpose(N_kM__) - A00_at_new_Mk__;
plot(sqrt(sum(abs(N_diff_Mk__).^2,2)).*abs(sin(theta_Ap0_vs_Np0_new_at_new_Mk__)) , sqrt(sum(abs(N00_new_at_new_vs_N00_new_at_two_Mk__).^2,2)),'.')
disp(sprintf(' %% predicted N00_new_at_new_vs_N00_new_at_two_Mk__: %0.9f',mean(abs(N00_new_at_new_vs_N00_new_at_two_Mk__))));
disp(sprintf(' %% predicted distance between gamma_new_M_ and gamma_two_M_: %0.9f',mean(sin(theta_Ap0_vs_Np0_new_w_))))
disp(sprintf(' %% typical sqrt(sum(abs(A00_at_new_vs_N00_new_at_new_Mk__).^2,2): %0.9f',sqrt(mean(sum(abs(A00_at_new_vs_N00_new_at_new_Mk__).^2,2)))));
disp(sprintf(' %% typical sqrt(sum(abs(N00_new_at_new_vs_N00_new_at_two_Mk__).^2,2): %0.9f',sqrt(mean(sum(abs(N00_new_at_new_vs_N00_new_at_two_Mk__).^2,2)))));


end;%if flag_check;

disp('returning'); return;
											
flag_check=1;
if flag_check;
%%%%%%%%;
% Start at ground-truth and take one step. ;
%%%%%%%%;
sigma_ = exp(transpose(log(1e-12):0.5:0)); n_sigma = numel(sigma_);
L2_A00_s_ = zeros(n_sigma,1);
%%%%%%%%;
for nsigma=0:n_sigma-1;
sigma = sigma_(1+nsigma);
M_kM__ = zeros(n_k_p_r,n_M);
Rp0_M_ = zeros(n_M,1);
gamma_new_M_ = zeros(n_M,1);
tmp_gamma_M_ = gamma_ini_M_;
for nM=0:n_M-1;
tmp_gamma_M = tmp_gamma_M_(1+nM); gamma_pre = tmp_gamma_M;
M_k_ = f_A00_w(gamma_pre) + sigma*eps_Mk__(1+nM,:);
tmp_R_w_ = sum(abs(bsxfun(@minus,A00_wk__,M_k_)).^2,2);
[~,ij_gamma_upd] = min(tmp_R_w_); index_gamma_upd = ij_gamma_upd - 1;
gamma_pos = gamma_(1+index_gamma_upd);
gamma_upd = gamma_pos + dgamma*quadratic_1d_interpolation_0(-tmp_R_w_(1+periodize(index_gamma_upd+[-1:+1],0,n_gamma)));
gamma_new = gamma_upd;
for nnewt=0:n_newt-1;
tmp_Rpp = f_Rpp_w_(gamma_new,M_k_); if abs(tmp_Rpp)<1e-16; tmp_Rpp = 1e-16*sign(tmp_Rpp); end;
gamma_new = gamma_new - f_Rp0_w_(gamma_new,M_k_)/tmp_Rpp;
end;%for nnewt=0:n_newt-1;
M_kM__(:,1+nM) = M_k_;
Rp0_M_(1+nM) = f_Rp0_w_(gamma_new,M_k_);
gamma_new_M_(1+nM) = gamma_new;
end;%for nM=0:n_M-1;
%%%%;
tmp_F_Mq__ = exp(+i*gamma_new_M_*transpose(q_));
tmp_mat__ = pinv(tmp_F_Mq__);
C_zero_qk__ = tmp_mat__*transpose(M_kM__);
C_zero_wk__ = F_wq__*C_zero_qk__;
%%%%;
flag_disp=0;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 1; p_col = 3; np=0;
subplot(p_row,p_col,1+np);np=np+1;
plot(1:n_M,Rp0_M_,'.'); xlabel('nM'); ylabel('Rp0');
subplot(p_row,p_col,1+np);np=np+1;
plot(tmp_gamma_M_,gamma_new_M_,'.'); xlabel('gamma old'); ylabel('gamma new');
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot(A00_wk__,'-','LineWidth',3);
plot(C_zero_wk__,'-','LineWidth',1);
hold off;
end;%if flag_disp;
%%%%;
L2_A00 = L2_X00_wk__(A00_wk__ - C_zero_wk__);
L2_A00_s_(1+nsigma) = L2_A00;
if (flag_verbose);
disp(sprintf(' %% nsigma %d/%d sigma = exp(%+6.2f) <-- L2_A00 %0.16f',nsigma,n_sigma,log(sigma),L2_A00));
end;%if (flag_verbose);
end;%for nsigma=0:n_sigma-1;
%%%%%%%%;
end;%if flag_check;

flag_check=1;
if flag_check;
%%%%%%%%;
% Now take a few more steps. ;
%%%%%%%%;
n_iteration = 64;
C_zero_qki___ = zeros(n_q,n_k_p_r,n_iteration);
C_zero_wki___ = zeros(n_gamma,n_k_p_r,n_iteration);
gamma_new_Mi__ = zeros(n_M,n_iteration);
Rp0_Mi__ = zeros(n_M,n_iteration);
kappa_F_inv_i_ = zeros(n_iteration,1);
kappa_F_i_ = zeros(n_iteration,1);
tmp_C_zero_qk__ = A00_qk__;
for niteration=0:n_iteration-1;
if (flag_verbose); if (mod(niteration,16)==0); disp(sprintf(' %% niteration %d/%d',niteration,n_iteration)); end; end;
parameter = [];
[ ...
 parameter ...
,tmp_C_zero_qk__ ...
,tmp_C_zero_wk__ ...
,tmp_gamma_new_M_ ...
,tmp_Rp0_M_ ...
,tmp_F_inv_qM__ ...
,tmp_F_Mq__ ...
] = ...
test_msa_update_3( ...
 parameter ...	   
,tmp_C_zero_qk__ ...
,M_kM__ ...
);
C_zero_qki___(:,:,1+niteration) = tmp_C_zero_qk__;
C_zero_wki___(:,:,1+niteration) = tmp_C_zero_wk__;
gamma_new_Mi__(:,1+niteration) = tmp_gamma_new_M_;
Rp0_Mi__(:,1+niteration) = tmp_Rp0_M_;
kappa_F_inv_i_(1+niteration) = norm(tmp_F_inv_qM__,2);
kappa_F_i_(1+niteration) = norm(tmp_F_Mq__,2);
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
end;%if flag_check;

flag_check=1;
if flag_check;
%%%%%%%%;
% Now take multiple steps for each sigma. ;
%%%%%%%%;
n_iteration = 64;
sigma_ = exp(transpose(log(1e-12):0.5:0)); n_sigma = numel(sigma_);
L2_A00_is__ = zeros(n_iteration,n_sigma);
C_zero_qkis____ = zeros(n_q,n_k_p_r,n_iteration,n_sigma);
C_zero_wkis____ = zeros(n_gamma,n_k_p_r,n_iteration,n_sigma);
gamma_new_Mis___ = zeros(n_M,n_iteration,n_sigma);
Rp0_Mis___ = zeros(n_M,n_iteration,n_sigma);
kappa_F_inv_is__ = zeros(n_iteration,n_sigma);
kappa_F_is__ = zeros(n_iteration,n_sigma);
%%%%%%%%;
for nsigma=0:n_sigma-1;
sigma = sigma_(1+nsigma);
%%%%;
M_kM__ = zeros(n_k_p_r,n_M);
tmp_gamma_M_ = gamma_ini_M_;
for nM=0:n_M-1;
tmp_gamma_M = tmp_gamma_M_(1+nM); gamma_pre = tmp_gamma_M;
M_k_ = f_A00_w(gamma_pre) + sigma*eps_Mk__(1+nM,:);
M_kM__(:,1+nM) = M_k_;
end;%for nM=0:n_M-1;
%%%%;
C_zero_qki___ = zeros(n_q,n_k_p_r,n_iteration);
C_zero_wki___ = zeros(n_gamma,n_k_p_r,n_iteration);
gamma_new_Mi__ = zeros(n_M,n_iteration);
Rp0_Mi__ = zeros(n_M,n_iteration);
kappa_F_inv_i_ = zeros(n_iteration,1);
kappa_F_i_ = zeros(n_iteration,1);
tmp_C_zero_qk__ = A00_qk__;
for niteration=0:n_iteration-1;
if (flag_verbose); if (mod(niteration,16)==0); disp(sprintf(' %% nsigma %d/%d, niteration %d/%d',nsigma,n_sigma,niteration,n_iteration)); end; end;
parameter = [];
[ ...
 parameter ...
,tmp_C_zero_qk__ ...
,tmp_C_zero_wk__ ...
,tmp_gamma_new_M_ ...
,tmp_Rp0_M_ ...
,tmp_F_inv_qM__ ...
,tmp_F_Mq__ ...
] = ...
test_msa_update_3( ...
 parameter ...	   
,tmp_C_zero_qk__ ...
,M_kM__ ...
);
C_zero_qki___(:,:,1+niteration) = tmp_C_zero_qk__;
C_zero_wki___(:,:,1+niteration) = tmp_C_zero_wk__;
gamma_new_Mi__(:,1+niteration) = tmp_gamma_new_M_;
Rp0_Mi__(:,1+niteration) = tmp_Rp0_M_;
kappa_F_inv_i_(1+niteration) = norm(tmp_F_inv_qM__,2);
kappa_F_i_(1+niteration) = norm(tmp_F_Mq__,2);
end;%for niteration=0:n_iteration-1;
%%%%;
L2_A00_is__(:,1+nsigma) = L2_A00_qk__vs_X00_qki___(A00_qk__,C_zero_qki___);
C_zero_qkis____(:,:,:,1+nsigma) = C_zero_qki___;
C_zero_wkis____(:,:,:,1+nsigma) = C_zero_wki___;
gamma_new_Mis___(:,:,1+nsigma) = gamma_new_Mi__;
Rp0_Mis___(:,:,1+nsigma) = Rp0_Mi__;
kappa_F_inv_is__(:,1+nsigma) = kappa_F_inv_i_;
kappa_F_is__(:,1+nsigma) = kappa_F_i_;
%%%%%%%%;
end;%for nsigma=0:n_sigma-1;
%%%%%%%%;
end;%if flag_check;


%%%%;
plot(log(sigma_),log(L2_A00_is__(end,:)),'ro',log(sigma_),log(sigma_),'k-')
%%%%;

gamma_std_is__ = zeros(n_iteration,n_sigma);
for nsigma=0:n_sigma-1;
gamma_new_Mi__ = gamma_new_Mis___(:,:,1+nsigma);
gamma_new_M_ = gamma_new_Mi__(:,1);
gamma_err_iM__ = transpose(periodize(bsxfun(@minus,gamma_new_M_,gamma_new_Mi__),-pi,+pi));
gamma_std_i_ = std(gamma_err_iM__,1,2);
gamma_std_is__(:,1+nsigma) = gamma_std_i_;
end;%for nsigma=0:n_sigma-1;

%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
plot(gamma_err_iM__,'k-');
xlabel('iteration'); xlim([0,n_iteration+1]); ylabel('gamma error');
subplot(1,2,2);
hist(gamma_err_iM__(end,:),128);
title('hist');
%%%%;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

fname_pre = sprintf('./dir_pm_mat/test_MSA_sheres_%s',str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
rseed_ = 0:n_rseed-1;
Error_sher_l2_ir__ = zeros(n_iteration,n_rseed);
Error_zero_l2_ir__ = zeros(n_iteration,n_rseed);
Error_emax_l2_ir__ = zeros(n_iteration,n_rseed);
Error_alte_l2_ir__ = zeros(n_iteration,n_rseed);
Error_sher_k1_ir__ = zeros(n_iteration,n_rseed);
Error_zero_k1_ir__ = zeros(n_iteration,n_rseed);
Error_emax_k1_ir__ = zeros(n_iteration,n_rseed);
Error_alte_k1_ir__ = zeros(n_iteration,n_rseed);
A_qr__ = zeros(n_q,n_rseed);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nrseed=0:n_rseed-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
rseed = rseed_(1+nrseed);
%%%%;
rng(rseed);rseed=rseed+1;
A_q_ = randn(n_q,1) + i*randn(n_q,1);
A_qr__(:,1+nrseed) = A_q_;
A_w_ = F_wq__*A_q_;

%%%%;
rng(rseed);rseed=rseed+1;
gamma_true_M_ = 2*pi*rand(n_M,1);
tmp_index_ = knnsearch(gamma_,gamma_true_M_,'K',1)-1;
rng(rseed);rseed=rseed+1;
M_M_ = A_w_(1+tmp_index_) + sigma_true*(randn(n_M,1) + i*randn(n_M,1));
%%%%;
if flag_disp>1;
figure(1+nf);nf=nf+1;clf;figsml;
c_use__ = colormap('hsv'); n_c_use = size(c_use__,1);
linewidth_use = 4;
markersize_use = 4;
hold on;
plot(real(A_w_),imag(A_w_),'k-','LineWidth',linewidth_use);
for nM=0:n_M-1;
gamma_true = gamma_true_M_(1+nM);
nc_use = max(0,min(n_c_use-1,floor(n_c_use*gamma_true/(2*pi))));
M = M_M_(1+nM);
plot(real(M),imag(M),'ko','MarkerFaceColor',c_use__(1+nc_use,:),'MarkerSize',markersize_use);
end;%for nM=0:n_M-1;
hold off;
xlabel('real'); ylabel('imag'); grid on;
end;%if flag_disp>1;
%%%%;
sigma_sher = 1*sigma_true;
xlim_ = prctile(real(A_w_),[  0,100]); xlim_ = mean(xlim_) + 1.125*0.5*diff(xlim_)*[-1,+1];
ylim_ = prctile(imag(A_w_),[  0,100]); ylim_ = mean(ylim_) + 1.125*0.5*diff(ylim_)*[-1,+1];
rng(rseed);rseed=rseed+1;
lambda = 0.0;
rand_q_ = randn(n_q,1) + i*randn(n_q,1);
B_sher_q_ = lambda*A_q_ + (1-lambda)*fnorm(A_q_)*rand_q_;
B_sher_w_ = F_wq__*B_sher_q_;
B_zero_w_ = B_sher_w_;
B_emax_w_ = B_sher_w_;
B_alte_w_ = B_sher_w_;
%%%%;
if flag_disp>0;
nf_base = nf; nf=nf+1;
figure(1+nf_base);clf;set(gcf,'Position',1+[0,0,1024*2,768*2]);
c_use__ = colormap_80s; n_c_use = size(c_use__,1);
linewidth_big = 6;
linewidth_sml = 4;
markersize_use = 4;
for ns=0:4-1;
subplot(2,4,1+ns);
hold on;
plot(real(A_w_),imag(A_w_),'k-','LineWidth',linewidth_big);
xlabel('real'); ylabel('imag'); grid on;
end;%for ns=0:4-1;
end;%if flag_disp>0;
%%%%;
n_M_factor_use = ceil(n_M_factor);
Error_sher_l2_i_ = zeros(n_iteration,1);
Error_zero_l2_i_ = zeros(n_iteration,1);
Error_emax_l2_i_ = zeros(n_iteration,1);
Error_alte_l2_i_ = zeros(n_iteration,1);
Error_sher_k1_i_ = zeros(n_iteration,1);
Error_zero_k1_i_ = zeros(n_iteration,1);
Error_emax_k1_i_ = zeros(n_iteration,1);
Error_alte_k1_i_ = zeros(n_iteration,1);
for niteration=0:n_iteration-1;
%%%%;
R2_sher_wM__ = abs(bsxfun(@minus,reshape(B_sher_w_,[n_w,1]),reshape(M_M_,[1,n_M]))).^2;
R2_sher_wM__ = bsxfun(@minus,R2_sher_wM__,min(R2_sher_wM__,[],1));
p_sher_wM__ = exp(-R2_sher_wM__/(2*sigma_sher^2));
p_sher_wM__ = bsxfun(@rdivide,p_sher_wM__,max(1e-24,sum(p_sher_wM__,1)));
p_sher_wM__ = bsxfun(@rdivide,p_sher_wM__,max(1e-24,sum(p_sher_wM__,2)));
tmp_index_ = efind(sum(p_sher_wM__,2)<1e-24);
p_sher_wM__(1+tmp_index_,:) = 1/max(1,n_M);
p_sher_wM__ = bsxfun(@rdivide,p_sher_wM__,max(1e-24,sum(p_sher_wM__,2)));
C_sher_w_ = p_sher_wM__*M_M_;
C_sher_q_ = F_inv_qw__*C_sher_w_;
C_sher_w_ = F_wq__*C_sher_q_;
B_sher_q_ = C_sher_q_;
B_sher_w_ = C_sher_w_;
%%%%;
R2_zero_wM__ = abs(bsxfun(@minus,reshape(B_zero_w_,[n_w,1]),reshape(M_M_,[1,n_M]))).^2;
R2_zero_wM__ = bsxfun(@minus,R2_zero_wM__,min(R2_zero_wM__,[],1));
[~,tmp_ij_] = min(R2_zero_wM__,[],1); tmp_index_ = tmp_ij_-1;
p_zero_wM__ = sparse(1+tmp_index_,1:n_M,1,n_w,n_M);
p_zero_sum_w_ = sum(p_zero_wM__,2);
tmp_index_ = efind(p_zero_sum_w_>0);
%C_zero_q_ = F_inv_qw__(:,1+tmp_index_)*(bsxfun(@rdivide,p_zero_wM__(1+tmp_index_,:),p_zero_sum_w_(1+tmp_index_))*M_M_);
%C_zero_q_ = (transpose(p_zero_wM__(1+tmp_index_,:))*F_wq__(1+tmp_index_,:)) \ M_M_;
tmp_mat__ = pinv(transpose(p_zero_wM__(1+tmp_index_,:))*F_wq__(1+tmp_index_,:)); C_zero_q_ = tmp_mat__*M_M_;
C_zero_w_ = F_wq__*C_zero_q_;
B_zero_q_ = C_zero_q_;
B_zero_w_ = C_zero_w_;
%%%%;
R2_emax_wM__ = abs(bsxfun(@minus,reshape(B_emax_w_,[n_w,1]),reshape(M_M_,[1,n_M]))).^2;
R2_emax_wM__ = bsxfun(@minus,R2_emax_wM__,min(R2_emax_wM__,[],1));
[~,tmp_ij_wM__] = sort(R2_emax_wM__,2);
p_emax_wM__ = sparse([],[],[],n_w,n_M,0);
for nM_factor_use=0:n_M_factor_use-1;
tmp_index_ = tmp_ij_wM__(:,1+nM_factor_use)-1;
p_emax_wM__ = p_emax_wM__ + sparse(1:n_w,1+tmp_index_,1/n_M_factor_use,n_w,n_M);
end;%for nM_factor_use=0:n_M_factor_use-1;
C_emax_q_ = F_inv_qw__*(p_emax_wM__*M_M_);
B_emax_q_ = C_emax_q_;
B_emax_w_ = F_wq__*B_emax_q_;
%%%%;
R2_alte_wM__ = abs(bsxfun(@minus,reshape(B_alte_w_,[n_w,1]),reshape(M_M_,[1,n_M]))).^2;
R2_alte_wM__ = bsxfun(@minus,R2_alte_wM__,min(R2_alte_wM__,[],1));
if mod(niteration,2)==0;
[~,tmp_ij_] = min(R2_alte_wM__,[],1); tmp_index_ = tmp_ij_-1;
p_alte_wM__ = sparse(1+tmp_index_,1:n_M,1,n_w,n_M);
p_alte_sum_w_ = sum(p_alte_wM__,2);
tmp_index_ = efind(p_alte_sum_w_>0);
%C_alte_q_ = F_inv_qw__(:,1+tmp_index_)*(bsxfun(@rdivide,p_alte_wM__(1+tmp_index_,:),p_alte_sum_w_(1+tmp_index_))*M_M_);
%C_alte_q_ = (transpose(p_alte_wM__(1+tmp_index_,:))*F_wq__(1+tmp_index_,:)) \ M_M_;
tmp_mat__ = pinv(transpose(p_alte_wM__(1+tmp_index_,:))*F_wq__(1+tmp_index_,:)); C_alte_q_ = tmp_mat__*M_M_;
C_alte_w_ = F_wq__*C_alte_q_;
end;%if mod(niteration,2)==0;
if mod(niteration,2)==1;
[~,tmp_ij_wM__] = sort(R2_alte_wM__,2);
p_alte_wM__ = sparse([],[],[],n_w,n_M,0);
for nM_factor_use=0:n_M_factor_use-1;
tmp_index_ = tmp_ij_wM__(:,1+nM_factor_use)-1;
p_alte_wM__ = p_alte_wM__ + sparse(1:n_w,1+tmp_index_,1/n_M_factor_use,n_w,n_M);
end;%for nM_factor_use=0:n_M_factor_use-1;
C_alte_q_ = F_inv_qw__*(p_alte_wM__*M_M_);
end;%if mod(niteration,2)==1;
B_alte_q_ = C_alte_q_;
B_alte_w_ = F_wq__*B_alte_q_;
%%%%;
[Error_sher_l2,Error_sher_k1] = test_MSA_sheres_error_0(q_max,n_q,q_,F_wq__,n_w,gamma_,[],A_w_,A_q_,B_sher_w_,B_sher_q_);
Error_sher_l2_i_(1+niteration) = Error_sher_l2;
Error_sher_k1_i_(1+niteration) = Error_sher_k1;
[Error_zero_l2,Error_zero_k1] = test_MSA_sheres_error_0(q_max,n_q,q_,F_wq__,n_w,gamma_,[],A_w_,A_q_,B_zero_w_,B_zero_q_);
Error_zero_l2_i_(1+niteration) = Error_zero_l2;
Error_zero_k1_i_(1+niteration) = Error_zero_k1;
[Error_emax_l2,Error_emax_k1] = test_MSA_sheres_error_0(q_max,n_q,q_,F_wq__,n_w,gamma_,[],A_w_,A_q_,B_emax_w_,B_emax_q_);
Error_emax_l2_i_(1+niteration) = Error_emax_l2;
Error_emax_k1_i_(1+niteration) = Error_emax_k1;
[Error_alte_l2,Error_alte_k1] = test_MSA_sheres_error_0(q_max,n_q,q_,F_wq__,n_w,gamma_,[],A_w_,A_q_,B_alte_w_,B_alte_q_);
Error_alte_l2_i_(1+niteration) = Error_alte_l2;
Error_alte_k1_i_(1+niteration) = Error_alte_k1;
%%%%;
if (niteration==n_iteration-1) | (mod(niteration,32)==0);
if flag_disp>0;
figure(1+nf_base); hold on;
nc_use = max(0,min(n_c_use-1,floor(n_c_use*niteration/n_iteration)));
subplot(2,4,1);
plot(real(B_sher_w_),imag(B_sher_w_),'-','LineWidth',linewidth_sml,'Color',c_use__(1+nc_use,:));
xlim(xlim_); ylim(ylim_); grid on;
title(sprintf('sher'));
subplot(2,4,2);
plot(real(B_zero_w_),imag(B_zero_w_),'-','LineWidth',linewidth_sml,'Color',c_use__(1+nc_use,:));
xlim(xlim_); ylim(ylim_); grid on;
title(sprintf('zero'));
subplot(2,4,3);
plot(real(B_emax_w_),imag(B_emax_w_),'-','LineWidth',linewidth_sml,'Color',c_use__(1+nc_use,:));
xlim(xlim_); ylim(ylim_); grid on;
title(sprintf('emax'));
subplot(2,4,4);
plot(real(B_alte_w_),imag(B_alte_w_),'-','LineWidth',linewidth_sml,'Color',c_use__(1+nc_use,:));
xlim(xlim_); ylim(ylim_); grid on;
title(sprintf('alte'));
subplot(2,4,[5:6]); cla;
hold on;
plot(1+[0:niteration],log10(Error_sher_l2_i_(1+[0:niteration])),'ko');
plot(1+[0:niteration],log10(Error_zero_l2_i_(1+[0:niteration])),'bo');
plot(1+[0:niteration],log10(Error_emax_l2_i_(1+[0:niteration])),'go');
plot(1+[0:niteration],log10(Error_alte_l2_i_(1+[0:niteration])),'co');
hold off;
title('log10(Error_l2)','Interpreter','none');
xlim([1,n_iteration]);
grid on;
subplot(2,4,[7:8]); cla;
hold on;
plot(1+[0:niteration],log10(Error_sher_k1_i_(1+[0:niteration])),'kx');
plot(1+[0:niteration],log10(Error_zero_k1_i_(1+[0:niteration])),'bx');
plot(1+[0:niteration],log10(Error_emax_k1_i_(1+[0:niteration])),'gx');
plot(1+[0:niteration],log10(Error_alte_k1_i_(1+[0:niteration])),'cx');
hold off;
title('log10(Error_k1)','Interpreter','none');
xlim([1,n_iteration]);
grid on;
%%%%;
drawnow();
end;%if flag_disp>0;
end;%if (niteration==n_iteration-1) | (mod(niteration,32)==0);
end;%for niteration=0:n_iteration-1;
%%%%;
Error_sher_l2_ir__(:,1+nrseed) = Error_sher_l2_i_;
Error_zero_l2_ir__(:,1+nrseed) = Error_zero_l2_i_;
Error_emax_l2_ir__(:,1+nrseed) = Error_emax_l2_i_;
Error_alte_l2_ir__(:,1+nrseed) = Error_alte_l2_i_;
Error_sher_k1_ir__(:,1+nrseed) = Error_sher_k1_i_;
Error_zero_k1_ir__(:,1+nrseed) = Error_zero_k1_i_;
Error_emax_k1_ir__(:,1+nrseed) = Error_emax_k1_i_;
Error_alte_k1_ir__(:,1+nrseed) = Error_alte_k1_i_;
%%%%;
if flag_disp>0;
close(gcf);
end;%if flag_disp>0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nrseed=0:n_rseed-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
save(fname_mat ...
     ,'q_max' ...
     ,'n_M' ...
     ,'sigma_true' ...
     ,'n_iteration' ...
     ,'n_rseed' ...
     ,'str_xfix' ...
     ,'A_qr__' ...
     ,'Error_sher_l2_ir__' ...
     ,'Error_zero_l2_ir__' ...
     ,'Error_emax_l2_ir__' ...
     ,'Error_alte_l2_ir__' ...
     ,'Error_sher_k1_ir__' ...
     ,'Error_zero_k1_ir__' ...
     ,'Error_emax_k1_ir__' ...
     ,'Error_alte_k1_ir__' ...
     );
close_fname_tmp(fname_pre);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if  exist(fname_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
load(fname_mat);
fname_fig = sprintf('./dir_pm_jpg/test_MSA_sheres_%s',str_xfix);
if ( flag_replot | ~exist(sprintf('%s.jpg',fname_fig)) );
Error_sher_l2_r_ = Error_sher_l2_ir__(end,:);
Error_zero_l2_r_ = Error_zero_l2_ir__(end,:);
Error_emax_l2_r_ = Error_emax_l2_ir__(end,:);
Error_alte_l2_r_ = Error_alte_l2_ir__(end,:);
Error_sher_k1_r_ = Error_sher_k1_ir__(end,:);
Error_zero_k1_r_ = Error_zero_k1_ir__(end,:);
Error_emax_k1_r_ = Error_emax_k1_ir__(end,:);
Error_alte_k1_r_ = Error_alte_k1_ir__(end,:);
emax_l2 = max(max(Error_sher_l2_r_),max(Error_alte_l2_r_));
emax_k1 = max(max(Error_sher_k1_r_),max(Error_alte_k1_r_));
Error_zero_l2_r_ = min(emax_l2-1e-12,Error_zero_l2_r_);
Error_emax_l2_r_ = min(emax_l2-1e-12,Error_emax_l2_r_);
Error_alte_l2_r_ = min(emax_l2-1e-12,Error_alte_l2_r_);
Error_zero_k1_r_ = min(emax_k1-1e-12,Error_zero_k1_r_);
Error_emax_k1_r_ = min(emax_k1-1e-12,Error_emax_k1_r_);
Error_alte_k1_r_ = min(emax_k1-1e-12,Error_alte_k1_r_);
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 3; np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
plot([0,emax_l2],[0,emax_l2],'k-',Error_sher_l2_r_,Error_zero_l2_r_,'ko','MarkerFaceColor','c');
xlim([0,emax_l2]);ylim([0,emax_l2]);
axis square; grid on;
xlabel('Error_sher','Interpreter','none');
ylabel('Error_zero','Interpreter','none');
title('sher vs zero l2');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
plot([0,emax_l2],[0,emax_l2],'k-',Error_sher_l2_r_,Error_emax_l2_r_,'ko','MarkerFaceColor','c');
xlim([0,emax_l2]);ylim([0,emax_l2]);
axis square; grid on;
xlabel('Error_sher','Interpreter','none');
ylabel('Error_emax','Interpreter','none');
title('sher vs emax l2');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
plot([0,emax_l2],[0,emax_l2],'k-',Error_sher_l2_r_,Error_alte_l2_r_,'ko','MarkerFaceColor','c');
xlim([0,emax_l2]);ylim([0,emax_l2]);
axis square; grid on;
xlabel('Error_sher','Interpreter','none');
ylabel('Error_alte','Interpreter','none');
title('sher vs alte l2');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
plot([0,emax_k1],[0,emax_k1],'k-',Error_sher_k1_r_,Error_zero_k1_r_,'ko','MarkerFaceColor','c');
xlim([0,emax_k1]);ylim([0,emax_k1]);
axis square; grid on;
xlabel('Error_sher','Interpreter','none');
ylabel('Error_zero','Interpreter','none');
title('sher vs zero k1');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
plot([0,emax_k1],[0,emax_k1],'k-',Error_sher_k1_r_,Error_emax_k1_r_,'ko','MarkerFaceColor','c');
xlim([0,emax_k1]);ylim([0,emax_k1]);
axis square; grid on;
xlabel('Error_sher','Interpreter','none');
ylabel('Error_emax','Interpreter','none');
title('sher vs emax k1');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
plot([0,emax_k1],[0,emax_k1],'k-',Error_sher_k1_r_,Error_alte_k1_r_,'ko','MarkerFaceColor','c');
xlim([0,emax_k1]);ylim([0,emax_k1]);
axis square; grid on;
xlabel('Error_sher','Interpreter','none');
ylabel('Error_alte','Interpreter','none');
title('sher vs alte k1');
%%%%%%%%;
sgtitle(fname_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
end;%if ( flag_replot | ~exist(sprintf('%s.jpg',fname_fig)) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if  exist(fname_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;



%{

  flag_replot=1;
  fname_fig = sprintf('./dir_pm_jpg/test_MSA_sheres_callout_FIGA');
  if ( flag_replot | ~exist(sprintf('%s.jpg',fname_fig)) );
  figure(1);clf;set(gcf,'Position',1+[0,0,1024*2,768]);
  p_row = 1; p_col = 3; np=0;
  %%%%%%%%;
  q_max = 3; str_q_max = sprintf('q%d',q_max);
  n_M = 1024; str_n_M = sprintf('M%.4d',round(n_M));
  sigma_true = 2/8; str_sigma_true = sprintf('s%.3d',round(1000*sigma_true));
  n_iteration = 128*2; str_n_iteration = sprintf('i%.3d',n_iteration);
  n_rseed = 32; str_n_rseed = sprintf('r%.2d',n_rseed);
  str_xfix = sprintf('%s%s%s%s%s%s',str_q_max,str_n_M,str_sigma_true,str_n_iteration,str_n_rseed);
  fname_mat = sprintf('./dir_pm_mat/test_MSA_sheres_%s.mat',str_xfix);
  load(fname_mat);
  Error_sher_l2_r_ = Error_sher_l2_ir__(end,:);
  Error_zero_l2_r_ = Error_zero_l2_ir__(end,:);
  Error_emax_l2_r_ = Error_emax_l2_ir__(end,:);
  Error_alte_l2_r_ = Error_alte_l2_ir__(end,:);
  Error_sher_k1_r_ = Error_sher_k1_ir__(end,:);
  Error_zero_k1_r_ = Error_zero_k1_ir__(end,:);
  Error_emax_k1_r_ = Error_emax_k1_ir__(end,:);
  Error_alte_k1_r_ = Error_alte_k1_ir__(end,:);
  emax_l2 = max(max(Error_sher_l2_r_),max(Error_alte_l2_r_));
  emax_k1 = max(max(Error_sher_k1_r_),max(Error_alte_k1_r_));
  Error_zero_l2_r_ = min(emax_l2-1e-12,Error_zero_l2_r_);
  Error_emax_l2_r_ = min(emax_l2-1e-12,Error_emax_l2_r_);
  Error_alte_l2_r_ = min(emax_l2-1e-12,Error_alte_l2_r_);
  Error_zero_k1_r_ = min(emax_k1-1e-12,Error_zero_k1_r_);
  Error_emax_k1_r_ = min(emax_k1-1e-12,Error_emax_k1_r_);
  Error_alte_k1_r_ = min(emax_k1-1e-12,Error_alte_k1_r_);
  subplot(p_row,p_col,1+np);np=np+1;
  hold on;
  plot([0,emax_k1],[0,emax_k1],'k-',Error_sher_k1_r_,Error_emax_k1_r_,'ko','MarkerFaceColor','c');
  plot([0,emax_k1],[0,emax_k1],'k-',Error_sher_k1_r_,Error_alte_k1_r_,'ko','MarkerFaceColor','m');
  hold off;
  xlim([0,emax_k1]);ylim([0,emax_k1]);
  axis square; grid on;
  xlabel('Error_sher','Interpreter','none');
  ylabel('Error_alte','Interpreter','none');
  title(sprintf('%s',str_xfix),'Interpreter','none');
  %%%%%%%%;
  q_max = 4; str_q_max = sprintf('q%d',q_max);
  n_M = 128; str_n_M = sprintf('M%.4d',round(n_M));
  sigma_true = 1/8; str_sigma_true = sprintf('s%.3d',round(1000*sigma_true));
  n_iteration = 128*2; str_n_iteration = sprintf('i%.3d',n_iteration);
  n_rseed = 32; str_n_rseed = sprintf('r%.2d',n_rseed);
  str_xfix = sprintf('%s%s%s%s%s%s',str_q_max,str_n_M,str_sigma_true,str_n_iteration,str_n_rseed);
  fname_mat = sprintf('./dir_pm_mat/test_MSA_sheres_%s.mat',str_xfix);
  load(fname_mat);
  Error_sher_l2_r_ = Error_sher_l2_ir__(end,:);
  Error_zero_l2_r_ = Error_zero_l2_ir__(end,:);
  Error_emax_l2_r_ = Error_emax_l2_ir__(end,:);
  Error_alte_l2_r_ = Error_alte_l2_ir__(end,:);
  Error_sher_k1_r_ = Error_sher_k1_ir__(end,:);
  Error_zero_k1_r_ = Error_zero_k1_ir__(end,:);
  Error_emax_k1_r_ = Error_emax_k1_ir__(end,:);
  Error_alte_k1_r_ = Error_alte_k1_ir__(end,:);
  emax_l2 = max(max(Error_sher_l2_r_),max(Error_alte_l2_r_));
  emax_k1 = max(max(Error_sher_k1_r_),max(Error_alte_k1_r_));
  Error_zero_l2_r_ = min(emax_l2-1e-12,Error_zero_l2_r_);
  Error_emax_l2_r_ = min(emax_l2-1e-12,Error_emax_l2_r_);
  Error_alte_l2_r_ = min(emax_l2-1e-12,Error_alte_l2_r_);
  Error_zero_k1_r_ = min(emax_k1-1e-12,Error_zero_k1_r_);
  Error_emax_k1_r_ = min(emax_k1-1e-12,Error_emax_k1_r_);
  Error_alte_k1_r_ = min(emax_k1-1e-12,Error_alte_k1_r_);
  subplot(p_row,p_col,1+np);np=np+1;
  hold on;
  plot([0,emax_k1],[0,emax_k1],'k-',Error_sher_k1_r_,Error_emax_k1_r_,'ko','MarkerFaceColor','c');
  plot([0,emax_k1],[0,emax_k1],'k-',Error_sher_k1_r_,Error_alte_k1_r_,'ko','MarkerFaceColor','m');
  hold off;
  xlim([0,emax_k1]);ylim([0,emax_k1]);
  axis square; grid on;
  xlabel('Error_sher','Interpreter','none');
  ylabel('Error_alte','Interpreter','none');
  title(sprintf('%s',str_xfix),'Interpreter','none');
  %%%%%%%%;
  q_max = 8; str_q_max = sprintf('q%d',q_max);
  n_M = 32; str_n_M = sprintf('M%.4d',round(n_M));
  sigma_true = 1/10; str_sigma_true = sprintf('s%.3d',round(1000*sigma_true));
  n_iteration = 128*4; str_n_iteration = sprintf('i%.3d',n_iteration);
  n_rseed = 64; str_n_rseed = sprintf('r%.2d',n_rseed);
  str_xfix = sprintf('%s%s%s%s%s%s',str_q_max,str_n_M,str_sigma_true,str_n_iteration,str_n_rseed);
  fname_mat = sprintf('./dir_pm_mat/test_MSA_sheres_%s.mat',str_xfix);
  load(fname_mat);
  Error_sher_l2_r_ = Error_sher_l2_ir__(end,:);
  Error_zero_l2_r_ = Error_zero_l2_ir__(end,:);
  Error_emax_l2_r_ = Error_emax_l2_ir__(end,:);
  Error_alte_l2_r_ = Error_alte_l2_ir__(end,:);
  Error_sher_k1_r_ = Error_sher_k1_ir__(end,:);
  Error_zero_k1_r_ = Error_zero_k1_ir__(end,:);
  Error_emax_k1_r_ = Error_emax_k1_ir__(end,:);
  Error_alte_k1_r_ = Error_alte_k1_ir__(end,:);
  emax_l2 = max(max(Error_sher_l2_r_),max(Error_alte_l2_r_));
  emax_k1 = max(max(Error_sher_k1_r_),max(Error_alte_k1_r_));
  Error_zero_l2_r_ = min(emax_l2-1e-12,Error_zero_l2_r_);
  Error_emax_l2_r_ = min(emax_l2-1e-12,Error_emax_l2_r_);
  Error_alte_l2_r_ = min(emax_l2-1e-12,Error_alte_l2_r_);
  Error_zero_k1_r_ = min(emax_k1-1e-12,Error_zero_k1_r_);
  Error_emax_k1_r_ = min(emax_k1-1e-12,Error_emax_k1_r_);
  Error_alte_k1_r_ = min(emax_k1-1e-12,Error_alte_k1_r_);
  subplot(p_row,p_col,1+np);np=np+1;
  hold on;
  plot([0,emax_k1],[0,emax_k1],'k-',Error_sher_k1_r_,Error_emax_k1_r_,'ko','MarkerFaceColor','c');
  plot([0,emax_k1],[0,emax_k1],'k-',Error_sher_k1_r_,Error_alte_k1_r_,'ko','MarkerFaceColor','m');
  hold off;
  xlim([0,emax_k1]);ylim([0,emax_k1]);
  axis square; grid on;
  xlabel('Error_sher','Interpreter','none');
  ylabel('Error_alte','Interpreter','none');
  title(sprintf('%s',str_xfix),'Interpreter','none');
  %%%%%%%%;
  sgtitle('sher vs alte k1:','Interpreter','none');
  disp(sprintf(' %% writing %s',fname_fig));
  print('-depsc',sprintf('%s.eps',fname_fig));
  print('-djpeg',sprintf('%s.jpg',fname_fig));
  end;%if ( flag_replot | ~exist(sprintf('%s.jpg',fname_fig)) );

  %}
