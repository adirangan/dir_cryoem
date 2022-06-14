function test_msa_error_5(n_k_p_r,q_max,n_M);

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

%%%%%%%%%%%%%%%%;
if nargin<1;
%%%%%%%%%%%%%%%%;

dir_base = sprintf('/data/rangan/dir_cryoem/dir_rangan_playroom');
dir_pm_mat = sprintf('%s/dir_pm_mat',dir_base);
dir_pm_jpg = sprintf('%s/dir_pm_jpg',dir_base);
n_k_p_r_ = [2:2:8]; n_n_k_p_r = numel(n_k_p_r_);
q_max_ = 16:4:32; n_q_max = numel(q_max_);
n_M_ = 1024*[1,2,4,8]; n_n_M = numel(n_M_);
sigma_ = exp(floor(log(1e-4)):0.5:0); n_sigma = numel(sigma_);
L2_Ap0_skqM____ = zeros(n_sigma,n_n_k_p_r,n_q_max,n_n_M);
L2_App_skqM____ = zeros(n_sigma,n_n_k_p_r,n_q_max,n_n_M);
sigma_skqM____ = zeros(n_sigma,n_n_k_p_r,n_q_max,n_n_M);
n_k_p_r_skqM____ = zeros(n_sigma,n_n_k_p_r,n_q_max,n_n_M);
q_max_skqM____ = zeros(n_sigma,n_n_k_p_r,n_q_max,n_n_M);
n_M_skqM____ = zeros(n_sigma,n_n_k_p_r,n_q_max,n_n_M);
gamma_dif_fin_std_skqM____ = zeros(n_sigma,n_n_k_p_r,n_q_max,n_n_M);
gamma_dif_fin_std_qua_skqM____ = zeros(n_sigma,n_n_k_p_r,n_q_max,n_n_M);
%%%%%%%%;
for nn_k_p_r=0:n_n_k_p_r-1;
n_k_p_r = n_k_p_r_(1+nn_k_p_r);
for nq_max=0:n_q_max-1;
q_max = q_max_(1+nq_max);
for nn_M=0:n_n_M-1;
n_M = n_M_(1+nn_M);
%%%%%%%%;
str_k_p_r = sprintf('k%d',n_k_p_r);
str_q_max = sprintf('q%d',q_max);
str_n_M = sprintf('M%.4d',round(n_M));
str_xfix = sprintf('%s%s%s',str_k_p_r,str_q_max,str_n_M);
fname_pre = sprintf('%s/test_msa_error_%s',dir_pm_mat,str_xfix);
fname_mat = sprintf('%s.mat',fname_pre);
if ~exist(fname_mat,'file');
disp(sprintf(' %% %s not found, creating',fname_mat));
test_msa_error_5(n_k_p_r,q_max,n_M);
end;%if ~exist(fname_mat,'file');
if  exist(fname_mat,'file');
disp(sprintf(' %% %s found, not creating',fname_mat));
tmp_ = load(fname_mat);
L2_Ap0_skqM____(:,1+nn_k_p_r,1+nq_max,1+nn_M) = tmp_.L2_Ap0;
L2_App_skqM____(:,1+nn_k_p_r,1+nq_max,1+nn_M) = tmp_.L2_App;
sigma_skqM____(:,1+nn_k_p_r,1+nq_max,1+nn_M) = tmp_.sigma_;
n_k_p_r_skqM____(:,1+nn_k_p_r,1+nq_max,1+nn_M) =  tmp_.n_k_p_r;
q_max_skqM____(:,1+nn_k_p_r,1+nq_max,1+nn_M) = tmp_.q_max;
n_M_skqM____(:,1+nn_k_p_r,1+nq_max,1+nn_M) = tmp_.n_M;
gamma_dif_fin_std_skqM____(:,1+nn_k_p_r,1+nq_max,1+nn_M) = tmp_.gamma_dif_fin_std_s_;
gamma_dif_fin_std_qua_skqM____(:,1+nn_k_p_r,1+nq_max,1+nn_M) = tmp_.gamma_dif_fin_std_qua_s_;
clear tmp_;
end;%if  exist(fname_mat,'file');
%%%%%%%%;
end;%for nn_M=0:n_n_M-1;
end;%for nq_max=0:n_q_max-1;
end;%for nn_k_p_r=0:n_n_k_p_r-1;
%%%%%%%%;

fname_fig_pre = sprintf('%s/test_msa_error_5_FIGA',dir_pm_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
figure(1);clf;figbig;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
linewidth_use = 3;
markersize_sml = 4;
markersize_use = 6;
p_row = n_n_k_p_r; p_col = n_q_max; np=0;
hold on;
%%%%%%%%;
for nn_k_p_r=0:n_n_k_p_r-1;
n_k_p_r = n_k_p_r_(1+nn_k_p_r);
for nq_max=0:n_q_max-1;
q_max = q_max_(1+nq_max);
subplot(p_row,p_col,1+np);np=np+1;
%%%%;
hold on;
for nn_M=0:n_n_M-1;
n_M = n_M_(1+nn_M);
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*nn_M/n_n_M)));
tmp_x_ = log(sigma_skqM____(:,1+nn_k_p_r,1+nq_max,1+nn_M));
%tmp_f_ = log(L2_Ap0_skqM____(:,1+nn_k_p_r,1+nq_max,1+nn_M)) - 0.5*log(n_M_skqM____(:,1+nn_k_p_r,1+nq_max,1+nn_M));
tmp_y0_ = log(gamma_dif_fin_std_skqM____(:,1+nn_k_p_r,1+nq_max,1+nn_M));
tmp_n_M = n_M;
tmp_gamma_dif_fin_std_qua_s_ = zeros(n_sigma,1);
for nsigma=0:n_sigma-1;
tmp_sigma = sigma_(1+nsigma);
tmp_L2_Ap0 = L2_Ap0_skqM____(1+nsigma,1+nn_k_p_r,1+nq_max,1+nn_M);
tmp_L2_App = L2_App_skqM____(1+nsigma,1+nn_k_p_r,1+nq_max,1+nn_M);
tmp_b = -1; tmp_c = tmp_sigma/tmp_L2_Ap0*tmp_L2_Ap0/sqrt(tmp_n_M); tmp_a = (2*pi)*tmp_L2_App/tmp_L2_Ap0;
tmp_d = tmp_b.^2 - 4*tmp_a*tmp_c; if tmp_d< 0; tmp_r = 1; end; if (tmp_d>=0); tmp_r = (-tmp_b - sqrt(tmp_d))/(2*tmp_a); end;
tmp_gamma_dif_fin_std_qua_s_(1+nsigma) = tmp_r;
end;%for nsigma=0:n_sigma-1;
%tmp_y1_ = log(gamma_dif_fin_std_qua_skqM____(:,1+nn_k_p_r,1+nq_max,1+nn_M));
tmp_y1_ = log(tmp_gamma_dif_fin_std_qua_s_);
%plot(tmp_x_,tmp_x_,'k-','LineWidth',linewidth_use);
plot(tmp_x_,tmp_y0_,'ko-','Color',c_80s__(1+nc_80s,:),'MarkerFaceColor',c_80s__(1+nc_80s,:),'MarkerSize',markersize_sml);
plot(tmp_x_,tmp_y1_,'rd','MarkerFaceColor',c_80s__(1+nc_80s,:),'MarkerSize',markersize_use);
end;%for nn_M=0:n_n_M-1;
hold off;
xlabel('sigma','Interpreter','none');
ylabel('gamma_dif_std','Interpreter','none');
title(sprintf('n_k_p_r %d q_max %d',n_k_p_r,q_max),'Interpreter','none');
%%%%;
end;%for nq_max=0:n_q_max-1;
end;%for nn_k_p_r=0:n_n_k_p_r-1;
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
%%%%%%%%;
end;%if ~exist(fname_fig_jpg,'file');

disp('returning'); return;
%%%%%%%%%%%%%%%%;
end;%if nargin<1;
%%%%%%%%%%%%%%%%;

na=0;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); q_max=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;

flag_check=0;
if flag_check;
%%%%%%%%;
% First testing asymptotics for bayesian inference on circle. ;
%%%%%%%%;
clear;
%h_ = @(r,sigma) sqrt(2*pi) * sigma/sqrt(1+r) * exp(-r.^2./(2*sigma^2));
h_ = @(r,sigma) sqrt(2*pi) * sigma/sqrt(1+r) * exp(-r.^2./(2*sigma^2));
sigma = 1e-3;
r_ = 1e-3*transpose([-10:0.01:+10]); n_r = numel(r_);
lI_h_ = zeros(n_r,1);
lI_g_ = zeros(n_r,1);
for nr=0:n_r-1;
r = r_(1+nr);
f_ = @(theta) exp(-( 2 + 2*r + r.^2 - 2*cos(theta).*(1+r))./(2*sigma^2));
%f_ = @(theta) exp(-( theta.^2.*(1+r) + r.^2 )./(2*sigma^2));
g = integral(f_,-pi,+pi);
h = h_(r,sigma);
lI_h_(1+nr) = log(h);
lI_g_(1+nr) = log(g);
end;%for nr=0:n_r-1;
plot(r_,lI_h_,'k-',r_,lI_g_,'r-');
%%%%%%%%;
r_ = -1:1e-2:+1;
plot(r_,sqrt(1+r_),'k-',r_,1 + 0.5*r_ - 0.125*r_.^2,'r-');
plot(r_,1./sqrt(1+r_),'k-',r_,1 - 0.5*r_ + 0.375*r_.^2,'r-');
plot(r_,(1+r_).^(3/2),'k-',r_,1 + (3/2)*r_ + 0.375*r_.^2,'r-');
%%%%%%%%;
sigma = 1e-1;
f_ = @(r) sqrt(2*pi)*sigma.*sqrt(1+r).*exp(-r.^2./(2*sigma.^2));
lI_f = log(integral(f_,-Inf,+Inf));
lI_h = log(2*pi*sigma.^2*(1 - sigma.^2/8));
disp(sprintf(' lI_f vs lI_h: %0.16f',fnorm(lI_f - lI_h)/fnorm(lI_h)));
%%%%%%%%;
end;%if flag_check;

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

if isempty(n_k_p_r); n_k_p_r = 4; end; str_k_p_r = sprintf('k%d',n_k_p_r);
if isempty(q_max); q_max = 24; end; str_q_max = sprintf('q%d',q_max);
if isempty(n_M); n_M = 1024*2; end; str_n_M = sprintf('M%.4d',round(n_M));
n_iteration = 32; n_trial = 2;
str_xfix = sprintf('%s%s%s',str_k_p_r,str_q_max,str_n_M);
%%%%%%%%;

dir_base = sprintf('/data/rangan/dir_cryoem/dir_rangan_playroom');
dir_pm_mat = sprintf('%s/dir_pm_mat',dir_base);
fname_pre = sprintf('%s/test_msa_error_%s',dir_pm_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

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

rseed = 0; rng(rseed);
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

%%%%%%%%;
% Check the terms involved in a single step. ;
%%%%%%%%;
F_ini_Mq__ = exp(+i*gamma_ini_M_*transpose(q_));
F_ini_inv_qM__ = pinv(F_ini_Mq__);
N_ini_kM__ = transpose(f_A00_w_(gamma_ini_M_));
%%%%%%%%;
sigma_ = exp(floor(log(1e-4)):0.5:0); n_sigma = numel(sigma_);
%%%%%%%%;
gamma_Mits____ = zeros(n_M,n_iteration,n_trial,n_sigma);
F_inv_norm_its___ = zeros(n_iteration,n_trial,n_sigma);
L2_A00_vs_N00_predict_its___ = zeros(n_iteration,n_trial,n_sigma);
L2_A00_vs_N00_its___ = zeros(n_iteration,n_trial,n_sigma);
gamma_dif_fin_Mts___ = zeros(n_M,n_trial,n_sigma);
%%%%%%%%;
for nsigma=0:n_sigma-1;
sigma = sigma_(1+nsigma);
L2_A00_vs_N00_it__ = zeros(n_iteration,n_trial);
F_inv_norm_it__ = zeros(n_iteration,n_trial);
L2_A00_vs_N00_predict_it__ = zeros(n_iteration,n_trial);
gamma_Mit___ = zeros(n_M,n_iteration,n_trial);
%%%%%%%%;
for ntrial=0:n_trial-1;
rng(64*ntrial)
N_kM__ = N_ini_kM__ + sigma*crandn(n_k_p_r,n_M);
%%%%;
N00_qk__ = A00_qk__;
for niteration=0:n_iteration-1;
if (flag_verbose); if (mod(niteration,16)==0); disp(sprintf(' %% nsigma %d/%d ntrial %d/%d niteration %d/%d',nsigma,n_sigma,ntrial,n_trial,niteration,n_iteration)); end; end;
[ ...
 ~ ...
,N00_qk__ ...
,N00_wk__ ...
,gamma_M_ ...
,Rp0_M_ ...
,F_inv_qM__ ...
,F_Mq__ ...
] = ...
test_msa_update_3( ...
 [] ...	   
,N00_qk__ ...
,N_kM__ ...
);
gamma_Mit___(:,1+niteration,1+ntrial) = gamma_M_;
F_inv_norm = norm(F_inv_qM__,2);
F_inv_norm_it__(1+niteration,1+ntrial) = F_inv_norm;
L2_A00_vs_N00_predict = sqrt(2*pi)*F_inv_norm*sqrt(n_q)*sqrt(n_k_p_r)*(sigma);
L2_A00_vs_N00_predict_it__(1+niteration,1+ntrial) = L2_A00_vs_N00_predict;
L2_A00_vs_N00 = L2_A00_qk__vs_X00_qki___(A00_qk__,N00_qk__);
if (flag_verbose>1); disp(sprintf(' %% ntrial %d/%d: A00_qk__ vs N00_qk__: %0.16f',ntrial,n_trial,L2_A00_vs_N00)); end;
L2_A00_vs_N00_it__(1+niteration,1+ntrial) = L2_A00_vs_N00;
end;%for niteration=0:n_iteration-1;
end;%for ntrial=0:n_trial-1;
%%%%%%%%;
A00_i0_Mkt___ = zeros(n_M,n_k_p_r,n_trial);
Ap0_i0_Mkt___ = zeros(n_M,n_k_p_r,n_trial);
App_i0_Mkt___ = zeros(n_M,n_k_p_r,n_trial);
A00_fro_i0_Mt__ = zeros(n_M,n_trial);
Ap0_fro_i0_Mt__ = zeros(n_M,n_trial);
App_fro_i0_Mt__ = zeros(n_M,n_trial);
for ntrial=0:n_trial-1;
A00_i0_Mkt___(:,:,1+ntrial) = f_A00_w_(gamma_Mit___(:,1+0,1+ntrial));
A00_fro_i0_Mt__(:,1+ntrial) = sqrt(sum(abs(A00_i0_Mkt___(:,:,1+ntrial)).^2,2));
Ap0_i0_Mkt___(:,:,1+ntrial) = f_Ap0_w_(gamma_Mit___(:,1+0,1+ntrial));
Ap0_fro_i0_Mt__(:,1+ntrial) = sqrt(sum(abs(Ap0_i0_Mkt___(:,:,1+ntrial)).^2,2));
App_i0_Mkt___(:,:,1+ntrial) = f_App_w_(gamma_Mit___(:,1+0,1+ntrial));
App_fro_i0_Mt__(:,1+ntrial) = sqrt(sum(abs(App_i0_Mkt___(:,:,1+ntrial)).^2,2));
end;%for ntrial=0:n_trial-1;
gamma_dif_Mit___ = periodize(bsxfun(@minus,gamma_Mit___,gamma_Mit___(:,1+0,:)),-pi,+pi);
gamma_dif_Mit___ = bsxfun(@minus,gamma_dif_Mit___,mean(gamma_dif_Mit___,1));
gamma_dif_fin_Mt__ = squeeze(gamma_dif_Mit___(:,end,:));
%%%%%%%%;
gamma_Mits____(:,:,:,1+nsigma) = gamma_Mit___;
F_inv_norm_its___(:,:,1+nsigma) = F_inv_norm_it__;
L2_A00_vs_N00_predict_its___(:,:,1+nsigma) = L2_A00_vs_N00_predict_it__;
L2_A00_vs_N00_its___(:,:,1+nsigma) = L2_A00_vs_N00_it__;
gamma_dif_fin_Mts___(:,:,1+nsigma) = gamma_dif_fin_Mt__;
%%%%%%%%;
end;%for nsigma=0:n_sigma-1;

Ap0_ini_Mk__ = f_Ap0_w_(gamma_ini_M_);
Ap0_fro_ini_M_ = sqrt(sum(abs(Ap0_ini_Mk__).^2,2));
Ap0_fro_ini_avg = mean(Ap0_fro_ini_M_);
L2_Ap0 = L2_X00_qk__(Ap0_qk__);
disp(sprintf(' %% L2_Ap0.^2 - mean(2*pi*Ap0_fro_ini_M_.^2): %0.16f',L2_Ap0.^2 - mean(2*pi*Ap0_fro_ini_M_.^2)));
App_ini_Mk__ = f_App_w_(gamma_ini_M_);
App_fro_ini_M_ = sqrt(sum(abs(App_ini_Mk__).^2,2));
App_fro_ini_avg = mean(App_fro_ini_M_);
L2_App = L2_X00_qk__(App_qk__);
disp(sprintf(' %% L2_App.^2 - mean(2*pi*App_fro_ini_M_.^2): %0.16f',L2_App.^2 - mean(2*pi*App_fro_ini_M_.^2)));
%%%%;
gamma_dif_fin_std_s_ = zeros(n_sigma,1);
F_inv_norm_ini_s_ = zeros(n_sigma,1);
for nsigma=0:n_sigma-1;
gamma_dif_fin_std_s_(1+nsigma) = std(gamma_dif_fin_Mts___(:,:,1+nsigma),1,'all');
F_inv_norm_ini_s_(1+nsigma) = mean(F_inv_norm_its___(1+0,:,1+nsigma),'all');
end;%for nsigma=0:n_sigma-1;
%%%%;

gamma_dif_fin_std_qua_s_ = zeros(n_sigma,1);
for nsigma=0:n_sigma-1;
sigma = sigma_(1+nsigma);
%tmp_b = -1; tmp_c = sigma/L2_Ap0; tmp_a = L2_App/L2_Ap0;
tmp_b = -1; tmp_c = sigma/L2_Ap0; tmp_a = (2*pi)*L2_App/L2_Ap0;
tmp_d = tmp_b.^2 - 4*tmp_a*tmp_c; if tmp_d< 0; tmp_r = 1; end; if (tmp_d>=0); tmp_r = (-tmp_b - sqrt(tmp_d))/(2*tmp_a); end;
gamma_dif_fin_std_qua_s_(1+nsigma) = tmp_r;
end;%for nsigma=0:n_sigma-1;
%%%%;
plot(log(sigma_/L2_Ap0),log(gamma_dif_fin_std_s_),'ko',log(sigma_/L2_Ap0),log(sigma_/L2_Ap0),'k-',log(sigma_/L2_Ap0),log(gamma_dif_fin_std_qua_s_),'r-');

save(fname_mat ...
     ,'n_k_p_r','q_max','n_M','n_iteration','n_trial' ...
     ,'sigma_','n_sigma' ...
     ,'A00_qk__' ...
     ,'gamma_Mits____' ...
     ,'F_inv_norm_its___' ...
     ,'L2_A00_vs_N00_predict_its___' ...
     ,'L2_A00_vs_N00_its___' ...
     ,'gamma_dif_fin_Mts___' ...
     ,'Ap0_ini_Mk__' ...
     ,'Ap0_fro_ini_M_' ...
     ,'Ap0_fro_ini_avg' ...
     ,'L2_Ap0' ...
     ,'App_ini_Mk__' ...
     ,'App_fro_ini_M_' ...
     ,'App_fro_ini_avg' ...
     ,'L2_App' ...
     ,'gamma_dif_fin_std_s_' ...
     ,'F_inv_norm_ini_s_' ...
     ,'gamma_dif_fin_std_qua_s_' ...
     );

close_fname_tmp(fname_pre);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

disp('returning'); return;
											
