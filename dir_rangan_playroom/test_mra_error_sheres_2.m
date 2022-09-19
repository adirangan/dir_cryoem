% simple test of multi-reference-alignment error. ;
% also including a test of sheres-style maximum-likelihood. ;
% as well as a test of ampm. ;
function parameter = test_mra_error_sheres_2(parameter);

if nargin<1;
for q_max=[8:2:16];
for n_M = 2.^(6:12);
parameter = struct('type','parameter');
parameter.q_max = q_max;
parameter.n_M = n_M;
test_mra_error_sheres_2(parameter);
end;%for n_M = 2.^(6:12);
end;%for q_max=[8:2:16];
disp('returning'); return;
end;%if nargin<1;

na=0;
if nargin<1+na; parameter=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'q_max'); parameter.q_max = 16; end;
q_max = parameter.q_max;
if ~isfield(parameter,'n_M'); parameter.n_M = 2048; end;
n_M = parameter.n_M;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

flag_verbose = 1;
flag_disp = 1; nf=0;
flag_replot = 1;

dir_pm_mat = sprintf('/%s/rangan/dir_cryoem/dir_rangan_playroom/dir_pm_mat',string_root);
dir_pm_jpg = sprintf('/%s/rangan/dir_cryoem/dir_rangan_playroom/dir_pm_jpg',string_root);

rseed = 0;
str_rseed = sprintf('r%d',rseed);
%q_max = 16; %<-- largest bessel-order. ;
str_q_max = sprintf('q%d',q_max);
%n_M = 1024*2.0; %<-- number of images. ;
str_n_M = sprintf('M%.4d',round(n_M));
n_iteration = 64*1; str_n_iteration = sprintf('i%.3d',n_iteration);
n_sigma = 32+1; str_n_sigma = sprintf('s%.2d',n_sigma);
str_xfix = sprintf('%s%s%s%s%s',str_q_max,str_n_M,str_n_sigma,str_n_iteration,str_rseed);

fname_pre = sprintf('%s/test_mra_error_sheres_%s',dir_pm_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if 0 & ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% define signal. ;
%%%%%%%%;
n_q = 1+2*q_max;
q_ = transpose(-q_max:+q_max);
n_psi = 1024;
psi_ = transpose(linspace(0,2*pi,n_psi+1)); psi_ = psi_(1:end-1);
dpsi = mean(diff(psi_));
exp_iqpsi__ = zeros(n_psi,n_q);
for nq=0:n_q-1;
q_val = q_(1+nq);
exp_iqpsi__(:,1+nq) = exp(+i*q_val*psi_);
end;%for nq=0:n_q-1;
sigma_q = q_max/3;
rng(rseed); rseed = rseed+1;
F_q_tru_ = (randn(n_q,1) + i*randn(n_q,1)).*exp(-q_.^2/(2*sigma_q^2)); %<-- attenuated for smoothness. ;
F_q_tru_ = F_q_tru_/fnorm(F_q_tru_); %<-- normalize ;
l2_signal = sqrt(mean(abs(F_q_tru_).^2));
%%%%%%%%;
% define images. ;
%%%%%%%%;
gamma_tru_ = transpose(linspace(0,2*pi,n_M+1)); gamma_tru_ = gamma_tru_(1:end-1);
dgamma = mean(diff(gamma_tru_));
exp_iqdgamma_ = exp(+i*q_*dgamma);
rng(rseed); rseed = rseed + 1;
eps__ = (randn(n_q,n_M) + i*randn(n_q,n_M));
rng(rseed); rseed = rseed + 1;
gamma_est_ini_ = 2*pi*rand(n_M,1);
sigma_eps_use_ = exp(transpose(linspace(-1,0,n_sigma)));
snr_ = l2_signal./sigma_eps_use_; snr_(find(~isfinite(snr_)))=0; n_snr = numel(snr_);
sigma_set_ = sigma_eps_use_; n_sigma_set = numel(sigma_set_);
FvF_l2_zero_ = zeros(n_snr,1);
sigma_gam_emp_zero_ = zeros(n_snr,1);
sigma_eps_use_zero_ = zeros(n_snr,1);
FvF_l2_sher__ = zeros(n_sigma_set,n_snr);
sigma_gam_emp_sher__ = zeros(n_sigma_set,n_snr);
sigma_eps_use_sher__ = zeros(n_sigma_set,n_snr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nsnr=0:n_snr-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
snr = snr_(1+nsnr);
if (snr==0); sigma_eps_use=0; end;
if (snr> 0); sigma_eps_use = l2_signal/snr; end;
M_q__ = zeros(n_q,n_M);
for nM=0:n_M-1;
M_q__(:,1+nM) = F_q_tru_.*exp(-i*q_*gamma_tru_(1+nM)) + sigma_eps_use * eps__(:,1+nM);
end;%for nM=0:n_M-1;
M_l2_ = sqrt(sum(abs(M_q__).^2,1))*sqrt(2*pi);
M_l2_ = sqrt(sum(abs(exp_iqpsi__*M_q__).^2,1)*dpsi);
M_p__ = exp_iqpsi__*M_q__;
%%%%%%%%;
% plot each M_q__. ;
%%%%%%%%;
if flag_disp>1;
figure(1+nf);nf=nf+1;
imagesc(real(exp_iqpsi__*M_q__));
colormap(colormap_beach());
xlabel('image'); ylabel('psi');
end;%if flag_disp>1;
%%%%%%%%;
% Compare F_q_tru_ with F_q_avg_ ;
%%%%%%%%;
F_p_tru_ = exp_iqpsi__*F_q_tru_;
F_q_avg_ = mean(M_q__.*exp(+i*q_*transpose(gamma_tru_)),2);
if flag_disp>1;
figure(1+nf);nf=nf+1;
plot(psi_,real(exp_iqpsi__*F_q_tru_),'k-',psi_,real(exp_iqpsi__*F_q_avg_),'r-');
xlim([0,2*pi]);
xlabel('psi');ylabel('F'); title('tru vs avg');
end;%if flag_disp>1;

%%%%%%%%%%%%%%%%;
% Now run AMPM. ;
%%%%%%%%%%%%%%%%;
gamma_est_ = gamma_est_ini_;
MvM_ = sum(conj(M_p__).*M_p__,1)*dpsi/(2*pi);
FvF_tru = sum(conj(F_p_tru_).*F_p_tru_,1)*dpsi/(2*pi);
tmp_0in__ = zeros(n_psi,n_M);
tmp_mid__ = zeros(n_psi,n_M);
tmp_out__ = zeros(n_psi,n_M);
if flag_disp>2;
figure(1+nf);nf_base=nf;nf=nf+1;clf;figmed;
c_ = colormap_beach(); n_c = size(c_,1);
end;%if flag_disp>2;
for niteration=0:n_iteration-1;
F_q_est_ = mean(M_q__.*exp(+i*q_*transpose(gamma_est_)),2); %<-- back-propagation. ;
tmp_0in__(:) = 0;
tmp_0in__(1:n_q,:) = repmat(conj(F_q_est_),[1,n_M]).*M_q__;
tmp_mid__ = circshift(tmp_0in__,-q_max);
tmp_out__ = ifft(tmp_mid__)*n_psi;
F_p_est_ = exp_iqpsi__*F_q_est_;
FvF_est = sum(conj(F_p_est_).*F_p_est_,1)*dpsi/(2*pi);
FvM_R2_est__ = FvF_est + repmat(MvM_,n_psi,1) - 2*real(tmp_out__);
[~,tmp_ij_] = min(FvM_R2_est__,[],1); gamma_est_ = psi_(tmp_ij_);
if flag_disp>2;
subplot(1,2,1);hold on;
nc = max(0,min(n_c-1,floor(n_c*niteration/n_iteration)));
gamma_shift = median(periodize(gamma_est_-gamma_tru_,-pi,+pi));
gamma_shift_ = periodize(gamma_est_ - gamma_shift,0,2*pi);
plot(gamma_tru_,gamma_shift_,'.','Color',c_(1+nc,:));
xlim([0,2*pi]); ylim([0,2*pi]);
axisnotick; axis square;
xlabel('tru'); ylabel('shift (est)');
title(sprintf('snr %0.2f sigma_eps %0.2f',snr,sigma_eps_use),'Interpreter','none');
if (mod(niteration,128)==0); drawnow(); end;
end;%if flag_plot;
if flag_disp>2;
figure(1+nf_base);
subplot(1,2,2); hold on;
nc = max(0,min(n_c-1,floor(n_c*niteration/n_iteration)));
plot(psi_,real(F_p_est_),'-','Color',c_(1+nc,:));
xlim([0,2*pi]);
axisnotick;
xlabel('psi'); ylabel('F_p_est','Interpreter','none');
if (mod(niteration,16)==0); drawnow(); end;
end;%if flag_disp>2;
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
% Now calculate the errors in gamma_est_. ;
%%%%%%%%;
gamma_shift = median(periodize(gamma_est_-gamma_tru_,-pi,+pi));
F_q_est_ = F_q_est_.*exp(-i*q_*gamma_shift);
F_p_est_ = exp_iqpsi__*F_q_est_;
FvF_l2_zero_(1+nsnr) = sum(abs(F_p_est_ - F_p_tru_).^2)*dpsi/FvF_tru;
gamma_shift = median(periodize(gamma_est_-gamma_tru_,-pi,+pi));
gamma_shift_ = periodize(gamma_est_ - gamma_shift,0,2*pi);
gamma_error_ = periodize(gamma_shift_ - gamma_tru_,-pi,pi);
sigma_gam_emp = std(gamma_error_);
sigma_eps_use_zero_(1+nsnr) = sigma_eps_use;
sigma_gam_emp_zero_(1+nsnr) = sigma_gam_emp;
if flag_disp>2;
figure(1+nf);nf=nf+1;clf;figsml;plot(psi_,real(F_p_est_),'k-',psi_,real(F_p_tru_),'r-');
end;%if flag_disp>2;
%%%%%%%%;
if flag_disp>2; close(gcf); end;%if flag_disp>2;
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%;
for nsigma_set=0:n_sigma_set-1;
sigma_set = sigma_set_(1+nsigma_set);
if (flag_verbose); disp(sprintf(' %% nsigma_set %d/%d nsnr %d/%d',nsigma_set,n_sigma_set,nsnr,n_snr)); end;
%%%%%%%%%%%%%%%%;
% Now run scheres-style EM. ;
%%%%%%%%%%%%%%%%;
gamma_est_ = gamma_est_ini_;
%gamma_est_ = gamma_tru_;
gamma_dist_wM__ = zeros(n_psi,n_M);
sigma_dist_M_ = sigma_set*ones(n_M,1);
for nM=0:n_M-1;
sigma_dist = sigma_dist_M_(1+nM);
gamma_est = gamma_est_(1+nM);
gamma_dist_ = 1/sqrt(2*pi)/sigma_dist .* exp(-periodize(psi_ - gamma_est,-pi,+pi).^2/(2*sigma_dist^2));
tmp_int = sum(gamma_dist_)*dpsi;
gamma_dist_wM__(:,1+nM) = gamma_dist_/tmp_int;
end;%for nM=0:n_M-1;
gamma_dist_qM__ = ctranspose(exp_iqpsi__)*gamma_dist_wM__/n_psi;
%%%%;
flag_check=0;
if flag_check;
nM = 48;
M_q_ = M_q__(:,1+nM);
M_p_ = exp_iqpsi__*M_q_;
gamma_dist_w_ = gamma_dist_wM__(:,1+nM);
gamma_dist_q_ = gamma_dist_qM__(:,1+nM);
F_q_0_ = zeros(n_q,1);
for npsi=0:n_psi-1;
gamma = psi_(1+npsi);
gamma_dist = gamma_dist_w_(1+npsi);
F_q_0_ = F_q_0_ + M_q_.*exp(+i*q_*gamma)*gamma_dist/n_psi*(2*pi);
end;%for npsi=0:n_psi-1;
F_q_1_ = M_q_.*conj(gamma_dist_q_)*(2*pi);
disp(sprintf(' %% F_q_0_ vs F_q_1_: %0.16f',fnorm(F_q_0_-F_q_1_)/fnorm(F_q_0_)));
F_p_0_ = exp_iqpsi__*F_q_0_; F_p_1_ = exp_iqpsi__*F_q_1_;
plot(psi_,real(F_p_0_),'r-',psi_,real(F_p_1_),'k-',psi_,real(F_p_tru_),'b-',psi_,real(M_p_),'go');
tmp_0in_ = zeros(n_psi,1);
tmp_0in_(1:n_q) = conj(F_q_0_).*M_q_;
tmp_mid_ = circshift(tmp_0in_,-q_max);
tmp_out_ = ifft(tmp_mid_)*n_psi;
tmp_FvM_ = zeros(n_psi,1);
for npsi=0:n_psi-1;
psi = psi_(1+npsi);
tmp_M_p_ = exp_iqpsi__*(M_q_.*exp(+i*q_*psi));
tmp_FvM_(1+npsi) = sum(conj(F_p_0_).*tmp_M_p_)*dpsi/(2*pi);
end;%for npsi=0:n_psi-1;
tmp_FvF = sum(conj(F_p_0_).*F_p_0_)*dpsi/(2*pi);
tmp_MvM = sum(conj(tmp_M_p_).*tmp_M_p_)*dpsi/(2*pi);
disp(sprintf(' %% tmp_FvM_ vs tmp_out_: %0.16f',fnorm(tmp_FvM_-tmp_out_)/fnorm(tmp_FvM_)));
tmp_FvM_R2_ = tmp_FvF + tmp_MvM - 2*real(tmp_FvM_);
tmp_weight_ = dpsi*ones(n_psi,1)/(2*pi);
tmp_sigma = sigma_me_0(tmp_FvM_R2_,tmp_weight_);
end;%if flag_check;
%%%%;
scheres_weight_w_ = dpsi*ones(n_psi,1)/(2*pi);
MvM_ = sum(conj(M_p__).*M_p__,1)*dpsi/(2*pi);
FvF_tru = sum(conj(F_p_tru_).*F_p_tru_,1)*dpsi/(2*pi);
tmp_0in__ = zeros(n_psi,n_M);
tmp_mid__ = zeros(n_psi,n_M);
tmp_out__ = zeros(n_psi,n_M);
if flag_disp>2;
figure(1+nf);nf_base=nf;nf=nf+1;clf;figmed;
c_ = colormap_beach(); n_c = size(c_,1);
end;%if flag_disp>2;
for niteration=0:n_iteration-1;
F_q_est_ = mean(M_q__.*conj(gamma_dist_qM__)*(2*pi),2); %<-- back-propagation. ;
tmp_0in__(:) = 0;
tmp_0in__(1:n_q,:) = repmat(conj(F_q_est_),[1,n_M]).*M_q__;
tmp_mid__ = circshift(tmp_0in__,-q_max);
tmp_out__ = ifft(tmp_mid__)*n_psi;
F_p_est_ = exp_iqpsi__*F_q_est_;
FvF_est = sum(conj(F_p_est_).*F_p_est_,1)*dpsi/(2*pi);
FvM_R2_est__ = FvF_est + repmat(MvM_,n_psi,1) - 2*real(tmp_out__);
for nM=0:n_M-1;
FvM_R2_est_ = FvM_R2_est__(:,1+nM);
%sigma_dist = sigma_dist_M_(1+nM);
%sigma_dist = sigma_me_0(FvM_R2_est_,scheres_weight_w_,sigma_dist);
%sigma_dist = sqrt(min(FvM_R2_est_));
sigma_dist = sigma_set;
sigma_dist_M_(1+nM) = sigma_dist;
gamma_dist_ = 1/sqrt(2*pi)/sigma_dist .* exp(-FvM_R2_est_/(2*sigma_dist^2));
tmp_int = sum(gamma_dist_)*dpsi;
gamma_dist_wM__(:,1+nM) = gamma_dist_/tmp_int;
end;%for nM=0:n_M-1;
gamma_dist_qM__ = ctranspose(exp_iqpsi__)*gamma_dist_wM__/n_psi;
[~,tmp_ij_] = max(gamma_dist_wM__,[],1); gamma_est_ = psi_(tmp_ij_);
if flag_disp>2;
subplot(1,2,1);hold on;
nc = max(0,min(n_c-1,floor(n_c*niteration/n_iteration)));
gamma_shift = median(periodize(gamma_est_-gamma_tru_,-pi,+pi));
gamma_shift_ = periodize(gamma_est_ - gamma_shift,0,2*pi);
plot(gamma_tru_,gamma_shift_,'.','Color',c_(1+nc,:));
xlim([0,2*pi]); ylim([0,2*pi]);
axisnotick; axis square;
xlabel('tru'); ylabel('shift (est)');
title(sprintf('snr %0.2f sigma_eps %0.2f',snr,sigma_eps_use),'Interpreter','none');
if (mod(niteration,128)==0); drawnow(); end;
end;%if flag_plot;
if flag_disp>2;
figure(1+nf_base);
subplot(1,2,2); hold on;
nc = max(0,min(n_c-1,floor(n_c*niteration/n_iteration)));
plot(psi_,real(F_p_est_),'-','Color',c_(1+nc,:));
xlim([0,2*pi]);
axisnotick;
xlabel('psi'); ylabel('F_p_est','Interpreter','none');
if (mod(niteration,16)==0); drawnow(); end;
end;%if flag_disp>2;
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
% Now calculate the errors in gamma_est_. ;
%%%%%%%%;
gamma_shift = median(periodize(gamma_est_-gamma_tru_,-pi,+pi));
F_q_est_ = F_q_est_.*exp(-i*q_*gamma_shift);
F_p_est_ = exp_iqpsi__*F_q_est_;
FvF_l2_sher__(1+nsigma_set,1+nsnr) = sum(abs(F_p_est_ - F_p_tru_).^2)*dpsi/FvF_tru;
gamma_shift = median(periodize(gamma_est_-gamma_tru_,-pi,+pi));
gamma_shift_ = periodize(gamma_est_ - gamma_shift,0,2*pi);
gamma_error_ = periodize(gamma_shift_ - gamma_tru_,-pi,pi);
sigma_gam_emp = std(gamma_error_);
sigma_eps_use_sher__(1+nsigma_set,1+nsnr) = sigma_eps_use;
sigma_gam_emp_sher__(1+nsigma_set,1+nsnr) = sigma_gam_emp;
if flag_disp>2;
figure(1+nf);nf=nf+1;clf;figsml;plot(psi_,real(F_p_est_),'k-',psi_,real(F_p_tru_),'r-');
end;%if flag_disp>2;
%%%%%%%%;
if flag_disp>2; close(gcf); end;%if flag_disp>2;
%%%%%%%%%%%%%%%%;
end;%for nsigma_set=0:n_sigma_set-1;
%%%%%%%%%%%%%%%%;

if (flag_verbose); disp(sprintf(' %% nsnr %d/%d, saving %s',nsnr,n_snr,fname_mat)); end;
save(fname_mat ...
     ,'rseed' ...
     ,'str_rseed' ...
     ,'q_max' ...
     ,'str_q_max' ...
     ,'n_M' ...
     ,'str_n_M' ...
     ,'n_iteration' ...
     ,'n_sigma' ...
     ,'str_xfix' ...
     ,'fname_pre' ...
     ,'n_q' ...
     ,'q_' ...
     ,'n_psi' ...
     ,'psi_' ...
     ,'dpsi' ...
     ,'F_q_tru_' ...
     ,'l2_signal' ...
     ,'gamma_tru_' ...
     ,'gamma_est_ini_' ...
     ,'sigma_eps_use_' ...
     ,'snr_' ...
     ,'sigma_set_' ...
     ,'FvF_l2_zero_' ...
     ,'sigma_gam_emp_zero_' ...
     ,'sigma_eps_use_zero_' ...
     ,'FvF_l2_sher__' ...
     ,'sigma_gam_emp_sher__' ...
     ,'sigma_eps_use_sher__' ...
     );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nsnr=0:n_snr-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

close_fname_tmp(fname_pre);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if  exist(fname_mat,'file');
load(fname_mat);
n_snr = numel(snr_);
n_sigma_set = numel(sigma_set_);
n_psi = 1024;
psi_ = transpose(linspace(0,2*pi,n_psi+1)); psi_ = psi_(1:end-1);
dpsi = mean(diff(psi_));
exp_iqpsi__ = zeros(n_psi,n_q);
for nq=0:n_q-1;
q_val = q_(1+nq);
exp_iqpsi__(:,1+nq) = exp(+i*q_val*psi_);
end;%for nq=0:n_q-1;
F_p_tru_ = exp_iqpsi__*F_q_tru_;
FvF_tru = sum(conj(F_p_tru_).*F_p_tru_,1)*dpsi/(2*pi);
F_p_one_ = mean(F_p_tru_)*ones(size(F_p_tru_));
FvF_l2_base = sum(abs(F_p_one_ - F_p_tru_).^2)*dpsi/FvF_tru;

flag_plot=0;
if flag_plot;
%%%%%%%%;
fname_fig_pre = sprintf('%s/test_mra_error_sheres_2_%s_FIGA',dir_pm_jpg,str_xfix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;
figure(1);clf;set(gcf,'Position',1+[0,0,1024+1024,768]);fig80s;
%%%%;
subplot(1,2,1);
hold on;
imagesc([transpose(sigma_gam_emp_zero_);sigma_gam_emp_sher__;],[0,pi]);axis image;colorbar;
set(gca,'Ydir','normal');
ylabel('sigma_set','Interpreter','none');
set(gca,'Ytick',1:n_sigma_set+1,'YTickLabel',num2str([0;sigma_set_(:)],2));
xlabel('sigma_true','Interpreter','none');
set(gca,'Xtick',1:n_snr,'XTickLabel',num2str(sigma_eps_use_(:),2)); xtickangle(90);
title('sigma_gam_emp__','Interpreter','none');
set(gca,'FontSize',18);
%%%%;
subplot(1,2,2);
hold on;
tmp_FvF_l2__ = FvF_l2_sher__;
tmp_FvF_l2__(find(~isfinite(tmp_FvF_l2__(:)))) = 1;
imagesc([transpose(FvF_l2_zero_);tmp_FvF_l2__],[0,FvF_l2_base]);axis image;colorbar;
set(gca,'Ydir','normal');
ylabel('sigma_set','Interpreter','none');
set(gca,'Ytick',1:n_sigma_set+1,'YTickLabel',num2str([0;sigma_set_(:)],2));
xlabel('sigma_true','Interpreter','none');
set(gca,'Xtick',1:n_snr,'XTickLabel',num2str(sigma_eps_use_(:),2)); xtickangle(90);
title('FvF_l2__','Interpreter','none');
set(gca,'FontSize',18);
%%%%%%%%;
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
%close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;
end;%if flag_plot;

flag_plot=1;
if flag_plot;
%%%%%%%%;
fname_fig_pre = sprintf('%s/test_mra_error_sheres_2_%s_FIGB',dir_pm_jpg,str_xfix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;
figure(1);clf;set(gcf,'Position',1+[0,0,1024,512]);figbeach;
%%%%;
subplot(1,2,1);
hold on;
plot(psi_,real(F_p_tru_),'r-','LineWidth',2);
plot(psi_,imag(F_p_tru_),'g-','LineWidth',2);
hold off;
xlim([0,2*pi]); xlabel('$\psi$','Interpreter','latex');
set(gca,'XTick',0:pi/2:2*pi,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'});
ylabel('F');
title('F true','Interpreter','none');
set(gca,'FontSize',18);
%%%%;
subplot(1,2,2);
hold on;
tmp_FvF_l2__ = FvF_l2_sher__;
tmp_FvF_l2__(find(~isfinite(tmp_FvF_l2__(:)))) = 1;
imagesc(tmp_FvF_l2__/FvF_l2_base,[0,1]);axis image;
tmp_c_ = colorbar; set(tmp_c_,'Ticks',[0,1]);
set(gca,'Ydir','normal');
ylabel('sigma_BI','Interpreter','none');
set(gca,'Ytick',1:8:n_sigma_set,'YTickLabel',num2str([sigma_set_(1:8:end)],2));
xlabel('sigma_true','Interpreter','none');
set(gca,'Xtick',1:8:n_snr,'XTickLabel',num2str(sigma_eps_use_(1:8:end),2)); xtickangle(90);
title('L2 error (scaled)','Interpreter','none');
set(gca,'FontSize',18);
%%%%%%%%;
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
sgtitle('');
tmp_dir = sprintf('/%s/rangan/dir_cryoem/dir_ampm_manuscript/dir_ampm_fig_misc',string_root);
fname_fig_jpg_strip = sprintf('%s/test_mra_error_sheres_2_%s_FIGB_strip.jpg',tmp_dir,str_xfix);
disp(sprintf(' %% writing %s',fname_fig_jpg_strip));
print('-djpeg',fname_fig_jpg_strip);
%close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;
end;%if flag_plot;

end;%if  exist(fname_mat,'file');

disp('returning'); return;

