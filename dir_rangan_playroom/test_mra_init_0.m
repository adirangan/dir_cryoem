% simple teset of multi-reference-alignment error. ;
% Here we compare an alternating-minimization paradigm with ;
% repeated alignment to a random-function, followed by averaging. ;
clear;
rng(0);
q_max = 16; %<-- largest bessel-order. ;
n_q = 1+2*q_max;
q_ = transpose(-q_max:+q_max);
n_psi = 1024;
psi_ = transpose(linspace(0,2*pi,n_psi+1)); psi_ = psi_(1:end-1);
exp_iqpsi__ = zeros(n_psi,n_q);
for nq=0:n_q-1;
q_val = q_(1+nq);
exp_iqpsi__(:,1+nq) = exp(+i*q_val*psi_);
end;%for nq=0:n_q-1;
sigma_q = q_max/3;
F_q_tru_ = crandn(n_q,1).*exp(-q_.^2/(2*sigma_q^2)); %<-- attenuated for smoothness. ;
F_q_tru_ = F_q_tru_/fnorm(F_q_tru_); %<-- normalize ;
%%%%%%%%;
% define images. ;
%%%%%%%%;
n_M = 1024*0.5; %<-- note that, as long as n_M is sufficiently large, the error is determined by \partial_gamma_F_ and the snr. ;
gamma_tru_ = transpose(linspace(0,2*pi,n_M+1)); gamma_tru_ = gamma_tru_(1:end-1);
dgamma = mean(diff(gamma_tru_));
exp_iqdgamma_ = exp(+i*q_*dgamma);
eps__ = crandn(n_q,n_M);
n_iteration = 128; 
gamma_est_ini__ = 2*pi*rand(n_M,n_iteration); %<-- first column will be used for standard alternating-minimization. ;
snr_ = [0;0.10;0.20;0.30;0.40;0.50;0.60;0.70;0.80;0.90;1.00;2.00;3.00;4.00;5.00;10.00]; n_snr = numel(snr_);
sigma_eps_use_ = zeros(n_snr,1);
sigma_am_gam_emp_ = zeros(n_snr,1);
sigma_am_gam_emp__ = zeros(n_iteration,n_snr);
sigma_ra_gam_emp_ = zeros(n_snr,1);
sigma_ra_gam_emp__ = zeros(n_iteration,n_snr);
for nsnr=0:n_snr-1;
snr = snr_(1+nsnr);
l2_signal = sqrt(mean(abs(F_q_tru_).^2));
if (snr==0); sigma_eps_use=0; end;
if (snr> 0); sigma_eps_use = l2_signal/snr; end;
sigma_eps_use_(1+nsnr) = sigma_eps_use;
M_q__ = zeros(n_q,n_M);
for nM=0:n_M-1;
M_q__(:,1+nM) = F_q_tru_.*exp(-i*q_*gamma_tru_(1+nM)) + sigma_eps_use * eps__(:,1+nM);
end;%for nM=0:n_M-1;
%%%%%%%%;
% plot each M_q__. ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
imagesc(real(exp_iqpsi__*M_q__));
colormap(colormap_beach());
xlabel('image'); ylabel('psi');
end;%if flag_plot;
%%%%%%%%;
% Compare F_q_tru_ with F_q_avg_ ;
%%%%%%%%;
F_q_avg_ = mean(M_q__.*exp(+i*q_*transpose(gamma_tru_)),2);
flag_plot=0;
if flag_plot;
plot(psi_,real(exp_iqpsi__*F_q_tru_),'k-',psi_,real(exp_iqpsi__*F_q_avg_),'r-');
xlim([0,2*pi]);
xlabel('psi');ylabel('F'); title('tru vs avg');
end;%if flag_plot;
%%%%%%%%;
% Now run standard alternating-minimization. ;
%%%%%%%%;
gamma_est_ = gamma_est_ini__(:,1); %<-- use first column for standard alternating-minimization. ;
stepsize = 1;
tmp_0in__ = zeros(n_psi,n_M);
tmp_mid__ = zeros(n_psi,n_M);
tmp_out__ = zeros(n_psi,n_M);
flag_plot=0;
if flag_plot;
c_ = colormap_beach(); n_c = size(c_,1);
figure(1);clf; hold on;
end;%if flag_plot;
for niteration=0:n_iteration-1;
F_q_est_ = mean(M_q__.*exp(+i*q_*transpose(gamma_est_)),2); %<-- lsq-step. ;
tmp_0in__(:) = 0;
tmp_0in__(1:n_q,:) = repmat(conj(F_q_est_),[1,n_M]).*M_q__;
tmp_mid__ = circshift(tmp_0in__,-q_max);
tmp_out__ = ifft(tmp_mid__);
[~,index_gamma_upd_] = max(tmp_out__,[],1); index_gamma_upd_ = index_gamma_upd_ - 1;
gamma_upd_ = psi_(1+index_gamma_upd_); %<-- alignment. ;
if flag_plot;
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
gamma_est_ = gamma_est_ + stepsize*(gamma_upd_ - gamma_est_);
gamma_shift = median(periodize(gamma_est_-gamma_tru_,-pi,+pi));
gamma_shift_ = periodize(gamma_est_ - gamma_shift,0,2*pi);
gamma_error_ = periodize(gamma_shift_ - gamma_tru_,-pi,pi);
sigma_gam_emp = std(gamma_error_);
sigma_am_gam_emp__(1+niteration,1+nsnr) = sigma_gam_emp;
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
% Now calculate the errors in gamma_est_. ;
%%%%%%%%;
gamma_shift = median(periodize(gamma_est_-gamma_tru_,-pi,+pi));
gamma_shift_ = periodize(gamma_est_ - gamma_shift,0,2*pi);
gamma_error_ = periodize(gamma_shift_ - gamma_tru_,-pi,pi);
sigma_gam_emp = std(gamma_error_);
sigma_am_gam_emp_(1+nsnr) = sigma_gam_emp;
%%%%%%%%;
% Now run trial-averaged random-alignment. ;
%%%%%%%%;
ra_F_q_est__ = zeros(n_q,n_iteration);
ra_gamma_est__ = zeros(n_M,n_iteration);
stepsize = 1;
for niteration=0:n_iteration-1;
gamma_est_ = gamma_est_ini__(:,1+niteration); %<-- use different columns for random-alignment. ;
tmp_0in__ = zeros(n_psi,n_M);
tmp_mid__ = zeros(n_psi,n_M);
tmp_out__ = zeros(n_psi,n_M);
F_q_est_ = mean(M_q__.*exp(+i*q_*transpose(gamma_est_)),2); %<-- lsq-step. ;
tmp_0in__(:) = 0;
tmp_0in__(1:n_q,:) = repmat(conj(F_q_est_),[1,n_M]).*M_q__;
tmp_mid__ = circshift(tmp_0in__,-q_max);
tmp_out__ = ifft(tmp_mid__);
[~,index_gamma_upd_] = max(tmp_out__,[],1); index_gamma_upd_ = index_gamma_upd_ - 1;
gamma_upd_ = psi_(1+index_gamma_upd_); %<-- alignment. ;
gamma_est_ = gamma_est_ + stepsize*(gamma_upd_ - gamma_est_);
ra_F_q_est__(:,1+niteration) = mean(M_q__.*exp(+i*q_*transpose(gamma_est_)),2); %<-- lsq-step. ;
ra_gamma_est__(:,1+niteration) = gamma_est_;
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
% align the various ra_F_q_est_ against the first. ;
%%%%%%%%;
ra_F_q_fix_ = ra_F_q_est__(:,1+0);
tmp_0in__ = zeros(n_psi,n_iteration);
tmp_mid__ = zeros(n_psi,n_iteration);
tmp_out__ = zeros(n_psi,n_iteration);
tmp_0in__(:) = 0;
tmp_0in__(1:n_q,:) = repmat(conj(ra_F_q_fix_),[1,n_iteration]).*ra_F_q_est__;
tmp_mid__ = circshift(tmp_0in__,-q_max);
tmp_out__ = ifft(tmp_mid__);
[~,index_gamma_upd_] = max(tmp_out__,[],1); index_gamma_upd_ = index_gamma_upd_ - 1;
gamma_upd_ = psi_(1+index_gamma_upd_); %<-- alignment. ;
ra_gamma_upd__ = periodize(ra_gamma_est__ + repmat(transpose(gamma_upd_),[n_M,1]),0,2*pi);;
%%%%%%%%;
% Now we can average each of the gamma over the circle. ;
%%%%%%%%;
ra_cosg_upd__ = cos(ra_gamma_upd__);
ra_sing_upd__ = sin(ra_gamma_upd__);
for niteration=0:n_iteration-1;
ra_cosg_upd_med__(:,1+niteration) = median(ra_cosg_upd__(:,1:1+niteration),2);
ra_sing_upd_med__(:,1+niteration) = median(ra_sing_upd__(:,1:1+niteration),2);
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
% Now calculate the errors in gamma_upd_. ;
%%%%%%%%;
for niteration=0:n_iteration-1;
ra_gamma_upd_ = atan2(ra_sing_upd_med__(:,1+niteration),ra_cosg_upd_med__(:,1+niteration));
gamma_shift = median(periodize(ra_gamma_upd_-gamma_tru_,-pi,+pi));
gamma_shift_ = periodize(ra_gamma_upd_ - gamma_shift,0,2*pi);
gamma_error_ = periodize(gamma_shift_ - gamma_tru_,-pi,pi);
sigma_gam_emp = std(gamma_error_);
sigma_ra_gam_emp__(1+niteration,1+nsnr) = sigma_gam_emp;
end;%for niteration=0:n_iteration-1;
sigma_ra_gam_emp_(1+nsnr) = sigma_gam_emp;
%%%%%%%%;
end;%for nsnr=0:n_snr-1;

flag_plot=0;
if flag_plot;
figure(1);clf;
subplot(1,2,1);
hold on;
plot(snr_,sigma_am_gam_emp_,'ro-',snr_,sigma_ra_gam_emp_,'bx-');
plot(snr_,sigma_am_gam_emp__(1+ 4,:),'ro-',snr_,sigma_ra_gam_emp__(1+ 4,:),'bx-');
plot(snr_,sigma_am_gam_emp__(1+19,:),'ro-',snr_,sigma_ra_gam_emp__(1+19,:),'bx-');
hold off;
legend({'am lim','ra lim','am ni4','ra ni4','am ni19','ra ni19'});
xlabel('snr'); ylabel('\sigma_\gamma');
subplot(2,2,2);
figbeach(); 
imagesc(transpose(sigma_am_gam_emp__),[0,2]);
xlabel('iteration'); ylabel('nsnr'); xlim([1,20]);
title('am');
subplot(2,2,4);
figbeach(); 
imagesc(transpose(sigma_ra_gam_emp__),[0,2]);
xlabel('iteration'); ylabel('nsnr'); xlim([1,20]);
title('ra');
set(gcf,'Position',1+[0,0,1024+512,512]);
fname_fig = sprintf('./test_mra_init_N%d_0.jpg',n_M);
print('-djpeg',fname_fig);
fname_fig = sprintf('./test_mra_init_N%d_0.eps',n_M);
print('-depsc',fname_fig);
end;%if flag_plot;
