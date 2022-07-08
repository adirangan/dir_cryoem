function test_mra_error_analysis_2(q_max,n_M,flag_recalc,flag_replot);
% simple test of multi-reference-alignment error. ;
% trying to verify each term in local expansion. ;

if nargin<1;
for q_max = 2:16;
for n_M = [256,512,1024,2048];
test_mra_error_analysis_2(q_max,n_M,1,1);
end;%for n_M = [256,512,1024,2048];
end;%for q_max = 2:16;
disp('returning');
end;%if nargin<1;

na=0;
if (nargin<1+na); q_max = []; end; na=na+1;
if (nargin<1+na); n_M = []; end; na=na+1;
if (nargin<1+na); flag_recalc = []; end; na=na+1;
if (nargin<1+na); flag_replot = []; end; na=na+1;

verbose=1;
if isempty(flag_recalc); flag_recalc=0; end;
if isempty(flag_replot); flag_replot=1; end;
rng(1);

%%%%%%%%;
% define volume. ;
%%%%%%%%;
if isempty(q_max); q_max = 3; end; %<-- largest bessel-order. ;
n_q = 1+2*q_max;
q_ = transpose(-q_max:+q_max);
n_psi = max(256,n_q*32);
psi_ = transpose(linspace(0,2*pi,n_psi+1)); psi_ = psi_(1:end-1); dpsi = mean(diff(psi_));
exp_iqpsi__ = zeros(n_psi,n_q);
for nq=0:n_q-1;
q_val = q_(1+nq);
exp_iqpsi__(:,1+nq) = exp(+i*q_val*psi_);
end;%for nq=0:n_q-1;
sigma_q = q_max/3;
F_q_tru_ = crandn(n_q,1).*exp(-q_.^2/(2*sigma_q^2)); %<-- attenuated for smoothness. ;
F_q_tru_ = F_q_tru_/fnorm(F_q_tru_); %<-- normalize ;
l2_signal = sqrt(mean(abs(F_q_tru_).^2));
F_w_tru_ = exp_iqpsi__*F_q_tru_;
if (verbose); disp(sprintf(' %% norm-squard: %0.6f <-- %0.6f',sum(abs(F_q_tru_).^2),sum(abs(F_w_tru_).^2)*dpsi/(2*pi))); end;
if (verbose); disp(sprintf(' %% l2_signal: %0.6f <-- sqrt(mean(abs(F_q_tru_).^2))',l2_signal)); end;
%%%%%%%%;

%%%%%%%%;
% define noise. ;
%%%%%%%%;
sigma_ = exp(-10:0.125:0); n_sigma = numel(sigma_);
n_pass = 2;

flag_check=1;
if flag_check;
gamma_tru_ini = sqrt(2);
M_q_ = F_q_tru_.*exp(-i*q_*gamma_tru_ini);
tmp_0in_q_ = zeros(n_psi,1);
tmp_0in_q_(1:n_q) = conj(F_q_tru_).*M_q_;
tmp_mid_q_ = circshift(tmp_0in_q_,-q_max);
tmp_out_q_ = ifft(tmp_mid_q_);
[~,index_gamma_upd] = max(real(tmp_out_q_),[],1); index_gamma_upd = index_gamma_upd - 1;
gamma_upd = psi_(1+index_gamma_upd); %<-- alignment. ;
if (verbose); disp(sprintf(' %% dpsi %0.16f',dpsi)); end;
if (verbose); disp(sprintf(' %% npass x/x: gamma_tru_ini vs gamma_upd %0.16f',gamma_tru_ini - gamma_upd)); end;
gamma_upd = gamma_upd + dpsi*quadratic_1d_interpolation_0(real(tmp_out_q_(1+periodize(index_gamma_upd+[-1:+1],0,n_psi))));
if (verbose); disp(sprintf(' %% npass x/x: gamma_tru_upd vs gamma_upd %0.16f',gamma_tru_ini - gamma_upd)); end;
for npass=0:n_pass-1;
tmp_M_q_ = M_q_.*exp(+i*q_*gamma_upd); %<-- check alignment. ;
tmp_0in_q_ = zeros(n_psi,1);
tmp_0in_q_(1:n_q) = conj(F_q_tru_).*tmp_M_q_;
tmp_mid_q_ = circshift(tmp_0in_q_,-q_max);
tmp_out_q_ = ifft(tmp_mid_q_);
[~,index_gamma_set] = max(real(tmp_out_q_),[],1); index_gamma_set = index_gamma_set - 1;
gamma_set = psi_(1+index_gamma_set); %<-- alignment. ;
gamma_set = gamma_set + dpsi*quadratic_1d_interpolation_0(real(tmp_out_q_(1+periodize(index_gamma_set+[-1:+1],0,n_psi))));
gamma_upd = gamma_upd + gamma_set;
if (verbose); disp(sprintf(' %% npass %d/%d: gamma_tru_upd vs gamma_upd %0.16f',npass,n_pass,gamma_tru_ini - gamma_upd)); end;
end;%for npass=0:n_pass-1;
end;%if flag_check;

%%%%%%%%;
% define images. ;
%%%%%%%%;
if isempty(n_M); n_M = n_q*256; end; %<-- number of images. ;

dir_trunk = '/data/rangan/dir_cryoem/dir_rangan_playroom';
dir_pm_mat = sprintf('%s/dir_pm_mat',dir_trunk);
dir_pm_jpg = sprintf('%s/dir_pm_jpg',dir_trunk);
fname_mat_pre = sprintf('%s/test_mra_error_analysis_2_Q%dN%d_',dir_pm_mat,q_max,n_M);
fname_mat = sprintf('%s.mat',fname_mat_pre);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (verbose); disp(sprintf(' %% %s not found, creating',fname_mat)); end;
gamma_tru_ini_ = transpose(linspace(0,2*pi,n_M+1)); gamma_tru_ini_ = gamma_tru_ini_(1:end-1);
gamma_tru_fin_ = zeros(n_M,1);
eps_qM__ = crandn(n_q,n_M);
gamma_est_ini_ = 2*pi*rand(n_M,1);
sigma_eps_use_ = sigma_;
snr_ = l2_signal./sigma_eps_use_; snr_(find(~isfinite(snr_)))=0; n_snr = numel(snr_);
sigma_gam_emp_ = zeros(n_snr,1);
gamma_set_error_sM__ = zeros(n_snr,n_M);

%%%%%%%%%%%%%%%%;
for nsnr=0:n_snr-1;
%%%%%%%%%%%%%%%%;

snr = snr_(1+nsnr);
if (snr==0); sigma_eps_use=0; end;
if (snr> 0); sigma_eps_use = l2_signal/snr; end;
M_qM__ = zeros(n_q,n_M);
for nM=0:n_M-1;
eps_q_ = eps_qM__(:,1+nM);
gamma_tru_ini = gamma_tru_ini_(1+nM);
M_q_ = (F_q_tru_ + sigma_eps_use * eps_q_).*exp(-i*q_*gamma_tru_ini); %<-- add noise. ;
tmp_0in_q_ = zeros(n_psi,1);
tmp_0in_q_(1:n_q) = conj(F_q_tru_).*M_q_;
tmp_mid_q_ = circshift(tmp_0in_q_,-q_max);
tmp_out_q_ = ifft(tmp_mid_q_);
[~,index_gamma_upd] = max(real(tmp_out_q_),[],1); index_gamma_upd = index_gamma_upd - 1;
gamma_upd = psi_(1+index_gamma_upd); %<-- alignment. ;
gamma_upd = gamma_upd + dpsi*quadratic_1d_interpolation_0(real(tmp_out_q_(1+periodize(index_gamma_upd+[-1:+1],0,n_psi))));
%%%%;
for npass=0:n_pass-1;
tmp_M_q_ = M_q_.*exp(+i*q_*gamma_upd); %<-- check alignment. ;
tmp_0in_q_ = zeros(n_psi,1);
tmp_0in_q_(1:n_q) = conj(F_q_tru_).*tmp_M_q_;
tmp_mid_q_ = circshift(tmp_0in_q_,-q_max);
tmp_out_q_ = ifft(tmp_mid_q_);
[~,index_gamma_set] = max(real(tmp_out_q_),[],1); index_gamma_set = index_gamma_set - 1;
gamma_set = psi_(1+index_gamma_set); %<-- alignment. ;
gamma_set = gamma_set + dpsi*quadratic_1d_interpolation_0(real(tmp_out_q_(1+periodize(index_gamma_set+[-1:+1],0,n_psi))));
gamma_upd = gamma_upd + gamma_set;
end;%for npass=0:n_pass-1;
%%%%;
gamma_tru_fin = +gamma_upd; %<-- recenter. ;
gamma_tru_fin_(1+nM) = periodize(gamma_tru_fin,0,2*pi);
%%%%;
flag_check=0;
if flag_check;
tmp_M_q_ = M_q_.*exp(+i*q_*gamma_tru_fin); %<-- check alignment. ;
tmp_0in_q_ = zeros(n_psi,1);
tmp_0in_q_(1:n_q) = conj(F_q_tru_).*tmp_M_q_;
tmp_mid_q_ = circshift(tmp_0in_q_,-q_max);
tmp_out_q_ = ifft(tmp_mid_q_);
[~,index_gamma_set] = max(real(tmp_out_q_),[],1); index_gamma_set = index_gamma_set - 1;
gamma_set = psi_(1+index_gamma_set); %<-- alignment. ;
gamma_set = gamma_set + dpsi*quadratic_1d_interpolation_0(real(tmp_out_q_(1+periodize(index_gamma_set+[-1:+1],0,n_psi))));
gamma_set_error_sM__(1+nsnr,1+nM) = fnorm(gamma_set);
end;%if flag_check;
%%%%;
M_qM__(:,1+nM) = M_q_;
end;%for nM=0:n_M-1;

%%%%%%%%;
% plot each M_qM__. ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
imagesc(real(exp_iqpsi__*M_qM__));
colormap(colormap_beach());
xlabel('image'); ylabel('psi');
end;%if flag_plot;
%%%%%%%%;
% Compare F_q_tru_ with F_q_avg_ ;
%%%%%%%%;
F_q_avg_ = mean(M_qM__.*exp(+i*q_*transpose(gamma_tru_fin_)),2);
flag_plot=0;
if flag_plot;
plot(psi_,real(exp_iqpsi__*F_q_tru_),'k-',psi_,real(exp_iqpsi__*F_q_avg_),'r-');
xlim([0,2*pi]);
xlabel('psi');ylabel('F'); title('tru vs avg');
end;%if flag_plot;
%%%%%%%%;
% assign random viewing-angles. ;
%%%%%%%%;
gamma_est_ = gamma_est_ini_;
%%%%%%%%;
% Now run alternating minimization. ;
%%%%%%%%;
n_iteration = 128; stepsize = 1;
tmp_0in_wM__ = zeros(n_psi,n_M);
tmp_mid_wM__ = zeros(n_psi,n_M);
tmp_out_wM__ = zeros(n_psi,n_M);
flag_plot=1;
if flag_plot;
c_ = colormap_beach(); n_c = size(c_,1);
figure(1);clf; hold on;
end;%if flag_plot;
for niteration=0:n_iteration-1;
F_q_est_ = mean(M_qM__.*exp(+i*q_*transpose(gamma_est_)),2); %<-- lsq-step. ;
tmp_0in_wM__(:) = 0;
tmp_0in_wM__(1:n_q,:) = repmat(conj(F_q_est_),[1,n_M]).*M_qM__;
tmp_mid_wM__ = circshift(tmp_0in_wM__,-q_max);
tmp_out_wM__ = ifft(tmp_mid_wM__);
[~,index_gamma_upd_] = max(tmp_out_wM__,[],1); index_gamma_upd_ = index_gamma_upd_ - 1;
gamma_upd_ = psi_(1+index_gamma_upd_); %<-- alignment. ;
tmp_y_3M__ = zeros(3,n_M);
for nM=0:n_M-1;
tmp_y_3M__(:,1+nM) = real(tmp_out_wM__(1+periodize(index_gamma_upd_(1+nM)+[-1:+1],0,n_psi),1+nM));
end;%for nM=0:n_M-1;
gamma_upd_ = gamma_upd_ + dpsi*transpose(quadratic_1d_interpolation_0(tmp_y_3M__));
%%%%;
for npass=0:n_pass-1;
tmp_M_qM__ = M_qM__.*exp(+i*q_*transpose(gamma_upd_)); %<-- check alignment. ;
tmp_0in_wM__ = zeros(n_psi,n_M);
tmp_0in_(1:n_q,:) = conj(F_q_est_).*tmp_M_qM__;
tmp_mid_wM__ = circshift(tmp_0in_wM__,-q_max);
tmp_out_wM__ = ifft(tmp_mid_wM__);
[~,index_gamma_set_] = max(real(tmp_out_wM__),[],1); index_gamma_set_ = index_gamma_set_ - 1;
gamma_set_ = psi_(1+index_gamma_set_); %<-- alignment. ;
for nM=0:n_M-1;
tmp_y_3M__(:,1+nM) = real(tmp_out_wM__(1+periodize(index_gamma_set_(1+nM)+[-1:+1],0,n_psi),1+nM));
end;%for nM=0:n_M-1;
gamma_set_ = gamma_set_ + dpsi*transpose(quadratic_1d_interpolation_0(tmp_y_3M__));
gamma_upd_ = gamma_upd_ + gamma_set_;
end;%for npass=0:n_pass-1;
%%%%;
if flag_plot;
nc = max(0,min(n_c-1,floor(n_c*niteration/n_iteration)));
gamma_shift = median(periodize(gamma_est_-gamma_tru_fin_,-pi,+pi));
gamma_shift_ = periodize(gamma_est_ - gamma_shift,0,2*pi);
plot(gamma_tru_fin_,gamma_shift_,'.','Color',c_(1+nc,:));
xlim([0,2*pi]); ylim([0,2*pi]);
axisnotick; axis square;
xlabel('tru'); ylabel('shift (est)');
title(sprintf('snr=exp(%+0.2f) sigma_eps=exp(%+0.2f)',log(snr),log(sigma_eps_use)),'Interpreter','none');
if (mod(niteration,128)==0); drawnow(); end;
end;%if flag_plot;
gamma_est_ = gamma_est_ + stepsize*(gamma_upd_ - gamma_est_);
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
% Now calculate the errors in gamma_est_. ;
%%%%%%%%;
gamma_shift = median(periodize(gamma_est_-gamma_tru_fin_,-pi,+pi));
gamma_shift_ = periodize(gamma_est_ - gamma_shift,0,2*pi);
gamma_error_ = periodize(gamma_shift_ - gamma_tru_fin_,-pi,pi);
sigma_gam_emp = std(gamma_error_);
sigma_eps_use_(1+nsnr) = sigma_eps_use;
sigma_gam_emp_(1+nsnr) = sigma_gam_emp;

%%%%%%%%%%%%%%%%;
end;%for nsnr=0:n_snr-1;
%%%%%%%%%%%%%%%%;

%%%%%%%%;
if flag_plot; close(gcf); end;%if flag_plot;
%%%%%%%%;
save(fname_mat ...
     ,'n_pass','sigma_','n_sigma','q_max','n_psi','n_M','q_','F_q_tru_','l2_signal' ...
     ,'sigma_eps_use_','sigma_gam_emp_' ...
     );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%;
if exist(fname_mat,'file');
%%%%%%%%%%%%%%%%;
fname_fig_pre = sprintf('%s/test_mra_error_analysis_2_Q%dN%d_FIGA',dir_pm_jpg,q_max,n_M);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));

load(fname_mat);

alpha = 0;
%alpha = dpsi;
%alpha = 1;
%%%%%%%%;
flsigma = floor(log(min(sigma_)));
sigma_uni_ = transpose(exp(linspace(flsigma,+1,1024))); n_sigma_uni = numel(sigma_uni_);
sigma_gam_equi_upb_ = zeros(n_sigma_uni,1);
sigma_gam_equi_lob_ = zeros(n_sigma_uni,1);
partial_gamma_F_ = -i*q_.*F_q_tru_ ;
inv_partial_gamma_F_ = inv(ctranspose(partial_gamma_F_)*partial_gamma_F_) * ctranspose(partial_gamma_F_) ;
tmp_a = fnorm(inv_partial_gamma_F_ * (q_.^2 .* F_q_tru_));
tmp_b = -1.0d0;
tmp_c_ = fnorm(inv_partial_gamma_F_) * sigma_uni_ * sqrt(alpha^2 + 1/n_M);
for nsigma_uni=0:n_sigma_uni-1;
sigma_uni = sigma_uni_(1+nsigma_uni);
tmp_c = tmp_c_(1+nsigma_uni);
tmp_d = tmp_b^2 - 4*tmp_a*tmp_c;
if (tmp_d< 0); sigma_gam_equi_upb = +Inf; sigma_gam_equi_lob = +Inf; end;
if (tmp_d>=0); sigma_gam_equi_upb = (-tmp_b + sqrt(tmp_d))/(2*tmp_a); sigma_gam_equi_lob = (-tmp_b - sqrt(tmp_d))/(2*tmp_a); end;
sigma_gam_equi_upb_(1+nsigma_uni) = sigma_gam_equi_upb;
sigma_gam_equi_lob_(1+nsigma_uni) = sigma_gam_equi_lob;
end;%for nsigma_uni=0:n_sigma_uni-1;
%%%%%%%%;

figure(2); clf;
set(gcf,'Position',1+[0,0,512,512]);
subplot(1,1,1);
hold on;
%plot(sigma_uni_,sigma_gam_equi_upb_,'k-');
plot(log(sigma_uni_),log( tmp_c_ ),'-','Color',0.85*[1,1,1],'LineWidth',3);
plot(log(sigma_uni_),log(sigma_gam_equi_lob_),'k-','LineWidth',3);
tmp_index = max(efind(isfinite(sigma_gam_equi_lob_)));
plot(log(sigma_uni_(1+tmp_index))*[1,1],[log(sigma_gam_equi_lob_(1+tmp_index)),10],'k-','LineWidth',3);
plot(log(sigma_eps_use_),log(sigma_gam_emp_),'ro','MarkerFaceColor','r');
xlabel('$\log(\sigma_{\epsilon})$','Interpreter','latex');
ylabel('$\log(\sigma_{\gamma})$','Interpreter','latex');
title('error (log-scale)');
xlim([-flsigma,+1]);
ylim([-flsigma,+1]);
grid on;
set(gca,'FontSize',18);
%%%%%%%%;
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);

end;%if (flag_replot | ~exist(fname_fig,'file'));

%%%%%%%%%%%%%%%%;
end;%if exist(fname_mat,'file');
%%%%%%%%%%%%%%%%;
