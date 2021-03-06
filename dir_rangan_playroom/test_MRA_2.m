%%%%%%%%;
% test out simple MRA example. ;
% assuming n_q modes, ranging from q=-n_q/2,..,+n_q/2 (skipping 0). ;
%%%%%%%%;

rseed = 0;
n_q = 2;
n_q = rup(n_q,2);
q_ = transpose(setdiff(-n_q/2:+n_q/2,0));
n_q = numel(q_);
n_gamma = 64;
gamma_ = linspace(0,2*pi,n_gamma+1);
gamma_ = transpose(gamma_(1:end-1));

rng(rseed);
z_true_q_ = zeros(n_q,1);
z_true_q_ = randn(n_q,1) + i*randn(n_q,1);
r_true_q_ = abs(z_true_q_);
w_true_q_ = atan2(imag(z_true_q_),real(z_true_q_));
t_true = diff(w_true_q_,1,1)/2;

sigma = 3.5;
rseed = 3;
n_M = 16;
z_samp_qM__ = zeros(n_q,n_M);
rng(rseed);
z_nois_qM__ = (randn(n_q,n_M) + i*randn(n_q,n_M));
z_samp_qM__ = repmat(z_true_q_,[1,n_M]) + sigma*z_nois_qM__;
r_samp_qM__ = abs(z_samp_qM__);
w_samp_qM__ = atan2(imag(z_samp_qM__),real(z_samp_qM__));
t_samp_M_ = diff(w_samp_qM__,1,1)/2;

flag_plot=0;
%%%%%%%%;
% visualize theta for true and samp. ;
%%%%%%%%;
if flag_plot;
figure(1);clf;figsml;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
hold on;
%%%%%%%%;
plot(real(z_true_q_),imag(z_true_q_),'kx',r_true_q_.*cos(w_true_q_),r_true_q_.*sin(w_true_q_),'ko');
for nq=0:n_q-1;
r_true_q = r_true_q_(1+nq);
w_true_q = w_true_q_(1+nq);
q = q_(1+nq);
plot(r_true_q.*[cos(w_true_q - q.*t_true)],r_true_q.*[sin(w_true_q - q.*t_true)],'ko');
plot(r_true_q.*[0,cos(w_true_q - q.*t_true)],r_true_q.*[0,sin(w_true_q - q.*t_true)],'k-');
plot(r_true_q.*cos(gamma_),r_true_q.*sin(gamma_),'k-');
end;%for nq=0:n_q-1;
%%%%%%%%;
for nM=0:n_M-1;
nc_beach = max(0,min(n_c_beach-1,floor(n_c_beach*nM/n_M)));
z_samp_q_ = z_samp_qM__(:,1+nM);
r_samp_q_ = r_samp_qM__(:,1+nM);
w_samp_q_ = w_samp_qM__(:,1+nM);
t_samp = t_samp_M_(1+nM);
plot(real(z_samp_q_),imag(z_samp_q_),'x',r_samp_q_.*cos(w_samp_q_),r_samp_q_.*sin(w_samp_q_),'o','MarkerFaceColor',c_beach__(1+nc_beach,:));
for nq=0:n_q-1;
r_samp_q = r_samp_q_(1+nq);
w_samp_q = w_samp_q_(1+nq);
q = q_(1+nq);
plot(r_samp_q.*[cos(w_samp_q - q.*t_samp)],r_samp_q.*[sin(w_samp_q - q.*t_samp)],'o','MarkerFaceColor',c_beach__(1+nc_beach,:));
plot(r_samp_q.*[0,cos(w_samp_q - q.*t_samp)],r_samp_q.*[0,sin(w_samp_q - q.*t_samp)],'-','Color',c_beach__(1+nc_beach,:));
plot(r_samp_q.*cos(gamma_),r_samp_q.*sin(gamma_),'-','Color',c_beach__(1+nc_beach,:));
end;%for nq=0:n_q-1;
end;%for nM=0:n_M-1;
%%%%%%%%;
hold off;
lim_ = [-5,+5];
xlim(lim_);ylim(lim_);
set(gca,'XTick',min(lim_):max(lim_));
set(gca,'YTick',min(lim_):max(lim_));
axis square;
grid on;
%%%%%%%%;
disp('returning'); return;
end;%if flag_plot;

flag_plot=0;
%%%%%%%%;
% visualize gamma alignment. ;
%%%%%%%%;
if flag_plot;
figure(1);clf;figbig;
markersize_use = 16;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
str_symbol_ = {'o','p'}; n_str_symbol=numel(str_symbol_);
rseed = 3;
rng(rseed);
z_nois_qM__ = (randn(n_q,n_M) + i*randn(n_q,n_M));
%%%%;
%%%%%%%%;
sigma_ = linspace(0.5,4,8); n_sigma = numel(sigma_);
for nsigma=0:n_sigma-1;
sigma = sigma_(1+nsigma);
disp(sprintf(' %% nsigma %d/%d <-- %0.2f',nsigma,n_sigma,sigma));
z_samp_qM__ = repmat(z_true_q_,[1,n_M]) + sigma*z_nois_qM__;
r_samp_qM__ = abs(z_samp_qM__);
w_samp_qM__ = atan2(imag(z_samp_qM__),real(z_samp_qM__));
%%%%%%%%;
gamma_true_M_ = zeros(n_M,1);
for nM=0:n_M-1;
r_samp_q_ = r_samp_qM__(:,1+nM);
w_samp_q_ = w_samp_qM__(:,1+nM);
gamma_true = test_MRA_dlossdgamma_0(q_,r_true_q_,w_true_q_,r_samp_q_,w_samp_q_);
gamma_true_M_(1+nM) = gamma_true;
end;%for nM=0:n_M-1;
%%%%%%%%;
subplot(2,ceil(n_sigma/2),1+nsigma);
title(sprintf('sigma %0.2f',sigma),'Interpreter','none');
hold on;
for nq=0:n_q-1;
str_symbol = str_symbol_{1+mod(nq,n_str_symbol)};
r_true_q = r_true_q_(1+nq);
w_true_q = w_true_q_(1+nq);
plot(r_true_q.*cos(gamma_),r_true_q.*sin(gamma_),'k-');
plot(r_true_q.*cos(w_true_q),r_true_q.*sin(w_true_q),str_symbol,'Color','k','MarkerFaceColor','k','MarkerSize',markersize_use);
end;%for nq=0:n_q-1;
%%%%%%%%;
for nM=0:n_M-1;
nc_beach = max(0,min(n_c_beach-1,floor(n_c_beach*nM/n_M)));
z_samp_q_ = z_samp_qM__(:,1+nM);
r_samp_q_ = r_samp_qM__(:,1+nM);
w_samp_q_ = w_samp_qM__(:,1+nM);
gamma_true = gamma_true_M_(1+nM);
for nq=0:n_q-1;
str_symbol = str_symbol_{1+mod(nq,n_str_symbol)};
r_samp_q = r_samp_q_(1+nq);
w_samp_q = w_samp_q_(1+nq);
q = q_(1+nq);
plot(r_samp_q.*cos(gamma_),r_samp_q.*sin(gamma_),'-','Color',c_beach__(1+nc_beach,:));
plot(r_samp_q.*[cos(w_samp_q + q.*gamma_true)],r_samp_q.*[sin(w_samp_q + q.*gamma_true)],str_symbol,'MarkerFaceColor',c_beach__(1+nc_beach,:),'MarkerSize',markersize_use);
end;%for nq=0:n_q-1;
end;%for nM=0:n_M-1;
hold off;
%%%%%%%%;
lim_ = [-8,+8];
xlim(lim_);ylim(lim_);
set(gca,'XTick',min(lim_):max(lim_));
set(gca,'YTick',min(lim_):max(lim_));
axis square;
grid on;
%%%%%%%%;
end;%for nsigma=0:n_sigma-1;
%%%%%%%%;
disp('returning'); return;
end;%if flag_plot;

dlossdgamma = @(gamma,q_,r_true_q_,w_true_q_,r_samp_q_,w_samp_q_) sum( r_samp_q_.*r_true_q_.*q_.*sin(w_samp_q_ + q_.*reshape(gamma,[1,numel(gamma)]) - w_true_q_) , 1 );
lossofgamma = @(gamma,q_,r_true_q_,w_true_q_,r_samp_q_,w_samp_q_) sum( (r_true_q_.*cos(w_true_q_) - r_samp_q_.*cos(w_samp_q_ + q_.*reshape(gamma,[1,numel(gamma)]))).^2 + (r_true_q_.*sin(w_true_q_) - r_samp_q_.*sin(w_samp_q_ + q_.*reshape(gamma,[1,numel(gamma)]))).^2 , 1);

%{
n_M = 16;
rseed_ = 0:7; n_rseed = numel(rseed_);
sigma_ = 2.^[0:1:8]; n_sigma = numel(sigma_);
loss_gsr___ = zeros(n_gamma,n_sigma,n_rseed);
%%%%;
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
disp(sprintf(' %% nrseed %d/%d <-- %0.2f',nrseed,n_rseed,rseed));
rng(rseed);
z_nois_qM__ = (randn(n_q,n_M) + i*randn(n_q,n_M));
%%%%;
for nsigma=0:n_sigma-1;
sigma = sigma_(1+nsigma);
disp(sprintf(' %% nsigma %d/%d <-- %0.2f',nsigma,n_sigma,sigma));
z_samp_qM__ = repmat(z_true_q_,[1,n_M]) + sigma*z_nois_qM__;
r_samp_qM__ = abs(z_samp_qM__);
w_samp_qM__ = atan2(imag(z_samp_qM__),real(z_samp_qM__));
%%%%;
for ngamma=0:n_gamma-1;
gamma = gamma_(1+ngamma);
r_base_q_ = r_true_q_;
w_base_q_ = w_true_q_ + gamma; %<-- offset by gamma. ;
gamma_base_M_ = zeros(n_M,1);
loss_base_M_ = zeros(n_M,1);
%%%%;
for nM=0:n_M-1;
r_samp_q_ = r_samp_qM__(:,1+nM);
w_samp_q_ = w_samp_qM__(:,1+nM);
[gamma_base,loss_base] = test_MRA_dlossdgamma_0(q_,r_base_q_,w_base_q_,r_samp_q_,w_samp_q_);
gamma_base_M_(1+nM) = gamma_base;
loss_base_M_(1+nM) = loss_base;
end;%for nM=0:n_M-1;
%%%%;
loss_gsr___(1+ngamma,1+nsigma,1+nrseed) = sum(loss_base_M_);
end;%for ngamma=0:n_gamma-1;
%%%%;
end;%for nsigma=0:n_sigma-1;
%%%%;
end;%for nrseed=0:n_rseed-1;
%%%%;
loss_nrm_gsr___ = loss_gsr___;
loss_min_gsr___ = loss_gsr___;
for nrseed=0:n_rseed-1;
for nsigma=0:n_sigma-1;
loss_nrm_g_= loss_nrm_gsr___(:,1+nsigma,1+nrseed);
loss_nrm_g_ = (loss_nrm_g_ - mean(loss_nrm_g_)) / max(1e-12,std(loss_nrm_g_,1));
loss_min_g_ = (diff([loss_nrm_g_ ; loss_nrm_g_(1)])>=0) & (diff([loss_nrm_g_(end) ; loss_nrm_g_])<=0);
loss_nrm_gsr___(:,1+nsigma,1+nrseed) = loss_nrm_g_;
loss_min_gsr___(:,1+nsigma,1+nrseed) = loss_min_g_;
end;%for nsigma=0:n_sigma-1;
end;%for nrseed=0:n_rseed-1;
flag_plot=1;
%%%%%%%%;
% visualize gamma alignment. ;
%%%%%%%%;
if flag_plot;
figure(1);clf;figbig;
%subplot(1,3,1);plot(gamma_,loss_gs__,'k.-');
%xlim([0,2*pi]);xlabel('gamma');ylabel('loss');
subplot(1,2,1); imagesc(loss_nrm_gsr___(:,:,4)); fig80s;
set(gca,'XTick',1:n_sigma,'XTickLabel',log2(sigma_)); xlabel('log2(sigma)');
set(gca,'YTick',1:n_gamma,'YTickLabel',gamma_); xlabel('gamma');
subplot(1,2,2); imagesc(sum(loss_min_gsr___,3)); fig80s;
set(gca,'XTick',1:n_sigma,'XTickLabel',log2(sigma_)); xlabel('log2(sigma)');
set(gca,'YTick',1:n_gamma,'YTickLabel',gamma_); xlabel('gamma');
%%%%%%%%;
end;%if flag_plot;
 %}

n_M = 16;
n_gamma = 64;
gamma_ = linspace(0,2*pi,n_gamma+1);
gamma_ = transpose(gamma_(1:end-1));
rseed_ = 3:3; n_rseed = numel(rseed_);
sigma_ = linspace(3,5,6); n_sigma = numel(sigma_);
r_factor_ = 2.^linspace(-1,+1,1+n_gamma); n_r_factor = numel(r_factor_);
loss_gfsr____ = zeros(n_gamma,n_r_factor,n_sigma,n_rseed);
%%%%;
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
disp(sprintf(' %% nrseed %d/%d <-- %0.2f',nrseed,n_rseed,rseed));
rng(rseed);
z_nois_qM__ = (randn(n_q,n_M) + i*randn(n_q,n_M));
%%%%;
for nsigma=0:n_sigma-1;
sigma = sigma_(1+nsigma);
disp(sprintf(' %% nsigma %d/%d <-- %0.2f',nsigma,n_sigma,sigma));
z_samp_qM__ = repmat(z_true_q_,[1,n_M]) + sigma*z_nois_qM__;
r_samp_qM__ = abs(z_samp_qM__);
w_samp_qM__ = atan2(imag(z_samp_qM__),real(z_samp_qM__));
%%%%;
for nr_factor=0:n_r_factor-1;
r_factor = r_factor_(1+nr_factor);
r_base_q_ = r_true_q_;
r_base_q_(1+0) = r_base_q_(1+0) * r_factor;
r_base_q_(1+1) = r_base_q_(1+1) / r_factor;
for ngamma=0:n_gamma-1;
gamma = gamma_(1+ngamma);
w_base_q_ = w_true_q_ + gamma; %<-- offset by gamma. ;
gamma_base_M_ = zeros(n_M,1);
loss_base_M_ = zeros(n_M,1);
%%%%;
for nM=0:n_M-1;
r_samp_q_ = r_samp_qM__(:,1+nM);
w_samp_q_ = w_samp_qM__(:,1+nM);
[gamma_base,loss_base] = test_MRA_dlossdgamma_0(q_,r_base_q_,w_base_q_,r_samp_q_,w_samp_q_);
gamma_base_M_(1+nM) = gamma_base;
loss_base_M_(1+nM) = loss_base;
end;%for nM=0:n_M-1;
%%%%;
loss_gfsr____(1+ngamma,1+nr_factor,1+nsigma,1+nrseed) = sum(loss_base_M_);
end;%for ngamma=0:n_gamma-1;
end;%for nr_factor=0:n_r_factor-1;
%%%%;
end;%for nsigma=0:n_sigma-1;
%%%%;
end;%for nrseed=0:n_rseed-1;
%%%%;

loss_nrm_gfsr____ = loss_gfsr____;
%loss_min_gfsr____ = loss_gfsr____;
for nrseed=0:n_rseed-1;
for nsigma=0:n_sigma-1;
loss_nrm_gf__= loss_nrm_gfsr____(:,:,1+nsigma,1+nrseed);
loss_nrm_gf__ = (loss_nrm_gf__ - mean(loss_nrm_gf__,'all')) / max(1e-12,std(loss_nrm_gf__,1,'all'));
%loss_min_gf__ = (diff([loss_nrm_gf__ ; loss_nrm_gf__(1)])>=0) & (diff([loss_nrm_gf__(end) ; loss_nrm_gf__])<=0);
loss_nrm_gfsr____(:,:,1+nsigma,1+nrseed) = loss_nrm_gf__;
%loss_min_gfsr____(:,:,1+nsigma,1+nrseed) = loss_min_gf__;
end;%for nsigma=0:n_sigma-1;
end;%for nrseed=0:n_rseed-1;

 
flag_plot=1;
%%%%%%%%;
% visualize gamma alignment. ;
%%%%%%%%;
if flag_plot;
figure(1);clf;figbig;
for nsigma=0:n_sigma-1;
sigma = sigma_(1+nsigma);
subplot(2,ceil(n_sigma/2),1+nsigma);
loss_nrm_gf__ = squeeze(loss_nrm_gfsr____(:,:,1+nsigma,1+0));
%imagesc(loss_nrm_gf__); fig80s;
%contour(log2(r_factor_),gamma_,loss_nrm_gf__,transpose(prctile(loss_nrm_gf__,[0:1:50],'all'))); fig80s;
contour(loss_nrm_gf__,transpose(prctile(loss_nrm_gf__,[0:2.5:100],'all'))); fig80s;
set(gca,'YTick',1:n_gamma,'YTickLabel',gamma_); xlabel('gamma');
set(gca,'XTick',1:n_r_factor,'XTickLabel',log2(r_factor_)); xlabel('log2(r_factor)','Interpreter','none');
title(sprintf('sigma %0.2f',sigma),'Interpreter','none');
end;%for nsigma=0:n_sigma-1;
%%%%%%%%;
end;%if flag_plot;
