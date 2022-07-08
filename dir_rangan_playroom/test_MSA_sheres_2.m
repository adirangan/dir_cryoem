%%%%%%%%;
% setting up simple multi-slice-alignment for testing sheres-style qbp. ;
% assume a 2d-volume x, comprising a single ring. ;
% Assume that x is restricted to just 2*Q+1 bessel-coefficients. ;
% Allowing for a nonuniform viewing-angle distribution. ;
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

flag_recalc = 0;
flag_verbose = 1;
flag_disp = 0*flag_verbose;
flag_replot = 1;
nf = 0;

q_max = 8; str_q_max = sprintf('q%d',q_max);
n_M_factor = 1/32; n_gamma = 1024; n_M = ceil(n_gamma*n_M_factor);
str_n_M = sprintf('M%.4d',round(n_M));
n_sigma = 11; str_n_sigma = sprintf('s%.2d',n_sigma);
%n_lambda = 4+1; str_n_lambda = sprintf('l%.2d',n_lambda);
n_lambda = 1; str_n_lambda = sprintf('');
n_beta = 6+1; str_n_beta = sprintf('b%.2d',n_beta);
n_iteration = 128*1; str_n_iteration = sprintf('i%.3d',n_iteration);
n_rseed = 512; str_n_rseed = sprintf('r%.2d',n_rseed);
str_xfix = sprintf('%s%s%s%s%s%s',str_q_max,str_n_M,str_n_sigma,str_n_beta,str_n_iteration,str_n_rseed);
%%%%%%%%;
q_ = transpose(-q_max:1:+q_max);
n_q = numel(q_);
gamma_ = linspace(0,2*pi,1+n_gamma);
gamma_ = reshape(gamma_(1:n_gamma),[n_gamma,1]);
dgamma = mean(diff(gamma_));
sigma_ = transpose(linspace(0,1,n_sigma));
lambda_ = transpose(linspace(0,1,n_lambda)); if n_lambda==1; lambda_ = [0]; end;
beta_ = transpose(linspace(0,3,n_beta)); if n_beta==1; beta_ = [0]; end;
n_w = n_gamma;
F_wq__ = exp(+i*gamma_*transpose(q_));
F_inv_qw__ = ctranspose(F_wq__)/n_gamma;
%%%%%%%%;

fname_pre = sprintf('./dir_pm_mat/test_MSA_sheres_%s',str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre);
fname_tmp_mat = sprintf('%s_tmp.mat',fname_pre);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_recalc | ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
rseed_ = 0:n_rseed-1;
Error_sher_l2_ilbsr_____ = zeros(n_iteration,n_lambda,n_beta,n_sigma,n_rseed);
Error_zero_l2_ilbsr_____ = zeros(n_iteration,n_lambda,n_beta,n_sigma,n_rseed);
Error_emax_l2_ilbsr_____ = zeros(n_iteration,n_lambda,n_beta,n_sigma,n_rseed);
Error_alte_l2_ilbsr_____ = zeros(n_iteration,n_lambda,n_beta,n_sigma,n_rseed);
Error_sher_k1_ilbsr_____ = zeros(n_iteration,n_lambda,n_beta,n_sigma,n_rseed);
Error_zero_k1_ilbsr_____ = zeros(n_iteration,n_lambda,n_beta,n_sigma,n_rseed);
Error_emax_k1_ilbsr_____ = zeros(n_iteration,n_lambda,n_beta,n_sigma,n_rseed);
Error_alte_k1_ilbsr_____ = zeros(n_iteration,n_lambda,n_beta,n_sigma,n_rseed);
A_qr__ = zeros(n_q,n_rseed);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nrseed=0:n_rseed-1;
for nbeta=0:n_beta-1;
beta = beta_(1+nbeta);
I0_val = besseli(0,beta);
p_den_ = @(gamma) exp(beta*cos(gamma))/I0_val/(2*pi);
p_cdf_ = cumsum(p_den_(gamma_)*dgamma);
for nlambda=0:n_lambda-1;
for nsigma=0:n_sigma-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
rseed = rseed_(1+nrseed);
lambda = lambda_(1+nlambda);
sigma = sigma_(1+nsigma); sigma_true = sigma;
%%%%;
rng(rseed);rseed=rseed+1;
A_q_ = randn(n_q,1) + i*randn(n_q,1);
A_qr__(:,1+nrseed) = A_q_;
A_w_ = F_wq__*A_q_;
%%%%;
rng(rseed);rseed=rseed+1;
%gamma_true_M_ = 2*pi*rand(n_M,1); %<-- uniform distribution. ;
gamma_true_M_ = interp1(p_cdf_,gamma_,rand(n_M,1)); %<-- distribution exp(beta*cos(theta))/I0(beta)/(2*pi). ;
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
p_sher_wM__ = exp(-R2_sher_wM__/max(1e-12,2*sigma_sher^2));
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
sgtitle(sprintf('rseed %d lambda %0.2f beta %0.2f sigma %0.2f',rseed,lambda,beta,sigma),'Interpreter','none');
drawnow();
end;%if flag_disp>0;
end;%if (niteration==n_iteration-1) | (mod(niteration,32)==0);
end;%for niteration=0:n_iteration-1;
%%%%;
Error_sher_l2_ilbsr_____(:,1+nlambda,1+nbeta,1+nsigma,1+nrseed) = Error_sher_l2_i_;
Error_zero_l2_ilbsr_____(:,1+nlambda,1+nbeta,1+nsigma,1+nrseed) = Error_zero_l2_i_;
Error_emax_l2_ilbsr_____(:,1+nlambda,1+nbeta,1+nsigma,1+nrseed) = Error_emax_l2_i_;
Error_alte_l2_ilbsr_____(:,1+nlambda,1+nbeta,1+nsigma,1+nrseed) = Error_alte_l2_i_;
Error_sher_k1_ilbsr_____(:,1+nlambda,1+nbeta,1+nsigma,1+nrseed) = Error_sher_k1_i_;
Error_zero_k1_ilbsr_____(:,1+nlambda,1+nbeta,1+nsigma,1+nrseed) = Error_zero_k1_i_;
Error_emax_k1_ilbsr_____(:,1+nlambda,1+nbeta,1+nsigma,1+nrseed) = Error_emax_k1_i_;
Error_alte_k1_ilbsr_____(:,1+nlambda,1+nbeta,1+nsigma,1+nrseed) = Error_alte_k1_i_;
%%%%;
if flag_disp>0;
close(gcf);
end;%if flag_disp>0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nsigma=0:n_sigma-1;
end;%for nlambda=0:n_lambda-1;
end;%for nbeta=0:n_beta-1;
%%%%%%%%;
% temporary save. ;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% nrseed %d/%d" saving %s',nrseed,n_rseed,fname_tmp_mat)); end;
save(fname_tmp_mat ...
     ,'q_max' ...
     ,'n_M' ...
     ,'n_sigma' ...
     ,'n_lambda' ...
     ,'n_beta' ...
     ,'n_iteration' ...
     ,'n_rseed' ...
     ,'str_xfix' ...
     ,'A_qr__' ...
     ,'Error_sher_l2_ilbsr_____' ...
     ,'Error_zero_l2_ilbsr_____' ...
     ,'Error_emax_l2_ilbsr_____' ...
     ,'Error_alte_l2_ilbsr_____' ...
     ,'Error_sher_k1_ilbsr_____' ...
     ,'Error_zero_k1_ilbsr_____' ...
     ,'Error_emax_k1_ilbsr_____' ...
     ,'Error_alte_k1_ilbsr_____' ...
     );
end;%for nrseed=0:n_rseed-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
save(fname_mat ...
     ,'q_max' ...
     ,'n_M' ...
     ,'n_sigma' ...
     ,'n_lambda' ...
     ,'n_beta' ...
     ,'n_iteration' ...
     ,'n_rseed' ...
     ,'str_xfix' ...
     ,'A_qr__' ...
     ,'Error_sher_l2_ilbsr_____' ...
     ,'Error_zero_l2_ilbsr_____' ...
     ,'Error_emax_l2_ilbsr_____' ...
     ,'Error_alte_l2_ilbsr_____' ...
     ,'Error_sher_k1_ilbsr_____' ...
     ,'Error_zero_k1_ilbsr_____' ...
     ,'Error_emax_k1_ilbsr_____' ...
     ,'Error_alte_k1_ilbsr_____' ...
     );
close_fname_tmp(fname_pre);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_recalc | ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if  exist(fname_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
load(fname_mat);

fname_fig = sprintf('../dir_jpg/test_MSA_sheres_%s_FIGC',str_xfix);
if ( flag_replot | ~exist(sprintf('%s.jpg',fname_fig)) );
Error_sher_k1_bs__ = abs(squeeze( test_MSA_sheres_error_mean_0(Error_sher_k1_ilbsr_____(end,:,:,:,:),5) )).^2 / (n_gamma);
Error_zero_k1_bs__ = abs(squeeze( test_MSA_sheres_error_mean_0(Error_zero_k1_ilbsr_____(end,:,:,:,:),5) )).^2 / (n_gamma);
Error_emax_k1_bs__ = abs(squeeze( test_MSA_sheres_error_mean_0(Error_emax_k1_ilbsr_____(end,:,:,:,:),5) )).^2 / (n_gamma);
Error_alte_k1_bs__ = abs(squeeze( test_MSA_sheres_error_mean_0(Error_alte_k1_ilbsr_____(end,:,:,:,:),5) )).^2 / (n_gamma);
%Elim_ = prctile(Error_emax_k1_bs__(:),[  0,100]);
%Elim_ = mean(Elim_) + 1.125*0.5*diff(Elim_)*[-1,+1];
Elim_ = 1.0 + 1.0*[-1,+1];
%Elim_ = 2.0 + 2.0*[-1,+1];
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,3*768,512]);figbeach;
p_row = 1; p_col = 4; np=0;
for np=0:3-1;
if np==0; tmp_Error_k1_bs__ = min(Error_sher_k1_bs__,Error_zero_k1_bs__); tmp_str = 'min(sher,zero)'; end;
if np==1; tmp_Error_k1_bs__ = Error_emax_k1_bs__; tmp_str = 'emax'; end;
if np==2; tmp_Error_k1_bs__ = Error_alte_k1_bs__; tmp_str = 'alte'; end;
subplot(p_row,p_col,1+np);
imagesc(tmp_Error_k1_bs__,Elim_);
set(gca,'Ydir','normal');
set(gca,'XTick',1:n_sigma,'XTickLabel',num2str(sigma_,2)); xtickangle(90);
%set(gca,'YTick',1:n_lambda,'YTickLabel',num2str(lambda_,3));
set(gca,'YTick',1:n_beta,'YTickLabel',num2str(beta_,3));
set(gca,'TickLength',[0,0]);
xlabel('sigma');
ylabel('beta');%ylabel('lambda');
title(sprintf('Error_%s_k1_bs__',tmp_str),'Interpreter','none');
colorbar;
set(gca,'FontSize',16);
end;%for np=0:p_row*p_col-1;
subplot(p_row,p_col,p_row*p_col);
set(gca,'FontSize',16);
psi_ = linspace(-pi,+pi,1024+1);
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
linewidth_use = 3;
hold on;
for nbeta=0:n_beta-1;
beta = beta_(1+nbeta);
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*nbeta/n_beta)));
plot(psi_,exp(beta*cos(psi_))/besseli(0,beta),'-','LineWidth',linewidth_use,'Color',c_80s__(1+nc_80s,:));
end;%for nbeta=0:n_beta-1;
xlim([-pi,+pi]); xlabel('$\psi$','Interpreter','latex');
ylim([-0.25,5.00]);
title('Distribution');
%%%%;
sgtitle(fname_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
end;%if ( flag_replot | ~exist(sprintf('%s.jpg',fname_fig)) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if  exist(fname_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;




