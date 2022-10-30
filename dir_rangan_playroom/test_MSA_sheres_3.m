%%%%%%%%;
% A multi-slice-alignment toy-problem for testing sheres-style heterogeneity. ;
% Assumes a 2d-volume x, comprising a single ring. ;
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
flag_disp = 1*flag_verbose;
flag_replot = 1;
nf = 0;

q_max = 12; str_q_max = sprintf('q%d',q_max);
n_gamma = 1024; str_n_gamma = sprintf('g%.4d',n_gamma);
n_sigma = 11; str_n_sigma = sprintf('s%.2d',n_sigma);
n_iteration = 128*1; str_n_iteration = sprintf('i%.3d',n_iteration);
str_xfix = sprintf('%s%s%s%s%s%s',str_q_max,str_n_gamma,str_n_sigma,str_n_iteration);
%%%%%%%%;
q_ = transpose(-q_max:1:+q_max);
n_q = numel(q_);
gamma_w_ = linspace(0,2*pi,1+n_gamma);
gamma_w_ = reshape(gamma_w_(1:n_gamma),[n_gamma,1]);
dgamma = mean(diff(gamma_w_));
sigma_ = transpose(linspace(0,0.1,n_sigma));
n_w = n_gamma;
F_wq__ = exp(+i*gamma_w_*transpose(q_));
F_inv_qw__ = ctranspose(F_wq__)/n_gamma;
%%%%%%%%;

parameter = struct('type','parameter');
parameter.flag_verbose = 0;
parameter.flag_disp = 0;
parameter.str_shape = 'horseshoe';
[ ...
 parameter ...
,q_ ...
,gamma_w_ ...
,F_wq__ ...
,F_inv_qw__ ...
,A_q_ ...
,A_w_ ...
] = ...
MSA_shape_0( ...
 parameter ...
,q_max ...
,n_gamma ...
);

if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
linewidth_use = 4;
plot(real(A_w_),imag(A_w_),'k-','LineWidth',linewidth_use);
end;%if flag_disp;

Error_sher1_l2_is__ = zeros(n_iteration,n_sigma);
Error_sher1_k1_is__ = zeros(n_iteration,n_sigma);
for nsigma=4;%for nsigma=0:n_sigma-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
sigma = sigma_(1+nsigma); sigma_true = sigma;
parameter = struct('type','parameter');
parameter.flag_verbose = 0;
parameter.flag_disp = 0;
parameter.str_shape = 'horseshoe';
parameter.sigma = sigma_true;
[ ...
 parameter ...
,q_ ...
,gamma_w_ ...
,F_wq__ ...
,F_inv_qw__ ...
,A_q_ ...
,A_w_ ...
,x_ ...
,x_0__ ...
,x_1__ ...
,Ag__ ...
] = ...
MSA_shape_0( ...
 parameter ...
,q_max ...
,n_gamma ...
);
n_M = (parameter.n_x).^2;
z_M_ = reshape(x_0__ + i*x_1__,[n_M,1]); %<-- location (i.e., value) of each image-location. ;
p_M_ = reshape(Ag__,[n_M,1]); %<-- density (of images) at each image-location. ;
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
%%%%;
subplot_{1+0} = subplot(1,2,1);
linewidth_use = 8; x_max = parameter.x_max;
s = surfline_0(real(A_w_),imag(A_w_),gamma_w_); set(s,'LineWidth',linewidth_use); colormap('hsv');
xlabel('real(A_w_)','Interpreter','none');
ylabel('imag(A_w_)','Interpreter','none');
xlim(x_max*[-1,+1]);
ylim(x_max*[-1,+1]);
set(gca,'XTick',-x_max:0.1:+x_max);
set(gca,'YTick',-x_max:0.1:+x_max);
xtickangle(90); grid on;
axis square;
title('image ring');
%%%%;
Ilim_ = [0,max(Ag__,[],'all')];
subplot_{1+1} = subplot(1,2,2);
imagesc(Ag__,Ilim_); set(gca,'ydir','normal'); axis image; axisnotick;
xlabel('imag(A_w_)','Interpreter','none');
ylabel('real(A_w_)','Interpreter','none');
title(sprintf('sigma=%0.3f',sigma),'Interpreter','none');
%%%%;
colormap(subplot_{1+0},colormap('hsv'));
colormap(subplot_{1+1},colormap_80s);
end;%if flag_disp;
%%%%;
sigma_sher = 1*sigma_true;
xlim_ = prctile(real(A_w_),[  0,100]); xlim_ = mean(xlim_) + 1.125*0.5*diff(xlim_)*[-1,+1];
ylim_ = prctile(imag(A_w_),[  0,100]); ylim_ = mean(ylim_) + 1.125*0.5*diff(ylim_)*[-1,+1];
rng(0); dA_q_ = randn(n_q,1) + i*randn(n_q,1); dA0_q_ = randn(n_q,1) + i*randn(n_q,1); dA1_q_ = randn(n_q,1) + i*randn(n_q,1);
B_sher_q_ = A_q_ + sigma_true*dA_q_;
B_sher_w_ = F_wq__*B_sher_q_;
B0_sher_q_ = A_q_ + sigma_true*dA0_q_;
B0_sher_w_ = F_wq__*B0_sher_q_;
B1_sher_q_ = A_q_ + sigma_true*dA1_q_;
B1_sher_w_ = F_wq__*B1_sher_q_;
%%%%;
if flag_disp>0;
nf_base = nf; nf=nf+1;
figure(1+nf_base);clf;figmed;
c_use__ = colormap_80s; n_c_use = size(c_use__,1);
linewidth_big = 6;
linewidth_sml = 4;
markersize_use = 4;
subplot(1,3,1);
hold on;
plot(real(A_w_),imag(A_w_),'k-','LineWidth',linewidth_big);
xlabel('real'); ylabel('imag'); grid on;
end;%if flag_disp>0;
%%%%;
Error_sher1_l2_i_ = zeros(n_iteration,1);
Error_sher1_k1_i_ = zeros(n_iteration,1);
for niteration=0:n_iteration-1;
[ ...
 ~ ...
,~ ...
,B_sher_q_ ...
,B_sher_w_ ...
] = ...
MSA_sheres_update_1( ...
 [] ...
,n_w ...
,n_q ...
,F_wq__...
,F_inv_qw__...
,n_M...
,z_M_...
,p_M_...
,sigma_sher...
,1 ...
,[1] ...
,B_sher_q_...
,B_sher_w_...
);
%%%%;
[Error_sher1_l2,Error_sher1_k1] = test_MSA_sheres_error_0(q_max,n_q,q_,F_wq__,n_w,gamma_w_,[],A_w_,A_q_,B_sher_w_,B_sher_q_);
Error_sher1_l2_i_(1+niteration) = Error_sher1_l2;
Error_sher1_k1_i_(1+niteration) = Error_sher1_k1;
%%%%;
if (niteration==n_iteration-1) | (mod(niteration,32)==0);
if flag_disp>0;
figure(1+nf_base); hold on;
nc_use = max(0,min(n_c_use-1,floor(n_c_use*niteration/n_iteration)));
subplot(1,3,1);
plot(real(B_sher_w_),imag(B_sher_w_),'-','LineWidth',linewidth_sml,'Color',c_use__(1+nc_use,:));
xlim(xlim_); ylim(ylim_); grid on;
title(sprintf('sher'));
subplot(1,3,2); cla;
hold on;
plot(1+[0:niteration],log10(Error_sher1_l2_i_(1+[0:niteration])),'ko');
hold off;
title('log10(Error_l2)','Interpreter','none');
xlim([1,n_iteration]);
grid on;
subplot(1,3,3); cla;
hold on;
plot(1+[0:niteration],log10(Error_sher1_k1_i_(1+[0:niteration])),'kx');
hold off;
title('log10(Error_k1)','Interpreter','none');
xlim([1,n_iteration]);
grid on;
%%%%;
sgtitle(sprintf('sigma %0.2f',sigma),'Interpreter','none');
drawnow();
end;%if flag_disp>0;
end;%if (niteration==n_iteration-1) | (mod(niteration,32)==0);
end;%for niteration=0:n_iteration-1;
%%%%;
Error_sher1_l2_is__(:,1+nsigma) = Error_sher1_l2_i_;
Error_sher1_k1_is__(:,1+nsigma) = Error_sher1_k1_i_;
%%%%;
if flag_disp>0;
%close(gcf);
end;%if flag_disp>0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nsigma=0:n_sigma-1;


%%%%%%%%;
% temporary save. ;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% nrseed %d/%d" saving %s',nrseed,n_rseed,fname_tmp_mat)); end;
save(fname_tmp_mat ...
     ,'q_max' ...
     ,'n_M' ...
     ,'n_sigma' ...
     ,'n_iteration' ...
     ,'str_xfix' ...
     ,'Error_sher1_l2_is__' ...
     ,'Error_sher1_k1_is__' ...
     );
end;%for nrseed=0:n_rseed-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
save(fname_mat ...
     ,'q_max' ...
     ,'n_M' ...
     ,'n_sigma' ...
     ,'n_iteration' ...
     ,'str_xfix' ...
     ,'Error_sher1_l2_is__' ...
     ,'Error_sher1_k1_is__' ...
     );
close_fname_tmp(fname_pre);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_recalc | ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if  exist(fname_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
load(fname_mat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if  exist(fname_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;




