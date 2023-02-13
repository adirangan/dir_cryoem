function ...
[ ...
 parameter ...
,gamma_ ...
] = ...
test_heterogeneity_spurious_helper_0( ...
 parameter ...
,n_w_sum ...
,n_S ...
,gamma_true_S_ ...
,S_k_p_wkS__ ...
,weight_2d_k_all_ ...
);

str_thisfunction = 'test_heterogeneity_spurious_helper_0';

if nargin< 1;
disp(sprintf(' %% testing %s',str_thisfunction));
n_w_sum = 1;
n_S = 1024;
gamma_true_S_ = linspace(0,2*pi,1+n_S); gamma_true_S_ = transpose(gamma_true_S_(1:n_S));
S_k_p_wkS__ = reshape(sin(2*gamma_true_S_) + i*cos(gamma_true_S_),[n_w_sum,n_S]);
parameter = struct('type','parameter');
parameter.n_gamma = 128;
parameter.q_max = 16;
parameter.n_iteration = 128;
parameter.flag_verbose = 1;
parameter.flag_display = 1;
test_heterogeneity_spurious_helper_0( ...
 parameter ...
,n_w_sum ...
,n_S ...
,gamma_true_S_ ...
,S_k_p_wkS__ ...
);
disp('returning'); return;
end;%if nargin< 1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_w_sum=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); gamma_true_S_=[]; end; na=na+1;
if (nargin<1+na); S_k_p_wkS__=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_all_=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
tolerance_master = parameter.tolerance_master;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'flag_display'); parameter.flag_display = 0; end;
flag_display = parameter.flag_display;
if ~isfield(parameter,'n_iteration'); parameter.n_iteration = 32; end;
n_iteration = parameter.n_iteration;
if ~isfield(parameter,'q_max'); parameter.q_max = 48; end;
q_max = parameter.q_max;
if ~isfield(parameter,'n_gamma'); parameter.n_gamma = parameter.q_max; end;
n_gamma = parameter.n_gamma;
if ~isfield(parameter,'rseed'); parameter.rseed = 0; end;
rseed = parameter.rseed; 

if (flag_verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

rng(rseed);

if isempty(gamma_true_S_); gamma_true_S_ = linspace(0,2*pi,1+n_S); gamma_true_S_ = transpose(gamma_true_S_(1:n_S)); end;
if isempty(weight_2d_k_all_); weight_2d_k_all_ = ones(n_w_sum,1); end;

%%%%%%%%;
% Here we change notation: ;
% k <- wk (compressing all information within each template). ;
% w <-- S (indexing templates, ordered by azimuthal angle, by w). ;
% M <-- S (indexing images, also ordered by azimuthal angle w). ;
%%%%%%%%;
A_0in_Mk__ = transpose(bsxfun(@times,S_k_p_wkS__,sqrt(weight_2d_k_all_))); %<-- interpreted as images. ;
if flag_display;
[U_A_M1__,S_A_11__,V_A_1k__] = svds(A_0in_Mk__,1);
end;%if flag_display;
n_M = size(A_0in_Mk__,1);

%%%%%%%%;
q_ = transpose(-q_max:1:+q_max);
n_q = numel(q_);
%%%%%%%%;

%%%%%%%%;
gamma_ = linspace(0,2*pi,1+n_gamma); gamma_ = reshape(gamma_(1:n_gamma),[n_gamma,1]);
n_w = n_gamma;
F_wq__ = exp(+i*gamma_*transpose(q_));
F_inv_qw__ = ctranspose(F_wq__)/n_gamma;
%%%%%%%%;

index_gamma_tru_M_ = transpose(round(linspace(0,n_gamma-1,n_M)));

B_qk__ = zeros(n_q,n_w_sum);
B_qk__(1+q_max-1,1) = -1;
B_qk__(1+q_max+1,1) = +1;
B_wk__ = F_wq__*B_qk__;

if flag_display;
figure(1);clf;figmed;
subplot(1,2,1); local_display(A_0in_Mk__,V_A_1k__);
subplot(1,2,2); local_display(B_wk__,V_A_1k__);
end;%if flag_display;
%%%%%%%%;
for niteration=0:n_iteration-1;
%%%%%%%%;

R2_wM__ = calculate_R2_wM(B_wk__,A_0in_Mk__);
[~,tmp_ij_wM__] = sort(R2_wM__,2);
p_wM__ = sparse([],[],[],n_w,n_M,0);
n_M_factor_use = ceil(n_M/max(1,n_w));
for nM_factor_use=0:n_M_factor_use;
tmp_index_ = tmp_ij_wM__(:,1+nM_factor_use)-1;
p_wM__ = p_wM__ + sparse(1:n_w,1+tmp_index_,1/n_M_factor_use,n_w,n_M);
end;%for nM_factor_use=0:n_M_factor_use-1;
C_qk__ = F_inv_qw__*(p_wM__*A_0in_Mk__);
B_qk__ = C_qk__;
B_wk__ = F_wq__*B_qk__;

%{
R2_wM__ = calculate_R2_wM(B_wk__,A_0in_Mk__);
[~,ij_wM__] = sort(R2_wM__,2,'ascend'); index_wM__ = ij_wM__-1;
index_w_ = zeros(n_w,1);
flag_M_ = zeros(n_M,1);
nw_use_ = randperm(n_w)-1; nw_use_ = 0:n_w-1;
flag_val = 0;
for nw=0:n_w-1;
nw_use = nw_use_(1+nw);
index_M_ = index_wM__(1+nw_use,:);
tmp_flag_M_ = flag_M_(1+index_M_);
tmp_nM = index_M_(find(tmp_flag_M_==flag_val,1,'first'));
index_w_(1+nw_use) = tmp_nM;
flag_M_(1+tmp_nM) = flag_M_(1+tmp_nM) + 1;
if (nw==n_M-1); flag_val = flag_val+1; end;
end;%for nw=0:n_w-1;
index_gamma_est_w_ = index_w_;
%}

%{
subplot(1,3,1); cla;
plot(index_gamma_est_w_,'.'); 
subplot(1,3,2); cla;
tmp_B_w_ = B_wk__*transpose(V_A_1k__);
tmp_A_M_ = A_0in_Mk__*transpose(V_A_1k__);
hold on;
plot(real(tmp_B_w_),imag(tmp_B_w_),'ko');
plot(real(tmp_A_M_),imag(tmp_A_M_),'r.');
for nw=0:8:n_w-1;%for nw=0:n_w-1;
plot([real(tmp_B_w_(1+nw));real(tmp_A_M_(1+index_gamma_est_w_(1+nw)))],[imag(tmp_B_w_(1+nw));imag(tmp_A_M_(1+index_gamma_est_w_(1+nw)))],'k-');
end;%for nw=0:n_w-1;
hold off;
B_qk__ = F_inv_qw__*A_0in_Mk__(1+index_gamma_est_w_,:); %<-- maximum-entropy using quadrature. ;
B_wk__ = F_wq__*B_qk__;
subplot(1,3,3); cla;
tmp_B_w_ = B_wk__*transpose(V_A_1k__);
tmp_A_M_ = A_0in_Mk__*transpose(V_A_1k__);
hold on;
plot(real(tmp_B_w_),imag(tmp_B_w_),'ko');
plot(real(tmp_A_M_),imag(tmp_A_M_),'r.');
hold off;
error('stopping');
%}

%{
B_qk__ = F_inv_qw__*A_0in_Mk__(1+index_gamma_est_w_,:); %<-- maximum-entropy using quadrature. ;
B_wk__ = F_wq__*B_qk__;
%%%%%%%%;
R2_wM__ = calculate_R2_wM(B_wk__,A_0in_Mk__);
[~,ij_gamma_est_M_] = min(R2_wM__,[],1); %<-- match templates to images. ;
index_gamma_est_M_ = ij_gamma_est_M_-1;
pinv_qM__ = pinv(F_wq__(1+index_gamma_est_M_,:),tolerance_master);
B_qk__ = pinv_qM__*A_0in_Mk__; %<-- maximum likelihood using least-squares. ;
B_wk__ = F_wq__*B_qk__;
%%%%%%%%;
%}
if flag_display;
subplot(1,2,2); cla; local_display(B_wk__,V_A_1k__); drawnow();
end;%if flag_display;
%%%%%%%%;
end;%for niteration=0:n_iteration-1;

error('stopping');

if (flag_verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

A_w_ = cos(gamma_) + i*2*sin(2*gamma_);
R2_ww__ = bsxfun(@plus,reshape(sum(abs(A_w_).^2,2),[n_w,1]),reshape(sum(abs(A_w_).^2,2),[1,n_w])) - 2*real(A_w_*ctranspose(A_w_));

flag_disp = 2; nf=0;
%%%%%%%%%%%%%%%%;
q_max = 8;
n_M_factor = 1/32; n_gamma = 1024; n_M = ceil(n_gamma*n_M_factor);
n_iteration = 128*1;
sigma_true = 0.1;
lambda = 0.0;
%%%%%%%%;
q_ = transpose(-q_max:1:+q_max);
n_q = numel(q_);
gamma_ = linspace(0,2*pi,1+n_gamma);
gamma_ = reshape(gamma_(1:n_gamma),[n_gamma,1]);
dgamma = mean(diff(gamma_));
n_w = n_gamma;
F_wq__ = exp(+i*gamma_*transpose(q_));
F_inv_qw__ = ctranspose(F_wq__)/n_gamma;
%%%%%%%%;
rseed = 0;
rng(rseed);rseed=rseed+1;
%A_q_ = randn(n_q,1) + i*randn(n_q,1);
%A_w_ = F_wq__*A_q_;
A_w_ = cos(gamma_) + i*2*sin(2*gamma_);
A_q_ = F_inv_qw__*A_w_;
%%%%;
rng(rseed);rseed=rseed+1;
gamma_true_M_ = 2*pi*rand(n_M,1); %<-- uniform distribution. ;
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
[B_emax_w_] = MSA_shape_speed1_1(1,n_w,gamma_,B_emax_w_); B_emax_q_ = F_inv_qw__*B_emax_w_;
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
if flag_disp>0;
%close(gcf);
end;%if flag_disp>0;
%%%%%%%%%%%%%%%%;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

function R2_wM__ = calculate_R2_wM(B_wk__,A_Mk__);
%%%%%%%%;
% |B_w_-A_M_|^2 = |B_w_|.^2 + |A_M_|.^2 - 2<B_w_,A_M_> ;
%%%%%%%%;
n_w = size(B_wk__,1); n_M = size(A_Mk__,1); 
R2_wM__ = bsxfun(@plus,reshape(sum(abs(B_wk__).^2,2),[n_w,1]),reshape(sum(abs(A_Mk__).^2,2),[1,n_M])) - 2*real(B_wk__*ctranspose(A_Mk__));

function local_display(A_Mk__,V_1k__,c__);
na=0;
if (nargin<1+na); A_Mk__=[]; end; na=na+1;
if (nargin<1+na); V_1k__=[]; end; na=na+1;
if (nargin<1+na); c__=[]; end; na=na+1;
if isempty(V_1k__); [~,~,V_1k__] = svds(A_Mk__,1); end;
if isempty(c__); c__ = colormap('hsv'); end;
colormap(c__);
n_M = size(A_Mk__,1);
if ~isempty(A_Mk__);
tmp_A_M_ = A_Mk__*transpose(V_1k__); tmp_gamma_M_ = linspace(0,2*pi,1+n_M); tmp_gamma_M_ = transpose(tmp_gamma_M_(1:n_M));
hold on; 
surfline_0(real(tmp_A_M_),imag(tmp_A_M_),tmp_gamma_M_);
hold off;
end;%if ~isempty(A_Mk__);
