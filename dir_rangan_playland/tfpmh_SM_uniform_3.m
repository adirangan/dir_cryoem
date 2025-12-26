function  ...
[ ...
 parameter ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,image_delta_x_M_ ...
,image_delta_y_M_ ...
,image_I_value_M_ ...
,image_X_value_M_ ...
,image_S_index_M_ ...
] = ...
tfpmh_SM_uniform_3( ...
 parameter ...
,n_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,n_M ...
,X_SM__ ...
,delta_x_SM__ ...
,delta_y_SM__ ...
,gamma_z_SM__ ...
,I_value_SM__ ...
);

%%%%%%%%;
if nargin<1;
flag_verbose=1; flag_disp=0;nf=0;
parameter = struct('type','parameter');
parameter.flag_deterministic = 1;
parameter.f_rand = 0.15;
n_S=113;
viewing_azimu_b_S_=transpose(linspace(0,2*pi,n_S));
viewing_polar_a_S_=transpose(linspace(0,1*pi,n_S));
n_M=117;
%%%%;
n_SM = n_S*n_M;
X_SM__ = reshape(mod([0:n_SM-1],sqrt(19)),[n_S,n_M]);
delta_x_SM__ = reshape(mod([0:n_SM-1],sqrt(23)),[n_S,n_M]);
delta_y_SM__ = reshape(mod([0:n_SM-1],sqrt(29)),[n_S,n_M]);
gamma_z_SM__ = reshape(mod([0:n_SM-1],sqrt(31)),[n_S,n_M]);
I_value_SM__ = reshape(mod([0:n_SM-1],sqrt(37)),[n_S,n_M]);
tmp_t=tic();
[ ...
 parameter ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,image_delta_x_M_ ...
,image_delta_y_M_ ...
,image_I_value_M_ ...
,image_X_value_M_ ...
,image_S_index_M_ ...
] = ...
tfpmh_SM_uniform_3( ...
 parameter ...
,n_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,n_M ...
,X_SM__ ...
,delta_x_SM__ ...
,delta_y_SM__ ...
,gamma_z_SM__ ...
,I_value_SM__ ...
);
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% tfpmh_SM_uniform_3: time %0.6fs',tmp_t)); end;
%%%%;
dir_pymat = '/data/rangan/dir_cryoem/dir_rangan_python/dir_pymat' ;
fname_pymat = sprintf('%s/test_tfpmh_SM_uniform_3_SM.mat',dir_pymat);
if ~exist(fname_pymat,'file'); disp(sprintf(' %% %s not found',fname_pymat)); end;
if  exist(fname_pymat,'file');
tmp_ = load(fname_pymat);
fnorm_disp(flag_verbose,'viewing_azimu_b_S_',viewing_azimu_b_S_,'tmp_.viewing_azimu_b_S_',tmp_.viewing_azimu_b_S_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'viewing_polar_a_S_',viewing_polar_a_S_,'tmp_.viewing_polar_a_S_',tmp_.viewing_polar_a_S_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'X_SM__',X_SM__,'tmp_.X_SM__',tmp_.X_SM__,' %%<-- should be zero');
fnorm_disp(flag_verbose,'delta_x_SM__',delta_x_SM__,'tmp_.delta_x_SM__',tmp_.delta_x_SM__,' %%<-- should be zero');
fnorm_disp(flag_verbose,'delta_y_SM__',delta_y_SM__,'tmp_.delta_y_SM__',tmp_.delta_y_SM__,' %%<-- should be zero');
fnorm_disp(flag_verbose,'gamma_z_SM__',gamma_z_SM__,'tmp_.gamma_z_SM__',tmp_.gamma_z_SM__,' %%<-- should be zero');
fnorm_disp(flag_verbose,'I_value_SM__',I_value_SM__,'tmp_.I_value_SM__',tmp_.I_value_SM__,' %%<-- should be zero');
fnorm_disp(flag_verbose,'euler_polar_a_M_',euler_polar_a_M_,'tmp_.euler_polar_a_M_',tmp_.euler_polar_a_M_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'euler_azimu_b_M_',euler_azimu_b_M_,'tmp_.euler_azimu_b_M_',tmp_.euler_azimu_b_M_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'euler_gamma_z_M_',euler_gamma_z_M_,'tmp_.euler_gamma_z_M_',tmp_.euler_gamma_z_M_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'image_delta_x_M_',image_delta_x_M_,'tmp_.image_delta_x_M_',tmp_.image_delta_x_M_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'image_delta_y_M_',image_delta_y_M_,'tmp_.image_delta_y_M_',tmp_.image_delta_y_M_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'image_I_value_M_',image_I_value_M_,'tmp_.image_I_value_M_',tmp_.image_I_value_M_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'image_X_value_M_',image_X_value_M_,'tmp_.image_X_value_M_',tmp_.image_X_value_M_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'image_S_index_M_',image_S_index_M_,'double(tmp_.image_S_index_M_)',double(tmp_.image_S_index_M_),' %%<-- should be zero');
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 4; np=0;
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(euler_polar_a_M_-tmp_.euler_polar_a_M_,'o'); title(sprintf('euler_polar_a_M_'),'Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(euler_azimu_b_M_-tmp_.euler_azimu_b_M_,'o'); title(sprintf('euler_azimu_b_M_'),'Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(euler_gamma_z_M_-tmp_.euler_gamma_z_M_,'o'); title(sprintf('euler_gamma_z_M_'),'Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(image_delta_x_M_-tmp_.image_delta_x_M_,'o'); title(sprintf('image_delta_x_M_'),'Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(image_delta_y_M_-tmp_.image_delta_y_M_,'o'); title(sprintf('image_delta_y_M_'),'Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(image_I_value_M_-tmp_.image_I_value_M_,'o'); title(sprintf('image_I_value_M_'),'Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(image_X_value_M_-tmp_.image_X_value_M_,'o'); title(sprintf('image_X_value_M_'),'Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(image_S_index_M_-double(tmp_.image_S_index_M_),'o'); title(sprintf('image_S_index_M_'),'Interpreter','none');
end;%if flag_disp;
end;%if  exist(fname_pymat,'file');
%%%%;
if flag_disp;
n_test = 3;
for ntest=0:n_test-1;
figure(1+nf);nf=nf+1;clf;
p_row = 2; p_col = 3; np=0;
nM = max(0,min(n_M-1,floor(n_M*(ntest+0.5)/max(1,n_test))));
euler_polar_a = euler_polar_a_M_(1+nM);
euler_azimu_b = euler_azimu_b_M_(1+nM);
euler_gamma_z = euler_gamma_z_M_(1+nM);
image_delta_x = image_delta_x_M_(1+nM);
image_delta_y = image_delta_y_M_(1+nM);
image_I_value = image_I_value_M_(1+nM);
image_X_value = image_X_value_M_(1+nM);
image_S_index = image_S_index_M_(1+nM);
X_S_ = X_SM__(:,1+nM);
delta_x_S_ = delta_x_SM__(:,1+nM);
delta_y_S_ = delta_y_SM__(:,1+nM);
gamma_z_S_ = gamma_z_SM__(:,1+nM);
I_value_S_ = I_value_SM__(:,1+nM);
[X_srt_S_,ij_srt_S_] = sort(X_S_,'descend'); [~,ij_trs_S_] = sort(ij_srt_S_,'ascend');
delta_x_srt_S_ = delta_x_S_(ij_srt_S_);
delta_y_srt_S_ = delta_y_S_(ij_srt_S_);
gamma_z_srt_S_ = gamma_z_S_(ij_srt_S_);
I_value_srt_S_ = I_value_S_(ij_srt_S_);
azimu_b_srt_S_ = viewing_azimu_b_S_(ij_srt_S_);
polar_a_srt_S_ = viewing_polar_a_S_(ij_srt_S_);
ijS_srt = ij_trs_S_(1+image_S_index);
subplot(p_row,p_col,1+np);np=np+1; hold on; plot(1:n_S,X_srt_S_,'kx'); plot(ijS_srt,image_X_value,'ro'); xlim([0,n_S+1]); title('X');
subplot(p_row,p_col,1+np);np=np+1; hold on; plot(1:n_S,delta_x_srt_S_,'kx'); plot(ijS_srt,image_delta_x,'ro'); xlim([0,n_S+1]); title('delta_x');
subplot(p_row,p_col,1+np);np=np+1; hold on; plot(1:n_S,delta_y_srt_S_,'kx'); plot(ijS_srt,image_delta_y,'ro'); xlim([0,n_S+1]); title('delta_y');
subplot(p_row,p_col,1+np);np=np+1; hold on; plot(1:n_S,gamma_z_srt_S_,'kx'); plot(ijS_srt,euler_gamma_z,'ro'); xlim([0,n_S+1]); title('gamma_z');
subplot(p_row,p_col,1+np);np=np+1; hold on; plot(1:n_S,azimu_b_srt_S_,'kx'); plot(ijS_srt,euler_azimu_b,'ro'); xlim([0,n_S+1]); title('azimu_b');
subplot(p_row,p_col,1+np);np=np+1; hold on; plot(1:n_S,polar_a_srt_S_,'kx'); plot(ijS_srt,euler_polar_a,'ro'); xlim([0,n_S+1]); title('polar_a');
sgtitle(sprintf('nM %d/%d',nM,n_M),'Interpreter','none');
end;%for ntest=0:n_test-1;
end;%if flag_disp;
%%%%;
disp('returning'); return;
end;%if nargin<1;
%%%%%%%%;

%%%%%%%%;
na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); viewing_azimu_b_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_S_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); X_SM__=[]; end; na=na+1;
if (nargin<1+na); delta_x_SM__=[]; end; na=na+1;
if (nargin<1+na); delta_y_SM__=[]; end; na=na+1;
if (nargin<1+na); gamma_z_SM__=[]; end; na=na+1;
if (nargin<1+na); I_value_SM__=[]; end; na=na+1;
%%%%%%%%;
if isempty(parameter); parameter = struct('type','parameter'); end;%if isempty(parameter);
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'flag_deterministic'); parameter.flag_deterministic = 0; end;
flag_deterministic = parameter.flag_deterministic;
if isempty(I_value_SM__); I_value_SM__ = ones(size(X_SM__)); end;
%%%%%%%%;
if (~isfield(parameter,'f_rand')); parameter.f_rand = 0.05; end; %<-- parameter_bookmark. ;
f_rand = parameter.f_rand;
%%%%%%%%;
euler_polar_a_M_ = zeros(n_M,1);
euler_azimu_b_M_ = zeros(n_M,1);
euler_gamma_z_M_ = zeros(n_M,1);
image_delta_x_M_ = zeros(n_M,1);
image_delta_y_M_ = zeros(n_M,1);
image_I_value_M_ = zeros(n_M,1);
image_X_value_M_ = zeros(n_M,1);
image_S_index_M_ = zeros(n_M,1);
%%%%%%%%;
assert(ndims(X_SM__)==2);
%%%%%%%%;
[X_srt_SM__,ij_srt_SM__] = sort(X_SM__,1,'descend'); index_srt_SM__ = ij_srt_SM__-1;
if (f_rand<=0); index_nS_srt_M_ = zeros(n_M,1); end;
if (f_rand> 0);
index_nS_srt_M_ = floor(f_rand*0.5*n_S).*ones(n_M,1);
if flag_deterministic==0;
index_nS_srt_M_ = floor(f_rand*rand(n_M,1)*n_S).*ones(n_M,1);
end;%if flag_deterministic==0;
end;%if (f_rand> 0);
index_nS_M_ = index_srt_SM__(1+index_nS_srt_M_ + transpose(0:n_M-1)*n_S);
index_nSM_ = index_nS_M_ + transpose(0:n_M-1)*n_S;
euler_polar_a_M_ = viewing_polar_a_S_(1+index_nS_M_);
euler_azimu_b_M_ = viewing_azimu_b_S_(1+index_nS_M_);
euler_gamma_z_M_ = gamma_z_SM__(1+index_nSM_);
image_delta_x_M_ = delta_x_SM__(1+index_nSM_);
image_delta_y_M_ = delta_y_SM__(1+index_nSM_);
image_X_value_M_ = X_SM__(1+index_nSM_);
image_S_index_M_ = index_nS_M_;
image_I_value_M_ = I_value_SM__(1+index_nSM_);
