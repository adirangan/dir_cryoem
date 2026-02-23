function ...
[ ...
 parameter ...
,M_k_q_k_outp_ ...
,fA_form_orig_z_ ...
] = ...
fht_2( ...
 parameter ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,weight_2d_x_p_r_ ...
,M_x_q_x_orig_ ...
,q_val ...
,flag_sign ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_2d_k_p_r_ ...
,M_k_q_k_refe_ ...
,fA_form_orig_z_ ...
);

str_thisfunction = 'fht_2';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
disp(sprintf(' %% see test_bb_3.m'));
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_x_p_r=[]; end; na=na+1;
if (nargin<1+na); x_p_r_=[]; end; na=na+1;
if (nargin<1+na); x_p_r_max=[]; end; na=na+1;
if (nargin<1+na); weight_2d_x_p_r_=[]; end; na=na+1;
if (nargin<1+na); M_x_q_x_orig_=[]; end; na=na+1;
if (nargin<1+na); q_val=[]; end; na=na+1;
if (nargin<1+na); flag_sign=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); M_k_q_k_refe_=[]; end; na=na+1;
if (nargin<1+na); fA_form_orig_z_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master=1e-12; end;
tolerance_master = parameter.tolerance_master;
if ~isfield(parameter,'tolerance_nufft'); parameter.tolerance_nufft=tolerance_master; end;
tolerance_nufft = parameter.tolerance_nufft;

if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp=0; end;
flag_disp=parameter.flag_disp; nf=0;
if ~isfield(parameter,'alpha'); parameter.alpha=1.0; end;
alpha=parameter.alpha;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(q_val); q_val = 0; end;
if isempty(flag_sign); flag_sign = +1; end;
flag_sign_n = -1;

if (flag_verbose>0);
disp(sprintf(' %% n_x_p_r %d x_p_r_max %0.2f',n_x_p_r,x_p_r_max));
disp(sprintf(' %% n_k_p_r %d k_p_r_max %0.2f',n_k_p_r,k_p_r_max));
end;%if (flag_verbose>0);

%%%%%%%%;
% Bessel. ;
%%%%%%%%;
if abs(q_val)==0; n_val = flag_sign_n*0.75; end;
if abs(q_val)> 0; n_val = flag_sign_n*1.00; end;
if (flag_verbose>0); disp(sprintf(' %% q_val %0.2f n_val %0.2f',q_val,n_val)); end;

Q = @(q_val,mu_val) (-1).^(q_val .* (q_val<0)) .* exp( mu_val*log(2) - (1+mu_val)*log(2*pi) + gamma_godfrey_1(0.5*(1 + abs(q_val) + mu_val)) - gamma_godfrey_1(0.5*(1 + abs(q_val) - mu_val)) ) ;

%%%%%%%%;
% calculate M-integral. ;
%%%%%%%%;
lk_lim_ = [log(min(k_p_r_)),log(max(k_p_r_))]/alpha; lk_dia = diff(lk_lim_);
lx_lim_ = [log(min(x_p_r_)),log(max(x_p_r_))]/alpha; lx_dia = diff(lx_lim_);
nyquist_factor = 2.0;
n_gamma_z_pre = max(1024+1,n_x_p_r);
gamma_z_max = max(1,n_gamma_z_pre)/max(1e-12,diff(lx_lim_))*nyquist_factor; %<-- first ensure that gamma_z_max is sufficient to resolve lx_. ;
n_gamma_z_pos = round(gamma_z_max * max(1e-12,diff(lk_lim_)) / nyquist_factor) ; %<-- now ensure that n_gamma_z is sufficient to resolve lk_. ;
n_gamma_z = max(n_gamma_z_pre,n_gamma_z_pos) ;
if (flag_verbose>0); disp(sprintf(' %% n_gamma_z_pre %d gamma_z_max %0.2f n_gamma_z_pos %d ',n_gamma_z_pre,gamma_z_max,n_gamma_z_pos)); end;
%%%%;
if (gamma_z_max >= 1.05*n_gamma_z/max(1e-12,lx_dia)*nyquist_factor);
disp(sprintf(' %% Warning, lx_dia %0.2f vs n_gamma_z/max(1e-12,gamma_z_max)*nyquist_factor %0.2f in %s',lx_dia,n_gamma_z/max(1e-12,gamma_z_max)*nyquist_factor,str_thisfunction));
end;%if (gamma_z_max >= n_gamma_z/max(1e-12,lx_dia)*nyquist_factor);
%%%%;
if (gamma_z_max >= 1.05*n_gamma_z/max(1e-12,lk_dia)*nyquist_factor);
disp(sprintf(' %% Warning, lk_dia %0.2f vs n_gamma_z/max(1e-12,gamma_z_max)*nyquist_factor %0.2f in %s',lk_dia,n_gamma_z/max(1e-12,gamma_z_max)*nyquist_factor,str_thisfunction));
end;%if (gamma_z_max >= n_gamma_z/max(1e-12,lk_dia)*nyquist_factor);
%%%%%%%%;
tmp_t=tic();
[gamma_z_,weight_1d_gamma_z_] = chebpts(n_gamma_z,gamma_z_max*[-1,+1]);
tmp_t=toc(tmp_t);
parameter = parameter_timing_update(parameter,'chebpts',tmp_t);
gamma_z_ = reshape(gamma_z_,[n_gamma_z,1]);
weight_1d_gamma_z_ = reshape(weight_1d_gamma_z_,[n_gamma_z,1]);
%%%%%%%%;
fB_nuff_orig_z_ = [];
if flag_disp;
fB_nuff_orig_z_ = finufft1d3(log(max(tolerance_master,x_p_r_))/alpha,M_x_q_x_orig_.*x_p_r_.^(n_val).*weight_2d_x_p_r_,-1,tolerance_nufft,gamma_z_);
end;%if flag_disp;
fB_quad_orig_z_ = [];
if flag_disp;
fB_quad_orig_z_ = zeros(n_gamma_z,1);
for ngamma_z=0:n_gamma_z-1;
gamma_z = gamma_z_(1+ngamma_z);
fB = sum(M_x_q_x_orig_.*x_p_r_.^(n_val).*x_p_r_.^(-i*gamma_z/alpha).*weight_2d_x_p_r_);
fB_quad_orig_z_(1+ngamma_z) = fB;
end;%for ngamma_z=0:n_gamma_z-1;
fnorm_disp(flag_verbose,'fB_nuff_orig_z_',fB_nuff_orig_z_,'fB_quad_orig_z_',fB_quad_orig_z_,'%<-- should be zero');
end;%if flag_disp;
%%%%%%%%;
tmp_t=tic();
n_lx = n_gamma_z; [node_lx_,weight_lx_] = chebpts(n_lx,lx_lim_);
tmp_t=toc(tmp_t);
parameter = parameter_timing_update(parameter,'chebpts',tmp_t);
node_lx_ = reshape(node_lx_,[n_lx,1]); weight_lx_ = reshape(weight_lx_,[n_lx,1]);
tmp_t=tic();
Mxn_x_q_lx_ = interp1(log(x_p_r_)/alpha,M_x_q_x_orig_.*x_p_r_.^(2+n_val),node_lx_,'spline') * (2*pi);
tmp_t=toc(tmp_t);
parameter = parameter_timing_update(parameter,'interp1',tmp_t);
%%%%%%%%;
tmp_t=tic();
fB_nuff_orlx_z_ = finufft1d3(node_lx_,Mxn_x_q_lx_.*weight_lx_,-1,tolerance_nufft,gamma_z_);
tmp_t=toc(tmp_t);
parameter = parameter_timing_update(parameter,'finufft1d3',tmp_t);
%%%%%%%%;

if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figmed;
p_row=1;p_col=3;np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;hold on;
plot(log(max(tolerance_master,x_p_r_))/alpha,M_x_q_x_orig_.*x_p_r_.^(n_val).*weight_2d_x_p_r_,'ro-');
plot(node_lx_,Mxn_x_q_lx_.*weight_lx_,'mx');
xlim(lx_lim_); xlabel('lx_lim_','Interpreter','none');
title('M_x_q_x_orig_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;hold on;
plot(gamma_z_,real(fB_nuff_orig_z_),'ro-');
if ~isempty(fB_quad_orig_z_);
plot(gamma_z_,real(fB_quad_orig_z_),'r.-');
end;%if ~isempty(fB_quad_orig_z_);
plot(gamma_z_,real(fB_nuff_orlx_z_),'mx-');
xlim(gamma_z_max*[-1,+1]); xlabel('gamma_z_','Interpreter','none');
title('real(fB_nuff_orig_z_)','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;hold on;
plot(gamma_z_,log(abs(fB_nuff_orig_z_)),'ro-');
if ~isempty(fB_quad_orig_z_);
plot(gamma_z_,log(abs(fB_quad_orig_z_)),'r.-');
end;%if ~isempty(fB_quad_orig_z_);
plot(gamma_z_,log(abs(fB_nuff_orlx_z_)),'mx-');
xlim(gamma_z_max*[-1,+1]); xlabel('gamma_z_','Interpreter','none');
title('log(abs(fB_nuff_orig_z_))','Interpreter','none');
%%%%;
end;%if flag_disp;

%%%%%%%%;
% Calculate J-integral. ;
%%%%%%%%;
if isempty(fA_form_orig_z_);
tmp_t = tic();
fA_form_orig_z_ = (flag_sign*i).^(q_val) / alpha .* Q(q_val,-(1+n_val+i*gamma_z_/alpha));
tmp_t=toc(tmp_t);
parameter = parameter_timing_update(parameter,'fA_form_orig_z_',tmp_t);
end;%if isempty(fA_form_orig_z_);
%%%%%%%%;
if flag_disp;
fC_comb_orig_z_ = fA_form_orig_z_.*flipud(fB_nuff_orig_z_);
bC_comb_orig_y_ = finufft1d3(gamma_z_,fC_comb_orig_z_.*weight_1d_gamma_z_,+1,tolerance_nufft,log(max(tolerance_master,k_p_r_))/alpha);
M_k_q_k_orig_ = k_p_r_.^n_val .* bC_comb_orig_y_  / (2*pi) ;
end;%if flag_disp;
%%%%%%%%;
fA_form_orlx_z_ = fA_form_orig_z_;
%%%%%%%%;
tmp_t=tic();
fC_comb_orlx_z_ = fA_form_orlx_z_.*flipud(fB_nuff_orlx_z_);
bC_comb_orlx_y_ = finufft1d3(gamma_z_,fC_comb_orlx_z_.*weight_1d_gamma_z_,+1,tolerance_nufft,log(max(tolerance_master,k_p_r_))/alpha);
M_k_q_k_orlx_ = k_p_r_.^n_val .* bC_comb_orlx_y_  / (2*pi) ;
tmp_t=toc(tmp_t);
parameter = parameter_timing_update(parameter,'finufft1d3',tmp_t);
%%%%%%%%;

if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figmed;
p_row=1;p_col=3;np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;hold on;
plot(gamma_z_,real(fA_form_orig_z_),'kx-');
xlim(gamma_z_max*[-1,+1]); xlabel('gamma_z_','Interpreter','none');
title('real(fA_form_orig_z_)','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;hold on;
plot(gamma_z_,real(fC_comb_orig_z_),'ro-');
plot(gamma_z_,real(fC_comb_orlx_z_),'mx-');
xlim(gamma_z_max*[-1,+1]); xlabel('gamma_z_','Interpreter','none');
title('real(fC_comb_orig_z_)','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;hold on;
plot(gamma_z_,log(abs(fC_comb_orig_z_)),'ro-');
plot(gamma_z_,log(abs(fC_comb_orlx_z_)),'mx-');
xlim(gamma_z_max*[-1,+1]); xlabel('gamma_z_','Interpreter','none');
title('log(abs(fC_comb_orig_z_))','Interpreter','none');
%%%%;
end;%if flag_disp;

if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figsml;
p_row=1;p_col=1;np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;hold on;
plot(log(k_p_r_),real(M_k_q_k_orig_),'ro-');
plot(log(k_p_r_),real(M_k_q_k_orlx_),'mx-');
if ~isempty(M_k_q_k_refe_);
plot(log(k_p_r_),real(M_k_q_k_refe_),'k.-');
end;%if ~isempty(M_k_q_k_refe_);
xlim([log(min(k_p_r_)),log(max(k_p_r_))]); xlabel('log(k_p_r_)','Interpreter','none');
title('real(M_k_q_k_orig_)','Interpreter','none');
%%%%;
end;%if flag_disp;

M_k_q_k_outp_ = M_k_q_k_orlx_;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

