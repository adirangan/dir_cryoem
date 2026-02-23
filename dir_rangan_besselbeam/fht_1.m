function ...
[ ...
 parameter ...
,M_k_q_k_outp_ ...
,n_gamma_z ...
,gamma_z_ ...
,gamma_z_max ...
,weight_1d_gamma_z_ ...
,fA_form_orig_z_ ...
,fA_form_appr_z_ ...
,fA_form_orlx_z_ ...
,fA_form_aplx_z_ ...
] = ...
fht_1( ...
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
,n_gamma_z ...
,gamma_z_ ...
,gamma_z_max ...
,weight_1d_gamma_z_ ...
,fA_form_orig_z_ ...
,fA_form_appr_z_ ...
,fA_form_orlx_z_ ...
,fA_form_aplx_z_ ...
,M_k_q_k_refe_ ...
);

str_thisfunction = 'fht_1';

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
if (nargin<1+na); n_gamma_z=[]; end; na=na+1;
if (nargin<1+na); gamma_z_=[]; end; na=na+1;
if (nargin<1+na); gamma_z_max=[]; end; na=na+1;
if (nargin<1+na); weight_1d_gamma_z_=[]; end; na=na+1;
if (nargin<1+na); fA_form_orig_z_=[]; end; na=na+1;
if (nargin<1+na); fA_form_appr_z_=[]; end; na=na+1;
if (nargin<1+na); fA_form_orlx_z_=[]; end; na=na+1;
if (nargin<1+na); fA_form_aplx_z_=[]; end; na=na+1;
if (nargin<1+na); M_k_q_k_refe_=[]; end; na=na+1;

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
if isempty(n_gamma_z); n_gamma_z = 1024+1; end;
if isempty(gamma_z_max); gamma_z_max = n_x_p_r/2 ; end; %<-- nyquist. ;
flag_sign_n = -1;

%%%%%%%%;
% Bessel. ;
%%%%%%%%;
if abs(q_val)==0; n_val_orig = flag_sign_n*0.75; end;
if abs(q_val)> 0; n_val_orig = flag_sign_n*1.00; end;
Q = @(q_val,mu_val) (-1).^(q_val .* (q_val<0)) .* exp( mu_val*log(2) - (1+mu_val)*log(2*pi) + gamma_godfrey_1(0.5*(1 + abs(q_val) + mu_val)) - gamma_godfrey_1(0.5*(1 + abs(q_val) - mu_val)) ) ;
%%%%%%%%;
% calculate M-integral. ;
%%%%%%%%;
if  isempty(gamma_z_) |  isempty(weight_1d_gamma_z_);
[gamma_z_,weight_1d_gamma_z_] = chebpts(n_gamma_z,gamma_z_max*[-1,+1]);
end;%if  isempty(gamma_z_) |  isempty(weight_1d_gamma_z_);
gamma_z_ = reshape(gamma_z_,[n_gamma_z,1]);
weight_1d_gamma_z_ = reshape(weight_1d_gamma_z_,[n_gamma_z,1]);
%%%%%%%%;
fB_nuff_orig_z_ = finufft1d3(log(max(tolerance_master,x_p_r_))/alpha,M_x_q_x_orig_.*x_p_r_.^(n_val_orig).*weight_2d_x_p_r_,-1,tolerance_nufft,gamma_z_);
%%%%%%%%;
if flag_disp;
fB_quad_orig_z_ = zeros(n_gamma_z,1);
for ngamma_z=0:n_gamma_z-1;
gamma_z = gamma_z_(1+ngamma_z);
fB = sum(M_x_q_x_orig_.*x_p_r_.^(n_val_orig).*x_p_r_.^(-i*gamma_z/alpha).*weight_2d_x_p_r_);
fB_quad_orig_z_(1+ngamma_z) = fB;
end;%for ngamma_z=0:n_gamma_z-1;
fnorm_disp(flag_verbose,'fB_nuff_orig_z_',fB_nuff_orig_z_,'fB_quad_orig_z_',fB_quad_orig_z_,'%<-- should be zero');
end;%if flag_disp;
%%%%%%%%;
% Calculate J-integral. ;
%%%%%%%%;
if  isempty(fA_form_orig_z_);
fA_form_orig_z_ = (flag_sign*i).^(q_val) / alpha .* Q(q_val,-(1+n_val_orig+i*gamma_z_/alpha));
end;%if  isempty(fA_form_orig_z_);
%%%%%%%%;
fC_comb_orig_z_ = fA_form_orig_z_.*flipud(fB_nuff_orig_z_);
bC_comb_orig_y_ = finufft1d3(gamma_z_,fC_comb_orig_z_.*weight_1d_gamma_z_,+1,tolerance_nufft,log(max(tolerance_master,k_p_r_))/alpha);
M_k_q_k_orig_ = k_p_r_.^n_val_orig .* bC_comb_orig_y_  / (2*pi) ;
%%%%%%%%;

%%%%%%%%;
% resample in logarithmic-space. ;
%%%%%%%%;
if abs(q_val)==0; n_val_orlx = flag_sign_n*0.75; end;
if abs(q_val)> 0; n_val_orlx = flag_sign_n*1.00; end;
n_lx = n_gamma_z; [node_lx_,weight_lx_] = chebpts(n_lx,[log(min(x_p_r_)),log(max(x_p_r_))]/alpha);
node_lx_ = reshape(node_lx_,[n_lx,1]); weight_lx_ = reshape(weight_lx_,[n_lx,1]);
Mxn_x_q_lx_ = interp1(log(x_p_r_)/alpha,M_x_q_x_orig_.*x_p_r_.^(2+n_val_orlx),node_lx_,'spline') * (2*pi);
%%%%%%%%;
fB_nuff_orlx_z_ = finufft1d3(node_lx_,Mxn_x_q_lx_.*weight_lx_,-1,tolerance_nufft,gamma_z_);
%%%%%%%%;
if flag_disp;
fB_quad_orlx_z_ = zeros(n_gamma_z,1);
for ngamma_z=0:n_gamma_z-1;
gamma_z = gamma_z_(1+ngamma_z);
fB = sum(Mxn_x_q_lx_.*exp(-i*gamma_z.*node_lx_).*weight_lx_);
fB_quad_orlx_z_(1+ngamma_z) = fB;
end;%for ngamma_z=0:n_gamma_z-1;
fnorm_disp(flag_verbose,'fB_nuff_orlx_z_',fB_nuff_orlx_z_,'fB_quad_orlx_z_',fB_quad_orlx_z_,'%<-- should be zero');
end;%if flag_disp;
%%%%%%%%;
% Calculate J-integral. ;
%%%%%%%%;
if  isempty(fA_form_orlx_z_);
fA_form_orlx_z_ = (flag_sign*i).^(q_val) / alpha .* Q(q_val,-(1+n_val_orlx+i*gamma_z_/alpha));
end;%if  isempty(fA_form_orlx_z_);
%%%%%%%%;
fC_comb_orlx_z_ = fA_form_orlx_z_.*flipud(fB_nuff_orlx_z_);
bC_comb_orlx_y_ = finufft1d3(gamma_z_,fC_comb_orlx_z_.*weight_1d_gamma_z_,+1,tolerance_nufft,log(max(tolerance_master,k_p_r_))/alpha);
M_k_q_k_orlx_ = k_p_r_.^n_val_orlx .* bC_comb_orlx_y_  / (2*pi) ;
%%%%%%%%;

%%%%%%%%;
% Approximation. ;
% J_approx = @(x_p_r_) (-1).^(q_val .* (q_val<0)) .* (x_p_r_/2).^(abs(q_val)) .* (1/gamma(1+abs(q_val)) - x_p_r_.^2/(4*gamma(2+abs(q_val))) + x_p_r_.^4/(16*2*gamma(3+abs(q_val)))) .* exp(-b.*x_p_r_.^6) ;
%%%%%%%%;
tmp_b = 1.0;
if abs(q_val)==0; n_val_appr = flag_sign_n*(0.75); end;
if abs(q_val)> 0; n_val_appr = flag_sign_n*(1.00); end;
R = @(mu,nu) exp(gamma_godfrey_1((mu+1)./nu)) ./ max(1e-12,nu.*tmp_b.^((mu+1)./nu)) ./ (2*pi) ;
S = @(q_val,mu_val) ...
 (-1).^(q_val .* (q_val<0)) ...
.* (2*pi).^(-mu_val) .* (1/2).^(abs(q_val)) ...
.* ( ...
     +(1/( 1*exp(gamma_godfrey_1(1+abs(q_val)))) .* R(mu_val+abs(q_val)+0,6)) ...
     -(1/( 4*exp(gamma_godfrey_1(2+abs(q_val)))) .* R(mu_val+abs(q_val)+2,6)) ...
     +(1/(32*exp(gamma_godfrey_1(3+abs(q_val)))) .* R(mu_val+abs(q_val)+4,6)) ...
     ) ...
;
%%%%%%%%;
% calculate M-integral. ;
%%%%%%%%;
fB_nuff_appr_z_ = finufft1d3(log(max(tolerance_master,x_p_r_))/alpha,M_x_q_x_orig_.*x_p_r_.^(n_val_appr).*weight_2d_x_p_r_,-1,tolerance_nufft,gamma_z_);
%%%%%%%%;
if flag_disp;
fB_quad_appr_z_ = zeros(n_gamma_z,1);
for ngamma_z=0:n_gamma_z-1;
gamma_z = gamma_z_(1+ngamma_z);
fB = sum(M_x_q_x_orig_.*x_p_r_.^(n_val_appr).*x_p_r_.^(-i*gamma_z/alpha).*weight_2d_x_p_r_);
fB_quad_appr_z_(1+ngamma_z) = fB;
end;%for ngamma_z=0:n_gamma_z-1;
fnorm_disp(flag_verbose,'fB_nuff_appr_z_',fB_nuff_appr_z_,'fB_quad_appr_z_',fB_quad_appr_z_,'%<-- should be zero');
end;%if flag_disp;
%%%%%%%%;
% Calculate approximate J-integral. ;
%%%%%%%%;
if  isempty(fA_form_appr_z_);
fA_form_appr_z_ = (flag_sign*i).^(q_val) / alpha .* S(q_val,-(1+n_val_appr+i*gamma_z_/alpha));
end;%if  isempty(fA_form_appr_z_);
%%%%%%%%;
fC_comb_appr_z_ = fA_form_appr_z_.*flipud(fB_nuff_appr_z_);
bC_comb_appr_y_ = finufft1d3(gamma_z_,fC_comb_appr_z_.*weight_1d_gamma_z_,+1,tolerance_nufft,log(max(tolerance_master,k_p_r_))/alpha);
M_k_q_k_appr_ = k_p_r_.^n_val_appr .* bC_comb_appr_y_  / (2*pi) ;
%%%%%%%%;

%%%%%%%%;
% resample in logarithmic-space. ;
%%%%%%%%;
if abs(q_val)==0; n_val_aplx = flag_sign_n*0.75; end;
if abs(q_val)> 0; n_val_aplx = flag_sign_n*1.00; end;
n_lx = n_gamma_z; [node_lx_,weight_lx_] = chebpts(n_lx,[log(min(x_p_r_)),log(max(x_p_r_))]/alpha);
node_lx_ = reshape(node_lx_,[n_lx,1]); weight_lx_ = reshape(weight_lx_,[n_lx,1]);
Mxn_x_q_lx_ = interp1(log(x_p_r_)/alpha,M_x_q_x_orig_.*x_p_r_.^(2+n_val_aplx),node_lx_,'spline') * (2*pi);
%%%%%%%%;
fB_nuff_aplx_z_ = finufft1d3(node_lx_,Mxn_x_q_lx_.*weight_lx_,-1,tolerance_nufft,gamma_z_);
%%%%%%%%;
if flag_disp;
fB_quad_aplx_z_ = zeros(n_gamma_z,1);
for ngamma_z=0:n_gamma_z-1;
gamma_z = gamma_z_(1+ngamma_z);
fB = sum(Mxn_x_q_lx_.*exp(-i*gamma_z.*node_lx_).*weight_lx_);
fB_quad_aplx_z_(1+ngamma_z) = fB;
end;%for ngamma_z=0:n_gamma_z-1;
fnorm_disp(flag_verbose,'fB_nuff_aplx_z_',fB_nuff_aplx_z_,'fB_quad_aplx_z_',fB_quad_aplx_z_,'%<-- should be zero');
end;%if flag_disp;
%%%%%%%%;
% Calculate J-integral. ;
%%%%%%%%;
if  isempty(fA_form_aplx_z_);
fA_form_aplx_z_ = (flag_sign*i).^(q_val) / alpha .* S(q_val,-(1+n_val_aplx+i*gamma_z_/alpha));
end;%if  isempty(fA_form_aplx_z_);
%%%%%%%%;
fC_comb_aplx_z_ = fA_form_aplx_z_.*flipud(fB_nuff_aplx_z_);
bC_comb_aplx_y_ = finufft1d3(gamma_z_,fC_comb_aplx_z_.*weight_1d_gamma_z_,+1,tolerance_nufft,log(max(tolerance_master,k_p_r_))/alpha);
M_k_q_k_aplx_ = k_p_r_.^n_val_aplx .* bC_comb_aplx_y_  / (2*pi) ;
%%%%%%%%;

%%%%%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 3; p_col = 2; np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla; hold on;
if ~isempty(fB_nuff_orig_z_);
plot(gamma_z_,real(fB_quad_orig_z_),'ro-',gamma_z_,real(fB_nuff_orig_z_),'rx-'); title('real(fB_nuff_orig_z_)','Interpreter','none');
end;%if ~isempty(fB_nuff_orig_z_);
if ~isempty(fB_nuff_orlx_z_);
plot(gamma_z_,real(fB_quad_orlx_z_),'mo-',gamma_z_,real(fB_nuff_orlx_z_),'mx-'); title('real(fB_nuff_orlx_z_)','Interpreter','none');
end;%if ~isempty(fB_nuff_orlx_z_);
if ~isempty(fB_nuff_appr_z_);
plot(gamma_z_,real(fB_quad_appr_z_),'go-',gamma_z_,real(fB_nuff_appr_z_),'gx-'); title('real(fB_nuff_appr_z_)','Interpreter','none');
end;%if ~isempty(fB_nuff_appr_z_);
if ~isempty(fB_nuff_aplx_z_);
plot(gamma_z_,real(fB_quad_aplx_z_),'co-',gamma_z_,real(fB_nuff_aplx_z_),'cx-'); title('real(fB_nuff_aplx_z_)','Interpreter','none');
end;%if ~isempty(fB_nuff_aplx_z_);
drawnow();
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla; hold on;
if ~isempty(fB_nuff_orig_z_);
plot(gamma_z_,imag(fB_quad_orig_z_),'ro-',gamma_z_,imag(fB_nuff_orig_z_),'rx-'); title('imag(fB_nuff_orig_z_)','Interpreter','none');
end;%if ~isempty(fB_nuff_orig_z_);
if ~isempty(fB_nuff_orlx_z_);
plot(gamma_z_,imag(fB_quad_orlx_z_),'mo-',gamma_z_,imag(fB_nuff_orlx_z_),'mx-'); title('imag(fB_nuff_orlx_z_)','Interpreter','none');
end;%if ~isempty(fB_nuff_orlx_z_);
if ~isempty(fB_nuff_appr_z_);
plot(gamma_z_,imag(fB_quad_appr_z_),'go-',gamma_z_,imag(fB_nuff_appr_z_),'gx-'); title('imag(fB_nuff_appr_z_)','Interpreter','none');
end;%if ~isempty(fB_nuff_appr_z_);
if ~isempty(fB_nuff_aplx_z_);
plot(gamma_z_,imag(fB_quad_aplx_z_),'co-',gamma_z_,imag(fB_nuff_aplx_z_),'cx-'); title('imag(fB_nuff_aplx_z_)','Interpreter','none');
end;%if ~isempty(fB_nuff_aplx_z_);
drawnow();
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla; hold on;
if ~isempty(fC_comb_orig_z_);
plot(gamma_z_,real(fA_form_orig_z_),'rx-',gamma_z_,real(fB_nuff_orig_z_),'rx-',gamma_z_,real(fC_comb_orig_z_),'rs-'); title('real(fC_comb_orig_z_)','Interpreter','none');
end;%if ~isempty(fC_comb_orig_z_);
if ~isempty(fC_comb_orlx_z_);
plot(gamma_z_,real(fA_form_orlx_z_),'mx-',gamma_z_,real(fB_nuff_orlx_z_),'mx-',gamma_z_,real(fC_comb_orlx_z_),'ms-'); title('real(fC_comb_orlx_z_)','Interpreter','none');
end;%if ~isempty(fC_comb_orlx_z_);
if ~isempty(fC_comb_appr_z_);
plot(gamma_z_,real(fA_form_appr_z_),'gx-',gamma_z_,real(fB_nuff_appr_z_),'gx-',gamma_z_,real(fC_comb_appr_z_),'gs-'); title('real(fC_comb_appr_z_)','Interpreter','none');
end;%if ~isempty(fC_comb_appr_z_);
if ~isempty(fC_comb_aplx_z_);
plot(gamma_z_,real(fA_form_aplx_z_),'cx-',gamma_z_,real(fB_nuff_aplx_z_),'cx-',gamma_z_,real(fC_comb_aplx_z_),'cs-'); title('real(fC_comb_aplx_z_)','Interpreter','none');
end;%if ~isempty(fC_comb_aplx_z_);
drawnow();
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla; hold on;
if ~isempty(fC_comb_orig_z_);
plot(gamma_z_,imag(fA_form_orig_z_),'rx-',gamma_z_,imag(fB_nuff_orig_z_),'rx-',gamma_z_,imag(fC_comb_orig_z_),'rs-'); title('imag(fC_comb_orig_z_)','Interpreter','none');
end;%if ~isempty(fC_comb_orig_z_);
if ~isempty(fC_comb_orlx_z_);
plot(gamma_z_,imag(fA_form_orlx_z_),'mx-',gamma_z_,imag(fB_nuff_orlx_z_),'mx-',gamma_z_,imag(fC_comb_orlx_z_),'ms-'); title('imag(fC_comb_orlx_z_)','Interpreter','none');
end;%if ~isempty(fC_comb_orlx_z_);
if ~isempty(fC_comb_appr_z_);
plot(gamma_z_,imag(fA_form_appr_z_),'gx-',gamma_z_,imag(fB_nuff_appr_z_),'gx-',gamma_z_,imag(fC_comb_appr_z_),'gs-'); title('imag(fC_comb_appr_z_)','Interpreter','none');
end;%if ~isempty(fC_comb_appr_z_);
if ~isempty(fC_comb_aplx_z_);
plot(gamma_z_,imag(fA_form_aplx_z_),'cx-',gamma_z_,imag(fB_nuff_aplx_z_),'cx-',gamma_z_,imag(fC_comb_aplx_z_),'cs-'); title('imag(fC_comb_aplx_z_)','Interpreter','none');
end;%if ~isempty(fC_comb_aplx_z_);
drawnow();
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla; hold on;
if ~isempty(M_k_q_k_orig_);
plot(log(k_p_r_),real(M_k_q_k_orig_),'ro-'); title('real(M_k_q_k_orig_)','Interpreter','none');
end;%if ~isempty(M_k_q_k_orig_);
if ~isempty(M_k_q_k_orlx_);
plot(log(k_p_r_),real(M_k_q_k_orlx_),'mo-'); title('real(M_k_q_k_orlx_)','Interpreter','none');
end;%if ~isempty(M_k_q_k_orlx_);
if ~isempty(M_k_q_k_appr_);
plot(log(k_p_r_),real(M_k_q_k_appr_),'go-'); title('real(M_k_q_k_appr_)','Interpreter','none');
end;%if ~isempty(M_k_q_k_appr_);
if ~isempty(M_k_q_k_aplx_);
plot(log(k_p_r_),real(M_k_q_k_aplx_),'co-'); title('real(M_k_q_k_aplx_)','Interpreter','none');
end;%if ~isempty(M_k_q_k_aplx_);
if ~isempty(M_k_q_k_refe_);
plot(log(k_p_r_),real(M_k_q_k_refe_),'ko-'); title('real(M_k_q_k_refe_)','Interpreter','none');
end;%if ~isempty(M_k_q_k_refe_);
drawnow();
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla; hold on;
if ~isempty(M_k_q_k_orig_);
plot(log(k_p_r_),imag(M_k_q_k_orig_),'ro-'); title('imag(M_k_q_k_orig_)','Interpreter','none');
end;%if ~isempty(M_k_q_k_orig_);
if ~isempty(M_k_q_k_orlx_);
plot(log(k_p_r_),imag(M_k_q_k_orlx_),'mo-'); title('imag(M_k_q_k_orlx_)','Interpreter','none');
end;%if ~isempty(M_k_q_k_orlx_);
if ~isempty(M_k_q_k_appr_);
plot(log(k_p_r_),imag(M_k_q_k_appr_),'go-'); title('imag(M_k_q_k_appr_)','Interpreter','none');
end;%if ~isempty(M_k_q_k_appr_);
if ~isempty(M_k_q_k_aplx_);
plot(log(k_p_r_),imag(M_k_q_k_aplx_),'co-'); title('imag(M_k_q_k_aplx_)','Interpreter','none');
end;%if ~isempty(M_k_q_k_aplx_);
if ~isempty(M_k_q_k_refe_);
plot(log(k_p_r_),imag(M_k_q_k_refe_),'ko-'); title('imag(M_k_q_k_refe_)','Interpreter','none');
end;%if ~isempty(M_k_q_k_refe_);
drawnow();
%%%%%%%%;
end;%if flag_disp>0;
%%%%%%%%;

%%%%%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 3; p_col = 2; np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla; hold on;
if ~isempty(fB_nuff_orig_z_);
plot(gamma_z_,log(abs(real(fB_quad_orig_z_))),'ro-',gamma_z_,log(abs(real(fB_nuff_orig_z_))),'rx-'); title('log(abs(real(fB_nuff_orig_z_)))','Interpreter','none');
end;%if ~isempty(fB_nuff_orig_z_);
if ~isempty(fB_nuff_orlx_z_);
plot(gamma_z_,log(abs(real(fB_quad_orlx_z_))),'mo-',gamma_z_,log(abs(real(fB_nuff_orlx_z_))),'mx-'); title('log(abs(real(fB_nuff_orlx_z_)))','Interpreter','none');
end;%if ~isempty(fB_nuff_orlx_z_);
if ~isempty(fB_nuff_appr_z_);
plot(gamma_z_,log(abs(real(fB_quad_appr_z_))),'go-',gamma_z_,log(abs(real(fB_nuff_appr_z_))),'gx-'); title('log(abs(real(fB_nuff_appr_z_)))','Interpreter','none');
end;%if ~isempty(fB_nuff_appr_z_);
if ~isempty(fB_nuff_aplx_z_);
plot(gamma_z_,log(abs(real(fB_quad_aplx_z_))),'co-',gamma_z_,log(abs(real(fB_nuff_aplx_z_))),'cx-'); title('log(abs(real(fB_nuff_aplx_z_)))','Interpreter','none');
end;%if ~isempty(fB_nuff_aplx_z_);
drawnow();
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla; hold on;
if ~isempty(fB_nuff_orig_z_);
plot(gamma_z_,log(abs(imag(fB_quad_orig_z_))),'ro-',gamma_z_,log(abs(imag(fB_nuff_orig_z_))),'rx-'); title('log(abs(imag(fB_nuff_orig_z_)))','Interpreter','none');
end;%if ~isempty(fB_nuff_orig_z_);
if ~isempty(fB_nuff_orlx_z_);
plot(gamma_z_,log(abs(imag(fB_quad_orlx_z_))),'mo-',gamma_z_,log(abs(imag(fB_nuff_orlx_z_))),'mx-'); title('log(abs(imag(fB_nuff_orlx_z_)))','Interpreter','none');
end;%if ~isempty(fB_nuff_orlx_z_);
if ~isempty(fB_nuff_appr_z_);
plot(gamma_z_,log(abs(imag(fB_quad_appr_z_))),'go-',gamma_z_,log(abs(imag(fB_nuff_appr_z_))),'gx-'); title('log(abs(imag(fB_nuff_appr_z_)))','Interpreter','none');
end;%if ~isempty(fB_nuff_appr_z_);
if ~isempty(fB_nuff_aplx_z_);
plot(gamma_z_,log(abs(imag(fB_quad_aplx_z_))),'co-',gamma_z_,log(abs(imag(fB_nuff_aplx_z_))),'cx-'); title('log(abs(imag(fB_nuff_aplx_z_)))','Interpreter','none');
end;%if ~isempty(fB_nuff_aplx_z_);
drawnow();
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla; hold on;
if ~isempty(fC_comb_orig_z_);
plot(gamma_z_,log(abs(real(fA_form_orig_z_))),'rx-',gamma_z_,log(abs(real(fB_nuff_orig_z_))),'rx-',gamma_z_,log(abs(real(fC_comb_orig_z_))),'rs-'); title('log(abs(real(fC_comb_orig_z_)))','Interpreter','none');
end;%if ~isempty(fC_comb_orig_z_);
if ~isempty(fC_comb_orlx_z_);
plot(gamma_z_,log(abs(real(fA_form_orlx_z_))),'mx-',gamma_z_,log(abs(real(fB_nuff_orlx_z_))),'mx-',gamma_z_,log(abs(real(fC_comb_orlx_z_))),'ms-'); title('log(abs(real(fC_comb_orlx_z_)))','Interpreter','none');
end;%if ~isempty(fC_comb_orlx_z_);
if ~isempty(fC_comb_appr_z_);
plot(gamma_z_,log(abs(real(fA_form_appr_z_))),'gx-',gamma_z_,log(abs(real(fB_nuff_appr_z_))),'gx-',gamma_z_,log(abs(real(fC_comb_appr_z_))),'gs-'); title('log(abs(real(fC_comb_appr_z_)))','Interpreter','none');
end;%if ~isempty(fC_comb_appr_z_);
if ~isempty(fC_comb_aplx_z_);
plot(gamma_z_,log(abs(real(fA_form_aplx_z_))),'cx-',gamma_z_,log(abs(real(fB_nuff_aplx_z_))),'cx-',gamma_z_,log(abs(real(fC_comb_aplx_z_))),'cs-'); title('log(abs(real(fC_comb_aplx_z_)))','Interpreter','none');
end;%if ~isempty(fC_comb_aplx_z_);
drawnow();
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla; hold on;
if ~isempty(fC_comb_orig_z_);
plot(gamma_z_,log(abs(imag(fA_form_orig_z_))),'rx-',gamma_z_,log(abs(imag(fB_nuff_orig_z_))),'rx-',gamma_z_,log(abs(imag(fC_comb_orig_z_))),'rs-'); title('log(abs(imag(fC_comb_orig_z_)))','Interpreter','none');
end;%if ~isempty(fC_comb_orig_z_);
if ~isempty(fC_comb_orlx_z_);
plot(gamma_z_,log(abs(imag(fA_form_orlx_z_))),'mx-',gamma_z_,log(abs(imag(fB_nuff_orlx_z_))),'mx-',gamma_z_,log(abs(imag(fC_comb_orlx_z_))),'ms-'); title('log(abs(imag(fC_comb_orlx_z_)))','Interpreter','none');
end;%if ~isempty(fC_comb_orlx_z_);
if ~isempty(fC_comb_appr_z_);
plot(gamma_z_,log(abs(imag(fA_form_appr_z_))),'gx-',gamma_z_,log(abs(imag(fB_nuff_appr_z_))),'gx-',gamma_z_,log(abs(imag(fC_comb_appr_z_))),'gs-'); title('log(abs(imag(fC_comb_appr_z_)))','Interpreter','none');
end;%if ~isempty(fC_comb_appr_z_);
if ~isempty(fC_comb_aplx_z_);
plot(gamma_z_,log(abs(imag(fA_form_aplx_z_))),'cx-',gamma_z_,log(abs(imag(fB_nuff_aplx_z_))),'cx-',gamma_z_,log(abs(imag(fC_comb_aplx_z_))),'cs-'); title('log(abs(imag(fC_comb_aplx_z_)))','Interpreter','none');
end;%if ~isempty(fC_comb_aplx_z_);
drawnow();
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla; hold on;
if ~isempty(M_k_q_k_orig_);
plot(log(k_p_r_),log(abs(real(M_k_q_k_orig_))),'ro-'); title('log(abs(real(M_k_q_k_orig_)))','Interpreter','none');
end;%if ~isempty(M_k_q_k_orig_);
if ~isempty(M_k_q_k_orlx_);
plot(log(k_p_r_),log(abs(real(M_k_q_k_orlx_))),'mo-'); title('log(abs(real(M_k_q_k_orlx_)))','Interpreter','none');
end;%if ~isempty(M_k_q_k_orlx_);
if ~isempty(M_k_q_k_appr_);
plot(log(k_p_r_),log(abs(real(M_k_q_k_appr_))),'go-'); title('log(abs(real(M_k_q_k_appr_)))','Interpreter','none');
end;%if ~isempty(M_k_q_k_appr_);
if ~isempty(M_k_q_k_aplx_);
plot(log(k_p_r_),log(abs(real(M_k_q_k_aplx_))),'co-'); title('log(abs(real(M_k_q_k_aplx_)))','Interpreter','none');
end;%if ~isempty(M_k_q_k_aplx_);
if ~isempty(M_k_q_k_refe_);
plot(log(k_p_r_),log(abs(real(M_k_q_k_refe_))),'ko-'); title('log(abs(real(M_k_q_k_refe_)))','Interpreter','none');
end;%if ~isempty(M_k_q_k_refe_);
drawnow();
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla; hold on;
if ~isempty(M_k_q_k_orig_);
plot(log(k_p_r_),log(abs(imag(M_k_q_k_orig_))),'ro-'); title('log(abs(imag(M_k_q_k_orig_)))','Interpreter','none');
end;%if ~isempty(M_k_q_k_orig_);
if ~isempty(M_k_q_k_orlx_);
plot(log(k_p_r_),log(abs(imag(M_k_q_k_orlx_))),'mo-'); title('log(abs(imag(M_k_q_k_orlx_)))','Interpreter','none');
end;%if ~isempty(M_k_q_k_orlx_);
if ~isempty(M_k_q_k_appr_);
plot(log(k_p_r_),log(abs(imag(M_k_q_k_appr_))),'go-'); title('log(abs(imag(M_k_q_k_appr_)))','Interpreter','none');
end;%if ~isempty(M_k_q_k_appr_);
if ~isempty(M_k_q_k_aplx_);
plot(log(k_p_r_),log(abs(imag(M_k_q_k_aplx_))),'co-'); title('log(abs(imag(M_k_q_k_aplx_)))','Interpreter','none');
end;%if ~isempty(M_k_q_k_aplx_);
if ~isempty(M_k_q_k_refe_);
plot(log(k_p_r_),log(abs(imag(M_k_q_k_refe_))),'ko-'); title('log(abs(imag(M_k_q_k_refe_)))','Interpreter','none');
end;%if ~isempty(M_k_q_k_refe_);
drawnow();
%%%%%%%%;
end;%if flag_disp>0;
%%%%%%%%;

M_k_q_k_outp_ = M_k_q_k_orlx_ ;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

