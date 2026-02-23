function ...
[ ...
 parameter ...
,M_k_q_k_ ...
,n_gamma_z ...
,gamma_z_ ...
,gamma_z_max ...
,weight_1d_gamma_z_ ...
,fA_form_z_ ...
,fA_for2_z_ ...
] = ...
fht_0( ...
 parameter ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,weight_2d_x_p_r_ ...
,M_x_q_x_ ...
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
,fA_form_z_ ...
,fA_for2_z_ ...
,M_k_q_k_ref_ ...
);

str_thisfunction = 'fht_0';

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
if (nargin<1+na); M_x_q_x_=[]; end; na=na+1;
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
if (nargin<1+na); fA_form_z_=[]; end; na=na+1;
if (nargin<1+na); fA_for2_z_=[]; end; na=na+1;
if (nargin<1+na); M_k_q_k_ref_=[]; end; na=na+1;

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

%%%%%%%%;
% Bessel. ;
%%%%%%%%;
if abs(q_val)==0; n_val = flag_sign*0.75; end;
if abs(q_val)> 0; n_val = flag_sign*1.00; end;
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
fB_nuff_z_ = finufft1d3(log(max(tolerance_master,x_p_r_))/alpha,M_x_q_x_.*x_p_r_.^(n_val).*weight_2d_x_p_r_,-1,tolerance_nufft,gamma_z_);
%%%%%%%%;
if flag_disp;
fB_quad_z_ = zeros(n_gamma_z,1);
for ngamma_z=0:n_gamma_z-1;
gamma_z = gamma_z_(1+ngamma_z);
fB = sum(M_x_q_x_.*x_p_r_.^(n_val).*x_p_r_.^(-i*gamma_z/alpha).*weight_2d_x_p_r_);
fB_quad_z_(1+ngamma_z) = fB;
end;%for ngamma_z=0:n_gamma_z-1;
fnorm_disp(flag_verbose,'fB_nuff_z_',fB_nuff_z_,'fB_quad_z_',fB_quad_z_,'%<-- should be zero');
end;%if flag_disp;
%%%%%%%%;
% Calculate J-integral. ;
%%%%%%%%;
if  isempty(fA_form_z_);
fA_form_z_ = (flag_sign*i).^(q_val) / alpha .* Q(q_val,-(1+n_val+i*gamma_z_/alpha));
end;%if  isempty(fA_form_z_);
%%%%%%%%;
fC_comb_z_ = fA_form_z_.*flipud(fB_nuff_z_);
bC_comb_y_ = finufft1d3(gamma_z_,fC_comb_z_.*weight_1d_gamma_z_,+1,tolerance_nufft,log(max(tolerance_master,k_p_r_))/alpha);
M_k_q_k_ = k_p_r_.^n_val .* bC_comb_y_  / (2*pi) ;
%%%%%%%%;

%%%%%%%%;
% Approximation. ;
% J_approx = @(x_p_r_) (-1).^(q_val .* (q_val<0)) .* (x_p_r_/2).^(abs(q_val)) .* (1/gamma(1+abs(q_val)) - x_p_r_.^2/(4*gamma(2+abs(q_val))) + x_p_r_.^4/(16*2*gamma(3+abs(q_val)))) .* exp(-b.*x_p_r_.^6) ;
%%%%%%%%;
tmp_b = 1.0;
if abs(q_val)==0; n_va2 = flag_sign*(0.75); end;
if abs(q_val)> 0; n_va2 = flag_sign*(1.00); end;
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
fB_nuf2_z_ = finufft1d3(log(max(tolerance_master,x_p_r_))/alpha,M_x_q_x_.*x_p_r_.^(n_va2).*weight_2d_x_p_r_,-1,tolerance_nufft,gamma_z_);
%%%%%%%%;
if flag_disp;
fB_qua2_z_ = zeros(n_gamma_z,1);
for ngamma_z=0:n_gamma_z-1;
gamma_z = gamma_z_(1+ngamma_z);
fB = sum(M_x_q_x_.*x_p_r_.^(n_va2).*x_p_r_.^(-i*gamma_z/alpha).*weight_2d_x_p_r_);
fB_qua2_z_(1+ngamma_z) = fB;
end;%for ngamma_z=0:n_gamma_z-1;
fnorm_disp(flag_verbose,'fB_nuf2_z_',fB_nuf2_z_,'fB_qua2_z_',fB_qua2_z_,'%<-- should be zero');
end;%if flag_disp;
%%%%%%%%;
% Calculate approximate J-integral. ;
%%%%%%%%;
if  isempty(fA_for2_z_);
fA_for2_z_ = (flag_sign*i).^(q_val) / alpha .* S(q_val,-(1+n_va2+i*gamma_z_/alpha));
end;%if  isempty(fA_for2_z_);
%%%%%%%%;
fC_com2_z_ = fA_for2_z_.*flipud(fB_nuf2_z_);
bC_com2_y_ = finufft1d3(gamma_z_,fC_com2_z_.*weight_1d_gamma_z_,+1,tolerance_nufft,log(max(tolerance_master,k_p_r_))/alpha);
N_k_q_k_ = k_p_r_.^n_va2 .* bC_com2_y_  / (2*pi) ;
%%%%%%%%;

%%%%%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 6; p_col = 2; np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(gamma_z_,real(fB_quad_z_),'ko-',gamma_z_,real(fB_nuff_z_),'rx-'); title('real(fB_nuff_z_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(gamma_z_,imag(fB_quad_z_),'ko-',gamma_z_,imag(fB_nuff_z_),'rx-'); title('imag(fB_nuff_z_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(gamma_z_,real(fA_form_z_),'gx-',gamma_z_,real(fB_nuff_z_),'rx-',gamma_z_,real(fC_comb_z_),'bs-'); title('real(fC_comb_z_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(gamma_z_,imag(fA_form_z_),'gx-',gamma_z_,imag(fB_nuff_z_),'rx-',gamma_z_,imag(fC_comb_z_),'bs-'); title('imag(fC_comb_z_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
hold on;
plot(log(k_p_r_),real(M_k_q_k_),'o-'); title('real(M_k_q_k_)','Interpreter','none');
if ~isempty(M_k_q_k_ref_); plot(log(k_p_r_),real(M_k_q_k_ref_),'gx-'); end;
subplot(p_row,p_col,1+np);np=np+1;cla;
hold on;
plot(log(k_p_r_),imag(M_k_q_k_),'o-'); title('imag(M_k_q_k_)','Interpreter','none');
if ~isempty(M_k_q_k_ref_); plot(log(k_p_r_),imag(M_k_q_k_ref_),'gx-'); end;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(gamma_z_,real(fB_qua2_z_),'ko-',gamma_z_,real(fB_nuf2_z_),'rx-'); title('real(fB_nuf2_z_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(gamma_z_,imag(fB_qua2_z_),'ko-',gamma_z_,imag(fB_nuf2_z_),'rx-'); title('imag(fB_nuf2_z_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(gamma_z_,real(fA_for2_z_),'gx-',gamma_z_,real(fB_nuf2_z_),'rx-',gamma_z_,real(fC_com2_z_),'bs-'); title('real(fC_com2_z_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(gamma_z_,imag(fA_for2_z_),'gx-',gamma_z_,imag(fB_nuf2_z_),'rx-',gamma_z_,imag(fC_com2_z_),'bs-'); title('imag(fC_com2_z_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
hold on;
plot(log(k_p_r_),real(N_k_q_k_),'o-'); title('real(N_k_q_k_)','Interpreter','none');
if ~isempty(M_k_q_k_ref_); plot(log(k_p_r_),real(M_k_q_k_ref_),'gx-'); end;
subplot(p_row,p_col,1+np);np=np+1;cla;
hold on;
plot(log(k_p_r_),imag(N_k_q_k_),'o-'); title('imag(N_k_q_k_)','Interpreter','none');
if ~isempty(M_k_q_k_ref_); plot(log(k_p_r_),imag(M_k_q_k_ref_),'gx-'); end;
end;%if flag_disp>0;
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

