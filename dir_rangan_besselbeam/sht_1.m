function ...
[ ...
 parameter ...
,M_k_q_k_outp_ ...
] = ...
sht_1( ...
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
);

str_thisfunction = 'sht_1';

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

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master=1e-12; end;
tolerance_master = parameter.tolerance_master;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp=0; end;
flag_disp=parameter.flag_disp; nf=0;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(q_val); q_val = 0; end;
if isempty(flag_sign); flag_sign = +1; end;

%%%%%%%%;
% Bessel. ;
%%%%%%%%;
M_k_q_k_slow_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
M_k_q_k_slow_(1+nk_p_r) = sum( (flag_sign.*i).^q_val .* besselj(q_val,2*pi*k_p_r*x_p_r_) .* M_x_q_x_orig_ .* weight_2d_x_p_r_ ); %<-- note that weight_2d_x_p_r_ accounts for radial weighting. ;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;

%%%%%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 1; p_col = 2; np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla; hold on;
if ~isempty(M_k_q_k_slow_);
plot(log(k_p_r_),real(M_k_q_k_slow_),'ro-'); title('real(M_k_q_k_slow_)','Interpreter','none');
end;%if ~isempty(M_k_q_k_slow_);
if ~isempty(M_k_q_k_refe_);
plot(log(k_p_r_),real(M_k_q_k_refe_),'ko-'); title('real(M_k_q_k_refe_)','Interpreter','none');
end;%if ~isempty(M_k_q_k_refe_);
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla; hold on;
if ~isempty(M_k_q_k_slow_);
plot(log(k_p_r_),imag(M_k_q_k_slow_),'ro-'); title('imag(M_k_q_k_slow_)','Interpreter','none');
end;%if ~isempty(M_k_q_k_slow_);
if ~isempty(M_k_q_k_refe_);
plot(log(k_p_r_),imag(M_k_q_k_refe_),'ko-'); title('imag(M_k_q_k_refe_)','Interpreter','none');
end;%if ~isempty(M_k_q_k_refe_);
drawnow();
%%%%%%%%;
end;%if flag_disp>0;
%%%%%%%%;

M_k_q_k_outp_ = M_k_q_k_slow_ ;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

