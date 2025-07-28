function ...
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_p_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
,weight_3d_k_p_r_ ...
);

str_thisfunction = 'get_weight_2d_2';

if nargin<1;
%%%%%%%%;
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
flag_verbose=0; flag_disp=1; nf=0;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi); str_T_vs_L = 'L';
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 0*flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_T_vs_L ...
);
%%%%%%%%;
template_k_eq_d = -1;
n_w_max = 98; n_w_0in_ = n_w_max*ones(n_k_p_r,1);
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_p_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 0*flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
,weight_3d_k_p_r_ ...
);
n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 1; p_col = 3; np=0;
end;%if flag_disp>0;
sigma = k_p_r_max/2;
%%%%;
f = @(k,w) ones(size(k));
F_tru = pi*k_p_r_max.^2;
F_est = dot(f(k_p_r_wk_,k_p_w_wk_),weight_2d_k_p_wk_)*(2*pi)^2;
errrel_0 = abs(F_tru-F_est)/abs(F_tru);
disp(sprintf(' %% F_tru vs F_est: %0.16f',errrel_0));
if flag_disp>0;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(f(k_p_r_wk_,k_p_w_wk_)),[],colormap_beach());
axis image; axisnotick;
title(sprintf('log10(errrel_0) %+0.2f',log10(errrel_0)),'Interpreter','none');
end;%if flag_disp>0;
%%%%;
f = @(k,w) -1.*(1/sigma.^2).*exp(-k.^2/(2*sigma^2));
F_tru = 2*pi*(exp(-k_p_r_max.^2/(2*sigma^2)) - 1);
F_est = dot(f(k_p_r_wk_,k_p_w_wk_),weight_2d_k_p_wk_)*(2*pi)^2;
errrel_1 = abs(F_tru-F_est)/abs(F_tru);
disp(sprintf(' %% F_tru vs F_est: %0.16f',errrel_1));
if flag_disp>0;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(f(k_p_r_wk_,k_p_w_wk_)),[],colormap_beach());
axis image; axisnotick;
title(sprintf('log10(errrel_1) %+0.2f',log10(errrel_1)),'Interpreter','none');
end;%if flag_disp>0;
%%%%;
f = @(k,w) -(sin(w).*cos(w)+1).*(1/sigma.^2).*exp(-k.^2/(2*sigma^2));
F_tru = 2*pi*(exp(-k_p_r_max.^2/(2*sigma^2)) - 1);
F_est = dot(f(k_p_r_wk_,k_p_w_wk_),weight_2d_k_p_wk_)*(2*pi)^2;
errrel_2 = abs(F_tru-F_est)/abs(F_tru);
disp(sprintf(' %% F_tru vs F_est: %0.16f',errrel_2));
if flag_disp>0;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(f(k_p_r_wk_,k_p_w_wk_)),[],colormap_beach());
axis image; axisnotick;
title(sprintf('log10(errrel_2) %+0.2f',log10(errrel_2)),'Interpreter','none');
end;%if flag_disp>0;
%%%%%%%%;
disp('returning'); return;
%%%%%%%%;
end;%if nargin<1;

na=0;
if (nargin<1+na); flag_verbose=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); template_k_eq_d=[]; end; na=na+1;
if (nargin<1+na); n_w_0in_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;

if isempty(flag_verbose); flag_verbose = 0; end;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

flag_calc = 0; 
if isempty(n_k_p_r) | isempty(k_p_r_); flag_calc=1; end;
if flag_calc;
if (flag_verbose>0); disp(sprintf(' %% calculating weight_3d_k_p_r_')); end;
if isempty(k_p_r_max); k_p_r_max = 48.0/(2*pi); end;
k_eq_d = 1.0/(2*pi); str_T_vs_L = 'L';
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 max(1,flag_verbose-1) ...
,k_p_r_max ...
,k_eq_d ...
,str_T_vs_L ...
);
end;%if flag_calc;
if isempty(template_k_eq_d); template_k_eq_d=-1; end;
if isempty(n_w_0in_);
if (flag_verbose>0); disp(sprintf(' %% calculating n_w_0in_')); end;
l_max_upb = round(2*pi*k_p_r_max); %<-- typically sufficient for 2-3 digits of precision. ;
n_w_max = 2*(l_max_upb+1);
n_w_0in_ = n_w_max*ones(n_k_p_r,1);
end;%if isempty(n_w_0in_);
if  isempty(weight_3d_k_p_r_); flag_pinv=1; end;
if ~isempty(weight_3d_k_p_r_); flag_pinv=0; end;

n_w_ = zeros(n_k_p_r,1);
if (template_k_eq_d>0);
if (flag_verbose>0); disp(sprintf(' %% template_k_eq_d> 0')); end;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_equator = 3+round(2*pi*k_p_r/template_k_eq_d);
n_polar_a = 3+round(n_equator/2);
n_w_(1+nk_p_r) = 2*n_polar_a;
end;%for nk_p_r=0:n_k_p_r-1;
end;%if (template_k_eq_d>0);
if (template_k_eq_d<=0);
if (flag_verbose>0); disp(sprintf(' %% template_k_eq_d<=0')); end;
n_w_ = n_w_0in_;
assert(numel(n_w_)==n_k_p_r); assert(min(n_w_)>0);
end;%if (template_k_eq_d<=0);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
if (flag_verbose>0); disp(sprintf(' %% n_w_max %d n_w_sum %d',n_w_max,n_w_sum)); end;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% Set up integration weights for the templates.')); end;
%%%%%%%%;
if flag_pinv==1;
if (flag_verbose>0); disp(sprintf(' %% calling pinv')); end;
tmp_P_ = zeros(n_k_p_r,n_k_p_r); %<-- polynomials of order 0:n_k_p_r-1 evaluated on k_p_r_/k_p_r_max. ;
tmp_I_ = zeros(n_k_p_r,1); %<-- integrals of those polynomials on the 2d-disc of radius 1. ;
for nk_p_r=0:n_k_p_r-1;
tmp_x = @(x) x.^nk_p_r;
tmp_P_(1+nk_p_r,:) = tmp_x(k_p_r_/k_p_r_max);
tmp_I_(1+nk_p_r) = 2*pi*1/(nk_p_r+2);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_W_ = pinv(tmp_P_,1e-6)*tmp_I_;
if (flag_verbose>1); disp(sprintf(' %% weight error: %0.16f',fnorm(tmp_P_*tmp_W_ - tmp_I_)/fnorm(tmp_I_))); end;
weight_2d_k_p_r_ = tmp_W_*k_p_r_max^2;
end;%if flag_pinv==1;
if flag_pinv==0;
if (flag_verbose>0); disp(sprintf(' %% rescaling weight_3d_k_p_r_')); end;
weight_2d_k_p_r_ = 2*pi*reshape(weight_3d_k_p_r_,[n_k_p_r,1])./k_p_r_;
if flag_verbose>0;
tmp_P_ = zeros(n_k_p_r,n_k_p_r); %<-- polynomials of order 0:n_k_p_r-1 evaluated on k_p_r_/k_p_r_max. ;
tmp_I_ = zeros(n_k_p_r,1); %<-- integrals of those polynomials on the 2d-disc of radius 1. ;
for nk_p_r=0:n_k_p_r-1;
tmp_x = @(x) x.^nk_p_r;
tmp_P_(1+nk_p_r,:) = tmp_x(k_p_r_/k_p_r_max);
tmp_I_(1+nk_p_r) = 2*pi*1/(nk_p_r+2);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_PW_ = tmp_P_*weight_2d_k_p_r_/k_p_r_max^2;
disp(sprintf(' %% tmp_I_ vs tmp_PW_: %0.16f',fnorm(tmp_I_-tmp_PW_)/max(1e-12,fnorm(tmp_I_))));
end;%if flag_verbose>0;
end;%if flag_pinv==0;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% constructing weight_2d_k_p_wk_')); end;
weight_2d_k_p_wk_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
weight_2d_k_p_wk_(1+tmp_index_) = weight_2d_k_p_r_(1+nk_p_r) / max(1,n_w_(1+nk_p_r)) / (2*pi)^2;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% Set up integration nodes for the templates.')); end;
%%%%%%%%;
k_c_0_wk_ = zeros(n_w_sum,1);
k_c_1_wk_ = zeros(n_w_sum,1);
k_p_r_wk_ = zeros(n_w_sum,1);
k_p_w_wk_ = zeros(n_w_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
k_p_r = k_p_r_(1+nk_p_r);
for nw=0:n_w-1;
gamma_z = (2*pi)*(1.0*nw)/max(1,n_w); cc = cos(gamma_z); sc = sin(gamma_z);
k_c_0_wk_(1+na) = k_p_r*cc;
k_c_1_wk_(1+na) = k_p_r*sc;
k_p_r_wk_(1+na) = k_p_r;
k_p_w_wk_(1+na) = gamma_z;
na=na+1;
end;%for nw=0:n_w-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
