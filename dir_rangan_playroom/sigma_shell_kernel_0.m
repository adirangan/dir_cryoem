function ...
[ ...
 parameter ...
,sigma_crit ...
] = ...
sigma_shell_kernel_0( ...
 parameter ...
,l_max ...
,l_upb ...
);

str_thisfunction = 'sigma_shell_kernel_0';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));

disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); l_max=[]; end; na=na+1;
if (nargin<1+na); l_upb=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master=1e-2; end;
tolerance_master=parameter.tolerance_master;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'fminsearch_TolFun'); parameter.fminsearch_TolFun=1e-6; end;
fminsearch_TolFun=parameter.fminsearch_TolFun;
if ~isfield(parameter,'fminsearch_MaxIter'); parameter.fminsearch_MaxIter=1024; end;
fminsearch_MaxIter=parameter.fminsearch_MaxIter;
if ~isfield(parameter,'fminsearch_MaxFunEval'); parameter.fminsearch_MaxFunEval=1024; end;
fminsearch_MaxFunEval=parameter.fminsearch_MaxFunEval;

l_zer = 0;
if isempty(l_max); l_max = 49; end;
if isempty(l_upb); l_upb = min(1024,8*49); end;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

l_all_ = transpose([l_zer:l_upb]);
l_pos_ = transpose([l_max:l_upb]);

f = @(sigma) sum(abs(kernel_coefficient_l2_(l_pos_,sigma)./kernel_coefficient_l2_(l_all_,sigma) - tolerance_master).^2);
sigma_crit = fminsearch(f,0,optimset('TolFun',fminsearch_TolFun,'MaxIter',fminsearch_MaxIter,'MaxFunEval',fminsearch_MaxFunEval));

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

function c = kernel_coefficient(l_val,sigma);
c = exp(-l_val.*(1+l_val).*0.5.*sigma^2) .* sqrt(1+2*l_val) .* sqrt(4*pi) ;

function l2 = kernel_coefficient_l2(l_val_,sigma);
l2 = sum(abs(kernel_coefficient(l_val_,sigma)).^2);

function l2_ = kernel_coefficient_l2_(l_val_,sigma_);
l2_ = size(sigma_);
for ns=0:numel(sigma_)-1;
l2_(1+ns) = kernel_coefficient_l2(l_val_,sigma_(1+ns));
end;%for ns=0:numel(sigma_)-1;

