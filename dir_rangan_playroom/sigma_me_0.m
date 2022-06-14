function sigma_out = sigma_me_0(R2_,weight_,sigma_0in,n_iteration,tolerance_sigma);
na=0;
if (nargin<1+na); R2_=[]; end; na=na+1;
if (nargin<1+na); weight_=[]; end; na=na+1;
if (nargin<1+na); sigma_0in=[]; end; na=na+1;
if (nargin<1+na); tolerance_sigma=[]; end; na=na+1;
if (nargin<1+na); n_iteration=[]; end; na=na+1;
R2_ = R2_(:);
if isempty(weight_); weight_ = ones(size(R2_))/numel(R2_); end;
weight_ = weight_(:);
if isempty(sigma_0in); sigma_0in = 1.0; end;
if isempty(n_iteration); n_iteration = 32; end;
if isempty(tolerance_sigma); tolerance_sigma = 1e-6; end;

str_thisfunction = 'sigma_me_0';
verbose=0;
if (verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

tmp_t = tic();
tmp_sigma_pre = sigma_0in;
niteration=0;
flag_continue=1;
while flag_continue;
tmp_p_ = exp(-R2_/(2*max(1e-12,tmp_sigma_pre)^2))/sqrt(2*pi)/tmp_sigma_pre;
tmp_p_sum = dot(tmp_p_,weight_);
tmp_p_ = tmp_p_/tmp_p_sum;
tmp_sigma_pos = sqrt(0.5*sum(tmp_p_.*R2_.*weight_));
tmp_error = fnorm(tmp_sigma_pos-tmp_sigma_pre)/max(1e-12,fnorm(tmp_sigma_pre));
if (verbose>2); disp(sprintf(' %% niteration %d/%d, tmp_sigma_pre %0.6f tmp_error %0.6f',niteration,n_iteration,tmp_sigma_pre,tmp_error)); end;
flag_continue = (tmp_error> tolerance_sigma) & (niteration<n_iteration);
niteration = niteration+1;
tmp_sigma_pre = tmp_sigma_pos;
end;%while flag_continue;
if (verbose>1); disp(sprintf(' %% niteration %d/%d, tmp_sigma_pre %0.6f tmp_error %0.6f',niteration,n_iteration,tmp_sigma_pre,tmp_error)); end;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% sigma_me_0_: %0.3fs',tmp_t)); end;

sigma_out = tmp_sigma_pre;
if (verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
