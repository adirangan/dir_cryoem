function ...
[ ...
 parameter ...
,ssnll_Ms__ ...
,D2_SM__ ...
] = ...
ssnll_from_X_0( ...
 parameter ...
,n_S ...
,n_M ...
,X_SM__ ...
,S_l2_S_ ...
,M_l2_M_ ...
,viewing_weight_all_ ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
);

str_thisfunction = 'ssnll_from_X_0';

%%%%%%%%;
% calculates ssnll (sigma*sigma*log_unlikelihood). ;
%%%%%%%%;
% Inputs: ;
% n_S: integer number of templates. ;
% n_M: integer number of images. ;
% X_SM__: double array of size (n_S,nM). ;
%     X_SM__(1+nS,1+nM) = correlation between template nS and image nM. ;
% S_l2_S_: double array of size (n_S,1). ;
%    S_l2_S_(1+nS) = |S_k_p_wkS__(:,1+nS)|^{2}. ;
% M_l2_M_: double array of size (n_M,1). ;
%    M_l2_M_(1+nM) = |M_k_p_wkM__(:,1+nM)|^{2}. ;
% viewing_weight_all_: double array of size(n_S,1). ;
%    viewing_weight_all_(1+nS) = quadrature-weight for template nS. ;
% n_sigma_bayesian: integer number of sigma values. ;
% sigma_bayesian_: double array of size (n_sigma_bayesian,1). ;
%    sigma_bayesian := sigma_bayesian_(1+nsigma_bayesian) = temperature parameter for bayesian-integration. ;
%    Note that sigma_bayesian==0 corresponds to maximum-likelihood alignment. ;
%%%%%%%%;
% Outputs: ;
% ssnll_Ms__: double array of size (n_M,n_sigma_bayesian). ;
%     ssnll_Ms__(1+nM,1+nsigma_bayesian) = ssnll for image nM and sigma_bayesian. ;
% D2_SM__: double array of size (n_S,n_M). ;
%     D2_SM__(1+nS,1+nM) = distance-squared between template nS and image nM. ;
%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); X_SM__=[]; end; na=na+1;
if (nargin<1+na); S_l2_S_=[]; end; na=na+1;
if (nargin<1+na); M_l2_M_=[]; end; na=na+1;
if (nargin<1+na); viewing_weight_all_=[]; end; na=na+1;
if (nargin<1+na); n_sigma_bayesian=[]; end; na=na+1;
if (nargin<1+na); sigma_bayesian_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;

if isempty(viewing_weight_all_); viewing_weight_all_ = ones(n_S,1); end;
if isempty(n_sigma_bayesian); n_sigma_bayesian = 1; end;
if isempty(sigma_bayesian_); sigma_bayesian = zeros(n_sigma_bayesian,1); end;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

ssnll_Ms__ = zeros(n_M,n_sigma_bayesian);
D2_SM__ = bsxfun(@plus,reshape(S_l2_S_,[n_S,1]),reshape(M_l2_M_,[1,n_M])) ...
  - 2*X_SM__ .* (reshape(sqrt(S_l2_S_),[n_S,1])*reshape(sqrt(M_l2_M_),[1,n_M]));
for nsigma_bayesian=0:n_sigma_bayesian-1;
sigma_bayesian = sigma_bayesian_(1+nsigma_bayesian);
if (sigma_bayesian< 1e-12);
tmp_ssnll_M_ = 0.5*min(D2_SM__,[],1); %<-- sigma^2 * negative-log-likelihood for 1-molecule model (ignoring normalization). ;
end;%if (sigma_bayesian< 1e-12);
if (sigma_bayesian>=1e-12);
tmp_D2_min_M_ = min(D2_SM__,[],1);
tmp_l_SM__ = exp(-bsxfun(@minus,D2_SM__,reshape(tmp_D2_min_M_,[1,n_M]))/max(1e-12,2*sigma_bayesian^2)); %<-- at least one zero entry. ;
tmp_l_M_ = sum(bsxfun(@times,tmp_l_SM__,reshape(viewing_weight_all_,[n_S,1])),1)/max(1e-12,sum(viewing_weight_all_));
tmp_ssnll_M_ = sigma_bayesian^2 * log(max(1e-12,tmp_l_M_)) + tmp_D2_min_M_;
end;%if (sigma_bayesian>=1e-12);
ssnll_Ms__(:,1+nsigma_bayesian) = tmp_ssnll_M_;
clear tmp_D2_min_M_ tmp_l_SM__ tmp_l_M_ tmp_ssnll_M_ ;
end;%for nsigma_bayesian=0:n_sigma_bayesian-1;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
