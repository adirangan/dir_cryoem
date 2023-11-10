function ...
[ ...
 parameter ...
,ssnll_s_ ...
,ssnll_Ms__ ...
,inten_vs__ ...
] = ...
ssnll_from_X_1( ...
 parameter ...
,n_volume ... 
,n_S_v_ ...
,n_M ...
,X_vSM___ ...
,S_l2_vS__ ...
,M_l2_M_ ...
,viewing_weight_vS__ ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
,inten_0in_vs__ ... 
);

str_thisfunction = 'ssnll_from_X_1';

%%%%%%%%;
% Calculates ssnll (sigma*sigma*log_unlikelihood). ;
% Allows for multiple volumes, assuming hard-assignment of images to volumes. ;
% Optimizes over volume-intensity parameter (inten) below. ;
% Assumes a single inten for each volume. ;
%%%%%%%%;
% Inputs: ;
% n_volume: integer number of volumes. ;
% n_S_v_: integer array of size (n_volume,1). ;
%    n_S := n_S_v_(1+nvolume) = number of templates for volume nV. ;
% n_M: integer number of images. ;
% X_vSM___: cell array of size (n_volume,1). ;
%    X_SM__ := X_vSM___{1+nvolume} = double array of size (n_S,nM). ;
%    X_SM__(1+nS,1+nM) = correlation between template nS and image nM. ;
% S_l2_vS__: cell array of size (n_volume,1). ;
%    S_l2_S_ := S_l2_vS__{1+nvolume} = double array of size (n_S,1). ;
%    S_l2_S_(1+nS) = |S_k_p_wkS__(:,1+nS)|^{2}. ;
% M_l2_M_: double array of size (n_M,1). ;
%    M_l2_M_(1+nM) = |M_k_p_wkM__(:,1+nM)|^{2}. ;
% viewing_weight_vS__: cell array of size (n_volume,1). ;
%    viewing_weight_all_ := viewing_weight_vS__{1+nvolume}: double array of size(n_S,1). ;
%    viewing_weight_all_(1+nS) = quadrature-weight for template nS. ;
% n_sigma_bayesian: integer number of sigma values. ;
% sigma_bayesian_: double array of size (n_sigma_bayesian,1). ;
%    sigma_bayesian_(1+nsigma_bayesian) = temperature parameter for bayesian-integration. ;
%    Note that sigma_bayesian==0 corresponds to maximum-likelihood alignment. ;
% inten_0in_vs__: double array of size (n_volume,n_sigma_bayesian). ;
%     inten_0in_vs__(1+nvolume,1+nsigma_bayesian) = input inten for volume nvolume at sigma_bayesian. ;
%%%%%%%%;
% Outputs: ;
% ssnll_s_: double array of size (n_sigma_bayesian,1). ;
%     ssnll_s_(1+nsigma_bayesian) = ssnll for sigma_bayesian. ;
% ssnll_Ms__: double array of size (n_M,n_sigma_bayesian). ;
%     ssnll_Ms__(1+nM,1+nsigma_bayesian) = ssnll for image nM and sigma_bayesian. ;
% inten_vs__: double array of size (n_volume,n_sigma_bayesian). ;
%     inten_vs__(1+nvolume,1+nsigma_bayesian) = optimal inten for volume nvolume at sigma_bayesian. ;
%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_volume=[]; end; na=na+1;
if (nargin<1+na); n_S_v_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); X_vSM___=[]; end; na=na+1;
if (nargin<1+na); S_l2_vS__=[]; end; na=na+1;
if (nargin<1+na); M_l2_M_=[]; end; na=na+1;
if (nargin<1+na); viewing_weight__=[]; end; na=na+1;
if (nargin<1+na); n_sigma_bayesian=[]; end; na=na+1;
if (nargin<1+na); sigma_bayesian_=[]; end; na=na+1;
if (nargin<1+na); inten_0in_vs__=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'fminsearch_MaxFunEvals'); parameter.fminsearch_MaxFunEvals=512; end;
flag_verbose=parameter.fminsearch_MaxFunEvals;
if ~isfield(parameter,'fminsearch_MaxIter'); parameter.fminsearch_MaxIter=512; end;
flag_verbose=parameter.fminsearch_MaxIter;
if ~isfield(parameter,'fminsearch_Display'); parameter.fminsearch_Display='off'; end;
flag_verbose=parameter.fminsearch_Display;

if isempty(viewing_weight_vS__);
viewing_weight_vS__ = cell(n_volume,1);
for nvolume=0:n_volume-1;
viewing_weight_vS__{1+nvolume} = ones(n_S_v_(1+nvolume),1);
end;%for nvolume=0:n_volume-1;
end;%if isempty(viewing_weight_vS__);
if isempty(n_sigma_bayesian); n_sigma_bayesian = 1; end;
if isempty(sigma_bayesian_); sigma_bayesian = zeros(n_sigma_bayesian,1); end;
if isempty(inten_0in_vs__); inten_0in_vs__ = zeros(n_volume,n_sigma_bayesian); end;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

ssnll_s_ = zeros(n_sigma_bayesian,1);
ssnll_Ms__ = zeros(n_M,n_sigma_bayesian,1);
inten_vs__ = zeros(n_volume,n_sigma_bayesian);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nsigma_bayesian=0:n_sigma_bayesian-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
sigma_bayesian = sigma_bayesian_(1+nsigma_bayesian);
inten_0in_v_ = inten_0in_vs__(:,1+nsigma_bayesian);
f_ssnll_M_ = @(linten_v_) ...
 ssnll_from_X_excerpt_1( ...
 n_volume ...
,linten_v_ ...
,n_S_v_ ...
,n_M ...
,X_vSM___ ...
,S_l2_vS__ ...
,M_l2_M_ ...
,viewing_weight_vS__ ...
,1 ...
,sigma_bayesian ...
);
f_ssnll = @(linten_v_) sum(f_ssnll_M_(linten_v_),1);
tmp_options = optimset('Display',parameter.fminsearch_Display,'MaxFunEvals',parameter.fminsearch_MaxFunEvals,'MaxIter',parameter.fminsearch_MaxIter);
[linten_v_,ssnll] = fminsearch(f_ssnll,inten_0in_v_,tmp_options);
inten_v_ = exp(linten_v_);
inten_vs__(:,1+nsigma_bayesian) = inten_v_;
ssnll_s_(1+nsigma_bayesian) = ssnll;
ssnll_Ms__(:,1+nsigma_bayesian) = f_ssnll_M_(linten_v_);
clear sigma_bayesian f_ssnll_M_ f_ssnll tmp_options linten_v_ inten_v_ ssnll ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nsigma_bayesian=0:n_sigma_bayesian-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

function ...
[ ...
ssnll_Ms__ ...
] = ...
 ssnll_from_X_excerpt_1( ...
 n_volume ...
,linten_v_ ...
,n_S_v_ ...
,n_M ...
,X_vSM___ ...
,S_l2_vS__ ...
,M_l2_M_ ...
,viewing_weight_vS__ ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
);
ssnll_Msv___ = zeros(n_M,n_sigma_bayesian,n_volume);
for nvolume=0:n_volume-1;
ssnll_Msv___(:,:,1+nvolume) = ...
ssnll_from_X_excerpt_0( ...
 n_S_v_(1+nvolume) ...
,n_M ...
,X_vSM___{1+nvolume} ...
,S_l2_vS__{1+nvolume} .* exp(linten_v_(1+nvolume)).^2 ...
,M_l2_M_ ...
,viewing_weight_vS__{1+nvolume} ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
);
end;%for nvolume=0:n_volume-1;
ssnll_Ms__ = squeeze(min(ssnll_Msv___,[],3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

function ...
[ ...
ssnll_Ms__ ...
] = ...
 ssnll_from_X_excerpt_0( ...
 n_S ...
,n_M ...
,X_SM__ ...
,S_l2_S_ ...
,M_l2_M_ ...
,viewing_weight_all_ ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
);

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
