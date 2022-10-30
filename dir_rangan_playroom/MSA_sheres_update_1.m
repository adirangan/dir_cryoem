function ...
[ ...
 parameter ...
,p_sher_c_ ...
,B_sher_qc__ ...
,B_sher_wc__ ...
] = ...
MSA_sheres_update_1( ...
 parameter ...
,n_w...
,n_q...
,F_wq__...
,F_inv_qw__...
,n_M...
,z_M_...
,p_M_...
,sigma_sher...
,n_class ...
,p_sher_c_ ...
,B_sher_qc__...
,B_sher_wc__...
);

str_thisfunction = 'MSA_sheres_update_1';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_w=[]; end; na=na+1;
if (nargin<1+na); n_q=[]; end; na=na+1;
if (nargin<1+na); F_wq__=[]; end; na=na+1;
if (nargin<1+na); F_inv_qw__=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); z_M_=[]; end; na=na+1;
if (nargin<1+na); p_M_=[]; end; na=na+1;
if (nargin<1+na); sigma_sher=[]; end; na=na+1;
if (nargin<1+na); n_class=[]; end; na=na+1;
if (nargin<1+na); p_sher_c_=[]; end; na=na+1;
if (nargin<1+na); B_sher_qc__=[]; end; na=na+1;
if (nargin<1+na); B_sher_wc__=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp = 0; end;
flag_verbose = parameter.flag_verbose;
flag_disp = parameter.flag_disp;

p_sher_c_ = max(0,p_sher_c_); p_sher_c_ = p_sher_c_/max(1e-12,sum(p_sher_c_));

if flag_verbose; disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
R2_sher_wcM___ = zeros(n_w,n_class,n_M);
p_sher_wcM___ = zeros(n_w,n_class,n_M);
for nclass=0:n_class-1;
B_sher_w_ = B_sher_wc__(:,1+nclass);
R2_sher_wM__ = abs(bsxfun(@minus,reshape(B_sher_w_,[n_w,1]),reshape(z_M_,[1,n_M]))).^2;
R2_sher_wcM___(:,1+nclass,:) = R2_sher_wM__;
clear R2_sher_wM__;
end;%for nclass=0:n_class-1;
R2_sher_wcM___ = bsxfun(@minus,R2_sher_wcM___,min(R2_sher_wcM___,[],[1,2]));
p_sher_wcM___ = exp(-R2_sher_wcM___/max(1e-12,2*sigma_sher^2));
p_sher_wcM___ = bsxfun(@times,p_sher_wcM___,reshape(p_sher_c_,[1,n_class,1]));
p_sher_wcM___ = bsxfun(@rdivide,p_sher_wcM___,max(1e-24,sum(p_sher_wcM___,[1,2]))); %<-- ensure that each image-location has a normalized set of p(w,c|M). ;
p_sher_c_ = reshape(sum(p_sher_wcM___,[1,3]),[n_class,1]);
p_sher_c_ = max(0,p_sher_c_); p_sher_c_ = p_sher_c_/max(1e-12,sum(p_sher_c_));
tmp_denom_wc__ = sum(bsxfun(@times,p_sher_wcM___,reshape(p_M_,[1,1,n_M])),3);
p_sher_wcM___ = bsxfun(@rdivide,p_sher_wcM___,max(1e-12,tmp_denom_wc__)); %<-- ensure that each w,c has a normalized set of p(M|w,c). ;
tmp_index_ = efind(sum(p_sher_wcM___,3)<1e-24); %<-- find and correct w,c which are very far away from all image-locations. ;
p_sher_wcM__ = reshape(p_sher_wcM___,[n_w*n_class,n_M]);
p_sher_wcM__(1+tmp_index_,:) = 1/max(1,n_M);
p_sher_wcM___ = reshape(p_sher_wcM__,[n_w,n_class,n_M]);
tmp_denom_wc__ = sum(bsxfun(@times,p_sher_wcM___,reshape(p_M_,[1,1,n_M])),3);
p_sher_wcM___ = bsxfun(@rdivide,p_sher_wcM___,max(1e-24,tmp_denom_wc__)); %<-- ensure that each w,c has a normalized set of p(M|w,c). ;
for nclass=0:n_class-1;
p_sher_wM__ = squeeze(p_sher_wcM___(:,1+nclass,:));
tmp_denom_w_ = sum(bsxfun(@times,p_sher_wM__,reshape(p_M_,[1,n_M])),2);
p_sher_wM__ = bsxfun(@rdivide,p_sher_wM__,max(1e-12,tmp_denom_w_)); %<-- ensure that each w has a normalized set of p(M|w). ;
tmp_index_ = efind(sum(p_sher_wM__,2)<1e-24); %<-- find and correct w which are very far away from all image-locations. ;
p_sher_wM__(1+tmp_index_,:) = 1/max(1,n_M);
tmp_denom_w_ = sum(bsxfun(@times,p_sher_wM__,reshape(p_M_,[1,n_M])),2);
p_sher_wM__ = bsxfun(@rdivide,p_sher_wM__,max(1e-24,tmp_denom_w_)); %<-- ensure that each w has a normalized set of p(M|w). ;
C_sher_w_ = p_sher_wM__*(z_M_.*p_M_); %<-- now, for each w, accumulate \sum_{nM} p(M|w)*z_M_(1+nM)*p_M_(1+nM). ;
C_sher_q_ = F_inv_qw__*C_sher_w_;
C_sher_w_ = F_wq__*C_sher_q_;
B_sher_q_ = C_sher_q_;
B_sher_w_ = C_sher_w_;
B_sher_qc__(:,1+nclass) = B_sher_q_;
B_sher_wc__(:,1+nclass) = B_sher_w_;
clear tmp_denom_w_ p_sher_wM__ C_sher_w_ C_sher_q_ B_sher_q_ B_sher_w_ ;
end;%for nclass=0:n_class-1;
