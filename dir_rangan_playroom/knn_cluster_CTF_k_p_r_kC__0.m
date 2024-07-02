function ...
[ ...
 parameter ...
,index_ncluster_from_nCTF_ ...
] = ...
knn_cluster_CTF_k_p_r_kC__0( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_CTF ...
,CTF_k_p_r_kC__ ...
);

verbose=0;
if (verbose); disp(sprintf('[entering knn_cluster_CTF_k_p_r_kC__0]')); end;

if isempty(parameter); parameter = struct('type','parameter'); end;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'tolerance_cluster')); parameter.tolerance_cluster = parameter.tolerance_master; end; %<-- parameter_bookmark. ;
tolerance_master = parameter.tolerance_master;
tolerance_cluster = parameter.tolerance_cluster;

CTF_k_p_l2_C_ = zeros(n_CTF,1);
CTF_k_p_l2_C_ = reshape(sqrt(sum(bsxfun(@times,abs(CTF_k_p_r_kC__).^2,reshape(weight_2d_k_p_r_,[n_k_p_r,1])),1)),[n_CTF,1]);

index_ncluster_from_nCTF_ = zeros(n_CTF,1);
niteration=0;ncluster=0;
index_unassigned_ = 0:n_CTF-1;
flag_continue = numel(index_unassigned_)>0;
while flag_continue;
n_unassigned = numel(index_unassigned_);
if (verbose>1); disp(sprintf(' %% niteration %d: n_unassigned %d',niteration,n_unassigned)); end;
[tmp_l2_min,tmp_index_min] = min(CTF_k_p_l2_C_(1+index_unassigned_)); tmp_index_min = tmp_index_min-1;
if (verbose>1); disp(sprintf(' %% min at index %d (%d) --> %0.6f',tmp_index_min,index_unassigned_(1+tmp_index_min),tmp_l2_min)); end;
CTF_k_p_r_k_ = CTF_k_p_r_kC__(:,1+index_unassigned_(1+tmp_index_min));
tmp_l2_ = reshape(sqrt(sum(bsxfun(@times,abs(bsxfun(@minus,CTF_k_p_r_kC__(:,1+index_unassigned_),reshape(CTF_k_p_r_k_,[n_k_p_r,1]))).^2,reshape(weight_2d_k_p_r_,[n_k_p_r,1])),1)),[n_unassigned,1]);
tmp_index_nearby_ = efind(tmp_l2_/max(1e-12,tmp_l2_min)<= tolerance_cluster);
n_nearby = numel(tmp_index_nearby_);
if (n_nearby<=0); disp(sprintf(' %% Warning, n_nearby<=0 in knn_cluster_CTF_k_p_r_kC__0')); end;
tmp_index_nearby_ = union(tmp_index_nearby_,tmp_index_min);
if (verbose>1); disp(sprintf(' %% niteration %d: n_nearby = %d',niteration,n_nearby)); end;
index_ncluster_from_nCTF_(1+index_unassigned_(1+tmp_index_nearby_)) = ncluster;
ncluster = ncluster+1;
index_unassigned_ = index_unassigned_(1+setdiff([0:n_unassigned-1],tmp_index_nearby_));
flag_continue = numel(index_unassigned_)> 0;
end;%while flag_continue;

if (verbose); disp(sprintf('[finished knn_cluster_CTF_k_p_r_kC__0]')); end;
