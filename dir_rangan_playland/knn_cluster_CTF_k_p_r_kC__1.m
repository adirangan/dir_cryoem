function ...
[ ...
 parameter ...
,index_ncluster_from_nCTF_ ...
] = ...
knn_cluster_CTF_k_p_r_kC__1( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_CTF ...
,CTF_k_p_r_kC__ ...
);
%%%%%%%%;
% Groups isotropic CTF_k_p_r_kC__ into clusters. ;
% The clustering proceeds iteratively, ;
% at each step choosing the CTF_k_p_r_k_ with the lowest norm, ;
% and then assigning to it the other CTFs with relative l2-difference less than tolerance_cluster. ;
%%%%%%%%;

str_thisfunction = 'knn_cluster_CTF_k_p_r_kC__1';

%%%%%%%%;
if nargin<1;
%%%%%%%%;
flag_verbose=1;
if (flag_verbose>0); disp(sprintf(' %% testing %s',str_thisfunction)); end;
parameter = struct('type','parameter');
parameter.flag_verbose = flag_verbose;
parameter.tolerance_cluster = 1e-1;
n_k_p_r = 49;
k_p_r_ = transpose(linspace(1,7,n_k_p_r));
weight_2d_k_p_r_ = cos(k_p_r_).^2;
n_CTF = 128;
CTF_k_p_r_kC__ = zeros(n_k_p_r,n_CTF);
CTF_phi_C_ = transpose(linspace(0,pi,n_CTF));
for nCTF=0:n_CTF-1;
CTF_k_p_r_kC__(:,1+nCTF) = sin(k_p_r_ - CTF_phi_C_(1+nCTF));
end;%for nCTF=0:n_CTF-1;
[ ...
 parameter ...
,index_ncluster_from_nCTF_ ...
] = ...
knn_cluster_CTF_k_p_r_kC__1( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_CTF ...
,CTF_k_p_r_kC__ ...
);
%%%%;
% compare to ../dir_rangan_python/dir_pymat/test_knn_cluster_CTF_k_p_r_kC__1.mat. ;
%%%%;
dir_pymat = '../dir_rangan_python/dir_pymat';
fname_pymat = sprintf('%s/test_knn_cluster_CTF_k_p_r_kC__1.mat',dir_pymat);
if ~exist(fname_pymat,'file'); disp(sprintf(' %% Warning, %s not found',fname_pymat)); end;
if  exist(fname_pymat,'file');
disp(sprintf(' %% %s found, loading',fname_pymat));
tmp_ = load(fname_pymat);
fnorm_disp(flag_verbose,'index_ncluster_from_nCTF_',index_ncluster_from_nCTF_,'double(tmp_.index_ncluster_from_nCTF_)',double(tmp_.index_ncluster_from_nCTF_),' %%<-- should be <1e-6');
figure(1);clf;figsml;plot(index_ncluster_from_nCTF_,'o'); xlabel('nCTF'); ylabel('ncluster');
end;%if  exist(fname_pymat,'file');
%%%%%%%%;
disp('returning'); return;
end;%if nargin<1;
%%%%%%%%;

if isempty(parameter); parameter = struct('type','parameter'); end;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
flag_verbose = parameter.flag_verbose;
if (~isfield(parameter,'flag_disp')); parameter.flag_disp = 0; end; %<-- parameter_bookmark. ;
flag_disp = parameter.flag_disp; nf=0;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
tolerance_master = parameter.tolerance_master;
if (~isfield(parameter,'tolerance_cluster')); parameter.tolerance_cluster = parameter.tolerance_master; end; %<-- parameter_bookmark. ;
tolerance_cluster = parameter.tolerance_cluster;

if (flag_verbose>0); disp(sprintf('[entering %s]',str_thisfunction)); end;

CTF_k_p_l2_C_ = zeros(n_CTF,1);
CTF_k_p_l2_C_ = reshape(sqrt(sum(bsxfun(@times,abs(CTF_k_p_r_kC__).^2,reshape(weight_2d_k_p_r_,[n_k_p_r,1])),1+0)),[n_CTF,1]);

index_ncluster_from_nCTF_ = zeros(n_CTF,1);
niteration=0;ncluster=0;
index_unassigned_ = transpose(0:n_CTF-1);
flag_continue = numel(index_unassigned_)>0;
while flag_continue;
n_unassigned = numel(index_unassigned_);
if (flag_verbose>1); disp(sprintf(' %% niteration %d: n_unassigned %d',niteration,n_unassigned)); end;
[tmp_l2_min,tmp_index_min] = min(CTF_k_p_l2_C_(1+index_unassigned_)); tmp_index_min = tmp_index_min-1;
if (flag_verbose>1); disp(sprintf(' %% min at index %d (%d) --> %0.6f',tmp_index_min,index_unassigned_(1+tmp_index_min),tmp_l2_min)); end;
CTF_cen_k_p_r_k_ = CTF_k_p_r_kC__(:,1+index_unassigned_(1+tmp_index_min));
l2_cen = CTF_k_p_l2_C_(1+index_unassigned_(1+tmp_index_min));
assert(abs(l2_cen-tmp_l2_min)<1e-12);
tmp_l2_ = reshape(sqrt(sum(bsxfun(@times,abs(bsxfun(@minus,CTF_k_p_r_kC__(:,1+index_unassigned_),reshape(CTF_cen_k_p_r_k_,[n_k_p_r,1]))).^2,reshape(weight_2d_k_p_r_,[n_k_p_r,1])),1)),[n_unassigned,1]);
tmp_index_nearby_ = efind(tmp_l2_/max(1e-12,tmp_l2_min)<= tolerance_cluster); %<-- relative difference. ;
n_nearby = numel(tmp_index_nearby_); %<-- should be at least 1, including CTF_cen_k_p_r_k_ itself. ;
if (n_nearby<=0); disp(sprintf(' %% Warning, n_nearby<=0 in knn_cluster_CTF_k_p_r_kC__0')); end;
tmp_index_nearby_ = union(tmp_index_nearby_,tmp_index_min);
if (flag_verbose>1); disp(sprintf(' %% niteration %d: n_nearby = %d',niteration,n_nearby)); end;
index_ncluster_from_nCTF_(1+index_unassigned_(1+tmp_index_nearby_)) = ncluster;
ncluster = ncluster+1;
index_unassigned_ = index_unassigned_(1+setdiff([0:n_unassigned-1],tmp_index_nearby_));
flag_continue = numel(index_unassigned_)> 0;
end;%while flag_continue;

if flag_disp;
CTF_k_p_l2_CC__ = reshape(sqrt(sum((reshape(CTF_k_p_r_kC__,[n_k_p_r,n_CTF,1]) - reshape(CTF_k_p_r_kC__,[n_k_p_r,1,n_CTF])).^2.*reshape(weight_2d_k_p_r_,[n_k_p_r,1,1]),1+0)),[n_CTF,n_CTF]);
figure(1+nf);nf=nf+1;clf;figmed;figbeach();
[~,ij_srt_] = sort(index_ncluster_from_nCTF_,'ascend'); index_srt_ = ij_srt_-1;
subplot_{1} = subplot(1,3,1); imagesc(CTF_k_p_r_kC__(:,1+index_srt_));
axisnotick; xlabel('nCTF','Interpreter','none'); ylabel('n_k_p_r','Interpreter','none');
subplot_{2} = subplot(1,3,2); imagesc(CTF_k_p_l2_CC__(1+index_srt_,1+index_srt_));
axisnotick; xlabel('nCTF','Interpreter','none'); ylabel('nCTF','Interpreter','none');
subplot(1,3,3); cla;
c__ = colormap('lines'); n_c = size(c__,1);
hold on;
for nCTF=0:n_CTF-1;
ncluster = index_ncluster_from_nCTF_(1+index_srt_(1+nCTF));
nc = max(0,min(n_c-1,mod(ncluster,n_c)));
plot(nCTF,CTF_k_p_l2_C_(1+index_srt_(1+nCTF)),'ko','MarkerFaceColor',c__(1+nc,:));
end;%for nCTF=0:n_CTF-1;
axisnotick; xlabel('nCTF','Interpreter','none'); ylabel('CTF_k_p_l2_C_','Interpreter','none');
colormap(subplot_{1},colormap_beach());
colormap(subplot_{2},colormap_beach());
end;%if flag_disp;

if (flag_verbose>0); disp(sprintf('[finished %s]',str_thisfunction)); end;
