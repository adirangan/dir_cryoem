function ...
[ ...
 index_quad_from_data_wMj__ ...
,n_index_quad_from_data_wM_ ...
,d2_quad_from_data_wMj__ ...
] = ...
box3d_S2_search_0( ...
 quad_n_all ...
,quad_k_c_qd__ ...
,n_w ...
,n_M ...
,data_k_c_wMd__ ...
,d_req ...
,flag_knn ...
);
%%%%%%%%;
% This returns a list of indices and squared-distances from quad_k_c_qd__ which are allocated to each data_k_c_wMd__. ;
% We assume n_d==3. ;
% We also assume that the points all lie on the sphere of radius 1. ;
%%%%%%%%;
% Inputs. ;
% quad_n_all: integer number of quadrature-points. ;
% quad_k_c_qd__: double array of size (quad_n_all,n_d) containing quadrature-points. ;
% n_w = integer number of points per image-ring. ;
% n_M = integer number of images. ;
% data_k_c_wMd__: double array of size (n_w*n_M,n_d) containing image points. ;
% d_req: double storing the requested distance (i.e., lower-bound) for inclusion. ;
% flag_knn: integer maximum number of nearest-neighbors requested (ignored if flag_knn<=0). ;
%%%%%%%%;
% Outputs: ;
% n_j: integer storing maximum number of quadrature-points allocated to each image-point. ;
% index_quad_from_data_wMj__: integer array of size (n_w*n_M,n_j). ;
%                             tmp_index_ = index_quad_from_data_wMj__(1+tab,:) ;
%                             stores the indices from quad_k_c_qd__ allocated to image-point data_, ;
%                             with: tab := nw+nM*n_w and data_ := data_k_c_wMd__(1+tab,:). ;
%                             dummy-values are stored as -1. ;
% n_index_quad_froM_data_wM_: integer array of size(n_w*n_M,1). ;
%                             stores number of non-dummy values in index_quad_from_data_wMj__(tab,:), ;
%                             with: tab := nw+nM*n_w. ;
% d2_quad_from_data_wMj__: double array of size (n_w*n_M,n_j). ;
%                          tmp_d2_ = d2_quad_from_data_wMj__(1+tab,:) ;
%                          stores the squared-distances from quad_k_c_qd__ allocated to image-point data_, ;
%                          with: tab := nw+nM*n_w and data_ := data_k_c_wMd__(1+tab,:). ;
%                          dummy-values are stored as +Inf. ;
%%%%%%%%;

str_thisfunction = 'box3d_S2_search_0';

if nargin<1;
flag_verbose=1; flag_disp = 0; nf=0;
if (flag_verbose>0); disp(sprintf(' %% testing %s',str_thisfunction)); end;
rng(0);
n_d = 3;
quad_n_all = 10828;
quad_k_c_qd__ = randn(quad_n_all,n_d);
quad_k_c_qd__ = bsxfun(@rdivide,quad_k_c_qd__,sqrt(sum(quad_k_c_qd__.^2,2)));
n_w = 128;
n_M = 671692/n_w;
data_n_all = n_w*n_M;
data_k_c_ad__ = randn(data_n_all,n_d);
data_k_c_ad__ = bsxfun(@rdivide,data_k_c_ad__,sqrt(sum(data_k_c_ad__.^2,2)));
d_req = 0.125;
%%%%;
tmp_t = tic();
flag_knn = -1;
[ ...
 index_quad_from_data_aj__ ...
,n_index_quad_from_data_a_ ...
,d2_quad_from_data_aj__ ...
] = ...
box3d_S2_search_0( ...
 quad_n_all ...
,quad_k_c_qd__ ...
,n_w ...
,n_M ...
,data_k_c_ad__ ...
,d_req ...
,flag_knn ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% box3d_S2_search_0 (not sorting): %0.6fs',tmp_t)); end;
%%%%;
tmp_t = tic();
flag_knn = 0;
[ ...
 index_quad_from_data_aj__ ...
,n_index_quad_from_data_a_ ...
,d2_quad_from_data_aj__ ...
] = ...
box3d_S2_search_0( ...
 quad_n_all ...
,quad_k_c_qd__ ...
,n_w ...
,n_M ...
,data_k_c_ad__ ...
,d_req ...
,flag_knn ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% box3d_S2_search_0 (yes sorting): %0.6fs',tmp_t)); end;
%%%%;
tmp_t = tic();
n_knn = ceil(quad_n_all/8);
[ij_knn_aj__,d1_knn_aj__] = knnsearch(quad_k_c_qd__,data_k_c_ad__,'K',n_knn); index_knn_aj__ = ij_knn_aj__ - 1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% knnsearch: %0.6fs',tmp_t)); end;
%%%%;
n_q = quad_n_all;
quad_qd__ = quad_k_c_qd__;
n_a = data_n_all;
data_ad__ = data_k_c_ad__;
flag_check=1;
if flag_check;
n_check_sum = 0; n_check_min = +Inf; n_check_max = -Inf;
for na=0:n_a-1;
n_j = min(n_index_quad_from_data_a_(1+na),n_knn);
tmp_ij_diff = min(efind(index_quad_from_data_aj__(1+na,1+[0:n_j-1]) - index_knn_aj__(1+na,1+[0:n_j-1]))); %<-- first difference. ;
assert(fnorm(d2_quad_from_data_aj__(1+na,1+[0:tmp_ij_diff-1-1]) - d1_knn_aj__(1+na,1+[0:tmp_ij_diff-1-1]).^2)<=1e-12); %<-- same before. ;
assert(d1_knn_aj__(1+na,tmp_ij_diff)>=d_req); %<-- missed neighbor is farther than d_req away. ;
n_check_sum = n_check_sum + (tmp_ij_diff-1);
n_check_min = min(n_check_min,(tmp_ij_diff-1));
n_check_max = max(n_check_max,(tmp_ij_diff-1));
end;%for na=0:n_a-1;
if (flag_verbose>0); disp(sprintf(' %% n_check_min %d n_check_max %d n_check_sum %d <-- %0.2f neighbors/point',n_check_min,n_check_max,n_check_sum,n_check_sum/n_a)); end;
end;%if flag_check;
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
markersize_use = 8;
markersize_big = 8;
linewidth_use = 0.5;
hold on;
plot3(quad_qd__(:,1+0),quad_qd__(:,1+1),quad_qd__(:,1+2),'.','MarkerSize',markersize_use,'Color','b');
plot3(data_ad__(:,1+0),data_ad__(:,1+1),data_ad__(:,1+2),'.','MarkerSize',markersize_use,'Color','r');
for na=0:min(32,n_a);%for na=0:n_a-1;
data_1d__ = data_ad__(1+na,:);
n_index_quad_from_data = n_index_quad_from_data_a_(1+na);
index_quad_from_data_j_ = index_quad_from_data_aj__(1+na,1+[0:n_index_quad_from_data-1]);
assert(sum(index_quad_from_data_j_>=0)==n_index_quad_from_data);
disp(sprintf(' %% na %d n_index_quad_from_data %d',na,n_index_quad_from_data));
quad_jd__ = quad_qd__(1+index_quad_from_data_j_,:);
line_jd2__ = cat(3,repmat(data_1d__,[n_index_quad_from_data,1]),quad_jd__);
line_0_2j__ = permute(reshape(line_jd2__(:,1+0,:),[n_index_quad_from_data,2]),[2,1]);
line_1_2j__ = permute(reshape(line_jd2__(:,1+1,:),[n_index_quad_from_data,2]),[2,1]);
line_2_2j__ = permute(reshape(line_jd2__(:,1+2,:),[n_index_quad_from_data,2]),[2,1]);
plot3(data_1d__(1+0),data_1d__(1+1),data_1d__(1+2),'o','MarkerSize',markersize_big,'Color','m');
plot3(line_0_2j__,line_1_2j__,line_2_2j__,'-','LineWidth',linewidth_use,'Color','g');
end;%for na=0:n_a-1;
hold off;
axis equal; axis vis3d; set(gca,'XTick',[],'YTick',[],'ZTick',[]);
end;%if flag_disp;
%%%%;
if (flag_verbose>0); disp('returning'); end; return;
%%%%%%%%;
end;%if nargin<1;
%%%%%%%%;

na=0;
if (nargin<1+na); quad_n_all=[]; end; na=na+1;
if (nargin<1+na); quad_k_c_qd__=[]; end; na=na+1;
if (nargin<1+na); n_w=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); data_k_c_wMd__=[]; end; na=na+1;
if (nargin<1+na); d_req=[]; end; na=na+1;
if (nargin<1+na); flag_knn=[]; end; na=na+1;
if isempty(d_req); d_req = 0.125; end;
d_req = max(1e-12,d_req);
if isempty(flag_knn); flag_knn = 0; end;

flag_verbose=1; flag_disp=0; nf=0;
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_d = 3; %<-- assume 3-dimensional. ;
if size(quad_k_c_qd__,2)~=3; disp(sprintf(' %% Warning, improper size in %s',str_thisfunction)); end;
if size(data_k_c_wMd__,2)~=3; disp(sprintf(' %% Warning, improper size in %s',str_thisfunction)); end;

quad_max_d_ = max(quad_k_c_qd__,[],1); quad_min_d_ = min(quad_k_c_qd__,[],1); quad_abs = max(max(abs(quad_max_d_)),max(abs(quad_min_d_)));
if quad_abs> 1+1e-12; disp(sprintf(' %% Warning, quad_k_c_qd__ not on unit-sphere in %s',str_thisfunction)); end;
data_max_d_ = max(data_k_c_wMd__,[],1); data_min_d_ = min(data_k_c_wMd__,[],1); data_abs = max(max(abs(data_max_d_)),max(abs(data_min_d_)));
if data_abs> 1+1e-12; disp(sprintf(' %% Warning, data_k_c_wMd__ not on unit-sphere in %s',str_thisfunction)); end;

n_q = quad_n_all; quad_qd__ = quad_k_c_qd__; %<-- relabel. ;
n_a = n_w*n_M; data_ad__ = data_k_c_wMd__; %<-- relabel. ;
%%%%%%%%%%%%%%%%;
% bin ranges from [0:n_bin-1]. ;
% tab = xbin + ybin*n_bin + zbin*n_bin*n_bin and ranges from [0,n_tab-1], with n_tab:=n_bin^3. ;
%%%%%%%%%%%%%%%%;
n_bin = 1+floor(1.0d0/d_req); n_tab = n_bin^3;
%%%%%%%%%%%%%%%%;
% restrict attention to those boxes that intersect S2. ;
%%%%%%%%%%%%%%%%;
tmp_t_1 = tic();
bin_lhs_b_ = transpose(2.0d0*d_req*(0 + [0:n_bin-1])) - 1.0d0;
[bin_lhs_0___,bin_lhs_1___,bin_lhs_2___] = ndgrid(bin_lhs_b_,bin_lhs_b_,bin_lhs_b_);
bin_lhs_02___ = bin_lhs_0___.^2; bin_lhs_12___ = bin_lhs_1___.^2; bin_lhs_22___ = bin_lhs_2___.^2;
bin_rhs_b_ = transpose(2.0d0*d_req*(1 + [0:n_bin-1])) - 1.0d0;
[bin_rhs_0___,bin_rhs_1___,bin_rhs_2___] = ndgrid(bin_rhs_b_,bin_rhs_b_,bin_rhs_b_);
bin_rhs_02___ = bin_rhs_0___.^2; bin_rhs_12___ = bin_rhs_1___.^2; bin_rhs_22___ = bin_rhs_2___.^2;
bin_max_r___ = sqrt( max(bin_rhs_02___,bin_lhs_02___) + max(bin_rhs_12___,bin_lhs_12___) + max(bin_rhs_22___,bin_lhs_22___) );
bin_min_r___ = sqrt( min(bin_rhs_02___,bin_lhs_02___) + min(bin_rhs_12___,bin_lhs_12___) + min(bin_rhs_22___,bin_lhs_22___) );
index_S2_t_ = efind( bin_min_r___<=1.0d0 & bin_max_r___>=1.0d0 );
tmp_t_1 = toc(tmp_t_1); if (flag_verbose>0); disp(sprintf(' %% index_S2_t_: %0.6fs',tmp_t_1)); end;
%%%%%%%%;
% set tab_shift. ;
%%%%%%%%;
tab_shift_ = reshape(bsxfun(@plus,reshape([-1:+1]*n_bin^0,[3,1,1]),bsxfun(@plus,reshape([-1:+1]*n_bin^1,[1,3,1]),reshape([-1:+1]*n_bin^2,[1,1,3]))),[1,27]);

tmp_t_0 = tic();
%%%%%%%%%%%%%%%%;
tmp_t_1 = tic();
nbin_quad_qd__ = max(0,min(n_bin-1,floor( (1.0d0+quad_qd__)/2.0d0 / d_req )));
ntab_quad_q_ = nbin_quad_qd__*(n_bin.^[0;1;2]);
ntab_from_quad_tq__ = sparse(1+ntab_quad_q_,1+[0:n_q-1],1,n_tab,n_q);
%ntab_from_quad_tq__ = transpose(sparse(1+[0:n_q-1],1+ntab_quad_q_,1,n_q,n_tab)); %<-- unclear if this is faster/slower than the previous line. ;
tmp_t_1 = toc(tmp_t_1); if (flag_verbose>0); disp(sprintf(' %% ntab_from_quad_tq__: %0.6fs',tmp_t_1)); end;
%%%%%%%%;
% This full-loop is slower than the pruned-loop below. ;
%%%%%%%%;
flag_calc=0;
if flag_calc;
tmp_t_1 = tic();
index_quad_nq_from_ntab__ = cell(n_tab,1);
ntab_from_quad_sum_t_ = zeros(n_tab,1);
for ntab=0:n_tab-1;
index_quad_nq_from_ntab__{1+ntab} = [];
tmp_index_ = efind(ntab_from_quad_tq__(1+ntab,:));
if ~isempty(tmp_index_);
index_quad_nq_from_ntab__{1+ntab} = tmp_index_;
ntab_from_quad_sum_t_(1+ntab) = numel(tmp_index_);
end;%if ~isempty(tmp_index_);
clear tmp_index_ ;
end;%for ntab=0:n_tab-1;
assert(fnorm(ntab_from_quad_sum_t_ - sum(ntab_from_quad_tq__,2))==0);
tmp_t_1 = toc(tmp_t_1); if (flag_verbose>0); disp(sprintf(' %% index_quad_nq_from_ntab__: %0.6fs',tmp_t_1)); end;
end;%flag_calc=0;
%%%%%%%%;
% This pruned-loop is faster than the full-loop above. ;
%%%%%%%%;
tmp_t_2 = tic();
index_quad_nq_from_ntab__ = cell(n_tab,1);
ntab_from_quad_sum_t_ = zeros(n_tab,1);
for nindex_S2 = 0:numel(index_S2_t_)-1;
ntab = index_S2_t_(1+nindex_S2);
index_quad_nq_from_ntab__{1+ntab} = [];
tmp_index_ = efind(ntab_from_quad_tq__(1+ntab,:));
if ~isempty(tmp_index_);
index_quad_nq_from_ntab__{1+ntab} = tmp_index_;
ntab_from_quad_sum_t_(1+ntab) = numel(tmp_index_);
end;%if ~isempty(tmp_index_);
end;%for nindex_S2 = 0:numel(index_S2_t_)-1;
assert(fnorm(ntab_from_quad_sum_t_ - sum(ntab_from_quad_tq__,2))==0);
tmp_t_2 = toc(tmp_t_2); if (flag_verbose>0); disp(sprintf(' %% index_quad_nq_from_ntab__ (opt): %0.6fs',tmp_t_2)); end;
%%%%%%%%;
n_j = min(n_q, 27 * max(ntab_from_quad_sum_t_));
%%%%%%%%;
tmp_t_1 = tic();
nbin_data_ad__ = max(0,min(n_bin-1,floor( (1.0d0+data_ad__)/2.0d0 / d_req )));
ntab_data_a_ = nbin_data_ad__*(n_bin.^[0;1;2]);
ntab_from_data_ta__ = sparse(1+ntab_data_a_,1+[0:n_a-1],1,n_tab,n_a);
%ntab_from_data_ta__ = transpose(sparse(1+[0:n_a-1],1+ntab_data_a_,1,n_a,n_tab)); %<-- unclear if this is faster/slower than the previous line. ;
tmp_t_1 = toc(tmp_t_1); if (flag_verbose>0); disp(sprintf(' %% ntab_from_data_ta__: %0.6fs',tmp_t_1)); end;
%%%%%%%%;
% use a pruned-loop again. ;
%%%%%%%%;
tmp_t_2 = tic();
index_data_na_from_ntab__ = cell(n_tab,1);
ntab_from_data_sum_t_ = zeros(n_tab,1);
for nindex_S2 = 0:numel(index_S2_t_)-1;
ntab = index_S2_t_(1+nindex_S2);
index_data_na_from_ntab__{1+ntab} = [];
tmp_index_ = efind(ntab_from_data_ta__(1+ntab,:));
if ~isempty(tmp_index_);
index_data_na_from_ntab__{1+ntab} = tmp_index_;
ntab_from_data_sum_t_(1+ntab) = numel(tmp_index_);
end;%if ~isempty(tmp_index_);
end;%for nindex_S2 = 0:numel(index_S2_t_)-1;
assert(fnorm(ntab_from_data_sum_t_ - sum(ntab_from_data_ta__,2))==0);
tmp_t_2 = toc(tmp_t_2); if (flag_verbose>0); disp(sprintf(' %% index_data_na_from_ntab__: %0.6fs',tmp_t_2)); end;
%%%%%%%%;
tmp_t_0 = toc(tmp_t_0); if (flag_verbose>0); disp(sprintf(' %% index_quad_nq_from_ntab__ and index_data_na_from_ntab__: %0.6fs',tmp_t_0)); end;

%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
markersize_use = 8;
markersize_big = 8;
linewidth_use = 0.5;
subplot(1,2,1);
hold on;
plot(0:n_tab-1,ntab_from_quad_sum_t_,'bo');
plot(0:n_tab-1,ntab_from_data_sum_t_,'rx');
xlabel('ntab');
ylabel('sum_t_','Interpreter','none');
subplot(1,2,2);
hold on;
plot3(quad_qd__(:,1+0),quad_qd__(:,1+1),quad_qd__(:,1+2),'.','MarkerSize',markersize_use,'Color','b');
plot3(data_ad__(:,1+0),data_ad__(:,1+1),data_ad__(:,1+2),'.','MarkerSize',markersize_use,'Color','r');
[~,ij_ntab] = max(ntab_from_quad_sum_t_); ntab = ij_ntab-1;
disp(sprintf(' %% ntab %d, ntab_from_quad_sum_t_(1+ntab) %d ntab_from_data_sum_t_(1+ntab) %d',ntab,ntab_from_quad_sum_t_(1+ntab),ntab_from_data_sum_t_(1+ntab)));
index_quad_nq_from_ntab_ = index_quad_nq_from_ntab__{1+ntab};
quad_sub_qd__ = quad_qd__(1+index_quad_nq_from_ntab_,:);
plot3(quad_sub_qd__(:,1+0),quad_sub_qd__(:,1+1),quad_sub_qd__(:,1+2),'o','MarkerSize',markersize_use,'Color','b');
index_data_na_from_ntab_ = index_data_na_from_ntab__{1+ntab};
data_sub_ad__ = data_ad__(1+index_data_na_from_ntab_,:);
plot3(data_sub_ad__(:,1+0),data_sub_ad__(:,1+1),data_sub_ad__(:,1+2),'o','MarkerSize',markersize_use,'Color','r');
hold off;
axis equal; axis vis3d; set(gca,'XTick',[],'YTick',[],'ZTick',[]);
disp('returning'); return;
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% Note that we do not clear temporary variables inside the loop, ;
% as this slows things down. ;
%%%%%%%%;
flag_rearrange_pos_vs_cur=0; %<-- rearranging within the loop seems faster. ;
tmp_t = tic();
%%%%%%%%;
neighbor_aj__ = -ones(n_a,n_j);
n_neighbor_a_ = zeros(n_a,1);
d2_neighbor_aj__ = +Inf*ones(n_a,n_j);
flag_touch_data_a_ = zeros(n_a,1);
index_data_na_from_ntab_lookup_ = zeros(n_a,1);
n_a_sub=0;
for nindex_S2 = 0:numel(index_S2_t_)-1;
ntab = index_S2_t_(1+nindex_S2);
ntab_from_data_sum = ntab_from_data_sum_t_(1+ntab);
if (ntab_from_data_sum> 0);
data_tab = ntab;
index_data_na_from_ntab_ = index_data_na_from_ntab__{1+data_tab};
if flag_rearrange_pos_vs_cur==0; tmp_index_insert_ = index_data_na_from_ntab_; end;
if flag_rearrange_pos_vs_cur==1; tmp_index_insert_ = n_a_sub+[0:ntab_from_data_sum-1]; end;
index_data_na_from_ntab_lookup_(1+tmp_index_insert_) = index_data_na_from_ntab_ ;
assert(numel(index_data_na_from_ntab_)==ntab_from_data_sum);
flag_touch_data_a_(1+tmp_index_insert_) = flag_touch_data_a_(1+tmp_index_insert_) + 1; %<-- rearrange later. ;
if (flag_verbose>2); disp(sprintf(' %% data_tab %d %%<-- ntab_from_data_sum %d',data_tab,ntab_from_data_sum)); end;
%%%%;
tmp_neighbor_q_ = -ones(n_j,1); n_neighbor_q = 0;
for quad_tab = data_tab + tab_shift_;
if (quad_tab>=0 & quad_tab< n_tab);
ntab_from_quad_sum = ntab_from_quad_sum_t_(1+quad_tab);
if (ntab_from_quad_sum> 0);
index_quad_nq_from_ntab_ = index_quad_nq_from_ntab__{1+quad_tab};
assert(numel(index_quad_nq_from_ntab_)==ntab_from_quad_sum);
tmp_neighbor_q_(1 + n_neighbor_q + [0:ntab_from_quad_sum-1]) = index_quad_nq_from_ntab_;
n_neighbor_q = n_neighbor_q + ntab_from_quad_sum;
if (flag_verbose>2); disp(sprintf(' %% %% quad_tab %d ntab_from_quad_sum %d %%<-- n_neighbor_q %d',quad_tab,ntab_from_quad_sum,n_neighbor_q)); end;
end;%if (ntab_from_quad_sum> 0);
end;%if (quad_tab>=0 & quad_tab< n_tab);
end;%for quad_tab = data_tab + tab_shift_;
if (n_neighbor_q>n_j); disp(sprintf(' %% n_neighbor_q %d n_j %d n_q %d sum(ntab_from_quad_sum_t_) %d',n_neighbor_q,n_j,n_q,sum(ntab_from_quad_sum_t_))); end;
assert(n_neighbor_q<=n_j);
%%%%;
tmp_neighbor_q_ = tmp_neighbor_q_(1+[0:n_neighbor_q-1]);
neighbor_aj__(1+tmp_index_insert_,1+[0:n_neighbor_q-1]) = repmat(reshape(tmp_neighbor_q_,[1,n_neighbor_q]),[ntab_from_data_sum,1]);
n_neighbor_a_(1+tmp_index_insert_) = n_neighbor_q;
d2_neighbor_aj__(1+tmp_index_insert_,1+[0:n_neighbor_q-1]) = reshape(sum(bsxfun(@minus,reshape(data_ad__(1+index_data_na_from_ntab_,:),[ntab_from_data_sum,1,n_d]),reshape(quad_qd__(1+tmp_neighbor_q_,:),[1,n_neighbor_q,n_d])).^2,3),[ntab_from_data_sum,n_neighbor_q]);
%%%%;
n_a_sub = n_a_sub + ntab_from_data_sum;
end;%if (ntab_from_data_sum> 0);
end;%for nindex_S2 = 0:numel(index_S2_t_)-1;
%%%%%%%%;
if flag_rearrange_pos_vs_cur==1; %<-- a-posteriori rearrangement seems slower. ;
[~,ij_data_na_from_ntab_refile_] = sort(index_data_na_from_ntab_lookup_,'ascend'); index_data_na_from_ntab_refile_ = ij_data_na_from_ntab_refile_ - 1;
index_data_na_from_ntab_lookup_ = index_data_na_from_ntab_lookup_(1+index_data_na_from_ntab_refile_);
flag_touch_data_a_ = flag_touch_data_a_(1+index_data_na_from_ntab_refile_);
neighbor_aj__ = neighbor_aj__(1+index_data_na_from_ntab_refile_,:);
n_neighbor_a_ = n_neighbor_a_(1+index_data_na_from_ntab_refile_);
d2_neighbor_aj__ = d2_neighbor_aj__(1+index_data_na_from_ntab_refile_,:);
end;%if flag_rearrange_pos_vs_cur==1;
%%%%%%%%;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% neighbor_aj__, n_neighbor_a_ and d2_neighbor_aj__: %0.6fs',tmp_t)); end;

%%%%%%%%;
n_j_max = max(n_neighbor_a_);
index_quad_from_data_wMj__ = neighbor_aj__(:,1+[0:n_j_max-1]) ; %<-- relabel. ;
n_index_quad_from_data_wM_ = n_neighbor_a_ ; %<-- relabel. ;
d2_quad_from_data_wMj__ = d2_neighbor_aj__(:,1+[0:n_j_max-1]) ; %<-- relabel. ;
%%%%%%%%;
assert(sum(flag_touch_data_a_==1)==n_a);

if flag_knn>=0; %<-- sort. ;
tmp_t = tic();
%%%%%%%%;
n_knn = n_j_max; if (flag_knn> 0); n_knn = min(flag_knn,n_knn); end;
index_quad_from_data_wMj__ = -ones(n_a,n_knn);
n_index_quad_from_data_wM_ = zeros(n_a,1);
d2_quad_from_data_wMj__ = +Inf*ones(n_a,n_knn);
for na=0:n_a-1;
n_neighbor = n_neighbor_a_(1+na);
neighbor_j_ = neighbor_aj__(1+na,1+[0:n_neighbor-1]);
d2_neighbor_j_ = d2_neighbor_aj__(1+na,1+[0:n_neighbor-1]);
[~,tmp_ij_] = sort(d2_neighbor_j_,'ascend');
tmp_n_j = min(n_neighbor,n_knn);
n_index_quad_from_data_wM_(1+na) = tmp_n_j;
index_quad_from_data_wMj__(1+na,1+[0:tmp_n_j-1]) = neighbor_j_(tmp_ij_(1+[0:tmp_n_j-1]));
d2_quad_from_data_wMj__(1+na,1+[0:tmp_n_j-1]) = d2_neighbor_j_(tmp_ij_(1+[0:tmp_n_j-1]));
end;%for na=0:n_a-1;
%%%%%%%%;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% sorting: %0.6fs',tmp_t)); end;
end;%if flag_knn>=0; %<-- sort. ;


if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;


