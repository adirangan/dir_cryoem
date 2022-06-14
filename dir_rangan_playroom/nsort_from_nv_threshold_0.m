function ...
[ ...
 yrank_avg_from_nx_ ...
] = ...
nsort_from_nv_threshold_0( ...
 flag_ij_vs_val ...
,val_x_ ...
,val_y_ ...
,str_dir_x ...
,str_dir_y ...
);

%%%%%%%%;
% determines the average rank (1-based) of y_, ;
% assuming that we threshold by the rank of x_. ;
% If flag_ij_vs_val==1, ;
% we assume val_x_ and val_y_ are already ij_nsort_from_nx_ and ij_nsort_from_ny_. ;
% str_dir_x and str_dir_y list the directions to sort val_x_ and val_y_ (default ascend). ;
%%%%%%%%;

if nargin<1;
val_x_ = [102 101 104 103 205];
val_y_ = [ 21  22  23  35  24];
yrank_avg_from_nx_ = nsort_from_nv_threshold_0(0,val_x_,val_y_);
disp(yrank_avg_from_nx_);
disp('returning');return;
end;%if nargin<1;

na=0;
if (nargin<1+na); flag_ij_vs_val=[]; end; na=na+1;
if (nargin<1+na); val_x_=[]; end; na=na+1;
if (nargin<1+na); val_y_=[]; end; na=na+1;
if (nargin<1+na); str_dir_x=[]; end; na=na+1;
if (nargin<1+na); str_dir_y=[]; end; na=na+1;
if isempty(flag_ij_vs_val); flag_ij_vs_fal = 0; end; %<-- assume given values. ;
if isempty(str_dir_x); str_dir_x = 'ascend'; end; %<-- direction to sort x. ;
if isempty(str_dir_y); str_dir_y = 'ascend'; end; %<-- direction to sort y. ;

n_v = numel(val_x_);
assert(n_v==numel(val_y_));

if flag_ij_vs_val==0;
[~,ij_nx_from_nsort_] = sort(val_x_(:),str_dir_x);
[~,ij_nsort_from_nx_] = sort(ij_nx_from_nsort_,'ascend');
[~,ij_ny_from_nsort_] = sort(val_y_(:),str_dir_y);
[~,ij_nsort_from_ny_] = sort(ij_ny_from_nsort_,'ascend');
end;%if flag_ij_vs_val==0;
if flag_ij_vs_val==1;
ij_nsort_from_nx_ = val_x_(:);
ij_nsort_from_ny_ = val_y_(:);
end;%if flag_ij_vs_val==1;

%%%%%%%%;
% Now at this point we construct the matrices: ;
% A__ = sparse(ij_nsort_from_nx_,ij_nsort_from_ny_,1,n_v,n_v);
% B__ = cumsum(A__,1,'reverse');
% C__ = bsxfun(@times,B__,reshape(1:n_v,[1,n_v]));
%%%%%%%%;
% For example, if: ;
% val_x_ = [102 101 104 103 205];
% val_y_ = [ 21  22  23  35  24];
% ij_nsort_from_nx_ = [2 1 4 3 5];
% ij_nsort_from_ny_ = [1 2 3 5 4];
% then: ;
% A__ = ...
%  [ 0 1 0 0 0 ...
%    1 0 0 0 0 ...
%    0 0 0 0 1 ...
%    0 0 1 0 0 ...
%    0 0 0 1 0 ] ;
% will be invariant under permutation, ;
% so long as the val_x_ and val_y_ are each unique. ;
% I.e., given a permutation p_ = randperm(n_v), ;
% the matrix A__ calculated from val_x_(p_) and val_y_(p_) ;
% will equal the original A__. ;
% For this example, ;
% B__ = ...
%  [ 1 1 1 1 1 ...
%    1 0 1 1 1 ...
%    0 0 1 1 1 ...
%    0 0 1 1 0 ...
%    0 0 0 1 0 ] ;
% and in general, given an x-rank threshold nx, the entries of B__(1+nx,:) ;
% indicate which y-ranks are included in the data-set ;
% after excluding all entries with an x-rank lower than nx. ;
% For this example, ;
% C__ = ...
%  [ 1 2 3 4 5 ...
%    1 0 3 4 5 ...
%    0 0 3 4 5 ...
%    0 0 3 4 0 ...
%    0 0 0 4 0 ] ;
% such that sum(C__,2) can be used to determine the average y-rank ;
% associated with each x-rank threshold. ;
% For example, if the ij_nsort_from_nx_ and ij_nsort_from_ny_ are unique, ;
% we can simply write:
% yrank_avg_from_nx_ = sum(C__,2)./transpose(n_v:-1:1);
% However, for generality, we use the expression below instead. ;
%%%%%%%%;
A__ = sparse(ij_nsort_from_nx_,ij_nsort_from_ny_,1,n_v,n_v);
B__ = cumsum(A__,1,'reverse');
C__ = bsxfun(@times,B__,reshape(1:n_v,[1,n_v]));
yrank_avg_from_nx_ = sum(C__,2)./sum(B__,2);



