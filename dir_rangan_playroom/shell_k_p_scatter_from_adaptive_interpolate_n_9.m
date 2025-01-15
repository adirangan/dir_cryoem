function ...
[ ...
 scatter_from_tensor_sba__ ...
,dscatterda_from_tensor_sba__ ...
,dscatterdb_from_tensor_sba__ ...
,ddscatterdaa_from_tensor_sba__ ...
,ddscatterdab_from_tensor_sba__ ...
,ddscatterdbb_from_tensor_sba__ ...
] = ...
shell_k_p_scatter_from_adaptive_interpolate_n_9( ...
 n_order ...
,n_k_p_a ...
,n_k_p_b_ ...
,k_p_a_lim_ ...
,n_scatter ...
,k_s_a_ ...
,k_s_b_ ...
);
%%%%%%%%;
% We interpret the k_p_a_ array as a standard tensor array (e.g., height). ;
% We interpret the k_p_b_ array as periodic (e.g., angle). ; 
% However, the number of k_p_b_ points depends on nk_p_a. ;
% For this we store: ;
% n_k_p_b :=n_k_p_b_(1+nk_p_a) ;
% we assume that k_p_b_lim_ := [0,2*pi]. ;
%%%%%%%%;
% returns sparse matrix encoding the n_order interpolation operator between a tensor_grid ;
% using k_p_?_. ;
%%%%%%%%;
% We presume these points are associated with an unrolled array: ;
% a_k_p_ba_. ;
%%%%;
% The n_scatter points to be interpolated have coordinates stored in ;
% k_s_a_, ;
% k_s_b_, ;
% ;
% The (n_scatter)-by-(n_unrolled) scatter_from_tensor_sba__ matrix stores the ;
% interpolation weights linking a_k_p_ba_ to the unrolled a_k_s_. ;
%%%%;

str_thisfunction = 'shell_k_p_scatter_from_adaptive_interpolate_n_9';

if (nargin<7);
%%%%%%%%;
disp(sprintf(' %% testing %s',str_thisfunction));
disp(sprintf(' %% see test_%s',str_thisfunction));
disp(sprintf(' %% returning')); return;
%%%%%%%%;
end;%if (nargin<6);

flag_verbose = 0;
flag_disp = 0; nf=0;
flag_check = 0;
flag_d = (nargout>=2);
flag_dd = (nargout>=4);

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

na=0;
if (nargin<1+na); n_order=[]; end; na=na+1;
if (nargin<1+na); n_k_p_a=[]; end; na=na+1;
if (nargin<1+na); n_k_p_b_=[]; end; na=na+1;
if (nargin<1+na); k_p_a_lim_=[]; end; na=na+1;
if (nargin<1+na); n_scatter=[]; end; na=na+1;
if (nargin<1+na); k_s_a_=[]; end; na=na+1;
if (nargin<1+na); k_s_b_=[]; end; na=na+1;

if numel(n_k_p_b_)==1;
n_k_p_b_ = repmat(n_k_p_b_(1+0),[n_k_p_a,1]);
end;%if numel(n_k_p_b_)==1;


if n_k_p_b_(1+0)~=1; disp(sprintf(' %% Warning, n_k_p_b_(1+0) %d in %s',n_k_p_b_(1+0),str_thisfunction)); end;
if n_k_p_b_(end)~=1; disp(sprintf(' %% Warning, n_k_p_b_(end) %d in %s',n_k_p_b_(end),str_thisfunction)); end;
for nk_p_a=1:n_k_p_a-2;
if mod(n_k_p_b_(1+nk_p_a),2)~=0; disp(sprintf(' %% Warning, mod(n_k_p_b_(1+nk_p_a),2) %d in %s',mod(n_k_p_b_(1+nk_p_a),2),str_thisfunction)); end;
end;%for nk_p_a=1:n_k_p_a-2;

n_k_p_b_max = max(n_k_p_b_);
n_k_p_b_sum = sum(n_k_p_b_);
n_k_p_b_csum_ = cumsum([0;n_k_p_b_]);
if (flag_verbose>0); disp(sprintf(' %% n_k_p_b_max %d n_k_p_b_sum %d',n_k_p_b_max,n_k_p_b_sum)); end;
if (flag_verbose>0); disp(sprintf(' %% n_k_p_b_: %s',num2str(transpose(n_k_p_b_),'%d '))); end;
if (flag_verbose>0); disp(sprintf(' %% n_k_p_b_csum_: %s',num2str(transpose(n_k_p_b_csum_),'%d '))); end;

k_p_b_lim__ = [0,2*pi];
if numel(k_p_b_lim__)==2;
k_p_b_lim__ = repmat([k_p_b_lim__(1+0),k_p_b_lim__(1+1)],[n_k_p_a,1]);
end;%if numel(k_p_b_lim__)==2;
if (flag_verbose>1); disp(sprintf(' %% k_p_b_lim__:')); end;
if (flag_verbose>1); disp(transpose(k_p_b_lim__)); end;

flag_b_shift_ = [ ones(n_k_p_a-1,1) ; zeros(n_k_p_a,1) ; ones(n_k_p_a-1,1) ];
nk_p_a_extend_ = transpose([ n_k_p_a:-1:2 , 1:n_k_p_a , n_k_p_a-1:-1:1 ]) - 1;
n_k_p_a_extend = numel(nk_p_a_extend_);
n_k_p_b_extend_ = n_k_p_b_(1+nk_p_a_extend_);
n_k_p_b_extend_max = max(n_k_p_b_extend_);
n_k_p_b_extend_sum = sum(n_k_p_b_extend_);
n_k_p_b_extend_csum_ = cumsum([0;n_k_p_b_extend_]);
if (flag_verbose>0); disp(sprintf(' %% n_k_p_b_extend_max %d n_k_p_b_extend_sum %d',n_k_p_b_extend_max,n_k_p_b_extend_sum)); end;
if (flag_verbose>0); disp(sprintf(' %% n_k_p_b_extend_: %s',num2str(transpose(n_k_p_b_extend_),'%d '))); end;
if (flag_verbose>0); disp(sprintf(' %% n_k_p_b_extend_csum_: %s',num2str(transpose(n_k_p_b_extend_csum_),'%d '))); end;
k_p_a_extend_lim_ = [ k_p_a_lim_(1) - diff(k_p_a_lim_) , k_p_a_lim_(1+1) + diff(k_p_a_lim_) ];
k_p_b_extend_lim__ = repmat([0,2*pi],[n_k_p_a_extend,1]);

index_row_tensor_extend_from_tensor_ba_ = zeros(n_k_p_b_extend_sum,1);
index_col_tensor_extend_from_tensor_ba_ = zeros(n_k_p_b_extend_sum,1);
for nk_p_a_extend=0:n_k_p_a_extend-1;
nk_p_a = nk_p_a_extend_(1+nk_p_a_extend);
n_k_p_b = n_k_p_b_(1+nk_p_a);
n_k_p_b_extend = n_k_p_b_extend_(1+nk_p_a_extend);
assert(n_k_p_b==n_k_p_b_extend);
n_k_p_b_csum = n_k_p_b_csum_(1+nk_p_a);
n_k_p_b_extend_csum = n_k_p_b_extend_csum_(1+nk_p_a_extend);
flag_b_shift = flag_b_shift_(1+nk_p_a_extend);
index_nb_ = transpose(0:n_k_p_b-1);
index_nb_shift_ = circshift(index_nb_,round(n_k_p_b/2));
index_nba_ = index_nb_ + n_k_p_b_csum ;
index_nb_shifta_ = index_nb_shift_ + n_k_p_b_csum ;
index_nba_extend_ = index_nb_ + n_k_p_b_extend_csum ;
index_row_tensor_extend_from_tensor_ba_(1+index_nba_extend_) = index_nba_extend_ ;
if flag_b_shift==0; index_col_tensor_extend_from_tensor_ba_(1+index_nba_extend_) = index_nba_ ; end;
if flag_b_shift==1; index_col_tensor_extend_from_tensor_ba_(1+index_nba_extend_) = index_nb_shifta_ ; end;
end;%for nk_p_a_extend=0:n_k_p_a_extend-1;
tensor_extend_from_tensor_baba__ = sparse(1+index_row_tensor_extend_from_tensor_ba_,1+index_col_tensor_extend_from_tensor_ba_,1,n_k_p_b_extend_sum,n_k_p_b_sum);

[ ...
 scatter_from_tensor_extend_sba__ ...
 dscatterda_from_tensor_extend_sba__ ...
 dscatterdb_from_tensor_extend_sba__ ...
 ddscatterdaa_from_tensor_extend_sba__ ...
 ddscatterdab_from_tensor_extend_sba__ ...
 ddscatterdbb_from_tensor_extend_sba__ ...
] = ...
cylinder_k_c_scatter_from_tensor_adaptive_interpolate_n_9( ...
 n_order ...
,n_k_p_a_extend ...
,n_k_p_b_extend_ ...
,k_p_a_extend_lim_ ...
,k_p_b_extend_lim__ ...
,n_scatter ...
,k_s_a_ ...
,k_s_b_ ...
);

scatter_from_tensor_sba__ = scatter_from_tensor_extend_sba__ * tensor_extend_from_tensor_baba__ ;
if flag_d;
dscatterda_from_tensor_sba__ = dscatterda_from_tensor_extend_sba__ * tensor_extend_from_tensor_baba__ ;
dscatterdb_from_tensor_sba__ = dscatterdb_from_tensor_extend_sba__ * tensor_extend_from_tensor_baba__ ;
end;%if flag_d;
if flag_dd;
ddscatterdaa_from_tensor_sba__ = ddscatterdaa_from_tensor_extend_sba__ * tensor_extend_from_tensor_baba__ ;
ddscatterdab_from_tensor_sba__ = ddscatterdab_from_tensor_extend_sba__ * tensor_extend_from_tensor_baba__ ;
ddscatterdbb_from_tensor_sba__ = ddscatterdbb_from_tensor_extend_sba__ * tensor_extend_from_tensor_baba__ ;
end;%if flag_dd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% returning')); end; return;


