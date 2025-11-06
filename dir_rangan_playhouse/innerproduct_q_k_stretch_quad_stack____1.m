function M_k_q_rwlM____ = innerproduct_q_k_stretch_quad_stack____1(n_k_p_r,n_w_,n_M,M_k_q_wkM__,l_max,flag_each);

str_thisfunction = 'innerproduct_q_k_stretch_quad_stack____1';

na=0;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); M_k_q_wkM__=[]; end; na=na+1;
if (nargin<1+na); l_max=[]; end; na=na+1;
if (nargin<1+na); flag_each=[]; end; na=na+1;

if isempty(flag_each); flag_each = 0; end;

if  flag_each | numel(unique(n_w_))> 1;
n_w_max = max(n_w_);
M_k_q_rwlM____ = zeros(n_k_p_r,n_w_max,1+2*l_max,n_M);
for nM=0:n_M-1;
M_k_q_rwlM____(:,:,:,1+nM) = innerproduct_q_k_stretch_quad_stack___0(n_k_p_r,n_w_,M_k_q_wkM__(:,1+nM),l_max);
end;%for nM=0:n_M-1;
end;%if  flag_each | numel(unique(n_w_))> 1;

if ~flag_each & numel(unique(n_w_))==1;
n_w = n_w_(1+0); n_w_max = n_w;
if mod(n_w,2)==1; disp(sprintf(' %% Warning, odd n_w in %s',str_thisfunction)); end;
tmp_q_ = transpose([[0:+n_w/2-1],[-n_w/2:-1]]);
M_k_q_wklM____ = zeros(n_w,n_k_p_r,1+2*l_max,n_M);
M_k_q_kwM___ = permute(reshape(M_k_q_wkM__,[n_w,n_k_p_r,n_M]),1+[1,0,2]);
M_k_q_rwlM____(:,:,1+l_max+0,:) = M_k_q_kwM___;
for l_val=-l_max:+l_max;
l_use = -l_val;
M_k_q_rwlM____(:,:,1+l_max+l_val,:) = circshift(M_k_q_kwM___,l_use,1+1);
if l_use> 0; M_k_q_rwlM____(:,1+[n_w/2-0*l_use:n_w/2+1*l_use],1+l_max+l_val,:) = 0; end;
if l_use==0; M_k_q_rwlM____(:,1+[n_w/2-0*l_use:n_w/2+0*l_use],1+l_max+l_val,:) = 0; end;
if l_use< 0; M_k_q_rwlM____(:,1+[n_w/2+1*l_use:n_w/2+0*l_use],1+l_max+l_val,:) = 0; end;
end;%for l_val=-l_max:+l_max;
end;%if ~flag_each & numel(unique(n_w_))==1;

