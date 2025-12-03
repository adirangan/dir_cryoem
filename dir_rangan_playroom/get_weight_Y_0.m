function ...
[ ...
 l_max_upb ...
,l_max_ ...
,l_max_max ...
,n_m_max ...
,m_max_ ...
,n_ml_ ...
,n_ml_max ...
,n_ml_sum ...
,n_ml_csum_ ...
,Y_l_val_yk_ ...
,Y_m_val_yk_ ...
,Y_k_val_yk_ ...
,weight_Y_yk_ ...
] = ...
get_weight_Y_0( ...
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
);

str_thisfunction = 'get_weight_Y_0';

na=0;
if (nargin<1+na); flag_verbose=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;

if isempty(flag_verbose); flag_verbose = 0; end;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

l_max_upb = round(2*pi*k_p_r_max); %<-- typically sufficient for 2-3 digits of precision. ;
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
l_max_(1+nk_p_r) = max(0,min(l_max_upb,1+ceil(2*pi*k_p_r_(1+nk_p_r))));
end;%for nk_p_r=0:n_k_p_r-1;
n_ml_ = (l_max_+1).^2;
n_ml_max = max(n_ml_);
n_ml_sum = sum(n_ml_);
n_ml_csum_ = cumsum([0;n_ml_]);
l_max_max = max(l_max_);
m_max_ = -l_max_max : +l_max_max;
n_m_max = numel(m_max_);
Y_l_val_yk_ = zeros(n_ml_sum,1);
Y_m_val_yk_ = zeros(n_ml_sum,1);
Y_k_val_yk_ = zeros(n_ml_sum,1);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_ml = n_ml_(1+nk_p_r); 
tmp_l_val_y_ = zeros(n_ml,1);
tmp_m_val_y_ = zeros(n_ml,1);
nml=0; 
for l_val=0:l_max;
for m_val=-l_val:+l_val;
tmp_l_val_y_(1+nml) = l_val;
tmp_m_val_y_(1+nml) = m_val;
nml=nml+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
assert(nml==n_ml); assert(nml==(l_max+1)^2);
tmp_index_ = n_ml_csum_(1+nk_p_r) + (0:n_ml-1);
Y_l_val_yk_(1+tmp_index_) = tmp_l_val_y_;
Y_m_val_yk_(1+tmp_index_) = tmp_m_val_y_;
Y_k_val_yk_(1+tmp_index_) = k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
weight_Y_yk_ = zeros(n_ml_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_ml_csum_(1+nk_p_r) + (0:n_ml-1);
weight_Y_yk_(1+tmp_index_) = weight_3d_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;


 
