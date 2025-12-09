function ...
[ ...
 a_k_Y_yk_ ...
,n_quad_from_data_q_ ...
,a_k_p_qk_ ...
] = ...
qbp_uniform_over_n_k_p_r_10( ...
 qbp_eps ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,image_delta_x_M_ ...
,image_delta_y_M_ ...
,image_I_value_M_ ...
);
%%%%%%%%;
% Applies quadrature-back-propagation to solve for a_k_Y_yk_. ;
% Associates CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+nM)) with image M_k_p_wkM__(:,1+nM);
% ;
% Input: ;
% qbp_eps: pseudoinverse threshold used for solving local least-squares for ctf. ; %<-- Bug: this is actually ignored internally, and a value of 1e-12 is hard-coded. Not fixing because we are putting together a draft of the paper. ;
% n_k_p_r: integer number of shells. ;
% k_p_r_: real array of size n_k_p_r. k-values for each shell. ;
% l_max_: integer array of size n_k_p_r. l_max_(1+nk_p_r) is the order used for a_k_Y_yk_ on shell nk_p_r. ;
% n_w_: integer array of size n_k_p_r. n_w_(1+nk_p_r) is the number of inplane_gamma_z values recorded at that ring. ;
% n_M: integer number of images. ;
% M_k_p_wkM__: complex array of size (n_w_sum,n_M). stack of images in k_p_ format. ;
% index_nCTF_from_nM_: integer array of size n_M. index_nCTF_from_nM_(1+nM) is the (base 0) CTF_index used for image M_k_p_wkM__(:,1+nM). ;
%             This can be empty or set to 1, in which case the same CTF_k_p_ will be used for each image. ;
% CTF_k_p_wkC__: complex array of size(n_w_sum,n_CTF). stack of ctf-functions in k_p_ format. ;
%            If index_nCTF_from_nM_ is empty or set to 1, then we assume this contains only a single CTF_k_p_, ;
%            which will then be used for all images. ;
% euler_polar_a_M_: real array of size n_M. polar_a used for each image ;
% euler_azimu_b_M_: real array of size n_M. azimu_b used for each image ;
% euler_gamma_z_M_: real array of size n_M. gamma_z used for each image ;
% image_delta_x_M_: real array of size n_M. delta_x used for each image ;
% image_delta_y_M_: real array of size n_M. delta_y used for each image ;
% image_I_value_M_: real array of size n_M. I_value used for each image ;
% ;
% Output: ;
% a_k_Y_yk_: complex array of size n_lm_sum. output functions in k_Y_ format. ;
%%%%%%%%;

str_thisfunction = 'qbp_uniform_over_n_k_p_r_10';

na=0;
if (nargin<1+na); qbp_eps=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wkM__=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_wkC__=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_M_=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_M_=[]; end; na=na+1;
if (nargin<1+na); euler_gamma_z_M_=[]; end; na=na+1;
if (nargin<1+na); image_delta_x_M_=[]; end; na=na+1;
if (nargin<1+na); image_delta_y_M_=[]; end; na=na+1;
if (nargin<1+na); image_I_value_M_=[]; end; na=na+1;

flag_verbose=0;
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_w_ = n_w_(1:n_k_p_r); n_w_sum = sum(n_w_); n_w_max = max(n_w_); n_w_csum_ = cumsum([0;n_w_]);
if (std(diff(n_w_csum_),1)>1e-6); disp(sprintf(' %% Warning, flag_uniform_over_n_k_p_r==1 expected in %s',str_thisfunction)); end;
n_w_max = n_w_csum_(1+1);
assert(n_w_sum==n_w_max*n_k_p_r);

if isempty(qbp_eps); qbp_eps = 1e-3; end;
if isempty(image_delta_x_M_); image_delta_x_M_ = zeros(n_M,1); end;
if isempty(image_delta_y_M_); image_delta_y_M_ = zeros(n_M,1); end;
if isempty(image_I_value_M_); image_I_value_M_ = ones(n_M,1); end;
if isempty(index_nCTF_from_nM_); index_nCTF_from_nM_ = zeros(n_M,1); end;
if isempty(CTF_k_p_wkC__); CTF_k_p_wkC__ = ones(n_w_sum,1); end;
if (qbp_eps>1); qbp_eps = max(1e-12,0.1^(qbp_eps)); end; %<-- convert qbp_eps from nl10 scale to explicit value. ;

%%%%%%%%;
% construct CTF of same size as images. ;
%%%%%%%%
n_CTF = size(CTF_k_p_wkC__,1+1);
if (size(CTF_k_p_wkC__,1+0)==n_k_p_r); 
CTF_k_p_r_kC__ = CTF_k_p_wkC__;
CTF_k_p_wkC__ = reshape(repmat(reshape(CTF_k_p_r_kC__,[1,n_k_p_r,n_CTF]),[n_w_max,1,1]),[n_w_sum,n_CTF]);
end;%if (size(CTF_k_p_wkC__,1+0)==n_k_p_r);
%%%%%%%%;
CTF_k_p_wkM__ = CTF_k_p_wkC__(:,1+index_nCTF_from_nM_);
%%%%%%%%;

T_M_k_p_wkM__ = M_k_p_wkM__;
if ( (fnorm(image_delta_x_M_)>0) | (fnorm(image_delta_y_M_)>0) | (fnorm(image_I_value_M_-ones(n_M,1))>1e-6) );
T_M_k_p_wkM__ = bsxfun(@times,transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wkM__,+image_delta_x_M_,+image_delta_y_M_),reshape(image_I_value_M_,[1,n_M]));
end;%if ( (fnorm(image_delta_x_M_)>0) | (fnorm(image_delta_y_M_)>0) | (fnorm(image_I_value_M_-ones(n_M,1))>0) );

l_max_ = l_max_(1:n_k_p_r); l_max_max = max(l_max_);
n_lm_ = (1+l_max_(:)).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]);
quad_k_eq_d = sqrt(4*pi/max(1,n_lm_max));
if (flag_verbose>0); disp(sprintf(' %% quad_k_eq_d %0.6f',quad_k_eq_d)); end;

%%%%%%%%;
tmp_t = tic();
str_L = 'L'; flag_uniform_over_polar_a = 0;
[ ...
 n_q ...
,k_p_azimu_b_q_ ...
,k_p_polar_a_q_ ...
,weight_shell_q_ ...
,k_c_0_q_ ...
,k_c_1_q_ ...
,k_c_2_q_ ...
,n_polar_a ...
,polar_a_ ...
,n_azimu_b_ ...
] = ...
sample_shell_6( ...
 1.0 ...
,quad_k_eq_d ...
,str_L ...
,flag_uniform_over_polar_a ...
) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% sample_shell_6: %0.6fs',tmp_t)); end;
k_c_q3__ = [ k_c_0_q_ , k_c_1_q_ , k_c_2_q_ ];
[k_p_polar_a_unique_,ij_unique_,ij_return_] = unique(k_p_polar_a_q_); index_return_ = ij_return_ - 1;
n_polar_a_unique = numel(k_p_polar_a_unique_);
%%%%%%%%;

%%%%%%%%;
if ~exist('sqrt_2lp1_','var'); sqrt_2lp1_=[]; end;
if ~exist('sqrt_2mp1_','var'); sqrt_2mp1_=[]; end;
if ~exist('sqrt_rat0_m_','var'); sqrt_rat0_m_=[]; end;
if ~exist('sqrt_rat3_lm__','var'); sqrt_rat3_lm__=[]; end;
if ~exist('sqrt_rat4_lm__','var'); sqrt_rat4_lm__=[]; end;
tmp_t = tic();
parameter_ylgndr = struct('type','ylgndr');
parameter_ylgndr.flag_verbose = 0;
parameter_ylgndr.flag_d = 0;
parameter_ylgndr.flag_dd = 0;
[ ...
 parameter_ylgndr ...
,d0y_jlm___ ...
,~ ...
,~ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_m_ ...
,sqrt_rat3_lm__ ...
,sqrt_rat4_lm__ ...
] = ...
ylgndr_2( ...
 parameter_ylgndr ...
,l_max_max ...
,cos(k_p_polar_a_unique_) ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_m_ ...
,sqrt_rat3_lm__ ...
,sqrt_rat4_lm__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% ylgndr_2: %0.6fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% This version involves preliminary inflation of m_val_, ;
% as well as extraneous exp evaluations, ;
% all to avoid memory movement. ;
%%%%%%%%;

tmp_t = tic();
d0y_lmj___ = permute(d0y_jlm___,1+[1,2,0]); %<-- permutation before inflation is faster. ;
d0y_jml___ = permute(d0y_lmj___(:,1+abs(-l_max_max:+l_max_max),:),1+[2,1,0]); %<-- retain unique cos(polar_a_{j}) for now.; 
d0y_jy__ = zeros(n_polar_a_unique,n_lm_max);
for nl=0:numel(l_max_)-1;
l_max = l_max_(1+nl);
for l_val=0:l_max;
m_val_ = -l_val:+l_val;
d0y_jy__(:,1+l_val.^2+l_val+m_val_) = d0y_jml___(:,1+l_max_max+m_val_,1+l_val);
end;%for l_val=0:l_max;
end;%for nl=0:numel(l_max_)-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% d0y_jy__: %0.6fs',tmp_t)); end;

tmp_t = tic();
d0y_qy__ = d0y_jy__(1+index_return_,:);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% d0y_jy__: %0.6fs',tmp_t)); end;

tmp_t = tic();
m_val_y_ = zeros(n_lm_max,1);
for nl=0:numel(l_max_)-1;
l_max = l_max_(1+nl);
for l_val=0:l_max;
m_val_ = -l_val:+l_val;
m_val_y_(1+l_val.^2+l_val+m_val_) = m_val_;
end;%for l_val=0:l_max;
end;%for nl=0:numel(l_max_)-1;
expimb_qy__ = permute(reshape(exp(+i*bsxfun(@times,m_val_y_,reshape(k_p_azimu_b_q_,[1,n_q]))),[n_lm_max,n_q]),1+[1,0]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% expimb_qy__: %0.6fs',tmp_t)); end;

tmp_t = tic();
d0Y_qy__ = d0y_qy__.*expimb_qy__ / sqrt(4*pi);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% d0Y_qy__: %0.6fs',tmp_t)); end;

%%%%;
tmp_t = tic();
[ ...
 data_k_p_polar_a_wM__ ...
,data_k_p_azimu_b_wM__ ...
,data_k_c_0_wM__ ...
,data_k_c_1_wM__ ...
,data_k_c_2_wM__ ...
] = ...
cg_rhs_2( ...
 n_M ...
,n_w_max ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,+euler_gamma_z_M_ ...
);
data_k_c_wM3__ = [ data_k_c_0_wM__(:) , data_k_c_1_wM__(:) , data_k_c_2_wM__(:) ];
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% cg_rhs_2: %0.6fs',tmp_t)); end;
%%%%;
tmp_t = tic();
n_wM = n_w_max*n_M;
ij_nq_from_nwM_ = knnsearch(k_c_q3__,data_k_c_wM3__,'K',1); index_nq_from_nwM_ = ij_nq_from_nwM_ - 1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% knnsearch: %0.6fs',tmp_t)); end;
tmp_t = tic();
quad_from_data_qwM__ = sparse(1+index_nq_from_nwM_,1:n_wM,1,n_q,n_wM);
n_quad_from_data_q_ = quad_from_data_qwM__*ones(n_wM,1); %<-- number of data-points per quadrature-point. ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% quad_from_data_qwM__: %0.6fs',tmp_t)); end;
data_from_quad_wMq__ = bsxfun(@rdivide,transpose(quad_from_data_qwM__),max(1,transpose(n_quad_from_data_q_)));
tmp_t = tic();
CTF_k_p_wMk__ = reshape(permute(reshape(CTF_k_p_wkM__,[n_w_max,n_k_p_r,n_M]),1+[0,2,1]),[n_wM,n_k_p_r]);
CTF2_k_p_qk__ = quad_from_data_qwM__*abs(CTF_k_p_wMk__).^2;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% CTF2_k_p_qk__: %0.6fs',tmp_t)); end;
tmp_t = tic();
T_M_k_p_wMk__ = reshape(permute(reshape(T_M_k_p_wkM__,[n_w_max,n_k_p_r,n_M]),1+[0,2,1]),[n_wM,n_k_p_r]);
CTF_T_M_k_p_wMk__ = CTF_k_p_wMk__.*T_M_k_p_wMk__;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% CTF_T_M_k_p_wMk__: %0.6fs',tmp_t)); end;
%%%%;
tmp_t = tic();
a_k_p_qk__ = bsxfun(@rdivide,quad_from_data_qwM__*CTF_T_M_k_p_wMk__,max(1e-12,CTF2_k_p_qk__));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% a_k_p_qk__: %0.6fs',tmp_t)); end;
tmp_t = tic();
a_k_Y_yk__ = conj(transpose(bsxfun(@times,reshape(weight_shell_q_,[n_q,1]),d0Y_qy__)))*a_k_p_qk__;
a_k_Y_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,a_k_Y_yk__);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% a_k_Y_yk_: %0.6fs',tmp_t)); end;
%%%%;
if nargout>2; a_k_p_qk_ = a_k_p_qk__(:); end;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
