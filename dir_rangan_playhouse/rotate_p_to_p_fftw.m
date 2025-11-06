function M_p_ = rotate_p_to_p_fftw(n_r,n_w_,n_w_sum,S_p_,gamma_z,flag_each);

str_thisfunction = 'rotate_p_to_p_fftw';

%%%%%%%%;
if nargin<1;
%%%%%%%%;
flag_verbose=1;
n_w_int = 2;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi); str_L = 'L';
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_L ...
);
%%%%;
if (flag_verbose>0); disp(sprintf(' %% testing not adaptive')); end;
template_k_eq_d = -1;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
n_w_max = n_w_int*2*(l_max_max+1); n_w_0in_ = n_w_max*ones(n_k_p_r,1);
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_p_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 0*flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
,weight_3d_k_p_r_ ...
);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
if (flag_verbose>0); disp(sprintf(' %% n_w_ [ %d .. %d ]',n_w_(1),n_w_(end))); end;
n_S = 5;
S_k_p_wkS__ = randn(n_w_sum,n_S) + i*randn(n_w_sum,n_S);
gamma_z = +pi/3;
R_k_p_wkS__ = reshape(rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__,gamma_z),[n_w_sum,n_S]);
T_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
T_k_p_wkS__(:,1+nS) = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS),gamma_z,1);
end;%for nS=0:n_S-1;
fnorm_disp(flag_verbose,'T_k_p_wkS__',T_k_p_wkS__,'R_k_p_wkS__',R_k_p_wkS__,' %%<-- should be zero');
n_gamma_z = n_S; gamma_z_ = linspace(0,2*pi,n_gamma_z+1); gamma_z_ = transpose(gamma_z_(1:n_gamma_z));
R_k_p_wkS__ = reshape(rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__,gamma_z_),[n_w_sum,n_S]);
T_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
T_k_p_wkS__(:,1+nS) = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS),gamma_z_(1+nS));
end;%for nS=0:n_S-1;
fnorm_disp(flag_verbose,'T_k_p_wkS__',T_k_p_wkS__,'R_k_p_wkS__',R_k_p_wkS__,' %%<-- should be zero');
%%%%;
if (flag_verbose>0); disp(sprintf(' %% testing yes adaptive')); end;
template_k_eq_d = -1;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
l_max_ = round(2*pi*k_p_r_);
n_w_0in_ = n_w_int*2*(l_max_+1);
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_p_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 0*flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
,weight_3d_k_p_r_ ...
);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
if (flag_verbose>0); disp(sprintf(' %% n_w_ [ %d .. %d ]',n_w_(1),n_w_(end))); end;
n_S = 5;
S_k_p_wkS__ = randn(n_w_sum,n_S) + i*randn(n_w_sum,n_S);
gamma_z = +pi/3;
R_k_p_wkS__ = reshape(rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__,gamma_z),[n_w_sum,n_S]);
T_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
T_k_p_wkS__(:,1+nS) = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS),gamma_z,1);
end;%for nS=0:n_S-1;
fnorm_disp(flag_verbose,'T_k_p_wkS__',T_k_p_wkS__,'R_k_p_wkS__',R_k_p_wkS__,' %%<-- should be zero');
n_gamma_z = n_S; gamma_z_ = linspace(0,2*pi,n_gamma_z+1); gamma_z_ = transpose(gamma_z_(1:n_gamma_z));
R_k_p_wkS__ = reshape(rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__,gamma_z_),[n_w_sum,n_S]);
T_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
T_k_p_wkS__(:,1+nS) = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS),gamma_z_(1+nS));
end;%for nS=0:n_S-1;
fnorm_disp(flag_verbose,'T_k_p_wkS__',T_k_p_wkS__,'R_k_p_wkS__',R_k_p_wkS__,' %%<-- should be zero');
%%%%;
disp('returning'); return;
%%%%%%%%;
end;%if nargin<1;
%%%%%%%%;

na=0;
if nargin<1+na; n_r=[]; end; na=na+1;
if nargin<1+na; n_w_=[]; end; na=na+1;
if nargin<1+na; n_w_sum=[]; end; na=na+1;
if nargin<1+na; S_p_=[]; end; na=na+1;
if nargin<1+na; gamma_z=[]; end; na=na+1;
if nargin<1+na; flag_each=[]; end; na=na+1;

if isempty(flag_each); flag_each = 0; end;

n_S = idivide(cast(numel(S_p_),'int64'),cast(n_w_sum,'int64'));
if numel(S_p_)~=n_w_sum*n_S;
disp(sprintf(' %% Warning, n_w_sum %d n_S %d in %s',n_w_sum,n_S,str_thisfunction));
end;%if numel(S_p_)~=n_w_sum*n_S;

n_gamma_z = numel(gamma_z); gamma_z_ = reshape(gamma_z,[n_gamma_z,1]);
if (n_gamma_z> 1) & (n_S~=n_gamma_z);
disp(sprintf(' %% Warning, n_S %d n_gamma_z %d in %s',n_S,n_gamma_z,str_thisfunction));
end;%if (n_gamma_z> 1) & (n_S~=n_gamma_z);

%%%%%%%%%%%%%%%%;
if  flag_each | (numel(unique(n_w_))> 1);
M_p_ = zeros(n_w_sum*n_S,1);
%%%%%%%%;
if n_S==1;
%%%%;
% single transformation on adaptive grid. ;
%%%%;
M_p_ = zeros(size(S_p_));
n_w_max = n_w_(1+n_r-1);
C_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
q = nw - n_w_max/2;
C =  +cos(q*gamma_z) + i* -sin(q*gamma_z) ;
C_(1+nw) = C; %<-- C_(1+nw) = exp(i*(nw-n_w_max/2)*gamma_z). ;
end;%for nw=0:n_w_max-1;
nw_sum=0;
for nr=0:n_r-1;
n_w = n_w_(1+nr);
if (n_w>0);
fftw_out_ = fft(S_p_(1+nw_sum + (0:n_w-1)));
for nq=0:n_w-1;
q = nq;
if (q>n_w/2-1);
q = q - n_w;
end;%!if (q.ge.n_w/2-1);
C = C_(1+n_w_max/2 + q); %<-- C_(1+n_w_max/2+q) = exp(i*([n_w_max/2+q] - n_w_max/2)*gamma_z). ;
fftw_out_(1+nq) = fftw_out_(1+nq) * C;
end;%for nq=0:n_w-1;
fftw_0in_ = ifft(fftw_out_);
M_p_(1+nw_sum + (0:n_w-1)) = fftw_0in_;
nw_sum = nw_sum + n_w;
end;%if (n_w>0);
end;%for nr=0:n_r-1;
assert(nw_sum==n_w_sum);
%%%%;
end;%if n_S==1;
%%%%%%%%;
if n_S> 1;
for nS=0:n_S-1;
tmp_index_ = cast(nS*n_w_sum,'double') + [0:n_w_sum-1];
if n_gamma_z==1;
M_p_(1+tmp_index_) = rotate_p_to_p_fftw(n_r,n_w_,n_w_sum,S_p_(1+tmp_index_),gamma_z);
end;%if n_gamma_z==1;
if n_gamma_z> 1;
M_p_(1+tmp_index_) = rotate_p_to_p_fftw(n_r,n_w_,n_w_sum,S_p_(1+tmp_index_),gamma_z_(1+nS));
end;%if n_gamma_z> 1;
end;%for nS=0:n_S-1;
end;%if n_S> 1;
%%%%%%%%;
end;%if  flag_each | (numel(unique(n_w_))> 1);
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%;
if ~flag_each & (numel(unique(n_w_))==1);
n_w = n_w_(1+0); n_w_max = n_w;
tmp_q_ = transpose([[0:n_w_max/2-1],[-n_w_max/2:-1]]);
C_wz__ = exp(-i*tmp_q_*transpose(gamma_z_));
M_p_ = reshape(ifft(bsxfun(@times,reshape(C_wz__,[n_w,1,n_gamma_z]),fft(reshape(S_p_,[n_w,n_r,n_S]),[],1+0)),[],1+0),[numel(S_p_),1]);
end;%if ~flag_each & (numel(unique(n_w_))==1);
%%%%%%%%%%%%%%%%;

