function ...
[ ...
 d0W_ ...
,V_lmm___ ...
,L_lm__ ...
,d1W_ ...
,d2W_ ...
] = ...
wignerd_c( ...
 n_l ...
,beta ...
,V_lmm___ ...
,L_lm__ ...
) ;
% Generates wigner-d matrices up to n_l ;
% See Gumerov and Duraiswami 2014. ; %<-- I could not get this to work. ;
% Try Feng, Wang, Yang, Jin, 2015. ; %<-- This works. ;
% Test with: ;
%{
%%%%;
flag_verbose = 1; nf=0;
l_max = 16; n_l = l_max;
m_max_ = transpose(-l_max:+l_max);
n_m_max = numel(m_max_);
beta = 5*pi/12;
d0V_ = wignerd_b(l_max,beta);
tmp_V__ = d0V_{1+l_max};
%%%%;
l_val = l_max;
m_val_ = transpose([-l_val:+l_val]);
n_m_val = numel(m_val_);
X_ = sqrt( (l_val + m_val_) .* (1 + l_val - m_val_) ); 
S__ = spdiags([-X_ , +flip(X_)],[+1,-1],n_m_val,n_m_val);
[V__,L__] = eigs(S__,n_m_val); L_ = diag(L__);
tmp_W__ = zeros(n_m_val,n_m_val);
for m0_val=-l_val:+l_val;
for m1_val=-l_val:+l_val;
tmp_W = 0.0 + i*0.0;
for nl=0:n_m_val-1;
tmp_W = tmp_W + V__(1+l_val+m0_val,1+nl) * exp(+L_(1+nl)*beta/2) * conj(V__(1+l_val+m1_val,1+nl));
end;%for nl=0:n_m_val-1;
tmp_W__(1+l_val+m0_val,1+l_val+m1_val) = tmp_W;
end;%for m1_val=-l_val:+l_val;
end;%for m0_val=-l_val:+l_val;
tmp_U__ = V__ * bsxfun(@times,reshape(exp(+L_*beta/2),[n_m_val,1]),ctranspose(V__));
if (flag_verbose); disp(sprintf(' %% tmp_U__ vs tmp_W__: %0.16f',fnorm(tmp_U__-tmp_W__)/fnorm(tmp_U__))); end;
%%%%;
sgn_ = (-1).^m_val_.*(m_val_>=0) + (+1).*(m_val_< 0);
tmp_W__ = tmp_W__.*(sgn_*transpose(sgn_));
tmp_U__ = tmp_U__.*(sgn_*transpose(sgn_));
if (flag_verbose); disp(sprintf(' %% tmp_V__ vs tmp_W__: %0.16f',fnorm(tmp_V__-tmp_W__)/fnorm(tmp_V__))); end;
if (flag_verbose); disp(sprintf(' %% tmp_V__ vs tmp_U__: %0.16f',fnorm(tmp_V__-tmp_U__)/fnorm(tmp_V__))); end;
%%%%;
if (flag_verbose>1);
figure(1+nf);nf=nf+1;clf;figmed; fig80s;
p_row = 1; p_col = 3; np=0;
subplot(p_row,p_col,1+np);np=np+1;
imagesc(real(tmp_W__),[-1,+1]); axis image; axisnotick; title('real(W)'); colorbar;
subplot(p_row,p_col,1+np);np=np+1;
imagesc(real(tmp_V__),[-1,+1]); axis image; axisnotick; title('real(V)'); colorbar;
subplot(p_row,p_col,1+np);np=np+1;
imagesc(abs(tmp_V__)-abs(tmp_W__)); axis image; axisnotick; title('diff(abs)'); colorbar;
end;%if (flag_verbose>1);
%%%%;
%}

%%%%%%%%;
if nargin<1;
%%%%%%%%;
flag_verbose = 1; nf=0;
if (flag_verbose>0); disp(sprintf(' %% [testing wignerd_c]')); end;
beta=pi/6; 
%%%%;
% compare to ../dir_rangan_python/W3_??__.ascii. ;
% Note that we only compare d0W_, d1W_ and d2W_, ;
% as L_m_ and V_mm__ can be ordered and phased differently. ;
%%%%;
n_l = 32;
[d0W0_,d1W0_,d2W0_] = dwignerdda_b(n_l,beta);
[d0W_,V_lmm___,L_lm__,d1W_,d2W_] = wignerd_c(n_l,beta);
tmp_error_0 = 0; tmp_error_1 = 0; tmp_error_2 = 0;
for nl=0:n_l;
tmp_error_0 = tmp_error_0 + fnorm(d0W0_{1+nl}-d0W_{1+nl})/max(1e-12,fnorm(d0W0_{1+nl}));
tmp_error_1 = tmp_error_1 + fnorm(d1W0_{1+nl}-d1W_{1+nl})/max(1e-12,fnorm(d1W0_{1+nl}));
tmp_error_2 = tmp_error_2 + fnorm(d2W0_{1+nl}-d2W_{1+nl})/max(1e-12,fnorm(d2W0_{1+nl}));
end;%for nl=0:n_l;
if (flag_verbose>0); disp(sprintf(' %% tmp_error_0: %0.16f',tmp_error_0)); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_error_1: %0.16f',tmp_error_1)); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_error_2: %0.16f',tmp_error_2)); end;
%%%%;
dir_ascii = '../dir_rangan_python/dir_ascii';
n_d = 3;
for nd=0:n_d-1;
if nd==0; tmp_str_pre = 'd0W_'; tmp_str_pos = '__'; tmp_pre = d0W_; end;
if nd==1; tmp_str_pre = 'd1W_'; tmp_str_pos = '__'; tmp_pre = d1W_; end;
if nd==2; tmp_str_pre = 'd2W_'; tmp_str_pos = '__'; tmp_pre = d2W_; end;
tmp_error=0;
for nl=0:n_l;
tmp_str_ascii = sprintf('%s%d%s',tmp_str_pre,nl,tmp_str_pos);
tmp_pos = tmp_pre{1+nl};
fname_ascii = sprintf('%s/%s.ascii',dir_ascii,tmp_str_ascii);
if ~exist(fname_ascii,'file'); disp(sprintf(' %% Warning, %s not found',fname_ascii)); end;
if  exist(fname_ascii,'file');
if (flag_verbose>1); disp(sprintf(' %% %s found, loading',fname_ascii)); end;
fid = fopen(fname_ascii,'r');
if nl==0; tmp_ascii = reshape(cell2mat(textscan(fid,'%f')),size(tmp_pos)); end;
if nl> 0; tmp_ascii = reshape(cell2mat(textscan(fid,'(%f)')),size(tmp_pos)); end;
fclose(fid);
tmp_error = tmp_error + fnorm_disp(max(0,flag_verbose-1),tmp_str_pos,tmp_pos,'ascii',tmp_ascii,' %% <-- should be <1e-6');
end;%if  exist(fname_ascii,'file');
end;%for nl=0:n_l;
if (flag_verbose>0); disp(sprintf(' %% nd %d/%d: tmp_error: %0.16f; %%<-- should be <1e-6 per entry',nd,n_d,tmp_error)); end;
end;%for nd=0:n_d-1;
%%%%%%%%;
% Note that the discrepancy is larger than 1e-6 at nl==88. ;
%%%%%%%%;
n_l=88;
tic; W0_ = wignerd_b(n_l,beta); disp(sprintf(' %% wignerd_b : %0.2f seconds',toc));
tic; [W1_,V_lmm___,L_lm__] = wignerd_c(n_l,beta); disp(sprintf(' %% wignerd_c (not precomputation): %0.2f seconds',toc));
tic; [W1_] = wignerd_c(n_l,beta,V_lmm___,L_lm__); disp(sprintf(' %% wignerd_c (yes precomputation): %0.2f seconds',toc));
tic; W2_ = wignerd_lsq_b(n_l,beta); disp(sprintf(' %% wignerd_lsq_b: %0.2f seconds',toc));  
tmp_fnorm_numerator = 0;
tmp_fnorm_denomator = 0;
for nl=0:n_l;
tmp_fnorm_numerator = tmp_fnorm_numerator + fnorm(W1_{1+nl}-W2_{1+nl});
tmp_fnorm_denomator = tmp_fnorm_denomator + fnorm(W1_{1+nl});
if (flag_verbose>1); disp(sprintf(' %% nl %d/%d: error %0.16f',nl,n_l,fnorm(W1_{1+nl}-W2_{1+nl}))); end;
end;%for nl=0:n_l;
if (flag_verbose>0); disp(sprintf(' %% error %0.16f',tmp_fnorm_numerator/max(1e-12,tmp_fnorm_denomator))); end;
if (flag_verbose>1); 
% top plot (real) should show matrices of all ones ;
n_l=5;
W1_ = wignerd_c(n_l,beta); 
W2_ = wignerd_lsq_b(n_l,beta);
for nl = 1:n_l;
subplot(2,n_l,nl+0*n_l); imagesc(real(W2_{nl})./real(W1_{nl}),[-2,2]); title('real');
subplot(2,n_l,nl+1*n_l); imagesc(imag(W2_{nl})./imag(W1_{nl}),[-2,2]); title('imag');
if (flag_verbose>1); disp(sprintf(' %% nl %d/%d: error %0.16f',nl,n_l,fnorm(W1_{1+nl}-W2_{1+nl}))); end;
end;%for nl = 1:n_l;
end;%if (flag_verbose>1); 
%%%%%%%%;
% Now test derivatives. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% Now testing first-derivative: ')); end;
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
n_l = 88; beta = pi/6; da = pi*1e-6;
[d0W_mid_,V_lmm___,L_lm__,d1W_mid_,d2W_mid_] = wignerd_c(n_l,beta+0*da);
[d0W_pos_] = wignerd_c(n_l,beta+1*da,V_lmm___,L_lm__);
[d0W_neg_] = wignerd_c(n_l,beta-1*da,V_lmm___,L_lm__);
d1W_dif_ = cell(1+n_l,1);
tmp_fnorm_numerator = 0;
tmp_fnorm_denomator = 0;
for nl=0:n_l;
d1W_dif_{1+nl} = (d0W_pos_{1+nl} - d0W_neg_{1+nl})/max(1e-12,2*da);
tmp_fnorm_numerator = tmp_fnorm_numerator + fnorm(d1W_dif_{1+nl}-d1W_mid_{1+nl});
tmp_fnorm_denomator = tmp_fnorm_denomator + fnorm(d1W_dif_{1+nl});
if (flag_verbose>1); disp(sprintf(' %% nl %d/%d d1W_dif_{1+nl} vs d1W_mid_{1+nl}: %0.16f',nl,n_l,fnorm(d1W_dif_{1+nl}-d1W_mid_{1+nl})/max(1e-12,fnorm(d1W_dif_{1+nl})))); end;
end;%for nl=0:n_l;
if (flag_verbose>0); disp(sprintf(' %% error %0.16f',tmp_fnorm_numerator/max(1e-12,tmp_fnorm_denomator))); end;
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% Now testing second-derivative: ')); end;
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
d2W_dif_ = cell(1+n_l,1);
tmp_fnorm_numerator = 0;
tmp_fnorm_denomator = 0;
for nl=0:n_l;
d2W_dif_{1+nl} = (d0W_pos_{1+nl} - 2*d0W_mid_{1+nl} + d0W_neg_{1+nl})/max(1e-12,da^2);
tmp_fnorm_numerator = tmp_fnorm_numerator + fnorm(d2W_dif_{1+nl}-d2W_mid_{1+nl});
tmp_fnorm_denomator = tmp_fnorm_denomator + fnorm(d2W_dif_{1+nl});
if (flag_verbose>1); disp(sprintf(' %% nl %d/%d d2W_dif_{1+nl} vs d2W_mid_{1+nl}: %0.16f',nl,n_l,fnorm(d2W_dif_{1+nl}-d2W_mid_{1+nl})/max(1e-12,fnorm(d2W_dif_{1+nl})))); end;
end;%for nl=0:n_l;
if (flag_verbose>0); disp(sprintf(' %% error %0.16f',tmp_fnorm_numerator/max(1e-12,tmp_fnorm_denomator))); end;
%%%%%%%%;
disp('returning'); return;
%%%%%%%%;
end;%if nargin<1;
%%%%%%%%;

flag_verbose=0;

na=0;
if (nargin<1+na); n_l = []; end; na=na+1;
if (nargin<1+na); beta = []; end; na=na+1;
if (nargin<1+na); V_lmm___ = []; end; na=na+1;
if (nargin<1+na); L_lm__ = []; end; na=na+1;

flag_d0 = 1;
flag_d1 = (nargout>=4);
flag_d2 = (nargout>=5);

% n_l = 10; beta = pi/6; flag_d0 = 1; flag_d0=1; flag_d1 = 0; flag_d2 = 0; d0V_ = wignerd_b(n_l,beta);

l_max = n_l;
l_val_ = transpose(0:l_max);
m_val_ = transpose(0:l_max);
m_max_ = transpose(-l_max:+l_max);

flag_precomputation = isempty(V_lmm___) | isempty(V_lmm___);
if flag_precomputation;
if (flag_verbose>0); disp(sprintf(' %% flag_precomputation: %d l_max %d',flag_precomputation,l_max)); end;
V_lmm___ = cell(1+l_max,1);
L_lm__ = cell(1+l_max,1);
for l_val=1:l_max;
m_val_ = transpose([-l_val:+l_val]);
n_m_val = numel(m_val_);
X_ = sqrt( (l_val + m_val_) .* (1 + l_val - m_val_) ); 
S__ = spdiags([-X_ , +flip(X_)],[+1,-1],n_m_val,n_m_val);
[V_mm__,L__] = eigs(S__,n_m_val); L_m_ = diag(L__);
V_lmm___{1+l_val} = V_mm__;
L_lm__{1+l_val} = L_m_;
clear V_mm__ L_m_ L__ S__ X_ ;
end;%for l_val=1:l_max;
end;%if flag_precomputation;

if (flag_d0); d0W_ = cell(1+l_max,1); d0W_{1+0} = [1]; end;
if (flag_d1); d1W_ = cell(1+l_max,1); d1W_{1+0} = [0]; end;
if (flag_d2); d2W_ = cell(1+l_max,1); d2W_{1+0} = [0]; end;
%%%%%%%%;
for l_val=1:l_max;
m_val_ = transpose([-l_val:+l_val]);
n_m_val = numel(m_val_);
sgn_ = (-1).^m_val_.*(m_val_>=0) + (+1).*(m_val_< 0);
sgn__ = (sgn_*transpose(sgn_));
V_mm__ = V_lmm___{1+l_val};
L_m_ = L_lm__{1+l_val};
%%%%;
if (flag_d0);
tmp_d0W__ = V_mm__ * bsxfun(@times,reshape((+L_m_/2).^0.*exp(+L_m_*beta/2),[n_m_val,1]),ctranspose(V_mm__));
tmp_d0W__ = tmp_d0W__.*sgn__;
d0W_{1+l_val} = tmp_d0W__;
clear tmp_d0W__;
end;%if (flag_d0);
%%%%;
if (flag_d1);
tmp_d1W__ = V_mm__ * bsxfun(@times,reshape((+L_m_/2).^1.*exp(+L_m_*beta/2),[n_m_val,1]),ctranspose(V_mm__));
tmp_d1W__ = tmp_d1W__.*sgn__;
d1W_{1+l_val} = tmp_d1W__;
clear tmp_d1W__;
end;%if (flag_d1);
%%%%;
if (flag_d2);
tmp_d2W__ = V_mm__ * bsxfun(@times,reshape((+L_m_/2).^2.*exp(+L_m_*beta/2),[n_m_val,1]),ctranspose(V_mm__));
tmp_d2W__ = tmp_d2W__.*sgn__;
d2W_{1+l_val} = tmp_d2W__;
clear tmp_d2W__;
end;%if (flag_d2);
%%%%;
clear m_val n_m_val sgn_ V_mm__ L_m_ ;
end;%for l_val=1:l_max;
%%%%%%%%;


