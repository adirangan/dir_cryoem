function ...
[ ...
 d0W_ ...
,d1W_ ...
,d2W_ ...
] = ...
wignerd_c( ...
 n_l ...
,beta ...
) ;
% Generates wigner-d matrices up to n_l ;
% See Gumerov and Duraiswami 2014. ; %<-- I could not get this to work. ;
% Try Feng, Wang, Yang, Jin, 2015. ; %<-- This works. ;
% Test with: ;
%{
%%%%;
verbose = 1; nf=0;
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
if (verbose); disp(sprintf(' %% tmp_U__ vs tmp_W__: %0.16f',fnorm(tmp_U__-tmp_W__)/fnorm(tmp_U__))); end;
%%%%;
sgn_ = (-1).^m_val_.*(m_val_>=0) + (+1).*(m_val_< 0);
tmp_W__ = tmp_W__.*(sgn_*transpose(sgn_));
tmp_U__ = tmp_U__.*(sgn_*transpose(sgn_));
if (verbose); disp(sprintf(' %% tmp_V__ vs tmp_W__: %0.16f',fnorm(tmp_V__-tmp_W__)/fnorm(tmp_V__))); end;
if (verbose); disp(sprintf(' %% tmp_V__ vs tmp_U__: %0.16f',fnorm(tmp_V__-tmp_U__)/fnorm(tmp_V__))); end;
%%%%;
if (verbose>1);
figure(1+nf);nf=nf+1;clf;figmed; fig80s;
p_row = 1; p_col = 3; np=0;
subplot(p_row,p_col,1+np);np=np+1;
imagesc(real(tmp_W__),[-1,+1]); axis image; axisnotick; title('real(W)'); colorbar;
subplot(p_row,p_col,1+np);np=np+1;
imagesc(real(tmp_V__),[-1,+1]); axis image; axisnotick; title('real(V)'); colorbar;
subplot(p_row,p_col,1+np);np=np+1;
imagesc(abs(tmp_V__)-abs(tmp_W__)); axis image; axisnotick; title('diff(abs)'); colorbar;
end;%if (verbose>1);
%%%%;
%}

%%%%%%%%;
if nargin<1;
%%%%%%%%;
verbose = 1; nf=0;
if (verbose>0); disp(sprintf(' %% [testing wignerd_c]')); end;
beta=pi/6; 
%%%%%%%%;
% Note that the discrepancy is larger than 1e-6 at nl==88. ;
%%%%%%%%;
n_l=88;
tic; W1_ = wignerd_c(n_l,beta);     disp(sprintf(' %% wignerd_c    : %0.2f seconds',toc));
tic; W2_ = wignerd_lsq_b(n_l,beta); disp(sprintf(' %% wignerd_lsq_b: %0.2f seconds',toc));  
tmp_fnorm_numerator = 0;
tmp_fnorm_denomator = 0;
for nl=0:n_l;
tmp_fnorm_numerator = tmp_fnorm_numerator + fnorm(W1_{1+nl}-W2_{1+nl});
tmp_fnorm_denomator = tmp_fnorm_denomator + fnorm(W1_{1+nl});
if (verbose>1); disp(sprintf(' %% nl %d/%d: error %0.16f',nl,n_l,fnorm(W1_{1+nl}-W2_{1+nl}))); end;
end;%for nl=0:n_l;
if (verbose>0); disp(sprintf(' %% error %0.16f',tmp_fnorm_numerator/max(1e-12,tmp_fnorm_denomator))); end;
if (verbose>1); 
% top plot (real) should show matrices of all ones ;
n_l=5;
W1_ = wignerd_c(n_l,beta); 
W2_ = wignerd_lsq_b(n_l,beta);
for nl = 1:n_l;
subplot(2,n_l,nl+0*n_l); imagesc(real(W2_{nl})./real(W1_{nl}),[-2,2]); title('real');
subplot(2,n_l,nl+1*n_l); imagesc(imag(W2_{nl})./imag(W1_{nl}),[-2,2]); title('imag');
if (verbose>1); disp(sprintf(' %% nl %d/%d: error %0.16f',nl,n_l,fnorm(W1_{1+nl}-W2_{1+nl}))); end;
end;%for nl = 1:n_l;
end;%if (verbose>1); 
%%%%%%%%;
% Now test derivatives. ;
%%%%%%%%;
if (verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (verbose>0); disp(sprintf(' %% Now testing first-derivative: ')); end;
if (verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
n_l = 88; beta = pi/6; da = pi*1e-6;
[d0W_mid_,d1W_mid_,d2W_mid_] = wignerd_c(n_l,beta+0*da);
[d0W_pos_] = wignerd_c(n_l,beta+1*da);
[d0W_neg_] = wignerd_c(n_l,beta-1*da);
d1W_dif_ = cell(1+n_l,1);
tmp_fnorm_numerator = 0;
tmp_fnorm_denomator = 0;
for nl=0:n_l;
d1W_dif_{1+nl} = (d0W_pos_{1+nl} - d0W_neg_{1+nl})/max(1e-12,2*da);
tmp_fnorm_numerator = tmp_fnorm_numerator + fnorm(d1W_dif_{1+nl}-d1W_mid_{1+nl});
tmp_fnorm_denomator = tmp_fnorm_denomator + fnorm(d1W_dif_{1+nl});
if (verbose>1); disp(sprintf(' %% nl %d/%d d1W_dif_{1+nl} vs d1W_mid_{1+nl}: %0.16f',nl,n_l,fnorm(d1W_dif_{1+nl}-d1W_mid_{1+nl})/max(1e-12,fnorm(d1W_dif_{1+nl})))); end;
end;%for nl=0:n_l;
if (verbose>0); disp(sprintf(' %% error %0.16f',tmp_fnorm_numerator/max(1e-12,tmp_fnorm_denomator))); end;
if (verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (verbose>0); disp(sprintf(' %% Now testing second-derivative: ')); end;
if (verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
d2W_dif_ = cell(1+n_l,1);
tmp_fnorm_numerator = 0;
tmp_fnorm_denomator = 0;
for nl=0:n_l;
d2W_dif_{1+nl} = (d0W_pos_{1+nl} - 2*d0W_mid_{1+nl} + d0W_neg_{1+nl})/max(1e-12,da^2);
tmp_fnorm_numerator = tmp_fnorm_numerator + fnorm(d2W_dif_{1+nl}-d2W_mid_{1+nl});
tmp_fnorm_denomator = tmp_fnorm_denomator + fnorm(d2W_dif_{1+nl});
if (verbose>1); disp(sprintf(' %% nl %d/%d d2W_dif_{1+nl} vs d2W_mid_{1+nl}: %0.16f',nl,n_l,fnorm(d2W_dif_{1+nl}-d2W_mid_{1+nl})/max(1e-12,fnorm(d2W_dif_{1+nl})))); end;
end;%for nl=0:n_l;
if (verbose>0); disp(sprintf(' %% error %0.16f',tmp_fnorm_numerator/max(1e-12,tmp_fnorm_denomator))); end;
%%%%%%%%;
disp('returning'); return;
%%%%%%%%;
end;%if nargin<1;
%%%%%%%%;

na=0;
if (nargin<1+na); n_l = []; end; na=na+1;
if (nargin<1+na); beta = []; end; na=na+1;

flag_d0 = 1;
flag_d1 = (nargout>=2);
flag_d2 = (nargout>=3);

% n_l = 10; beta = pi/6; flag_d0 = 1; flag_d0=1; flag_d1 = 0; flag_d2 = 0; d0V_ = wignerd_b(n_l,beta);

l_max = n_l;
l_val_ = transpose(0:l_max);
m_val_ = transpose(0:l_max);
m_max_ = transpose(-l_max:+l_max);

if (flag_d0); d0W_ = cell(1+l_max,1); d0W_{1+0} = [1]; end;
if (flag_d1); d1W_ = cell(1+l_max,1); d1W_{1+0} = [0]; end;
if (flag_d2); d2W_ = cell(1+l_max,1); d2W_{1+0} = [0]; end;
%%%%%%%%;
for l_val=1:l_max;
m_val_ = transpose([-l_val:+l_val]);
n_m_val = numel(m_val_);
sgn_ = (-1).^m_val_.*(m_val_>=0) + (+1).*(m_val_< 0);
X_ = sqrt( (l_val + m_val_) .* (1 + l_val - m_val_) ); 
S__ = spdiags([-X_ , +flip(X_)],[+1,-1],n_m_val,n_m_val);
[V__,L__] = eigs(S__,n_m_val); L_ = diag(L__);
%%%%;
if (flag_d0);
tmp_d0W__ = V__ * bsxfun(@times,reshape((+L_/2).^0.*exp(+L_*beta/2),[n_m_val,1]),ctranspose(V__));
tmp_d0W__ = tmp_d0W__.*(sgn_*transpose(sgn_));
d0W_{1+l_val} = tmp_d0W__;
clear tmp_d0W__;
end;%if (flag_d0);
%%%%;
if (flag_d1);
tmp_d1W__ = V__ * bsxfun(@times,reshape((+L_/2).^1.*exp(+L_*beta/2),[n_m_val,1]),ctranspose(V__));
tmp_d1W__ = tmp_d1W__.*(sgn_*transpose(sgn_));
d1W_{1+l_val} = tmp_d1W__;
clear tmp_d1W__;
end;%if (flag_d1);
%%%%;
if (flag_d2);
tmp_d2W__ = V__ * bsxfun(@times,reshape((+L_/2).^2.*exp(+L_*beta/2),[n_m_val,1]),ctranspose(V__));
tmp_d2W__ = tmp_d2W__.*(sgn_*transpose(sgn_));
d2W_{1+l_val} = tmp_d2W__;
clear tmp_d2W__;
end;%if (flag_d2);
%%%%;
clear m_val n_m_val sgn_ X_ S__ V__ L__ L_ ;
end;%for l_val=1:l_max;
%%%%%%%%;


