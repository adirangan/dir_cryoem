function ...
[ ...
 d0W_ ...
,d1W_ ...
,d2W_ ...
] = ...
wignerd_c( ...
 n_l ...
,beta ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_ ...
,sqrt_rat3__ ...
,sqrt_rat4__ ...
) ;
% Generates wigner-d matrices up to n_l ;
% See Gumerov and Duraiswami 2014. ;
% Test with: ;
%{

  beta=pi/6; 
  %%%%%%%%;
  % Note that the discrepancy is larger than 1e-6 at nl==88. ;
  %%%%%%%%;
  n_l=88;
  tic; W1_ = wignerd_c(n_l,beta);     disp(sprintf(' %% wignerd    : %0.2f seconds',toc));
  tic; W2_ = wignerd_lsq_b(n_l,beta); disp(sprintf(' %% wignerd_lsq_b: %0.2f seconds',toc));  
  for nl=0:n_l;
  disp(sprintf(' %% nl %d/%d: error %0.16f',nl,n_l,fnorm(W1_{1+nl}-W2_{1+nl})));
  end;%for nl=0:n_l;
  % top plot (real) should show matrices of all ones ;
  n_l=5;
  W1_ = wignerd_c(n_l,beta); 
  W2_ = wignerd_lsq_b(n_l,beta);
  for nl = 1:n_l;
  subplot(2,n_l,nl+0*n_l); imagesc(real(W2_{nl})./real(W1_{nl}),[-2,2]); title('real');
  subplot(2,n_l,nl+1*n_l); imagesc(imag(W2_{nl})./imag(W1_{nl}),[-2,2]); title('imag');
  disp(sprintf(' %% nl %d/%d: error %0.16f',nl,n_l,fnorm(W1_{1+nl}-W2_{1+nl})));
  end;%for nl = 1:n_l;

  %}

na=0;
if (nargin<1+na); n_l = []; end; na=na+1;
if (nargin<1+na); beta = []; end; na=na+1;
if (nargin<1+na); sqrt_2lp1_ = []; end; na=na+1;
if (nargin<1+na); sqrt_2mp1_ = []; end; na=na+1;
if (nargin<1+na); sqrt_rat0_ = []; end; na=na+1;
if (nargin<1+na); sqrt_rat3__ = []; end; na=na+1;
if (nargin<1+na); sqrt_rat4__ = []; end; na=na+1;

flag_d0 = 1;
flag_d1 = (nargin>=2);
flag_d2 = (nargin>=3);

% n_l = 10; beta = pi/6; flag_d0 = 1; flag_d0=1; flag_d1 = 0; flag_d2 = 0; sqrt_2lp1_ = []; sqrt_2mp1_ = []; sqrt_rat0_ = []; sqrt_rat3__ = []; sqrt_rat4__ = []; d0V_ = wignerd_b(n_l,beta);

l_max = n_l+1;
l_val_ = transpose(0:l_max);
m_val_ = transpose(0:l_max);
m_max_ = transpose(-l_max:+l_max);
[l_lm__,m_lm__] = ndgrid(l_val_,m_max_);

flag_precomputation = ( isempty(sqrt_2lp1_) | isempty(sqrt_2mp1_) | isempty(sqrt_rat0_) | isempty(sqrt_rat3__) | isempty(sqrt_rat4__) );
if flag_precomputation;
l_val_ = transpose(0:l_max);
m_val_ = transpose(0:l_max);
base_2lm1_ = 2.0d0*l_val_-1.0d0;
sqrt_2lp1_ = sqrt(2.0d0*l_val_+1.0d0);
sqrt_2mp1_ = sqrt(2.0d0*m_val_+1.0d0);
sqrt_rat0_ = sqrt((2.0d0*m_val_-1.0d0)./(2.0d0*m_val_));
sqrt_rat1__ = ...
  sqrt( ...
	(bsxfun(@plus,reshape(l_val_,[1+l_max,1]),reshape(m_val_,[1,1+l_max]))-1.0d0) ...
	.* ...
	(bsxfun(@minus,reshape(l_val_,[1+l_max,1]),reshape(m_val_,[1,1+l_max]))-1.0d0) ...
	) ...
  ;
sqrt_rat2__ = ...
  sqrt( ...
	(bsxfun(@minus,reshape(l_val_,[1+l_max,1]),reshape(m_val_,[1,1+l_max]))+0.0d0) ...
	.* ...
	(bsxfun(@plus,reshape(l_val_,[1+l_max,1]),reshape(m_val_,[1,1+l_max]))+0.0d0) ...
	);
sqrt_rat3__ = sqrt_rat1__./sqrt_rat2__;
sqrt_rat4__ = bsxfun(@rdivide,reshape(base_2lm1_,[1+l_max,1]),sqrt_rat2__);
end;%if flag_precomputation;

%%%%%%%%;
% restriction. ;
%%%%%%%%;
r_lm__ = (abs(m_lm__)<=l_lm__);
index_r_ = efind(~r_lm__);

%%%%%%%%;
% \epsilon_{m} = +1 if m< 0, and (-1)^(m) if m>=0. ;
%%%%%%%%;
epsilon_lm__ = (+1).*(m_lm__< 0) + (-1).^(m_lm__).*(m_lm__>=0) ;
epsilon_lm__(1+index_r_) = 0;

%%%%%%%%;
% a_{l}^{m} = \sqrt{\frac{(l+1+m)(l+1-m)}{(2l+1)(2l+3)}} ;
%%%%%%%%;
a_lm__ = sqrt( ( (l_lm__ + 1 + m_lm__).*(l_lm__ + 1 - m_lm__) ) ./ ( (2*l_lm__ + 1).*(2*l_lm__ + 3) ) ) ;
a_lm__(1+index_r_) = 0;

%%%%%%%%;
% \sgn(m) = (+1)*(m>=0) + (-1)*(m< 0);
%%%%%%%%;
sgn_lm__ = (+1)*(m_lm__>=0) + (-1)*(m_lm__< 0) ;
sgn_lm__(1+index_r_) = 0;

%%%%%%%%;
% b_{l}^{m} = \sgn(m) \cdot \sqrt{\frac{(l-m-1)(l-m-0)}{(2l-1)(2l+1)}} ;
%%%%%%%%;
b_lm__ = sgn_lm__ .* sqrt( ( (l_lm__ - m_lm__ - 1).*(l_lm__ - m_lm__ - 0) ) ./ ( (2*l_lm__ - 1).*(2*l_lm__ + 1) ) ) ;
b_lm__(1+index_r_) = 0;

%%%%%%%%;
% c_{l}^{m} = 0.5 \cdot (-1)^(m) \cdot \sgn(m) \cdot \sqrt{(l-m)(l+m+1)} ;
%%%%%%%%;
c_lm__ = 0.5 .* (-1).^m_lm__ .* sgn_lm__ .* sqrt( (l_lm__ - m_lm__).*(l_lm__ + m_lm__ + 1) ) ;
c_lm__(1+index_r_) = 0;

%%%%%%%%;
% d_{l}^{m} = 0.5 \cdot \sgn(m) \cdot \sqrt{(l-m)(l+m+1)} ;
%%%%%%%%;
d_lm__ = 0.5 .* sgn_lm__ .* sqrt( (l_lm__ - m_lm__).*(l_lm__ + m_lm__ + 1) ) ;
d_lm__(1+index_r_) = 0;

%%%%%%%%;
% y_lm__(1+l_val,1+abs(m_val)) = scaled legendre-polynomial of degree l and order m. ;
% i.e.,: (l,m) --> \sqrt{\frac{l-|m|}{l+|m|}} \cdot P^{|m|}_{l}(\cos(\beta)) ;
%%%%%%%%;
x_ = cos(beta);
[ ...
 d0y_jlm___ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_ ...
,sqrt_rat3__ ...
,sqrt_rat4__ ...
,d1y_jlm___ ...
,d2y_jlm___ ...
] = ...
ylgndr_1( ...
 l_max ...
,x_ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_ ...
,sqrt_rat3__ ...
,sqrt_rat4__ ...
);
if flag_d0; d0y_lm__ = reshape(d0y_jlm___,[1+l_max,1+l_max]); end;
if flag_d1; d1y_lm__ = reshape(d1y_jlm___,[1+l_max,1+l_max]); end;
if flag_d2; d2y_lm__ = reshape(d2y_jlm___,[1+l_max,1+l_max]); end;

if flag_d0; d0W_ = cell(1+l_max,1); end;
if flag_d1; d1W_ = cell(1+l_max,1); end;
if flag_d2; d2W_ = cell(1+l_max,1); end;

%%%%%%%%;
% H^{0,0}_{0} = 1. ;
%%%%%%%%;
for l_val=0;
if flag_d0; d0W_{1} = [1]; end;
if flag_d1; d1W_{1} = [1]; end;
if flag_d2; d2W_{1} = [1]; end;
end;%for l_val=0;
%%%%%%%%;
% H^{m,0}_{l}(\beta) = H^{0,m}_{l}(\beta) = (-1).^m * \sqrt(\factorial(l-|m|)/\factorial(l+|m|)) * P^{|m|}_{l}(\cos(\beta)). ;
%%%%%%%%;
for l_val=1:l_max;
if flag_d0;
tmp_d0y_m_ = (-1).^(m_lm__(1+l_val,1+[0:+1:+l_val])) .* d0y_lm__(1+l_val,1+[0:+1:+l_val]);
d0W_{1+l_val} = zeros(1+2*l_val,1+2*l_val);
d0W_{1+l_val}(1+l_val+0,1+l_val+[0:+1:+l_val]) = tmp_d0y_m_(1+[0:+1:+l_val]);
%d0W_{1+l_val}(1+l_val+0,1+l_val+[0:-1:-l_val]) = tmp_d0y_m_(1+[0:+1:+l_val]);
%d0W_{1+l_val}(1+l_val+[0:+1:+l_val],1+l_val+0) = tmp_d0y_m_(1+[0:+1:+l_val]);
%d0W_{1+l_val}(1+l_val+[0:-1:-l_val],1+l_val+0) = tmp_d0y_m_(1+[0:+1:+l_val]);
end;%if flag_d0;
if flag_d1;
tmp_d1y_m_ = (-1).^(m_lm__(1+l_val,1+[0:+1:+l_val])) .* d1y_lm__(1+l_val,1+[0:+1:+l_val]);
d1W_{1+l_val} = zeros(1+2*l_val,1+2*l_val);
d1W_{1+l_val}(1+l_val+0,1+l_val+[0:+1:+l_val]) = tmp_d1y_m_(1+[0:+1:+l_val]);
%d1W_{1+l_val}(1+l_val+0,1+l_val+[0:-1:-l_val]) = tmp_d1y_m_(1+[0:+1:+l_val]);
%d1W_{1+l_val}(1+l_val+[0:+1:+l_val],1+l_val+0) = tmp_d1y_m_(1+[0:+1:+l_val]);
%d1W_{1+l_val}(1+l_val+[0:-1:-l_val],1+l_val+0) = tmp_d1y_m_(1+[0:+1:+l_val]);
end;%if flag_d1;
if flag_d2;
tmp_d2y_m_ = (-1).^(m_lm__(1+l_val,1+[0:+1:+l_val])) .* d2y_lm__(1+l_val,1+[0:+1:+l_val]);
d2W_{1+l_val} = zeros(1+2*l_val,1+2*l_val);
d2W_{1+l_val}(1+l_val+0,1+l_val+[0:+1:+l_val]) = tmp_d2y_m_(1+[0:+1:+l_val]);
d2W_{1+l_val}(1+l_val+0,1+l_val+[0:-1:-l_val]) = tmp_d2y_m_(1+[0:+1:+l_val]);
d2W_{1+l_val}(1+l_val+[0:+1:+l_val],1+l_val+0) = tmp_d2y_m_(1+[0:+1:+l_val]);
d2W_{1+l_val}(1+l_val+[0:-1:-l_val],1+l_val+0) = tmp_d2y_m_(1+[0:+1:+l_val]);
end;%if flag_d2;
end;%for l_val=1:l_max;

%%%%%%%%;
% b^{0}_{l+1} H^{1,m}_{l} = 0.5 * b^{-m-1}_{l+1}(1-\cos\beta)H^{0,m+1}_{l+1} - 0.5*b^{m-1}_{l+1}(1+\cos\beta)H^{0,m-1}_{l+1} - a^{m}_{l}\sin\beta H^{0,m}_{l+1} ;
%%%%%%%%;
for l_val=1:l_max;
lp0 = l_val+0;
lp1 = l_val+1;
for m0_val=1;
for m1_val=1:l_val;
if flag_d0;
tmp_0 = 0;
if (lp1<=l_max) & (abs(-m1_val-1)<=lp1) & (abs(+m1_val+1)<=lp1);
tmp_0 = b_lm__(1+lp1,1+l_max+(-m1_val-1))*(1-cos(beta))*d0W_{1+lp1}(1+lp1+0,1+lp1+(+m1_val+1));
end;%if;
tmp_1 = 0;
if (lp1<=l_max) & (abs(+m1_val-1)<=lp1) & (abs(+m1_val-1)<=lp1);
tmp_1 = b_lm__(1+lp1,1+l_max+(+m1_val-1))*(1+cos(beta))*d0W_{1+lp1}(1+lp1+0,1+lp1+(+m1_val-1));
end;%if;
tmp_2 = 0;
if (lp1<=l_max) & (abs(+m1_val+0)<=lp0) & (abs(+m1_val+0)<=lp1);
tmp_2 = a_lm__(1+lp0,1+l_max+(+m1_val+0))*(0+sin(beta))*d0W_{1+lp1}(1+lp1+0,1+lp1+(+m1_val+0));
end;%if;
if (lp1<=l_max) & (abs(+m0_val+0)<=lp0) & (abs(+m1_val+0)<=lp0);
d0W_{1+lp0}(1+lp0+m0_val,1+lp0+m1_val) = ( + 0.5*tmp_0 - 0.5*tmp_1 - 1.0*tmp_2 )/b_lm__(1+lp1,1+l_max+0) ;
end;%if;
end;%if flag_d0;
end;%for m1_val=1:l_val;
end;%for m0_val=1;
end;%for l_val=1:l_max;

%%%%%%%%;
% d^{m0}_{l}H^{m0+1,m1}_{l} = d^{m0-1}_{l}H^{m0-1,m1}_{l} - d^{m1-1}_{l}H^{m0,m1-1}_{l} + d^{m1}_{l}H^{m0,m1+1}_{l} ;
%%%%%%%%;
for l_val=0:l_max;
lp0 = l_val+0;
lp1 = l_val+1;
for m0_val=+1:+1:+l_val-1;
for m1_val=+m0_val:+1:+l_val;
if flag_d0;
tmp_0 = 0; if (m1_val<=lp0); tmp_0 = d_lm__(1+lp0,1+l_max+(+m0_val-1))*d0W_{1+lp0}(1+lp0+(+m0_val-1),1+lp0+(+m1_val+0)); end;
tmp_1 = 0; if (m1_val<=lp0); tmp_1 = d_lm__(1+lp0,1+l_max+(+m1_val-1))*d0W_{1+lp0}(1+lp0+(+m0_val+0),1+lp0+(+m1_val-1)); end;
tmp_2 = 0; if (m1_val< lp0); tmp_2 = d_lm__(1+lp0,1+l_max+(+m1_val+0))*d0W_{1+lp0}(1+lp0+(+m0_val+0),1+lp0+(+m1_val+1)); end;
d0W_{1+lp0}(1+lp0+(+m0_val+1),1+lp0+(+m1_val+0)) = ( + tmp_0 - tmp_1 + tmp_2 )/d_lm__(1+lp0,1+l_max+(+m0_val+0));
end;%if flag_d0;
end;%for m1_val=+m0_val:+1:+l_val;
end;%for m0_val=+1:+1:+l_val-1;
end;%for l_val=0:l_max;

%%%%%%%%;
% d^{m0-1}_{l}H^{m0-1,m1}_{l} = d^{m0}_{l}H^{m0+1,m1}_{l} - d^{m1-1}_{l}H^{m0,m1-1}_{l} + d^{m1}_{l}H^{m0,m1+1}_{l} ;
%%%%%%%%;
for l_val=0:l_max;
lp0 = l_val+0;
lp1 = l_val+1;
for m0_val=-1:-1:-l_val+1;
for m1_val=-m0_val:+1:+l_val;
if flag_d0;
tmp_0 = 0; if (m1_val<=lp0); tmp_0 = d_lm__(1+lp0,1+l_max+(+m0_val+0))*d0W_{1+lp0}(1+lp0+(+m0_val+1),1+lp0+(+m1_val+0)); end;
tmp_1 = 0; if (m1_val<=lp0); tmp_1 = d_lm__(1+lp0,1+l_max+(+m1_val-1))*d0W_{1+lp0}(1+lp0+(+m0_val+0),1+lp0+(+m1_val-1)); end;
tmp_2 = 0; if (m1_val< lp0); tmp_2 = d_lm__(1+lp0,1+l_max+(+m1_val+0))*d0W_{1+lp0}(1+lp0+(+m0_val+0),1+lp0+(+m1_val+1)); end;
d0W_{1+lp0}(1+lp0+(+m0_val-1),1+lp0+(+m1_val+0)) = ( + tmp_0 - tmp_1 + tmp_2 )/d_lm__(1+lp0,1+l_max+(+m0_val-1));
end;%if flag_d0;
end;%for m1_val=-m0_val:+1:+l_val;
end;%for m0_val=-0:-1:-l_val+1;
end;%for l_val=0:l_max;

%%%%%%%%;
% H^{+m0,+m1}_{l} = H^{+m1,+m0}_{l} ;
% H^{+m0,+m1}_{l} = H^{-m1,-m0}_{l} ;
%%%%%%%%;
for l_val=0:l_max;
lp0 = l_val+0;
lp1 = l_val+1;
for m0_val=-l_val+1:+l_val-1;
for m1_val=-l_val:-1;
if flag_d0;
d0W_{1+lp0}(1+lp0+(+m0_val+0),1+lp0+(+m1_val+0)) = d0W_{1+lp0}(1+lp0+(+m0_val+0),1+lp0+(-m1_val+0));
end;%if flag_d0;
end;%for m1_val=-l_val:-1;
end;%for m0_val=-l_val+1:+l_val-1;
end;%for l_val=0:l_max;
%%%%%%%%;
for l_val=0:l_max;
for m0_val=[-l_val,+l_val];
for m1_val=-l_val:+l_val;
if flag_d0;
d0W_{1+lp0}(1+lp0+(+m0_val+0),1+lp0+(+m1_val+0)) = d0W_{1+lp0}(1+lp0+(+m1_val+0),1+lp0+(-m0_val+0));
end;%if flag_d0;
end;%for m1_val=-l_val:+l_val;
end;%for m0_val=[-l_val,+l_val];
end;%for l_val=0:l_max;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Fix condon-shortley-phase ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nl=1:n_l;
m_ = -nl:+nl;
s=(-1).^((m_<0).*m_); % needed to preserve condon-shortley phase. ;
S = transpose(s)*s;
if flag_d0; d0W_{1+nl} = d0W_{1+nl}.*S; end;
end;%for nl=1:n_l;

