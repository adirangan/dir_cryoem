function ...
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

% Note that this returns: ;
% d0y_jlm___(1+nx,1+l_val,1+m_val) ; 
% such that ;
% y__(1+nx,1+l_val,1+m_val) == Y__(1+l_val,1+m_val)(x_(1+nx)), where: ;
% Y__ is constructed via: ;
% m_val_=transpose(0:l_val);
% tmp_a1 = ((1+2*l_val)/(4*pi));
% tmp_a2_ = factorial(l_val-abs(m_val_))./factorial(l_val+abs(m_val_));
% tmp_a3_ = sqrt(tmp_a1*tmp_a2_);
% Y__(1+l_val,:) = legendre(l_val,x,'unnorm').*tmp_a3_*sqrt(4*pi);

if (nargin<2);
verbose=2; nf=0;
if (verbose>0); disp(sprintf(' %% testing ylgndr_1')); end;
l_max = 4; n_x = 25;
x_ = transpose(2*rand(n_x,1)-1);
tmp_y_fortran___ = ylgndr_1(l_max,x_);
tmp_Y_matlab___ = zeros(n_x,1+l_max,1+l_max);
for l_val=0:l_max;
m_val_ = transpose(0:l_val);
tmp_a1 = ((1+2*l_val)/(4*pi));
tmp_a2_ = exp(0.5*lfactorial(l_val-abs(m_val_)) - 0.5*lfactorial(l_val+abs(m_val_)));
tmp_a3_ = sqrt(tmp_a1)*tmp_a2_;
for nx=0:n_x-1;
tmp_Y_matlab___(1+nx,1+l_val,1:(1+l_val)) = legendre(l_val,x_(1+nx),'unnorm').*tmp_a3_;
end;%for nx=0:n_x-1;
end;%for l_val=0:l_max;
if (verbose>0); disp(sprintf(' %% l_max %d: tmp_y_fortran___ vs sqrt(4*pi)*tmp_Y_matlab___: %0.16f',l_max,fnorm(tmp_y_fortran___-sqrt(4*pi)*tmp_Y_matlab___)/fnorm(tmp_y_fortran___))); end;
%%%%%%%%;
% testing first-derivative. ;
%%%%%%%%;
if (verbose>0); disp(sprintf(' %% testing first-derivative')); end;
l_max = 4; n_x = 25; x_ = transpose(2*rand(n_x,1)-1); x_ = cos(acos(x_)); da = pi*1e-4; dx_ = periodize(x_-cos(acos(x_)+da),-1,+1); plot(dx_,'.');
[ ...
 d0y_mid_jlm___ ...
, ~ ...
, ~ ...
, ~ ...
, ~ ...
, ~ ...
,d1y_mid_jlm___ ...
,d2y_mid_jlm___ ...
] = ...
ylgndr_1(l_max,x_);
d0y_pos_jlm___ = ylgndr_1(l_max,x_+dx_);
d0y_neg_jlm___ = ylgndr_1(l_max,x_-dx_);
d1y_dif_jlm___ = bsxfun(@rdivide,d0y_pos_jlm___-d0y_neg_jlm___,reshape(2*dx_,[n_x,1,1]));
if (verbose>0); disp(sprintf(' %% d1y_dif_jlm___ vs d1y_mid_jlm___: %0.16f',fnorm(d1y_dif_jlm___ - d1y_mid_jlm___)/fnorm(d1y_dif_jlm___))); end;
if (verbose>0);
figure(1+nf);nf=nf+1;clf;figsml;
plot(d1y_dif_jlm___(:),d1y_mid_jlm___(:),'.'); axis equal; grid on;
end;%if (verbose>0);
%%%%%%%%;
% testing second-derivative. ;
%%%%%%%%;
if (verbose>0); disp(sprintf(' %% testing second-derivative')); end;
d2y_dif_jlm___ = bsxfun(@rdivide,d0y_pos_jlm___ - 2*d0y_mid_jlm___ + d0y_neg_jlm___,reshape(dx_.^2,[n_x,1,1]));
if (verbose>0); disp(sprintf(' %% d2y_dif_jlm___ vs d2y_mid_jlm___: %0.16f',fnorm(d2y_dif_jlm___ - d2y_mid_jlm___)/fnorm(d2y_dif_jlm___))); end;
if (verbose>0);
figure(1+nf);nf=nf+1;clf;figsml;
plot(d2y_dif_jlm___(:),d2y_mid_jlm___(:),'.'); axis equal; grid on;
end;%if (verbose>0);

return;


%%%%%%%%;
% testing stability. ;
%%%%%%%%;
if (verbose>0); disp(sprintf(' %% testing stability')); end;
l_max = 96*2.0; n_x = 1024*2;
x_ = rand(n_x,1);
%%%%;
tmp_t = tic;
[ ...
 tmp_y_fortran___ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_ ...
,sqrt_rat3__ ...
,sqrt_rat4__ ...
] = ...
ylgndr_1( ...
 l_max ...
,x_ ...
);
tmp_t = toc(tmp_t); if (verbose>0); disp(sprintf(' %% y_fortran___ (not precomputation) : %0.2fs',tmp_t)); end;
%%%%;
tmp_t = tic;
[ ...
 bkp_y_fortran___ ...
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
tmp_t = toc(tmp_t); if (verbose>0); disp(sprintf(' %% y_fortran___ (yes precomputation) : %0.2fs',tmp_t)); end;
assert(fnorm(tmp_y_fortran___-bkp_y_fortran___)<1e-12);
%%%%;
tmp_t = tic;
tmp_Y_matlab___ = zeros(n_x,1+l_max,1+l_max);
for l_val=0:l_max;
m_val_ = transpose(0:l_val);
tmp_a1 = ((1+2*l_val)/(4*pi));
tmp_a2_ = exp(0.5*lfactorial(l_val-abs(m_val_)) - 0.5*lfactorial(l_val+abs(m_val_)));
tmp_a3_ = sqrt(tmp_a1)*tmp_a2_;
tmp_Y_matlab___(:,1+l_val,1:(1+l_val)) = reshape(transpose(legendre(l_val,x_,'unnorm').*tmp_a3_),[n_x,1,1+l_val]);
end;%for l_val=0:l_max;
tmp_t = toc(tmp_t); if (verbose>0); disp(sprintf(' %% Y_matlab___: %0.2fs',tmp_t)); end;
%%%%;
if (verbose>0); disp(sprintf(' %% l_max %d: tmp_y_fortran___ vs sqrt(4*pi)*tmp_Y_matlab___: %0.16f',l_max,fnorm(tmp_y_fortran___-sqrt(4*pi)*tmp_Y_matlab___)/fnorm(tmp_y_fortran___))); end;
if (verbose>1);
figure(1+nf);nf=nf+1;figsml;fig81s;
llim_ = [-16,0];
fontsize_use = 12;
imagesc(squeeze(mean(log10(abs(tmp_y_fortran___-sqrt(4*pi)*tmp_Y_matlab___)),1)),llim_);
xlabel('m_val','Interpreter','none'); xtickangle(90);
set(gca,'XTick',1:16:1+l_max,'XTickLabel',0:16:l_max);
ylabel('l_val','Interpreter','none');
set(gca,'YTick',1:16:1+l_max,'YTickLabel',0:16:l_max);
title('error=log10(abs(fortran-matlab))','Interpreter','none');
set(gca,'FontSize',fontsize_use);
c_ = colorbar;
end;%if (verbose>1);
disp('returning'); return;
end;%if (nargin<2);

na=0;
if (nargin<1+na); l_max = []; end; na=na+1;
if (nargin<1+na); x_ = []; end; na=na+1;
if (nargin<1+na); sqrt_2lp1_ = []; end; na=na+1;
if (nargin<1+na); sqrt_2mp1_ = []; end; na=na+1;
if (nargin<1+na); sqrt_rat0_ = []; end; na=na+1;
if (nargin<1+na); sqrt_rat3__ = []; end; na=na+1;
if (nargin<1+na); sqrt_rat4__ = []; end; na=na+1;

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

flag_d0 = 1;
flag_d1 = nargout>=7;
flag_d2 = nargout>=8;

x_ = x_(:);
n_x = numel(x_);
d0u_ = -sqrt((1-x_).*(1+x_));
d1u_ = x_./sqrt((1-x_).*(1+x_));
d2u_ = -1./d0u_ - x_.^2./d0u_.^3;
if flag_d0; d0y_jlm___ = zeros(n_x,l_max+1,l_max+1); end;
if flag_d1; d1y_jlm___ = zeros(n_x,l_max+1,l_max+1); end;
if flag_d2; d2y_jlm___ = zeros(n_x,l_max+1,l_max+1); end;
if flag_d0; d0y_jlm___(:,1+0,1+0) = 1.0d0; end;
if flag_d1; d1y_jlm___(:,1+0,1+0) = 0.0d0; end;
if flag_d2; d2y_jlm___(:,1+0,1+0) = 0.0d0; end;
for m_val=0:l_max;
%if (m_val>0); d0y_jlm___(:,1+m_val+0,1+m_val+0) = d0y_jlm___(:,1+m_val-1,1+m_val-1).*d0u_*sqrt((2.0d0*m_val-1.0d0)/(2.0d0*m_val)); end;
if (m_val>0); 
if flag_d0; d0y_jlm___(:,1+m_val+0,1+m_val+0) = d0y_jlm___(:,1+m_val-1,1+m_val-1).*d0u_*sqrt_rat0_(1+m_val); end;
if flag_d1; d1y_jlm___(:,1+m_val+0,1+m_val+0) = d0y_jlm___(:,1+m_val-0,1+m_val-0).*(-m_val).*x_./d0u_.^2; end;
%if flag_d2; d2y_jlm___(:,1+m_val+0,1+m_val+0) = d2y_jlm___(:,1+m_val-1,1+m_val-1).*d0u_*sqrt_rat0_(1+m_val) + 2*d1y_jlm___(:,1+m_val-1,1+m_val-1).*d1u_*sqrt_rat0_(1+m_val) + d0y_jlm___(:,1+m_val-1,1+m_val-1).*d2u_*sqrt_rat0_(1+m_val) ; end;
if flag_d2; d2y_jlm___(:,1+m_val+0,1+m_val+0) = d1y_jlm___(:,1+m_val-0,1+m_val-0).*(-m_val).*x_./d0u_.^2 + d0y_jlm___(:,1+m_val-0,1+m_val-0).*(-m_val).*(1 + x_.^2)./d0u_.^4; end;
end;%if (m_val>0); 
%if (m_val<l_max); d0y_jlm___(:,1+m_val+1,1+m_val+0) = x_.*d0y_jlm___(:,1+m_val+0,1+m_val+0)*sqrt(2.0d0*m_val+1.0d0); end;
if (m_val<l_max);
if flag_d0; d0y_jlm___(:,1+m_val+1,1+m_val+0) = x_.*d0y_jlm___(:,1+m_val+0,1+m_val+0)*sqrt_2mp1_(1+m_val); end;
if flag_d1; d1y_jlm___(:,1+m_val+1,1+m_val+0) = (d0y_jlm___(:,1+m_val+0,1+m_val+0) + x_.*d1y_jlm___(:,1+m_val+0,1+m_val+0))*sqrt_2mp1_(1+m_val); end;
if flag_d2; d2y_jlm___(:,1+m_val+1,1+m_val+0) = (2*d1y_jlm___(:,1+m_val+0,1+m_val+0) + x_.*d2y_jlm___(:,1+m_val+0,1+m_val+0))*sqrt_2mp1_(1+m_val); end;
end;%if (m_val<l_max);
for l_val=m_val+2:l_max;
%d0y_jlm___(:,1+l_val,1+m_val+0) = ((2.0d0*l_val-1.0d0)*x_.*d0y_jlm___(:,1+l_val-1,1+m_val+0) - sqrt((l_val+m_val-1.0d0)*(l_val-m_val-1.0d0))*d0y_jlm___(:,1+l_val-2,1+m_val+0))/sqrt((l_val-m_val+0.0d0)*(l_val+m_val+0.0d0));
if flag_d0; d0y_jlm___(:,1+l_val,1+m_val+0) = sqrt_rat4__(1+l_val,1+m_val)*x_.*d0y_jlm___(:,1+l_val-1,1+m_val+0) - sqrt_rat3__(1+l_val,1+m_val)*d0y_jlm___(:,1+l_val-2,1+m_val+0); end;
if flag_d1; d1y_jlm___(:,1+l_val,1+m_val+0) = sqrt_rat4__(1+l_val,1+m_val)*(d0y_jlm___(:,1+l_val-1,1+m_val+0) + x_.*d1y_jlm___(:,1+l_val-1,1+m_val+0)) - sqrt_rat3__(1+l_val,1+m_val)*d1y_jlm___(:,1+l_val-2,1+m_val+0); end;
if flag_d2; d2y_jlm___(:,1+l_val,1+m_val+0) = sqrt_rat4__(1+l_val,1+m_val)*(2*d1y_jlm___(:,1+l_val-1,1+m_val+0) + x_.*d2y_jlm___(:,1+l_val-1,1+m_val+0)) - sqrt_rat3__(1+l_val,1+m_val)*d2y_jlm___(:,1+l_val-2,1+m_val+0); end;
end;%for l_val=m_val+2:l_max;
end;%for m_val=0:l_max;

%%%%;
% for l_val=0:l_max;
% for m_val=0:l_val;
% d0y_jlm___(:,1+l_val,1+m_val+0) = d0y_jlm___(:,1+l_val,1+m_val+0)*sqrt(2.0d0*l_val+1.0d0);
% end;%for m_val=0:l_val;
% end;%for l_val=0:l_max;
%%%%;
if flag_d0; d0y_jlm___ = bsxfun(@times,d0y_jlm___,reshape(sqrt_2lp1_,[1,1+l_max,1])); end;
if flag_d1; d1y_jlm___ = bsxfun(@times,d1y_jlm___,reshape(sqrt_2lp1_,[1,1+l_max,1])); end;
if flag_d2; d2y_jlm___ = bsxfun(@times,d2y_jlm___,reshape(sqrt_2lp1_,[1,1+l_max,1])); end;
