function ...
[ ...
 y_jlm___ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_ ...
,sqrt_rat3__ ...
,sqrt_rat4__ ...
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
% y_jlm___(1+nx,1+l_val,1+m_val) ; 
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
tmp_a2_ = exp(lfactorial(l_val-abs(m_val_)) - lfactorial(l_val+abs(m_val_)));
tmp_a3_ = sqrt(tmp_a1*tmp_a2_);
for nx=0:n_x-1;
tmp_Y_matlab___(1+nx,1+l_val,1:(1+l_val)) = legendre(l_val,x_(1+nx),'unnorm').*tmp_a3_;
end;%for nx=0:n_x-1;
end;%for l_val=0:l_max;
if (verbose>0); disp(sprintf(' %% l_max %d: tmp_y_fortran___ vs sqrt(4*pi)*tmp_Y_matlab___: %0.16f',l_max,fnorm(tmp_y_fortran___-sqrt(4*pi)*tmp_Y_matlab___)/fnorm(tmp_y_fortran___))); end;
l_max = 96*1.5; n_x = 1024*2;
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
tmp_a2_ = exp(lfactorial(l_val-abs(m_val_)) - lfactorial(l_val+abs(m_val_)));
tmp_a3_ = sqrt(tmp_a1*tmp_a2_);
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

x_ = x_(:);
n_x = numel(x_);
u_ = -sqrt((1-x_).*(1+x_));
y_jlm___ = zeros(n_x,l_max+1,l_max+1);
y_jlm___(:,1+0,1+0) = 1.0d0;
for m_val=0:l_max;
%if (m_val>0); y_jlm___(:,1+m_val+0,1+m_val+0) = y_jlm___(:,1+m_val-1,1+m_val-1).*u_*sqrt((2.0d0*m_val-1.0d0)/(2.0d0*m_val)); end;
if (m_val>0); y_jlm___(:,1+m_val+0,1+m_val+0) = y_jlm___(:,1+m_val-1,1+m_val-1).*u_*sqrt_rat0_(1+m_val); end;
%if (m_val<l_max); y_jlm___(:,1+m_val+1,1+m_val+0) = x_.*y_jlm___(:,1+m_val+0,1+m_val+0)*sqrt(2.0d0*m_val+1.0d0); end;
if (m_val<l_max); y_jlm___(:,1+m_val+1,1+m_val+0) = x_.*y_jlm___(:,1+m_val+0,1+m_val+0)*sqrt_2mp1_(1+m_val); end;
for l_val=m_val+2:l_max;
%y_jlm___(:,1+l_val,1+m_val+0) = ((2.0d0*l_val-1.0d0)*x_.*y_jlm___(:,1+l_val-1,1+m_val+0) - sqrt((l_val+m_val-1.0d0)*(l_val-m_val-1.0d0))*y_jlm___(:,1+l_val-2,1+m_val+0))/sqrt((l_val-m_val+0.0d0)*(l_val+m_val+0.0d0));
y_jlm___(:,1+l_val,1+m_val+0) = sqrt_rat4__(1+l_val,1+m_val)*x_.*y_jlm___(:,1+l_val-1,1+m_val+0) - sqrt_rat3__(1+l_val,1+m_val)*y_jlm___(:,1+l_val-2,1+m_val+0);
end;%for l_val=m_val+2:l_max;
end;%for m_val=0:l_max;

%%%%;
% for l_val=0:l_max;
% for m_val=0:l_val;
% y_jlm___(:,1+l_val,1+m_val+0) = y_jlm___(:,1+l_val,1+m_val+0)*sqrt(2.0d0*l_val+1.0d0);
% end;%for m_val=0:l_val;
% end;%for l_val=0:l_max;
%%%%;
y_jlm___ = bsxfun(@times,y_jlm___,reshape(sqrt_2lp1_,[1,1+l_max,1]));
