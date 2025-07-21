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
%%%%%%%%;
% This returns: ;
% d0y_jlm___(1+nx,1+l_val,1+m_val) ; 
% such that ;
% y__(1+nx,1+l_val,1+m_val) == Y__(1+l_val,1+m_val)(x_(1+nx)), where: ;
% Y__ is constructed via: ;
% m_val_=transpose(0:l_val);
% tmp_a1 = ((1+2*l_val)/(4*pi));
% tmp_a2_ = factorial(l_val-abs(m_val_))./factorial(l_val+abs(m_val_));
% tmp_a3_ = sqrt(tmp_a1*tmp_a2_);
% Y__(1+l_val,:) = legendre(l_val,x,'unnorm').*tmp_a3_*sqrt(4*pi);
%%%%%%%%;
% Inputs: ;
% l_max: integer maximum l_value. ;
% x_: double array of size n_x. points to evaluate. ;
% sqrt_2lp1_: double array of size (1+l_max), precomputation. (optional) ;
% sqrt_2mp1_: double array of size (1+l_max), precomputation. (optional) ;
% sqrt_rat0_: double array of size (1+l_max), precomputation. (optional) ;
% sqrt_rat3__: double array of size (1+l_max,1+l_max), precomputation. (optional) ;
% sqrt_rat4__: double array of size (1+l_max,1+l_max), precomputation. (optional) ;
%%%%%%%%;
% Outputs: ;
% d0y_jlm___: double array of size (n_x,1+l_max,1+l_max). ;
% sqrt_2lp1_: double array of size (1+l_max), precomputation. (optional) ;
% sqrt_2mp1_: double array of size (1+l_max), precomputation. (optional) ;
% sqrt_rat0_: double array of size (1+l_max), precomputation. (optional) ;
% sqrt_rat3__: double array of size (1+l_max,1+l_max), precomputation. (optional) ;
% sqrt_rat4__: double array of size (1+l_max,1+l_max), precomputation. (optional) ;
% d1y_jlm___: double array of size (n_x,1+l_max,1+l_max). first-derivative with respect to x. ;
% d2y_jlm___: double array of size (n_x,1+l_max,1+l_max). second-derivative with respect to x. ;
%%%%%%%%;

if (nargin<2);
verbose=2; flag_disp=0; nf=0;
if (verbose>0); disp(sprintf(' %% testing ylgndr_1')); end;
l_max = 4; n_x = 25;
x_ = transpose(2*rand(n_x,1)-1);
tmp_y_ylgndr___ = ylgndr_1(l_max,x_);
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
if (verbose>0); disp(sprintf(' %% l_max %d: tmp_y_ylgndr___ vs sqrt(4*pi)*tmp_Y_matlab___: %0.16f',l_max,fnorm(tmp_y_ylgndr___-sqrt(4*pi)*tmp_Y_matlab___)/fnorm(tmp_y_ylgndr___))); end;
%%%%%%%%;
% Matching to legendre-polynomial. ;
%%%%%%%%;
l_max=16; n_z = 1+1024; z_ = linspace(-1,+1,n_z); dz = 2/max(1,n_z-1);
d0y_jlm___ = ylgndr_1(l_max,z_); d0y_jl0__ = d0y_jlm___(:,:,1+0);
leg_jl0__ = zeros(n_z,1+l_max);
for l_val=0:l_max;
tmp_j_ = legendreP(l_val,z_);
tmp_j2_ = abs(tmp_j_).^2;
tmp_s = 0.5 * dz*(sum(tmp_j2_(2:end-1)) + 0.5*(tmp_j2_(1+0)+tmp_j2_(end)));
tmp_j_ = tmp_j_./max(1e-12,sqrt(tmp_s));
leg_jl0__(:,1+l_val) = tmp_j_;
end;%for l_val=0:l_max;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);fig81s;
imagesc(d0y_jl0__); axisnotick; 
ylabel('z'); xlabel('degree'); title('d0y_jl0__','Interpreter','none');
subplot(1,2,2);fig81s;
imagesc(leg_jl0__); axisnotick; 
ylabel('z'); xlabel('degree'); title('leg_jl0__','Interpreter','none');
end;%if flag_disp;
if (verbose>0); disp(sprintf(' %% leg_jl0__ vs d0y_jl0__: %0.16f',fnorm(leg_jl0__ - d0y_jl0__)/fnorm(leg_jl0__))); end;
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
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
plot(d1y_dif_jlm___(:),d1y_mid_jlm___(:),'.'); axis equal; grid on;
end;%if flag_disp;
%%%%%%%%;
% testing second-derivative. ;
%%%%%%%%;
if (verbose>0); disp(sprintf(' %% testing second-derivative')); end;
d2y_dif_jlm___ = bsxfun(@rdivide,d0y_pos_jlm___ - 2*d0y_mid_jlm___ + d0y_neg_jlm___,reshape(dx_.^2,[n_x,1,1]));
if (verbose>0); disp(sprintf(' %% d2y_dif_jlm___ vs d2y_mid_jlm___: %0.16f',fnorm(d2y_dif_jlm___ - d2y_mid_jlm___)/fnorm(d2y_dif_jlm___))); end;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
plot(d2y_dif_jlm___(:),d2y_mid_jlm___(:),'.'); axis equal; grid on;
end;%if flag_disp;
%%%%%%%%;
% testing stability. ;
%%%%%%%%;
if (verbose>0); disp(sprintf(' %% testing stability')); end;
l_max = 96*2.0; n_x = 1024*2;
x_ = rand(n_x,1);
%%%%;
tmp_t = tic;
[ ...
 tmp_y_ylgndr___ ...
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
tmp_t = toc(tmp_t); if (verbose>0); disp(sprintf(' %% y_ylgndr___ (not precomputation) : %0.2fs',tmp_t)); end;
%%%%;
tmp_t = tic;
[ ...
 bkp_y_ylgndr___ ...
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
tmp_t = toc(tmp_t); if (verbose>0); disp(sprintf(' %% y_ylgndr___ (yes precomputation) : %0.2fs',tmp_t)); end;
assert(fnorm(tmp_y_ylgndr___-bkp_y_ylgndr___)<1e-12);
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
if (verbose>0); disp(sprintf(' %% l_max %d: tmp_y_ylgndr___ vs sqrt(4*pi)*tmp_Y_matlab___: %0.16f',l_max,fnorm(tmp_y_ylgndr___-sqrt(4*pi)*tmp_Y_matlab___)/fnorm(tmp_y_ylgndr___))); end;
if (flag_disp>-1);
figure(1+nf);nf=nf+1;figbig;fig81s;
p_row = 2; p_col = 3; np=0;
fontsize_use = 12;
for ntype=0:5-1;
subplot(p_row,p_col,1+np);np=np+1;cla;
if ntype==0;
llim_ = [-16,1];
tmp__ = squeeze(mean(log10(abs(tmp_y_ylgndr___-sqrt(4*pi)*tmp_Y_matlab___)),1));
str_title = 'error=log10(abs(ylgndr-matlab))';
end;%if ntype==0;
if ntype==1;
llim_ = [-16,1];
tmp__ = squeeze(mean(log10(abs(tmp_y_ylgndr___)),1));
str_title = 'magnitude=log10(abs(ylgnd))';
end;%if ntype==1;
if ntype==2;
llim_ = [-16,1];
tmp__ = squeeze(mean(log10(abs(sqrt(4*pi)*tmp_Y_matlab___)),1));
str_title = 'error=log10(abs(matlab))';
end;%if ntype==2;
if ntype==3;
llim_ = [0,1];
tmp__ = isfinite(squeeze(mean((abs(tmp_y_ylgndr___)),1)));
str_title = 'isfinite((abs(ylgnd)))';
end;%if ntype==3;
if ntype==4;
llim_ = [0,1];
tmp__ = isfinite(squeeze(mean((abs(sqrt(4*pi)*tmp_Y_matlab___)),1)));
str_title = 'isfinite((abs(matlab)))';
end;%if ntype==4;
imagesc(tmp__,llim_);
xlabel('m_val','Interpreter','none'); xtickangle(90);
set(gca,'XTick',1:16:1+l_max,'XTickLabel',0:16:l_max);
ylabel('l_val','Interpreter','none');
set(gca,'YTick',1:16:1+l_max,'YTickLabel',0:16:l_max);
title(str_title,'Interpreter','none');
set(gca,'FontSize',fontsize_use);
c_ = colorbar;
end;%for ntype=0:3-1;
end;%if (flag_disp>-1);
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
