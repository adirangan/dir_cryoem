function ...
[ ...
 parameter ...
,d0y_jlm___ ...
,d1y_jlm___ ...
,d2y_jlm___ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_m_ ...
,sqrt_rat3_lm__ ...
,sqrt_rat4_lm__ ...
] = ...
ylgndr_2( ...
 parameter ...
,l_max ...
,x_ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_m_ ...
,sqrt_rat3_lm__ ...
,sqrt_rat4_lm__ ...
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
% sqrt_rat0_m_: double array of size (1+l_max), precomputation. (optional) ;
% sqrt_rat3_lm__: double array of size (1+l_max,1+l_max), precomputation. (optional) ;
% sqrt_rat4_lm__: double array of size (1+l_max,1+l_max), precomputation. (optional) ;
%%%%%%%%;
% Outputs: ;
% d0y_jlm___: double array of size (n_x,1+l_max,1+l_max). ;
% sqrt_2lp1_: double array of size (1+l_max), precomputation. (optional) ;
% sqrt_2mp1_: double array of size (1+l_max), precomputation. (optional) ;
% sqrt_rat0_m_: double array of size (1+l_max), precomputation. (optional) ;
% sqrt_rat3_lm__: double array of size (1+l_max,1+l_max), precomputation. (optional) ;
% sqrt_rat4_lm__: double array of size (1+l_max,1+l_max), precomputation. (optional) ;
% d1y_jlm___: double array of size (n_x,1+l_max,1+l_max). first-derivative with respect to x. ;
% d2y_jlm___: double array of size (n_x,1+l_max,1+l_max). second-derivative with respect to x. ;
%%%%%%%%;

str_thisfunction = 'ylgndr_2';

if (nargin<2);
flag_verbose=2; flag_disp=0; nf=0;
if (flag_verbose>0); disp(sprintf(' %% testing ylgndr_2')); end;
l_max = 4; n_x = 25;
x_ = transpose(linspace(-0.9,+0.9,n_x));
[~,tmp_y_ylgndr___] = ylgndr_2([],l_max,x_);
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
fnorm_disp(flag_verbose,'tmp_y_ylgndr___',tmp_y_ylgndr___,'sqrt(4*pi)*tmp_Y_matlab___',sqrt(4*pi)*tmp_Y_matlab___,' %% <-- should be 0');
%%%%;
% compare to ../dir_rangan_python/tmp_y_ylgndr_jlm___.ascii. ;
%%%%;
dir_ascii = '../dir_rangan_python/dir_ascii';
fname_ascii = sprintf('%s/tmp_y_ylgndr_jlm___.ascii',dir_ascii);
if ~exist(fname_ascii,'file'); disp(sprintf(' %% Warning, %s not found',fname_ascii)); end;
if  exist(fname_ascii,'file');
disp(sprintf(' %% %s found, loading',fname_ascii));
tmp_y_ylgndr_ascii___ = reshape(textread(fname_ascii),size(tmp_y_ylgndr___));
fnorm_disp(flag_verbose,'tmp_y_ylgndr___',tmp_y_ylgndr___,'tmp_y_ylgndr_ascii___',tmp_y_ylgndr_ascii___,' %% <-- should be <1e-6');
end;%if  exist(fname_ascii,'file');
%%%%%%%%;
% Matching to legendre-polynomial. ;
%%%%%%%%;
l_max=16; n_z = 1+1024; z_ = linspace(-1,+1,n_z); dz = 2/max(1,n_z-1);
[~,d0y_jlm___] = ylgndr_2([],l_max,z_); d0y_jl0__ = d0y_jlm___(:,:,1+0);
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
fnorm_disp(flag_verbose,'leg_jl0__',leg_jl0__,'d0y_jl0__',d0y_jl0__,' %% <-- should be <1e-2');
%%%%%%%%;
% testing first-derivative. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% testing first-derivative')); end;
l_max = 4; n_x = 25; x_ = transpose(linspace(-0.9,+0.9,n_x)); x_ = cos(acos(x_)); da = pi*1e-4; dx_ = periodize(x_-cos(acos(x_)+da),-1,+1); if flag_disp>0; plot(dx_,'.'); end;
[ ...
 ~ ...
,d0y_mid_jlm___ ...
,d1y_mid_jlm___ ...
,d2y_mid_jlm___ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_m_ ...
,sqrt_rat3_lm__ ...
,sqrt_rat4_lm__ ...
] = ...
ylgndr_2( ...
 [] ...
,l_max ...
,x_ ...
);
[~,d0y_pos_jlm___] = ylgndr_2([],l_max,x_+dx_);
[~,d0y_neg_jlm___] = ylgndr_2([],l_max,x_-dx_);
d1y_dif_jlm___ = bsxfun(@rdivide,d0y_pos_jlm___-d0y_neg_jlm___,reshape(2*dx_,[n_x,1,1]));
fnorm_disp(flag_verbose,'d1y_dif_jlm___',d1y_dif_jlm___,'d1y_mid_jlm___',d1y_mid_jlm___,' %% <-- should be <1e-2');
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
plot(d1y_dif_jlm___(:),d1y_mid_jlm___(:),'.'); axis equal; grid on;
end;%if flag_disp;
%%%%;
% compare to ../dir_rangan_python/d?y_mid_jlm___.ascii. ;
%%%%;
for ntype=0:8-1;
na=0;
if ntype==na; tmp_str = 'd0y_mid_jlm___'; tmp = d0y_mid_jlm___; end; na=na+1;
if ntype==na; tmp_str = 'd1y_mid_jlm___'; tmp = d1y_mid_jlm___; end; na=na+1;
if ntype==na; tmp_str = 'd2y_mid_jlm___'; tmp = d2y_mid_jlm___; end; na=na+1;
if ntype==na; tmp_str = 'sqrt_2lp1_'; tmp = sqrt_2lp1_; end; na=na+1;
if ntype==na; tmp_str = 'sqrt_2mp1_'; tmp = sqrt_2mp1_; end; na=na+1;
if ntype==na; tmp_str = 'sqrt_rat0_m_'; tmp = sqrt_rat0_m_; end; na=na+1;
if ntype==na; tmp_str = 'sqrt_rat3_lm__'; tmp = sqrt_rat3_lm__; end; na=na+1;
if ntype==na; tmp_str = 'sqrt_rat4_lm__'; tmp = sqrt_rat4_lm__; end; na=na+1;
fname_ascii = sprintf('%s/%s.ascii',dir_ascii,tmp_str);
if ~exist(fname_ascii,'file'); disp(sprintf(' %% Warning, %s not found',fname_ascii)); end;
if  exist(fname_ascii,'file');
disp(sprintf(' %% %s found, loading',fname_ascii));
tmp_ascii = reshape(textread(fname_ascii),size(tmp));
fnorm_disp(flag_verbose,tmp_str,tmp,'ascii',tmp_ascii,' %% <-- should be <1e-6');
end;%if  exist(fname_ascii,'file');
end;%for ntype=0:8-1;
%%%%%%%%;
% testing second-derivative. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% testing second-derivative')); end;
d2y_dif_jlm___ = bsxfun(@rdivide,d0y_pos_jlm___ - 2*d0y_mid_jlm___ + d0y_neg_jlm___,reshape(dx_.^2,[n_x,1,1]));
fnorm_disp(flag_verbose,'d2y_dif_jlm___',d2y_dif_jlm___,'d2y_mid_jlm___',d2y_mid_jlm___,' %% <-- should be <1e-2');
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
plot(d2y_dif_jlm___(:),d2y_mid_jlm___(:),'.'); axis equal; grid on;
end;%if flag_disp;
%%%%%%%%;
% testing stability. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% testing stability')); end;
l_max = 96*2.0; n_x = 1024*2;
x_ = transpose(linspace(-0.9,+0.9,n_x));
%%%%;
tmp_t = tic;
[ ...
 ~ ...
,tmp_y_ylgndr___ ...
,~ ...
,~ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_m_ ...
,sqrt_rat3_lm__ ...
,sqrt_rat4_lm__ ...
] = ...
ylgndr_2( ...
 [] ...
,l_max ...
,x_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% y_ylgndr___ (not precomputation) : %0.2fs',tmp_t)); end;
%%%%;
tmp_t = tic;
[ ...
 ~ ...
,bkp_y_ylgndr___ ...
] = ...
ylgndr_2( ...
 [] ...
,l_max ...
,x_ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_m_ ...
,sqrt_rat3_lm__ ...
,sqrt_rat4_lm__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% y_ylgndr___ (yes precomputation) : %0.2fs',tmp_t)); end;
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
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% Y_matlab___: %0.2fs',tmp_t)); end;
%%%%;
fnorm_disp(flag_verbose,'tmp_y_ylgndr___',tmp_y_ylgndr___,'sqrt(4*pi)*tmp_Y_matlab___',sqrt(4*pi)*tmp_Y_matlab___,' %% <-- stability test');
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
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); l_max = []; end; na=na+1;
if (nargin<1+na); x_ = []; end; na=na+1;
if (nargin<1+na); sqrt_2lp1_ = []; end; na=na+1;
if (nargin<1+na); sqrt_2mp1_ = []; end; na=na+1;
if (nargin<1+na); sqrt_rat0_m_ = []; end; na=na+1;
if (nargin<1+na); sqrt_rat3_lm__ = []; end; na=na+1;
if (nargin<1+na); sqrt_rat4_lm__ = []; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_d'); parameter.flag_d=1; end;
flag_d=parameter.flag_d;
if ~isfield(parameter,'flag_dd'); parameter.flag_dd=1; end;
flag_dd=parameter.flag_dd;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

flag_precomputation = ( isempty(sqrt_2lp1_) | isempty(sqrt_2mp1_) | isempty(sqrt_rat0_m_) | isempty(sqrt_rat3_lm__) | isempty(sqrt_rat4_lm__) );
if flag_precomputation;
l_val_ = transpose(0:l_max);
m_val_ = transpose(0:l_max);
base_2lm1_ = 2.0d0*l_val_-1.0d0;
sqrt_2lp1_ = sqrt(2.0d0*l_val_+1.0d0);
sqrt_2mp1_ = sqrt(2.0d0*m_val_+1.0d0);
sqrt_rat0_m_ = sqrt((2.0d0*m_val_-1.0d0)./max(1,2.0d0*m_val_));
sqrt_rat1_lm__ = ...
  sqrt( ...
	max(0 ...
	    , (bsxfun(@plus,reshape(l_val_,[1+l_max,1]),reshape(m_val_,[1,1+l_max]))-1.0d0) ...
	      .* ...
	      (bsxfun(@minus,reshape(l_val_,[1+l_max,1]),reshape(m_val_,[1,1+l_max]))-1.0d0) ...
	   ) ...
      ) ...
  ;
sqrt_rat2_lm__ = ...
  sqrt( ...
	max(0 ...
	    ,(bsxfun(@minus,reshape(l_val_,[1+l_max,1]),reshape(m_val_,[1,1+l_max]))+0.0d0) ...
	     .* ...
	     (bsxfun(@plus,reshape(l_val_,[1+l_max,1]),reshape(m_val_,[1,1+l_max]))+0.0d0) ...
	   ) ...
      ) ...
  ;
sqrt_rat3_lm__ = sqrt_rat1_lm__./sqrt_rat2_lm__;
sqrt_rat4_lm__ = bsxfun(@rdivide,reshape(base_2lm1_,[1+l_max,1]),sqrt_rat2_lm__);
for m_val=0:l_max;
if m_val==0; sqrt_rat0_m_(1+m_val)=0; end;
for l_val=0:l_max;
if m_val>=l_val; sqrt_rat3_lm__(1+l_val,1+m_val)=0; sqrt_rat4_lm__(1+l_val,1+m_val)=0; end;
end;%for l_val=0:l_max;
end;%for m_val=0:l_max;
end;%if flag_precomputation;

flag_d0 = 1;
flag_d1 = flag_d;
flag_d2 = flag_dd;

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
if flag_d0; d0y_jlm___(:,1+m_val+0,1+m_val+0) = d0y_jlm___(:,1+m_val-1,1+m_val-1).*d0u_*sqrt_rat0_m_(1+m_val); end;
if flag_d1; d1y_jlm___(:,1+m_val+0,1+m_val+0) = d0y_jlm___(:,1+m_val-0,1+m_val-0).*(-m_val).*x_./d0u_.^2; end;
%if flag_d2; d2y_jlm___(:,1+m_val+0,1+m_val+0) = d2y_jlm___(:,1+m_val-1,1+m_val-1).*d0u_*sqrt_rat0_m_(1+m_val) + 2*d1y_jlm___(:,1+m_val-1,1+m_val-1).*d1u_*sqrt_rat0_m_(1+m_val) + d0y_jlm___(:,1+m_val-1,1+m_val-1).*d2u_*sqrt_rat0_(1+m_val) ; end;
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
if flag_d0; d0y_jlm___(:,1+l_val,1+m_val+0) = sqrt_rat4_lm__(1+l_val,1+m_val)*x_.*d0y_jlm___(:,1+l_val-1,1+m_val+0) - sqrt_rat3_lm__(1+l_val,1+m_val)*d0y_jlm___(:,1+l_val-2,1+m_val+0); end;
if flag_d1; d1y_jlm___(:,1+l_val,1+m_val+0) = sqrt_rat4_lm__(1+l_val,1+m_val)*(d0y_jlm___(:,1+l_val-1,1+m_val+0) + x_.*d1y_jlm___(:,1+l_val-1,1+m_val+0)) - sqrt_rat3_lm__(1+l_val,1+m_val)*d1y_jlm___(:,1+l_val-2,1+m_val+0); end;
if flag_d2; d2y_jlm___(:,1+l_val,1+m_val+0) = sqrt_rat4_lm__(1+l_val,1+m_val)*(2*d1y_jlm___(:,1+l_val-1,1+m_val+0) + x_.*d2y_jlm___(:,1+l_val-1,1+m_val+0)) - sqrt_rat3_lm__(1+l_val,1+m_val)*d2y_jlm___(:,1+l_val-2,1+m_val+0); end;
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
