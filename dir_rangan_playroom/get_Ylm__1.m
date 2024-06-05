function ...
[ ...
 d0Y_lmj___ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_ ...
,sqrt_rat3__ ...
,sqrt_rat4__ ...
,d1Y_lmj___ ...
,d2Y_lmj___ ...
] = ...
get_Ylm__1( ...
 n_l ...
,l_val_ ...
,n_all ...
,azimu_b_all_ ...
,polar_a_all_ ...
,flag_flip ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_ ...
,sqrt_rat3__ ...
,sqrt_rat4__ ...
);
%%%%%%%%;
% This returns a cell array d0Y_lmj___ such that: ;
% d0Y_lmj___{1+nl} corresponds to spherical harmonics associated with a particular l_val=l_val_(1+nl). ;
% The matrix d0Y_lmj___{1+nl} is of size (2*l_val+1,n_all), ;
% and contains the values d0Y_lmj___{l_val}(m_val,1+np) of the spherical harmonic: ;
% Y_{l_val}^{m_val}(azimu_b_all_(1+np),polar_a_all_(1+np));
% If flag_flip is set to 1, ;
% we calculate the values of Ylm for the antipodal point on the sphere instead. ;
%%%%%%%%;
% Inputs: ;
% n_l: integer number of l_val (typically 1+l_max). ;
% l_val_: integer array of size n_l storing l_val (typically 0:l_max). ;
% n_all: integer number of points to evaluate. ;
% azimu_b_all_ = double array of size n_all. azimu_b values. ;
% polar_a_all_ = double array of size n_all. polar_a values. ;
% flag_flip = integer. if 1 then evaluate at antipodal points. (default 0) ;
% sqrt_2lp1_: double array of size (1+l_max), precomputation. (optional) ;
% sqrt_2mp1_: double array of size (1+l_max), precomputation. (optional) ;
% sqrt_rat0_: double array of size (1+l_max), precomputation. (optional) ;
% sqrt_rat3__: double array of size (1+l_max,1+l_max), precomputation. (optional) ;
% sqrt_rat4__: double array of size (1+l_max,1+l_max), precomputation. (optional) ;
%%%%%%%%;
% Outputs: ;
% d0Y_lmj___: cell array of size (1+l_max). ;
%             d0Y_lmj___{1+l_val} is a complex array of size (1+2*l_val,n_all) ;
%                                 containing spherical-harmonic evaluations. ;
% d1Y_lmj___: cell array of size (1+l_max). contains first-derivative with respect to polar_a. ;
% d1Y_lmj___: cell array of size (1+l_max). contains second-derivative with respect to polar_a. ;
%%%%%%%%;

flag_verbose=0;
str_thisfunction = 'get_Ylm__1';

%%%%%%%%;
if (nargin<1);
%%%%%%%%;
disp(sprintf(' %% testing %s',str_thisfunction));
flag_verbose=2;
nf=0;
n_l = 48;
l_val_ = 3+[0:n_l-1];
%n_all = 1024*(1/8);
n_all = 2; %<-- reduce n_all if plotting first-derivative and second-derivative errors. ;
azimu_b_all_ = 2*pi*rand(n_all,1);
polar_a_all_ = 1*pi*rand(n_all,1);
polar_a_all_(1:n_all/2) = polar_a_all_(n_all/2+1:end); %<-- test nonunique polar_a_all_. ;
%%%%;
flag_flip = 0;
tmp_t = tic();
[Ylm_0__] = get_Ylm__(n_l,l_val_,n_all,azimu_b_all_,polar_a_all_);
tmp_t = toc(tmp_t); if (flag_verbose>-1); disp(sprintf(' %% get_Ylm__ : %0.6fs',tmp_t)); end;
tmp_t = tic();
[Ylm_1__,sqrt_2lp1_,sqrt_2mp1_,sqrt_rat0_,sqrt_rat3__,sqrt_rat4__] = get_Ylm__1(n_l,l_val_,n_all,azimu_b_all_,polar_a_all_);
tmp_t = toc(tmp_t); if (flag_verbose>-1); disp(sprintf(' %% get_Ylm__1 (not precomputation): %0.6fs',tmp_t)); end;
tmp_t = tic();
[Ylm_1__] = get_Ylm__1(n_l,l_val_,n_all,azimu_b_all_,polar_a_all_,[],sqrt_2lp1_,sqrt_2mp1_,sqrt_rat0_,sqrt_rat3__,sqrt_rat4__);
tmp_t = toc(tmp_t); if (flag_verbose>-1); disp(sprintf(' %% get_Ylm__1 (yes precomputation): %0.6fs',tmp_t)); end;
tmp_error_numerator=0;
tmp_error_denomator=0;
for nl=0:n_l-1;
tmp_error_numerator = tmp_error_numerator + fnorm(Ylm_0__{1+nl}-Ylm_1__{1+nl});
tmp_error_denomator = tmp_error_denomator + fnorm(Ylm_0__{1+nl});
end;%for nl=0:n_l-1;
disp(sprintf(' %% flag_flip==0: Ylm_0__ vs Ylm_1__: %0.16f',tmp_error_numerator/max(1e-12,tmp_error_denomator)));
%%%%;
flag_flip = 1;
tmp_t = tic();
[Ylm_0__] = get_Ylm__(n_l,l_val_,n_all,azimu_b_all_,polar_a_all_,flag_flip);
tmp_t = toc(tmp_t); if (flag_verbose>-1); disp(sprintf(' %% get_Ylm__ : %0.6fs',tmp_t)); end;
tmp_t = tic();
[Ylm_1__] = get_Ylm__1(n_l,l_val_,n_all,azimu_b_all_,polar_a_all_,flag_flip,sqrt_2lp1_,sqrt_2mp1_,sqrt_rat0_,sqrt_rat3__,sqrt_rat4__);
tmp_t = toc(tmp_t); if (flag_verbose>-1); disp(sprintf(' %% get_Ylm__1 (yes precomputation): %0.6fs',tmp_t)); end;
tmp_error_numerator=0;
tmp_error_denomator=0;
for nl=0:n_l-1;
tmp_error_numerator = tmp_error_numerator + fnorm(Ylm_0__{1+nl}-Ylm_1__{1+nl});
tmp_error_denomator = tmp_error_denomator + fnorm(Ylm_0__{1+nl});
end;%for nl=0:n_l-1;
disp(sprintf(' %% flag_flip==0: Ylm_0__ vs Ylm_1__: %0.16f',tmp_error_numerator/max(1e-12,tmp_error_denomator)));
%%%%;
disp(sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%; '));
disp(sprintf(' %% testing first-derivative'));
disp(sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%; '));
flag_flip = 1;
da = pi*1e-4;
[ ...
 d0Y_mid_lmj___ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_ ...
,sqrt_rat3__ ...
,sqrt_rat4__ ...
,d1Y_mid_lmj___ ...
,d2Y_mid_lmj___ ...
] = ...
get_Ylm__1( ...
 n_l ...
,l_val_ ...
,n_all ...
,azimu_b_all_ ...
,polar_a_all_ + 0*da ...
,flag_flip ...
);
[ ...
 d0Y_pos_lmj___ ...
] = ...
get_Ylm__1( ...
 n_l ...
,l_val_ ...
,n_all ...
,azimu_b_all_ ...
,polar_a_all_ + 1*da ...
,flag_flip ...
);
[ ...
 d0Y_neg_lmj___ ...
] = ...
get_Ylm__1( ...
 n_l ...
,l_val_ ...
,n_all ...
,azimu_b_all_ ...
,polar_a_all_ - 1*da ...
,flag_flip ...
);
d1Y_dif_lmj___ = cell(n_l,1); for nl=0:n_l-1; d1Y_dif_lmj___{1+nl} = (d0Y_pos_lmj___{1+nl} - d0Y_neg_lmj___{1+nl})/max(1e-12,2*da); end;
tmp_error_numerator = 0;
tmp_error_denomator = 0;
for nl=0:n_l-1;
tmp_error_numerator = tmp_error_numerator + fnorm(d1Y_dif_lmj___{1+nl} - d1Y_mid_lmj___{1+nl});
tmp_error_denomator = tmp_error_denomator + fnorm(d1Y_dif_lmj___{1+nl});
end;%for nl=0:n_l-1;
disp(sprintf(' %% tmp_error: %0.16f',tmp_error_numerator/max(1e-12,tmp_error_denomator)));
if (flag_verbose>1);
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
hold on;
for nl=0:n_l-1;
plot(real(d1Y_dif_lmj___{1+nl}(:)),real(d1Y_mid_lmj___{1+nl}(:)),'k.');
end;%for nl=0:n_l-1;
hold off;
xlabel('d1Y_dif_lmj___','Interpreter','none');
ylabel('d1Y_mid_lmj___','Interpreter','none');
subplot(1,2,2);
hold on;
for nl=0:n_l-1;
plot(l_val_(1+nl),log10(abs(d1Y_dif_lmj___{1+nl}(:)-d1Y_mid_lmj___{1+nl}(:))),'k.');
end;%for nl=0:n_l-1;
hold off;
xlabel('l_val_','Interpreter','none');
ylabel('log10(abs(d1Y_dif_lmj___-d1Y_mid_lmj___))','Interpreter','none');
end;%if (flag_verbose>1);
disp(sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%; '));
disp(sprintf(' %% testing second-derivative'));
disp(sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%; '));
d2Y_dif_lmj___ = cell(n_l,1); for nl=0:n_l-1; d2Y_dif_lmj___{1+nl} = (d0Y_pos_lmj___{1+nl} - 2*d0Y_mid_lmj___{1+nl} + d0Y_neg_lmj___{1+nl})/max(1e-12,da^2); end;
tmp_error_numerator = 0;
tmp_error_denomator = 0;
for nl=0:n_l-1;
tmp_error_numerator = tmp_error_numerator + fnorm(d2Y_dif_lmj___{1+nl} - d2Y_mid_lmj___{1+nl});
tmp_error_denomator = tmp_error_denomator + fnorm(d2Y_dif_lmj___{1+nl});
end;%for nl=0:n_l-1;
disp(sprintf(' %% tmp_error: %0.16f',tmp_error_numerator/max(1e-12,tmp_error_denomator)));
if (flag_verbose>1);
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
hold on;
for nl=0:n_l-1;
plot(real(d2Y_dif_lmj___{1+nl}(:)),real(d2Y_mid_lmj___{1+nl}(:)),'k.');
end;%for nl=0:n_l-1;
hold off;
xlabel('d2Y_dif_lmj___','Interpreter','none');
ylabel('d2Y_mid_lmj___','Interpreter','none');
subplot(1,2,2);
hold on;
for nl=0:n_l-1;
plot(l_val_(1+nl),log10(abs(d2Y_dif_lmj___{1+nl}(:)-d2Y_mid_lmj___{1+nl}(:))),'k.');
end;%for nl=0:n_l-1;
hold off;
xlabel('l_val_','Interpreter','none');
ylabel('log10(abs(d2Y_dif_lmj___-d2Y_mid_lmj___))','Interpreter','none');
end;%if (flag_verbose>1);
%%%%%%%%;
disp('returning');return;
end;%if (nargin<1);
%%%%%%%%;

na=0;
if (nargin<1+na); n_l=[]; end; na=na+1;
if (nargin<1+na); l_val_=[]; end; na=na+1;
if (nargin<1+na); n_all=[]; end; na=na+1;
if (nargin<1+na); azimu_b_all_=[]; end; na=na+1;
if (nargin<1+na); polar_a_all_=[]; end; na=na+1;
if (nargin<1+na); flag_flip=[]; end; na=na+1;
if (nargin<1+na); sqrt_2lp1_ = []; end; na=na+1;
if (nargin<1+na); sqrt_2mp1_ = []; end; na=na+1;
if (nargin<1+na); sqrt_rat0_ = []; end; na=na+1;
if (nargin<1+na); sqrt_rat3__ = []; end; na=na+1;
if (nargin<1+na); sqrt_rat4__ = []; end; na=na+1;

flag_precomputation = ( isempty(sqrt_2lp1_) | isempty(sqrt_2mp1_) | isempty(sqrt_rat0_) | isempty(sqrt_rat3__) | isempty(sqrt_rat4__) );
if isempty(flag_flip); flag_flip=0; end;
if flag_flip; polar_a_all_ = 1*pi - polar_a_all_; end; %<-- switch hemispheres. ;
if flag_flip; azimu_b_all_ = 1*pi + azimu_b_all_; end; %<-- add twelve hours. ;

flag_d0 = 1;
flag_d1 = nargout>=7;
flag_d2 = nargout>=8;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if flag_d0; d0Y_lmj___ = cell(n_l,1); end;
if flag_d1; d1Y_lmj___ = cell(n_l,1); end;
if flag_d2; d2Y_lmj___ = cell(n_l,1); end;
tmp_t = tic();
[polar_a_unique_,ij_unique_,ij_return_] = unique(polar_a_all_(1:n_all));
polar_a_unique_ = reshape(polar_a_unique_,[numel(polar_a_unique_),1]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% polar_a_unique_: %0.6fs',tmp_t)); end;
tmp_t = tic();
l_max = max(l_val_);
%%%%;
if ( flag_d0 & ~flag_d1 & ~flag_d2 );
[ ...
 d0y_jlm___ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_ ...
,sqrt_rat3__ ...
,sqrt_rat4__ ...
] = ...
ylgndr_1( ...
 l_max ...
,cos(polar_a_unique_) ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_ ...
,sqrt_rat3__ ...
,sqrt_rat4__ ...
);
end;%if ( flag_d0 & ~flag_d1 & ~flag_d2 );
%%%%;
if ( flag_d0 &  flag_d1 & ~flag_d2 );
[ ...
 d0y_jlm___ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_ ...
,sqrt_rat3__ ...
,sqrt_rat4__ ...
,d1y_jlm___ ...
] = ...
ylgndr_1( ...
 l_max ...
,cos(polar_a_unique_) ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_ ...
,sqrt_rat3__ ...
,sqrt_rat4__ ...
);
end;%if ( flag_d0 &  flag_d1 & ~flag_d2 );
%%%%;
if ( flag_d0 &  flag_d1 &  flag_d2 );
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
,cos(polar_a_unique_) ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_ ...
,sqrt_rat3__ ...
,sqrt_rat4__ ...
);
end;%if ( flag_d0 &  flag_d1 &  flag_d2 );
%%%%;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% ylgndr_1 (flag_precomputation %d): %0.6fs',flag_precomputation,tmp_t)); end;

tmp_t = tic();
expimb_mj__ = exp(+i*bsxfun(@times,reshape(-l_max:+l_max,[1+2*l_max,1]),reshape(azimu_b_all_(1:n_all),[1,n_all])));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% expimb_mj__: %0.6fs',tmp_t)); end;

flag_calculate = 0;
if flag_calculate;
%%%%%%%%;
% compare with loop below. ;
%%%%%%%%;
tmp_t = tic();
if flag_d0;
tmp_d0Y_lmj___ = permute(d0y_jlm___,[2,3,1]); %<-- permutation before inflation is faster. ;
tmp_d0Y_lmj___ = bsxfun(@times,tmp_d0Y_lmj___(:,1+abs(-l_max:+l_max),ij_return_),reshape(expimb_mj__,[1,1+2*l_max,n_all]))/sqrt(4*pi);
%tmp_d0Y_lmj___ = permute(d0y_jlm___(ij_return_,:,1+abs(-l_max:+l_max)),[2,3,1]); %<-- permutation after inflation is slower. ;
%tmp_d0Y_lmj___ = bsxfun(@times,tmp_d0Y_lmj___,reshape(expimb_mj__,[1,1+2*l_max,n_all]))/sqrt(4*pi);
for nl=0:n_l-1;
l_val = l_val_(1+nl);
d0Y_lmj___{1+nl} = reshape(tmp_d0Y_lmj___(1+l_val,1+l_max+[-l_val:+l_val],:),[1+2*l_val,n_all]);
end;%for nl=0:n_l-1;
end;%if flag_d0;
%%%%;
if flag_d1;
tmp_d1Y_lmj___ = ((-1)^flag_flip)^1 * permute(bsxfun(@times,d1y_jlm___,-sin(polar_a_unique_)),[2,3,1]);
tmp_d1Y_lmj___ = bsxfun(@times,tmp_d1Y_lmj___(:,1+abs(-l_max:+l_max),ij_return_),reshape(expimb_mj__,[1,1+2*l_max,n_all]))/sqrt(4*pi);
for nl=0:n_l-1;
l_val = l_val_(1+nl);
d1Y_lmj___{1+nl} = reshape(tmp_d1Y_lmj___(1+l_val,1+l_max+[-l_val:+l_val],:),[1+2*l_val,n_all]);
end;%for nl=0:n_l-1;
end;%if flag_d1;
%%%%;
if flag_d2;
tmp_d2Y_lmj___ = ((-1)^flag_flip)^2 * permute(bsxfun(@times,d2y_jlm___,(-sin(polar_a_unique_)).^2) + bsxfun(@times,d1y_jlm___,-cos(polar_a_unique_)),[2,3,1]);
tmp_d2Y_lmj___ = bsxfun(@times,tmp_d2Y_lmj___(:,1+abs(-l_max:+l_max),ij_return_),reshape(expimb_mj__,[1,1+2*l_max,n_all]))/sqrt(4*pi);
for nl=0:n_l-1;
l_val = l_val_(1+nl);
d2Y_lmj___{1+nl} = reshape(tmp_d2Y_lmj___(1+l_val,1+l_max+[-l_val:+l_val],:),[1+2*l_val,n_all]);
end;%for nl=0:n_l-1;
end;%if flag_d2;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% permutation before inflation: %0.6fs',tmp_t)); end;
%%%%%%%%;
end;%if flag_calculate;

flag_calculate = 1;
if flag_calculate;
%%%%%%%%;
% compare with vectorized version above. ;
%%%%%%%%;
if flag_d1; c1y_jlm___ = ((-1)^flag_flip)^1 * bsxfun(@times,d1y_jlm___,-sin(polar_a_unique_)); end;
if flag_d2; c2y_jlm___ = ((-1)^flag_flip)^2 * (bsxfun(@times,d2y_jlm___,(-sin(polar_a_unique_)).^2) + bsxfun(@times,d1y_jlm___,-cos(polar_a_unique_))); end;
tmp_t = tic();
for nl=0:n_l-1;
l_val = l_val_(1+nl);
if (flag_verbose>1); disp(sprintf(' %% nl %d/%d: l_val %d',nl,n_l,l_val)); end;
if flag_d0; d0Y_lmj___{1+nl} = zeros(2*l_val+1,n_all); end;
if flag_d1; d1Y_lmj___{1+nl} = zeros(2*l_val+1,n_all); end;
if flag_d2; d2Y_lmj___{1+nl} = zeros(2*l_val+1,n_all); end;
for m_val=-l_val:+l_val;
m_abs = abs(m_val);
if flag_d0; d0Y_lmj___{1+nl}(1+l_val+m_val,:) = reshape(d0y_jlm___(ij_return_,1+l_val,1+m_abs),[1,n_all]).*expimb_mj__(1+l_max+m_val,:)/sqrt(4*pi); end;
if flag_d1; d1Y_lmj___{1+nl}(1+l_val+m_val,:) = reshape(c1y_jlm___(ij_return_,1+l_val,1+m_abs),[1,n_all]).*expimb_mj__(1+l_max+m_val,:)/sqrt(4*pi); end;
if flag_d2; d2Y_lmj___{1+nl}(1+l_val+m_val,:) = reshape(c2y_jlm___(ij_return_,1+l_val,1+m_abs),[1,n_all]).*expimb_mj__(1+l_max+m_val,:)/sqrt(4*pi); end;
end;%for m_val=-l_val:+l_val;
end;%for nl=0:n_l-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% loop over nl and m_val: %0.6fs',tmp_t)); end;
%%%%%%%%;
end;%if flag_calculate;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
