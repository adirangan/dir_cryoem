function [Ylm__] = get_Ylm__(n_l,l_val_,n_all,azimu_b_all_,polar_a_all_,flag_flip);
% This returns a cell array Ylm__ such that: ;
% Ylm__{nl} corresponds to the spherical harmonics associated with a particular l_val=l_val_(nl). ;
% The matrix Ylm__{nl} is of size (2*l_val+1,n_all), ;
% and contains the values Ylm__{l_val}(m_val,np) of the spherical harmonic: ;
% Y_{l_val}^{m_val}(azimu_b_{np},polar_a_{np});
% If flag_flip is used, we calculate the values of Ylm for the polar-reflection, ;
% or antipodal point on the sphere instead. ;

%%%%%%%%;
if (nargin<1);
%%%%%%%%;
disp(sprintf(' %% testing get_Ylm__'));
%%%%%%%%;
verbose = 1; nf=0;
l_max = 48; k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(4*pi);
%%%%;
[ ...
 n_all ...
,azimu_b_all_ ...
,polar_a_all_ ...
,weight_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,n_polar_a ...
,polar_a_ ...
,n_azimu_b_ ...
] = ...
sample_shell_6( ...
 k_p_r_max ...
,k_eq_d ...
) ;
%%%%;
tmp_t = tic();
[Ylm__] = get_Ylm__(1+l_max,[0:l_max],n_all,azimu_b_all_,polar_a_all_);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% get_Ylm__: %0.2fs',tmp_t)); end;
tmp_t = tic();
Y_lmj___ = zeros((1+l_max)^2,n_all);
for l_val=0:l_max;
Y_lmj__(1+l_val^2 + [0:1+2*l_val-1],:) = Ylm__{1+l_val};
end;%for l_val=0:l_max;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% Y_lmj__: %0.2fs',tmp_t)); end;
tmp_t = tic();
YY_lmlm__ = (4*pi) * bsxfun(@times,conj(Y_lmj__),reshape(weight_all_/(4*pi*k_p_r_max^2),[1,n_all]))*transpose(Y_lmj__); %<-- should be identity. ;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% YY_lmlm__: %0.2fs',tmp_t)); end;
%%%%;
if (verbose>0);
figure(1+nf);nf=nf+1;clf;figsml;
plot(real(diag(YY_lmlm__)),'o'); xlabel('lm'); ylabel('should be 1'); title('diag(YY_lmlm__)','Interpreter','none');
end;%if (verbose>0);
%%%%;
if (verbose>1);
figure(1+nf);nf=nf+1;clf;figbig;fig81s; p_row = 2; p_col = 2; np=0; ylim_ = [0,1]; llim_ = [-12,0];
subplot(p_row,p_col,1+np);np=np+1; imagesc(real(YY_lmlm__),ylim_); axis image; axisnotick; title('real(YY_lmlm__)','Interpreter','none'); colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(imag(YY_lmlm__),ylim_); axis image; axisnotick; title('imag(YY_lmlm__)','Interpreter','none'); colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(log10(abs(real(YY_lmlm__))),llim_); axis image; axisnotick; title('log10(abs(real(YY_lmlm__)))','Interpreter','none'); colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(log10(abs(imag(YY_lmlm__))),llim_); axis image; axisnotick; title('log10(abs(imag(YY_lmlm__)))','Interpreter','none'); colorbar;
end;%if (verbose>1);
%%%%;
disp('returning');return;
end;%if (nargin<1);
%%%%%%%%;

str_thisfunction = 'get_Ylm__';
verbose=0;
if (verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
if (nargin<6); flag_flip=0; end;

flag_calc=0;
if flag_calc;
%%%%%%%%;
% Note that this claculation is unstable for l_max>=88. ;
%%%%%%%%;
Ylm__ = cell(n_l,1);
[polar_a_unique_,ij_unique_,ij_return_] = unique(polar_a_all_(1:n_all));
for nl=1:n_l;
l_val = l_val_(nl);
if (verbose>0); disp(sprintf(' %% nl %d/%d: l_val %d',nl,n_l,l_val)); end;
Ylm__{nl} = zeros(2*l_val+1,n_all);
a1=((2*l_val+1)/(4*pi));
if ( flag_flip); Llm__ = legendre(l_val,cos(1*pi-polar_a_unique_),'unnorm'); end;
if (~flag_flip); Llm__ = legendre(l_val,cos(0*pi+polar_a_unique_),'unnorm'); end;
for m_val=-l_val:+l_val;
if (l_val >0); Llm_ = Llm__(1+abs(m_val),:); end; 
if (l_val==0); Llm_ = Llm__; end;
a2=exp(0.5*lfactorial(l_val-abs(m_val)) - 0.5*lfactorial(l_val+abs(m_val)));
c=sqrt(a1)*a2;
%s=(-1)^((m_val<0)*m_val); % needed to preserve condon-shortley phase. ;
s=1; % original phase ;
Ylm_ = zeros(1,n_all);
if ( flag_flip); Ylm_ = s*c*reshape(Llm_(ij_return_),1,n_all).*reshape(exp(+i*m_val*(1*pi+azimu_b_all_(1:n_all))),1,n_all); end;
if (~flag_flip); Ylm_ = s*c*reshape(Llm_(ij_return_),1,n_all).*reshape(exp(+i*m_val*(0*pi+azimu_b_all_(1:n_all))),1,n_all); end;
Ylm__{nl}(1+l_val+m_val,:) = Ylm_;
end;%for m_val=-l_val:+l_val;
end;%for nl=1:n_l;
end;%if flag_calc;

flag_calc=1;
if flag_calc;
[Ylm__] = get_Ylm__1(n_l,l_val_,n_all,azimu_b_all_,polar_a_all_,flag_flip);
end;%if flag_calc;

if (verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
