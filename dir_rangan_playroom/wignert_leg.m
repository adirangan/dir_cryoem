function [W_,pp_,m2_,l2_] = wignert_leg(n_l,delta_,n_s) ;
% generates wigner-t matrices up to n_l;
% uses legendre quadrature to solve ;
% test with: ;
%{
  [W1_,pp_,m2_,l2_] = wignert_leg(10,[0,0,1],3);
  [W2_,pp_,m2_,l2_] = wignert_lsq(10,[0,0,1],3);
  imagesc(log10(abs(W1_-W2_)),[-2,0]); colorbar;

  % Note that an oversampling parameter of n_s=3 or so ;
  % is required to achieve 9 digits of precision ;
  % with respect to the 2-norm of wignert. ;
  n_l = 20; delta_ = [0,0,1];
  n_s_ = [0.25:0.25:3.0];
  E_ = zeros(length(n_s_),1);
  W1_ = cell(length(n_s_),1);
  for ns=1:length(n_s_);
  n_s = n_s_(ns);
  disp(sprintf(' %% n_s %d',n_s));
  W1_{ns} = wignert_leg(n_l,delta_,n_s);
  end;%for ns=1:length(n_s_);
  for ns=1:length(n_s_);
  E_(ns) = norm(W1_{ns}-W1_{length(n_s_)});
  end;%for ns=1:length(n_s_);
  plot(n_s_(1:length(n_s_)-1),log10(E_((1:length(n_s_)-1))),'.-');
  xlabel('n_s');ylabel('log10(E_s)');title('error');

  %}

verbose=0;

if nargin<3; n_s = 3; end; % n_s = oversampling parameter ;

n_lm = (n_l+1).^2;
k = 1.0; sample_d = 1.0/(n_s*max(1,n_l));
[length,theta_,phi_,weight_] = sample_shell_0(k,sample_d,'L');
l_ = []; m_ = []; for nl=0:n_l; l_ = [l_ , nl*ones(1,2*nl+1) ]; m_ = [m_ , [-nl:+nl] ]; end;%for nl=0:n_l;
A_ = zeros(n_lm,length);
for nl=0:n_l;
l_val = nl;
Llm__=legendre(l_val,cos(phi_),'unnorm');
for m_val = -l_val:+l_val;
ix = 1+l_val*(l_val+1)+m_val;
m_abs = abs(m_val);
if (length>0); 
if (l_val>0); Llm_ = squeeze(Llm__(1+m_abs,:,:)); end; 
if (l_val==0); Llm_ = Llm__(:,:); end; 
end; %if (length>0); 
a1=((2*l_val+1)/(4*pi));
a2=exp(lfactorial(l_val-m_abs) - lfactorial(l_val+m_abs));
c=sqrt(a1*a2);
s=(-1)^((m_val<0)*m_val); % needed to preserve condon-shortley phase. ;
%s=1; % original phase ;
if (length>0); 
Ylm_ = s*c*Llm_.*exp(+i*m_val*transpose(theta_)); 
A_(ix,:) = s*Ylm_;
end;%if (length>0); 
end;%for m_val = -l_val:+l_val;
end;%for nl=0:n_l;

kx_ = k .* cos(theta_) .* sin(phi_);
ky_ = k .* sin(theta_) .* sin(phi_);
kz_ = k .* cos(phi_);
if (verbose); disp(sprintf(' %% translating by [%0.2f %0.2f %0.2f]',delta_)); end;
kd_ = kx_*delta_(1) + ky_*delta_(2) + kz_*delta_(3);
%B_ = A_ * sparse(1:length,1:length,exp(i*kd_));
BW_ = A_ * sparse(1:length,1:length,exp(i*kd_).*weight_);

W_ = conj(A_)*transpose(BW_);
%hist(log10(abs(W_(:))),[-20:0]);

pp_ = [];
m2_ = [];
l2_ = [];
for nm=0:n_l;
l1_ = nm:n_l;
ll_ = l1_.*(l1_+1);
if (nm==0); pp_ = [pp_,ll_]; m2_ = [m2_,zeros(1,1+n_l)]; l2_ = [l2_,l1_]; end;
if (nm>0); 
%pp_tmp = [ll_-nm ; ll_+nm]; pp_tmp = reshape(pp_tmp,1,2*(n_l-nm+1));
%m2_tmp = [-nm*ones(1,n_l-nm+1) ; +nm*ones(1,n_l-nm+1)]; m2_tmp = reshape(m2_tmp,1,2*(n_l-nm+1));
%l2_tmp = [l1_ ; l1_]; l2_tmp = reshape(l2_tmp,1,2*(n_l-nm+1));
pp_tmp = [ll_-nm ; ll_+nm]; pp_tmp = reshape(transpose(pp_tmp),1,2*(n_l-nm+1));
m2_tmp = [-nm*ones(1,n_l-nm+1) ; +nm*ones(1,n_l-nm+1)]; m2_tmp = reshape(transpose(m2_tmp),1,2*(n_l-nm+1));
l2_tmp = [l1_ ; l1_]; l2_tmp = reshape(transpose(l2_tmp),1,2*(n_l-nm+1));
pp_ = [pp_,pp_tmp]; m2_ = [m2_,m2_tmp]; l2_ = [l2_,l2_tmp];
end;%if (nm>0); 
end;%for nm=0:n_l;

plot_flag=0;
if plot_flag;
figure;
subplot(1,2,1);
imagesc(abs(W_),[0,1]); colorbar;
set(gca,'XTick',1:n_lm,'XTickLabel',m_);
set(gca,'YTick',1:n_lm,'YTickLabel',m_);
title(sprintf('m,l'));
subplot(1,2,2);
imagesc(abs(W_(1+pp_,1+pp_)),[0,1]); colorbar;
set(gca,'XTick',1:n_lm,'XTickLabel',m2_);
set(gca,'YTick',1:n_lm,'YTickLabel',m2_);
title(sprintf('l,m'));
figure;
subplot(1,2,1);
imagesc(log10(abs(W_)),[-2,0]); colorbar;
set(gca,'XTick',1:n_lm,'XTickLabel',m_);
set(gca,'YTick',1:n_lm,'YTickLabel',m_);
title(sprintf('m,l'));
subplot(1,2,2);
imagesc(log10(abs(W_(1+pp_,1+pp_))),[-2,0]); colorbar;
set(gca,'XTick',1:n_lm,'XTickLabel',m2_);
set(gca,'YTick',1:n_lm,'YTickLabel',m2_);
title(sprintf('l,m'));
end;%if plot_flag;
