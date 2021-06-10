function [W_] = wignerz_leg(n_l,delta_,n_sample) ;
% generates wigner-t matrices up to n_l for a variety of delta;
% assumes each delta_vector_ is given by: [0,0,delta_(ndelta)];
% uses legendre quadrature to solve ;
% test with: ;
%{

  % Note that an oversampling parameter of n_sample==0.5 or so ;
  % is required to achieve high precision at n_l==25 ;
  % with respect to the 2-norm of wignerz. ;
  % Note that as n_l increases the oversampling parameter should decrease. ;
  % For example, at n_l==35 we only really need n_sample==0.4 or so. ;
  n_l = 15; n_m = 1+2*n_l; delta_ = linspace(0,4,9);
  n_sample_ = [0.25:0.05:1.0];
  W__ = cell(length(n_sample_),1);
  for nsample=1:length(n_sample_);
  n_sample = n_sample_(nsample);
  disp(sprintf(' %% n_sample %d',n_sample));
  W__{nsample} = wignerz_leg(n_l,delta_,n_sample);
  end;%for nsample=1:length(n_sample_);
  E_sum_ = zeros(length(n_sample_),1);
  E_sub_ = zeros(n_m,length(delta_),length(n_sample_));
  for nsample=1:length(n_sample_);
  E_sum_(nsample) = 0;
  for nm=1:n_m; for nd=1:length(delta_);
  E_sub_(nm,nd,nsample) = norm(W__{nsample}{nm,nd}-W__{length(n_sample_)}{nm,nd});
  E_sum_(nsample) = E_sum_(nsample) + E_sub_(nm,nd,nsample);
  end;end;%for nm=1:n_m; for nd=1:length(delta_);
  end;%for nsample=1:length(n_sample_);
  figure;
  prows=4;pcols=4;
  for nsample=1:length(n_sample_)-1;
  subplot(prows,pcols,nsample);
  imagesc(log10(squeeze(E_sub_(:,:,nsample))),[-12,0]); colorbar;
  xlabel('delta'); ylabel('nm'); title(sprintf(' ns(%d)-->%0.2f',nsample,n_sample_(nsample)));
  end;%for nsample=1:length(n_sample_)-1;
  figure;
  plot(n_sample_(1:length(n_sample_)-1),log10(E_sum_((1:length(n_sample_)-1))),'.-');
  xlabel('n_sample');ylabel('log10(E_sum_)');title('error');  

  %}

verbose=0;

if nargin<3; n_sample = 3; end; % n_sample = oversampling parameter ;

n_m = 1+2*n_l;
n_delta = length(delta_);
W_ = cell(n_m,n_delta);

k = 1.0; sample_d = 1.0/(n_sample*max(1,n_l));
[length_sample,theta_,phi_,weight_] = sample_shell_0(k,sample_d,'L');
kx_ = k .* cos(theta_) .* sin(phi_);
ky_ = k .* sin(theta_) .* sin(phi_);
kz_ = k .* cos(phi_);

[m1_,l1_,m2_,l2_,pp_,qq_] = permute_ml_to_lm(n_l);
m3_ = unique(m2_,'stable');

Llm_all__ = cell(1+n_l,1);
for nl=0:n_l;
l_val = nl;
Llm_all__{1+nl} = legendre(l_val,cos(phi_),'unnorm');
end;%for nl=0:n_l;

E2_ = cell(n_delta,1);
for ndelta=1:n_delta;
%kd_ = kx_*delta_vector_(1) + ky_*delta_vector_(2) + kz_*delta_vector_(3);
kd_ = kz_*delta_(ndelta);
E2_{ndelta} = sparse(1:length_sample,1:length_sample,exp(i*kd_).*weight_);
end;%for ndelta=1:n_delta;

for nm=1:n_m;
%if nm==1; m_val=0; end; if nm>1; m_val = (1+floor((nm-2)/2))*((-1)^(mod(nm-1,2))); end;
m_val = m3_(nm); m_abs = abs(m_val);
E1_ = exp(+i*m_val*transpose(theta_)); 
z_ij = find(m2_==m_val); n_z = length(z_ij);
A_ = zeros(n_z,length_sample);
for nz=1:n_z;
l_val = l2_(z_ij(nz));
Llm__ = Llm_all__{1+l_val};
if (l_val>0); Llm_ = squeeze(Llm__(1+m_abs,:,:)); end; 
if (l_val==0); Llm_ = Llm__(:,:); end; 
a1=((2*l_val+1)/(4*pi));
a2=exp(lfactorial(l_val-m_abs) - lfactorial(l_val+m_abs));
c=sqrt(a1*a2);
s=(-1)^((m_val<0)*m_val); % needed to preserve condon-shortley phase. ;
%s=1; % original phase ;
Ylm_ = s*c*Llm_.*E1_;
A_(nz,:) = s*Ylm_;
end;%for nz=1:n_z;
for ndelta=1:n_delta;
BW_ = A_ * E2_{ndelta};
W_{nm,ndelta} = conj(A_)*transpose(BW_);
end;%for ndelta=1:n_delta;
end;%for nm=1:n_m;

plot_flag=0;
if plot_flag;
prows = n_delta;
pcols = n_l+1;
for nr = 1:n_delta; for nc = 1:n_l+1;
subplot(prows,pcols,nc + (nr-1)*pcols);
m_val = nc-1; d_val = delta_(nr);
%imagesc(abs(W_{1+2*(nc-1),nr}),[0,1]); 
imagesc(log10(abs(W_{1+2*(nc-1),nr})),[-5,0]); 
%title(sprintf('m%d d%0.2f',m_val,d_val));
text(0.5,n_l+2.5-nc,sprintf('m%d d%0.2f',m_val,d_val));
set(gca,'XTick',[],'YTick',[]);
end;end;%for nr = 1:n_delta; for nc = 1:n_l+1;
end;%if plot_flag;
