function Ylm_ = ylm(l_val,THETA_,PHI_) ;
% test with: ;
%{
  theta_ = linspace( 0 , 2*pi , 48 );  % Azimuthal/Longitude/Circumferential ;
  phi_   = linspace( 0 ,   pi , 48 );  % Altitude /Latitude /Elevation ;
  [THETA_,PHI_] = meshgrid(theta_,phi_);
  Ylm_ = ylm(1,THETA_,PHI_);
  subplot(2,3,1); imagesc(real(Ylm_{1})); title(sprintf('Re Y%d%d',1,-1));
  subplot(2,3,4); imagesc(imag(Ylm_{1})); title(sprintf('Im Y%d%d',1,-1));
  subplot(2,3,2); imagesc(real(Ylm_{2})); title(sprintf('Re Y%d%d',1, 0));
  subplot(2,3,5); imagesc(imag(Ylm_{2})); title(sprintf('Im Y%d%d',1, 0));
  subplot(2,3,3); imagesc(real(Ylm_{3})); title(sprintf('Re Y%d%d',1,+1));
  subplot(2,3,6); imagesc(imag(Ylm_{3})); title(sprintf('Im Y%d%d',1,+1));
  %}

n_lm = 1+2*l_val;
Llm_ = cell(n_lm,1);
Ylm_ = cell(n_lm,1);

Llm_tmp=legendre(l_val,cos(PHI_),'unnorm');
for m_val = -l_val:+l_val;
ix = 1+l_val+m_val;
m_abs = abs(m_val);
if (l_val>0); Llm_{ix} = squeeze(Llm_tmp(1+m_abs,:,:)); end;
if (l_val==0); Llm_{ix} = Llm_tmp(:,:); end;
a1=((2*l_val+1)/(4*pi));
a2=exp(lfactorial(l_val-m_abs) - lfactorial(l_val+m_abs));
c=sqrt(a1*a2);
s=(-1)^((m_val<0)*m_val); % needed to preserve condon-shortley phase. ;
Ylm_{ix}=s*c*Llm_{ix}.*exp(+i*m_val*THETA_);
end;%for m_val = -l_val:+l_val;

