function Ylm_ = ylm_1(l_max,THETA_,PHI_) ;
% test with: ;
if (nargin<1);
verbose=1; nf=0;
if (verbose); disp(sprintf(' %% testing ylm_1')); end;
theta_ = linspace( 0 , 2*pi , 48 );  % Azimuthal/Longitude/Circumferential ;
phi_   = linspace( 0 ,   pi , 48 );  % Altitude /Latitude /Elevation ;
[THETA_,PHI_] = meshgrid(theta_,phi_);
l_max = 90;
Ylm_0_ =   ylm(l_max,THETA_,PHI_);
Ylm_1_ = ylm_1(l_max,THETA_,PHI_);
if (verbose>0);
for m_val=-l_max:+l_max;
nl=l_max+m_val;
disp(sprintf(' %% (l,m) (%d,%d): Ylm_0_ vs Ylm_1_: %0.16f',l_max,m_val,fnorm(Ylm_0_{1+nl}-Ylm_1_{1+nl})/fnorm(Ylm_0_{1+nl})));
end;%for m_val=-l_max:+l_max;
end;
if (verbose>1);
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,3,1); imagesc(real(Ylm_0_{1})); title(sprintf('Re Y%d%d',1,-1));
subplot(2,3,4); imagesc(imag(Ylm_0_{1})); title(sprintf('Im Y%d%d',1,-1));
subplot(2,3,2); imagesc(real(Ylm_0_{2})); title(sprintf('Re Y%d%d',1, 0));
subplot(2,3,5); imagesc(imag(Ylm_0_{2})); title(sprintf('Im Y%d%d',1, 0));
subplot(2,3,3); imagesc(real(Ylm_0_{3})); title(sprintf('Re Y%d%d',1,+1));
subplot(2,3,6); imagesc(imag(Ylm_0_{3})); title(sprintf('Im Y%d%d',1,+1));
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,3,1); imagesc(real(Ylm_1_{1})); title(sprintf('Re Y%d%d',1,-1));
subplot(2,3,4); imagesc(imag(Ylm_1_{1})); title(sprintf('Im Y%d%d',1,-1));
subplot(2,3,2); imagesc(real(Ylm_1_{2})); title(sprintf('Re Y%d%d',1, 0));
subplot(2,3,5); imagesc(imag(Ylm_1_{2})); title(sprintf('Im Y%d%d',1, 0));
subplot(2,3,3); imagesc(real(Ylm_1_{3})); title(sprintf('Re Y%d%d',1,+1));
subplot(2,3,6); imagesc(imag(Ylm_1_{3})); title(sprintf('Im Y%d%d',1,+1));
end;%if (verbose>1);
disp('returning'); return;
end;%if (nargin<1);

n_lm = 1+2*l_max;
Llm_ = cell(n_lm,1);
Ylm_ = cell(n_lm,1);

Llm_jlm___ = ylgndr_1(l_max,cos(PHI_(:)));
Llm_jm__ = squeeze(Llm_jlm___(:,1+l_max,:));
%l_val = l_max;
%m_val_ = transpose(0:l_val);
%tmp_a1 = ((1+2*l_val)/(4*pi));
%tmp_a2_ = exp(lfactorial(l_val-abs(m_val_)) - lfactorial(l_val+abs(m_val_)));
%tmp_a3_ = sqrt(tmp_a1*tmp_a2_);
%Llm_jm__ = bsxfun(@rdivide,Llm_jm__,reshape(sqrt(4*pi)*tmp_a3_,[1,1+l_max]));
Llm_mba___ = reshape(permute(Llm_jm__,[2,1]),[1+l_max,size(THETA_)])/sqrt(4*pi);

for m_val = -l_max:+l_max;
ix = 1+l_max+m_val;
m_abs = abs(m_val);
if (l_max>0); Llm_{ix} = squeeze(Llm_mba___(1+m_abs,:,:)); end;
if (l_max==0); Llm_{ix} = Llm_mba___(:,:); end;
c = 1;
%a1=((2*l_max+1)/(4*pi));
%a2=exp(lfactorial(l_max-m_abs) - lfactorial(l_max+m_abs));
%c=sqrt(a1*a2);
s=(-1)^((m_val<0)*m_val); % needed to preserve condon-shortley phase. ;
Ylm_{ix}=s*c*Llm_{ix}.*exp(+i*m_val*THETA_);
end;%for m_val = -l_max:+l_max;
