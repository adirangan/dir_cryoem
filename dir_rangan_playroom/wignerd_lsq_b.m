function W_ = wignerd_lsq_b(n_l,beta) ;
% generates wigner-d matrices up to n_l;
% uses least-squares solve ;
% test with: 
%{

  n_l=40; beta=pi/6; 
  tic; W1_ = wignerd(n_l,beta);     disp(sprintf(' %% wignerd    : %0.2f seconds',toc));
  tic; W2_ = wignerd_lsq_b(n_l,beta); disp(sprintf(' %% wignerd_lsq_b: %0.2f seconds',toc));  
  for nl=0:n_l;
  disp(sprintf(' %% nl %d/%d: error %0.16f',nl,n_l,norm(W1_{1+nl}-W2_{1+nl})));
  end;%for nl=0:n_l;

  %}

n_oversample = 5;
W_ = cell(1+n_l,1);
W_{1} = [1];
for nl=1:n_l;
n_m = 1+2*nl;
%theta_ = linspace( 0 , 2*pi , 2+2*nl );  % Azimuthal/Longitude/Circumferential ;
%phi_   = linspace( 0 ,   pi , 2+2*nl );  % Altitude /Latitude /Elevation ;
theta_ = sort(2*pi*rand(1,ceil(sqrt(n_oversample*n_m))));
phi_   = sort(1*pi*rand(1,ceil(sqrt(n_oversample*n_m))));
[THETA_,PHI_] = meshgrid(theta_,phi_);
Ylm_orig_ = ylm(nl,THETA_,PHI_);
cb = cos(+beta); sb = sin(+beta); sg = -1;
Xn_ = sin(PHI_).*cos(THETA_);
Yn_ = sin(PHI_).*sin(THETA_);
Zn_ = cos(PHI_);
Xt_ = +cb*Xn_ + sg*sb*Zn_;
Yt_ = Yn_;
Zt_ = -sg*sb*Xn_ + cb*Zn_;
THETA_Y_ = atan2(Yt_,Xt_);
PHI_Y_ = acos(Zt_);
Ylm_rota_ = ylm(nl,THETA_Y_,PHI_Y_);
n_x = length(theta_)*length(phi_);
Y_orig_ = zeros(n_m,n_x);
Y_rota_ = zeros(n_m,n_x);
for nm=1:n_m;
m_val = -nl-1+nm;
s=(-1)^((m_val<0)*m_val); % needed to preserve condon-shortley phase. ;
Y_orig_(nm,:) = s*reshape(Ylm_orig_{nm},1,n_x);
Y_rota_(nm,:) = s*reshape(Ylm_rota_{nm},1,n_x);
end;%for nm=1:n_m;
W_{1+nl} = transpose(Y_rota_ / Y_orig_) ;
end;% for nl=1:n_l;
