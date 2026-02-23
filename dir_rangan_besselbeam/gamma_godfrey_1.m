function [log_f_,ori_f_] = gamma_godfrey_1(z_)
% GAMMA  Gamma function valid in the entire complex plane.
%        Accuracy is 15 significant digits along the real axis
%        and 13 significant digits elsewhere.
%        This routine uses a superb Lanczos series
%        approximation for the complex Gamma function.
%
%        z may be complex and of any size.
%        Also  n! = prod(1:n) = gamma(n+1)
%
%usage: [log_f_,ori_f_] = gamma_godfrey_1(z_)
%       
%tested on versions 6.0 and 5.3.1 under Sun Solaris 5.5.1
%
%References: C. Lanczos, SIAM JNA  1, 1964. pp. 86-96
%            Y. Luke, "The Special ... approximations", 1969 pp. 29-31
%            Y. Luke, "Algorithms ... functions", 1977
%            J. Spouge,  SIAM JNA 31, 1994. pp. 931-944
%            W. Press,  "Numerical Recipes"
%            S. Chang, "Computation of special functions", 1996
%            W. J. Cody "An Overview of Software Development for Special
%            Functions", 1975
%
%see also:   GAMMA GAMMALN GAMMAINC PSI
%see also:   mhelp GAMMA
%
%Paul Godfrey
%pgodfrey@intersil.com
%http://winnie.fit.edu/~gabdo/gamma.txt
%Sept 11, 2001

str_thisfunction = 'gamma_godfrey_1';

%%%%%%%%;
if nargin<1;
%%%%%%%%;
disp(sprintf(' %% testing %s',str_thisfunction));
nf=0;
n_z = 128; n_r = 1+n_z; n_i = 0+n_z; %<-- to check dimension. ;
zr_ = linspace(-5,+5,n_r);
zi_ = linspace(-5,+5,n_i);
[zr__,zi__] = ndgrid(zr_,zi_);
z_ = zr__(:) + i*zi__(:);
[bkp_f_] = gamma_godfrey_0(z_);
bkp_f__ = reshape(bkp_f_,[n_r,n_i]);
[log_f_,ori_f_] = gamma_godfrey_1(z_);
ori_f__ = reshape(ori_f_,[n_r,n_i]);
log_f__ = reshape(log_f_,[n_r,n_i]);
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed; figbeach();
p_row = 1; p_col = 3; np=0;
alim_ = prctile(abs(bkp_f_),[ 1,99]); alim_ = mean(alim_) + 1.25*0.5*diff(alim_)*[-1,+1];
tlim_ = [0,2*pi];
subplot_{1+np} = subplot(p_row,p_col,1+np); np=np+1;
imagesc(transpose(abs(bkp_f__)),alim_); set(gca,'ydir','normal'); axisnotick; title('abs(bkp_f__)','Interpreter','none');
subplot_{1+np} = subplot(p_row,p_col,1+np); np=np+1;
imagesc(transpose(abs(ori_f__)),alim_); set(gca,'ydir','normal'); axisnotick; title('abs(ori_f__)','Interpreter','none');
subplot_{1+np} = subplot(p_row,p_col,1+np); np=np+1;
imagesc(transpose(abs(exp(log_f__))),alim_); set(gca,'ydir','normal'); axisnotick; title('abs(exp(log_f__))','Interpreter','none');
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed; colormap('hsv');
p_row = 1; p_col = 3; np=0;
alim_ = prctile(abs(bkp_f_),[ 1,99]); alim_ = mean(alim_) + 1.25*0.5*diff(alim_)*[-1,+1];
tlim_ = [-pi,+pi];
subplot_{1+np} = subplot(p_row,p_col,1+np); np=np+1;
imagesc(transpose(angle(bkp_f__)),tlim_); set(gca,'ydir','normal'); axisnotick; title('angle(bkp_f__)','Interpreter','none');
subplot_{1+np} = subplot(p_row,p_col,1+np); np=np+1;
imagesc(transpose(angle(ori_f__)),tlim_); set(gca,'ydir','normal'); axisnotick; title('angle(ori_f__)','Interpreter','none');
subplot_{1+np} = subplot(p_row,p_col,1+np); np=np+1;
imagesc(transpose(angle(exp(log_f__))),tlim_); set(gca,'ydir','normal'); axisnotick; title('angle(exp(log_f__))','Interpreter','none');
%%%%%%%%;
end;%if nargin<1;
%%%%%%%%;

size_ori_ = size(z_);
z_ = z_(:); n_z = numel(z_);
z_ori_ = z_;
log_f_ = 0.0.*z_; % reserve space in advance
ori_f_ = 1.0.*z_; % reserve space in advance
ij_neg_ = find(real(z_)<0);
if ~isempty(ij_neg_); z_(ij_neg_)=-z_(ij_neg_); end;
% 15 sig. digits for 0<=real(z)<=171
% coeffs should sum to about g*g/2+23/24
g = 607/128; % best results when 4<=g<=5
c_ = [ ...
 + 0.99999999999999709182e-0 ...
;+ 5.71562356658629235170e+1 ...
;- 5.95979603554754912480e+1 ...
;+ 1.41360979747417471740e+1 ...
;- 0.49191381609762019978e-0 ...
;+ 0.33994649984811888699e-4 ...
;+ 0.46523628927048575665e-4 ...
;- 0.98374475304879564677e-4 ...
;+ 0.15808870322491248884e-3 ...
;- 0.21026444172410488319e-3 ...
;+ 0.21743961811521264320e-3 ...
;- 0.16431810653676389022e-3 ...
;+ 0.84418223983852743293e-4 ...
;- 0.26190838401581408670e-4 ...
;+ 0.36899182659531622704e-5 ...
];
%Num Recipes used g=5 with 7 terms
%for a less effective approximation
z_ = z_ - 1;
zh_  = z_ + 0.5;
zgh_ = zh_ + g;
%trick for avoiding FP overflow above z=141
ori_zp_ = zgh_.^(zh_*0.5);
log_zp_ = (zh_*0.5).*log(zgh_);

ss_ = sum(bsxfun(@rdivide,reshape(c_(15:-1:2),[1,14]),bsxfun(@plus,reshape(z_,[n_z,1]),reshape(14:-1:1,[1,14]))),2);

%sqrt(2Pi)
sq2pi =  2.5066282746310005024157652848110;
ori_f_ = (sq2pi*(c_(1)+ss_)).*((ori_zp_.*exp(-zgh_)).*ori_zp_);
ori_f_(z_==0 | z_==1) = 1.0;
log_f_ = log(sq2pi) + log(c_(1)+ss_) + 2*log(ori_zp_) - zgh_ ;
log_f_(z_==0 | z_==1) = 0.0;
%adjust for negative real parts
if ~isempty(ij_neg_);
ori_f_(ij_neg_) = -pi./(z_ori_(ij_neg_).*ori_f_(ij_neg_).*sin(pi*z_ori_(ij_neg_)));
log_f_(ij_neg_) = log(-pi) - log(z_ori_(ij_neg_)) - log_f_(ij_neg_) - log(sin(pi*z_ori_(ij_neg_))) ;
end;
%adjust for negative poles
ij_pol_ = find(round(z_ori_)==z_ori_ & imag(z_ori_)==0 & real(z_ori_)<=0);
if ~isempty(ij_pol_);
ori_f_(ij_pol_)=+Inf;
log_f_(ij_pol_)=+Inf;
end;
ori_f_=reshape(ori_f_,size_ori_);
log_f_=reshape(log_f_,size_ori_);
return;

