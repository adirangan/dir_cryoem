function [B_wk__,s_w_] = MSA_shape_speed1_1(n_k,n_w,gamma_w_,A_wk__,weight_k_);
%%%%%%%%;
% Reparametrizes a single (complex) image-ring by arclength. ;
% Assumes that gamma_w_ is periodic. ;
%%%%%%%%;

if nargin<1;
n_k = 2;
n_w = 1024;
gamma_w_ = linspace(0,2*pi,1+n_w); gamma_w_ = transpose(gamma_w_(1:n_w));
A_wk__ = [ sin(gamma_w_) + 0.125*cos(3*gamma_w_) , 1.3*sin(2*gamma_w_) - 0.5*cos(7*gamma_w_) ];
[B_wk__,s_w_] = MSA_shape_speed1_1(n_k,n_w,gamma_w_,A_wk__);
if n_k==1;
figure(1);clf;figmed;
subplot(1,2,1); plot(gamma_w_,A_wk__,'b-'); xlabel('gamma'); title('orig');
subplot(1,2,2); plot(s_w_,A_wk__,'b-',gamma_w_,B_wk__,'ro'); xlabel('gamma'); title('speed1');
end;%if n_k==1;
if n_k==2;
figure(1);clf;figbig;colormap('hsv');
subplot(2,2,1); surfline_0(A_wk__(:,1),A_wk__(:,2),gamma_w_); title('orig');
subplot(2,2,2); surfline_0(B_wk__(:,1),B_wk__(:,2),gamma_w_); title('speed1');
dA_wk__ = [A_wk__(2:end,:);A_wk__(1+0,:)] - [A_wk__(end,:);A_wk__(1:end-1,:)];
dA_w_ = max(1e-12,sqrt(sum(abs(dA_wk__).^2,2)));
dB_wk__ = [B_wk__(2:end,:);B_wk__(1+0,:)] - [B_wk__(end,:);B_wk__(1:end-1,:)];
dB_w_ = max(1e-12,sqrt(sum(abs(dB_wk__).^2,2)));
dA_lim_ = [min(dA_w_);max(dA_w_)]; dA_lim_ = mean(dA_lim_) + 1.25*0.5*diff(dA_lim_)*[-1,+1];
subplot(2,2,3); plot(dA_w_,'k-'); ylim(dA_lim_); title(sprintf('dA_w_','Interpreter','none'));
subplot(2,2,4); plot(dB_w_,'k-'); ylim(dA_lim_); title(sprintf('dB_w_','Interpreter','none'));
end;%if n_k==2;
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); n_k=[]; end; na=na+1;
if (nargin<1+na); n_w=[]; end; na=na+1;
if (nargin<1+na); gamma_w_=[]; end; na=na+1;
if (nargin<1+na); A_wk__=[]; end; na=na+1;
if (nargin<1+na); weight_k_=[]; end; na=na+1;

if isempty(gamma_w_); gamma_w_ = linspace(0,2*pi,1+n_w); gamma_w_ = transpose(gamma_w_(1:n_w)); end;
if isempty(weight_k_); weight_k_ = ones(n_k,1); end;

dgamma_w_ = [gamma_w_(2:end);gamma_w_(1+0)] - [gamma_w_(end);gamma_w_(1:end-1)];
dA_wk__ = [A_wk__(2:end,:);A_wk__(1+0,:)] - [A_wk__(end,:);A_wk__(1:end-1,:)];
dA_w_ = max(1e-12,sqrt(sum(bsxfun(@times,abs(dA_wk__).^2,reshape(weight_k_,[1,n_k])),2)));
s_w_ = cumsum(dA_w_);
s_w_ = 2*pi*s_w_/max(1e-12,s_w_(end));
B_wk__ = interp1(s_w_,A_wk__,transpose(linspace(s_w_(1+0),s_w_(end),n_w)));



