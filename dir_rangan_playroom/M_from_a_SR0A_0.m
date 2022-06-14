function ...
[ M_from_a_pMq___ ] = ...
M_from_a_SR0A_0( ...
 n_w_max ...
,gamma_z_M_ ...
,n_pick ...
,dgamma_z_pM__ ...
);
na=0;
if (nargin<1+na); n_w_max=[]; end; na=na+1;
if (nargin<1+na); gamma_z_M_=[]; end; na=na+1;
if (nargin<1+na); n_pick=[]; end; na=na+1;
if (nargin<1+na); dgamma_z_pM__=[]; end; na=na+1;

if isempty(n_pick); n_pick=1; end;
if isempty(dgamma_z_pM__); dgamma_z_pM__ = zeros(1,numel(gamma_z_M_)); end;

n_M = numel(gamma_z_M_);
gamma_z_M_ = periodize(reshape(gamma_z_M_,[n_M,1]),0,2*pi);
for npick=0:n_pick-1;
gamma_z_m_ = gamma_z_M_ + reshape(dgamma_z_pM__(1+npick,:),[n_M,1]);
M_from_a_pMq___(1+npick,:,:) = m_from_a_SR0A_0(n_w_max,gamma_z_m_);
end;%for npick=0:n_pick-1;

%%%%%%%%;
% single npick. ;
%%%%%%%%;
function ...
[ m_from_a_mq__ ] = ...
m_from_a_SR0A_0( ...
 n_w_max ...
,gamma_z_m_ ...
);
n_m = numel(gamma_z_m_);
gamma_z_m_ = periodize(reshape(gamma_z_m_,[n_m,1]),0,2*pi);
q_ = periodize([0:n_w_max-1],-n_w_max/2,+n_w_max/2);
m_from_a_mq__ = exp(+i*(gamma_z_m_*q_)) / sqrt(n_w_max);





