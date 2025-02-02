function ...
[ ...
 parameter ...
,viewing_R_polar_a_S_ ...
,viewing_R_azimu_b_S_ ...
,viewing_R_gamma_z_S_ ...
] = ...
euler_to_euler_0( ...
 parameter ...
,n_viewing_S ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_gamma_z_S_ ...
,R_inp__ ...
);

%%%%%%%%;
% Applies a transformation R_use__ to the euler-angles. ;
% The convention is that euler_ is [-gamma_z,+polar_a,+azimu_b]. ;
% Note that, by default, the sign of viewing_gamma_z is flipped. ;
% Also, by default, R_use__ := transpose(R_inp__). ;
%%%%%%%%;

str_thisfunction = 'euler_to_euler_0';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));

disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_viewing_S=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_azimu_b_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_gamma_z_S_=[]; end; na=na+1;
if (nargin<1+na); R_inp__=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_transpose'); parameter.flag_transpose=1; end;
flag_transpose=parameter.flag_transpose;
if ~isfield(parameter,'sign_polar_a'); parameter.sign_polar_a=+1; end;
sign_polar_a=parameter.sign_polar_a;
if ~isfield(parameter,'sign_azimu_b'); parameter.sign_azimu_b=+1; end;
sign_azimu_b=parameter.sign_azimu_b;
if ~isfield(parameter,'sign_gamma_z'); parameter.sign_gamma_z=-1; end;
sign_gamma_z=parameter.sign_gamma_z;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(viewing_polar_a_S_); viewing_polar_a_S_=zeros(n_S,1); end;
if isempty(viewing_azimu_b_S_); viewing_azimu_b_S_=zeros(n_S,1); end;
if isempty(viewing_gamma_z_S_); viewing_gamma_z_S_=zeros(n_S,1); end;

R_use__ = R_inp__;
if flag_transpose==1; R_use__ = transpose(R_inp__); end;

n_S = n_viewing_S;
viewing_R_polar_a_S_ = zeros(n_viewing_S,1);
viewing_R_azimu_b_S_ = zeros(n_viewing_S,1);
viewing_R_gamma_z_S_ = zeros(n_viewing_S,1);
for nS=0:n_S-1;
viewing_polar_a_sub = viewing_polar_a_S_(1+nS);
viewing_azimu_b_sub = viewing_azimu_b_S_(1+nS);
viewing_gamma_z_sub = viewing_gamma_z_S_(1+nS);
euler_sub_ = [ ...
 sign_gamma_z*viewing_gamma_z_sub ...
,sign_polar_a*viewing_polar_a_sub ...
,sign_azimu_b*viewing_azimu_b_sub ...
]; %<-- note sign. ;
R_sub__ = euler_to_R_0(euler_sub_);
R_rot__ = R_use__*R_sub__;
euler_rot_ = R_to_euler_0(R_rot__);
viewing_gamma_z_rot = euler_rot_(1+0);
viewing_polar_a_rot = euler_rot_(1+1);
viewing_azimu_b_rot = euler_rot_(1+2);
viewing_R_polar_a_S_(1+nS) = sign_polar_a*viewing_polar_a_rot;
viewing_R_azimu_b_S_(1+nS) = sign_azimu_b*viewing_azimu_b_rot;
viewing_R_gamma_z_S_(1+nS) = sign_gamma_z*viewing_gamma_z_rot; %<-- note sign. ;
end;%for nS=0:n_S-1;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

