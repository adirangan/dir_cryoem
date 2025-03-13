function ...
[ ...
 parameter ...
,a_k_p_qk_ ...
] = ...
a_k_p_from_plane_wave_expansion_0( ...
 parameter ...
,n_qk ...
,k_c_0_qk_ ...
,k_c_1_qk_ ...
,k_c_2_qk_ ...
,n_source_a ...
,v_source_a_ ...
,delta_a_c__ ...
);

str_thisfunction = 'a_k_p_from_plane_wave_expansion_0';

if nargin<1;
disp(sprintf(' %% to test %s see test_ssnll_from_plane_wave_expansion_2',str_thisfunction));
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_qk=[]; end; na=na+1;
if (nargin<1+na); k_c_0_qk_=[]; end; na=na+1;
if (nargin<1+na); k_c_1_qk_=[]; end; na=na+1;
if (nargin<1+na); k_c_2_qk_=[]; end; na=na+1;
if (nargin<1+na); n_source_a=[]; end; na=na+1;
if (nargin<1+na); v_source_a_=[]; end; na=na+1;
if (nargin<1+na); delta_a_c__=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

a_k_p_qk_ = zeros(n_qk,1);
for nsource_a=0:n_source_a-1;
delta_a_c_ = delta_a_c__(:,1+nsource_a);
v_source_a = v_source_a_(1+nsource_a);
a_k_p_qk_ = a_k_p_qk_ + v_source_a*exp(+i*2*pi*(k_c_0_qk_*delta_a_c_(1+0) + k_c_1_qk_*delta_a_c_(1+1) + k_c_2_qk_*delta_a_c_(1+2)));
end;%for nsource_a=0:n_source_a-1;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

