function ...
[ ...
 parameter ...
,pinv_H___ ...
] = ...
pagepinv_0( ...
 parameter ...
,H___ ...
);

str_thisfunction = 'pagepinv_0';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
disp(sprintf(' %% (see test_interpolate_template_8.m)',str_thisfunction));
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); H___=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'tolerance_pinv'); parameter.tolerance_pinv=1e-6; end;
tolerance_pinv=parameter.tolerance_pinv;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

[U_n_H___,S_n_H___,V_n_H___] = pagesvd(H___);
U_t_H___ = pagetranspose(U_n_H___);
%V_t_H___ = pagetranspose(V_n_H___); %<-- unnecessary. ;
S_i_H___ = min(1./tolerance_pinv,pageinv(S_n_H___));
S_i_H___(isinf(S_i_H___)|isnan(S_i_H___))=0;
pinv_H___ = pagemtimes(V_n_H___,pagemtimes(S_i_H___,U_t_H___));

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;


