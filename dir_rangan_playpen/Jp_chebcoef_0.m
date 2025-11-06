function ...
[ ...
 parameter ...
,Jp_jk__ ...
] = ...
Jp_chebcoef_0( ...
 parameter ...
,n_k ...
);

%%%%%%%%;
% returns the n_k-by-n_k matrix Jp_jk__, ;
% with Jp_jk__(1+nj,1+nk) containing ;
% the tschebyscheff-coefficient of order nk ;
% for jacobi-(0,1)-polynomial of order nj. ;
% The normalization convention for both polynomials is: ;
% T_{k}(+1)=+1, J_{j}(+1)=+1. ;
% This should be stable up to high-order (e.g., past n_k==1024). ;
%%%%%%%%;

str_thisfunction = 'Jp_chebcoef_0';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
flag_verbose=1;
n_k = 1024;
[~,Jp_jk__] = Jp_chebcoef_0([],n_k);
Jq_jk__ = zeros(n_k,n_k);
for nk=0:n_k-1;
Jq_jk__(1+nk,1:nk+1) = reshape(chebcoeffs(jacpoly(nk,0,1)),[1,nk+1]);
end;%for nk=0:n_k-1;
if (flag_verbose>0); disp(sprintf(' %% Jp_jk__: ')); end;
if (flag_verbose>0); disp(Jp_jk__(1:8,1:8)); end;
if (flag_verbose>0); disp(sprintf(' %% Jq_jk__: ')); end;
if (flag_verbose>0); disp(Jq_jk__(1:8,1:8)); end;
for nk=0:n_k-1;
tmp_errrel = fnorm_disp(0*flag_verbose,'Jq_jk__(1:1+nk,1:1+nk)',Jq_jk__(1:1+nk,1:1+nk),'Jp_jk__(1:1+nk,1:1+nk)',Jp_jk__(1:1+nk,1:1+nk),' %%<-- should be small');
if (flag_verbose>0);
if mod(1+nk,8)==0;
disp(sprintf(' %% nk %.3d/%.3d errrel: %0.16f',nk,n_k,tmp_errrel));
end;%if mod(1+nk,8)==0;
end;%if (flag_verbose>0);
end;%for nk=0:n_k-1;
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
% Here we use the combination of two recurrences: ;
% x\cdot T_{k} = \frac{1}{2}\left( T_{k+1} + T_{k-1} \right), ;
% with T_{-1} = T_{1} = x. ;
% and: ;
% (j+1)(2j-1) J_{j} = (2j+1)(2j-1) x J_{j-1} - J_{j-1} - (j-1)(2j+1) J_{j-2}. ;
% from which we assume: ;
% J_{j} = \sum_{k} J_{j,k} T_{k}, ;
% yielding:
% (j+1)(2j-1) \sum_{k} J_{j,k} T_{k} = ... ;
% + (2j+1)(2j-1) \sum_{k} J_{j-1,k} (T_{k+1} + T_{k-1})/2 ... ;
% - \sum_{k} J_{j-1,k} T_{k} ... ;
% - (j-1)(2j+1) \sum_{k} J_{j-2,k}, ;
% from which we have: ;
% (j+1)(2j-1) \sum_{k} J_{j,k} T_{k} = ... ;
% + (2j+1)(2j-1)/2 \sum_{k} J_{j-1,k-1} T_{k} ...
% + (2j+1)(2j-1)/2 \sum_{k} J_{j-1,k+1} T_{k} ...
% - \sum_{k} J_{j-1,k} T_{k} ... ;
% - (j-1)(2j+1) \sum_{k} J_{j-2,k} T_{k}, ;
% from which we deduce: ;
% J_{j,k} = ...
% + (2j+1)/(2j+2) J_{j-1,k-1} ...
% + (2j+1)/(2j+2) J_{j-1,k+1} ...
% - 1/(j+1)/(2j-1) J_{j-1,k} ... ;
% - (j-1)(2j+1)/(j+1)/(2j-1) J_{j-2,k}, ;
% with any references to J_{j,-1} corresponding to J_{j,+1} instead. ;
% Note that this 'reflection' will double-up the contribution of the J_{j-1,0} term ;
% to the accumulation of J_{j,1} (for each j). ;
%%%%%%%%;
% These recurrences are initialized via: ;
% T_{0} = 1, J_{0} = 1;
% T_{1} = x, J_{1} = -1/2 + 3/2 x = -1/2 T_{0} + 3/2 T_{1}. ;
%%%%%%%%;

n_j = n_k;
Jp_jk__ = zeros(n_k,n_k);
for nj=0:n_j-1;
if nj==0; nk=0; Jp_jk__(1+nj,1+nk) = 1; end;
if nj==1; nk=0; Jp_jk__(1+nj,1+nk) = -0.5; nk=1; Jp_jk__(1+nj,1+nk) = +1.5; end;
if nj>=2;
nj1 = nj-1; nj2 = nj-2;
for nk=0:nj;
nkc = nk; nkp = nkc+1; nkn = nkc-1;
J = 0.0;
J = J - 1.0/(nj+1)/(2*nj-1) * Jp_jk__(1+nj1,1+nkc);
J = J - ((nj-1)/(nj+1))*((2*nj+1)/(2*nj-1)) * Jp_jk__(1+nj2,1+nkc);
if nkp<=n_k-1;
J = J + ((2*nj+1)/(2*nj+2)) * Jp_jk__(1+nj1,1+nkp);
end;%if nkp<=n_k-1;
if nkn>=0;
J = J + ((2*nj+1)/(2*nj+2)) * Jp_jk__(1+nj1,1+nkn);
end;%if nkn>=0;
if nk==1;
J = J + ((2*nj+1)/(2*nj+2)) * Jp_jk__(1+nj1,1+nkn); %<-- here is the doubling of the J_{j-1,0} term. ;
end;%if nk==1;
Jp_jk__(1+nj,1+nk) = J;
end;%for nk=0:nj;
end;%if nj>=2;
end;%for nj=0:n_j-1;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
