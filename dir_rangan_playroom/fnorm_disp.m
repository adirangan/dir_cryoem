function fnorm_disp(flag_verbose,str_v0,v0,str_v1,v1,str_postfix);
na=0;
if (nargin<1+na); flag_verbose=[]; end; na=na+1;
if (nargin<1+na); str_v0=[]; end; na=na+1;
if (nargin<1+na); v0=[]; end; na=na+1;
if (nargin<1+na); str_v1=[]; end; na=na+1;
if (nargin<1+na); v1=[]; end; na=na+1;
if (nargin<1+na); str_postfix=[]; end; na=na+1;
if isempty(str_postfix); str_postfix=''; end;

if ndims(v0)~=ndims(v1);
disp(sprintf(' %% Warning, %s ndims(v0) %d ~= %s ndims(v1) %d',str_v0,ndims(v0),str_v1,ndims(v1)));
end;%if ndims(v0)~=ndims(v1);
if numel(v0)~=numel(v1);
disp(sprintf(' %% Warning, %s numel(v0) %d ~= %s numel(v1) %d',str_v0,numel(v0),str_v1,numel(v1)));
end;%if numel(v0)~=numel(v1);
d0_ = size(v0); n_d0 = numel(d0_);
d1_ = size(v1); n_d1 = numel(d1_);
for nd0=0:n_d0-1;
nd1 = nd0;
d0 = d0_(1+nd0); d1 = d1_(1+nd1);
if (d0~=d1);
disp(sprintf(' %% Warning, size(%s,1+%d) %d ~= size(%s,1+%d) %d',str_v0,nd0,d0,str_v1,nd1,d1));
end;%if (d0~=d1);
end;%for nd0=0:n_d0-1;
if (flag_verbose>0);
disp(sprintf(' %% %s %+0.16f vs %s %+0.16f: r %0.16f%s',str_v0,fnorm(v0),str_v1,fnorm(v1),fnorm(v0-v1)/max(1e-12,fnorm(v0)),str_postfix));
end;%if (flag_verbose>0);
