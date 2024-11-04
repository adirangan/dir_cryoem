function fnorm_disp(flag_verbose,str_v0,v0,str_v1,v1);
if numel(v0)~=numel(v1);
disp(sprintf(' %% Warning, numel(v0) %d ~= numel(v1) %d',numel(v0),numel(v1)));
end;%if numel(v0)~=numel(v1);
if (flag_verbose>0);
disp(sprintf(' %% %s %+0.16f vs %s %+0.16f: r %0.16f',str_v0,fnorm(v0),str_v1,fnorm(v1),fnorm(v0-v1)/max(1e-12,fnorm(v0))));
end;%if (flag_verbose>0);
