function ...
[ ...
 parameter ...
,E_x_ ...
] = ...
max_filter_1d_0( ...
 parameter ...
,n_x ...
,A_x_ ...
,n_w ... 
);
%%%%%%%%;
% Just for testing purposes. ;
% Use imdilate instead. ;
%%%%%%%%;

str_thisfunction = 'max_filter_1d_0';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
rng(0);
n_x = 8; A_x_ = round(128*randn(n_x,1));
%%%%;
for flag_periodic=[0,1];
for str_strategy_={'loop+max0','loop+max1'};
str_strategy = str_strategy_{1};
disp(sprintf(' %% flag_periodic %d; str_strategy %s',flag_periodic,str_strategy));
tmp_parameter = struct('type','parameter','flag_periodic',flag_periodic,'str_strategy',str_strategy);
n_w = 0; disp(sprintf(' %% %% n_w %d A_x_: %s',n_w,num2str(transpose(A_x_(:)),' %+.3d')));
for n_w = [1,2];
[~,E_x_] = max_filter_1d_0(tmp_parameter,n_x,A_x_,n_w);
disp(sprintf(' %% %% n_w %d E_x_: %s',n_w,num2str(transpose(E_x_(:)),' %+.3d')));
end;%for n_w = [1,2];
end;%for str_strategy_={'loop+max0','loop+max1'};
end;%for flag_periodic=[0,1];
%%%%;
for str_strategy_={'loop+max0','loop+max1'};
str_strategy = str_strategy_{1};
disp(sprintf(' %% str_strategy: %s',str_strategy));
tmp_parameter = struct('type','parameter','str_strategy',str_strategy);
n_x = 1024*32; %A_x_ = round(128*randn(n_x,1));
A_x_ = 128*randn(n_x,1);
n_w = 128;%n_w = 512;
tmp_t = tic(); [~,E_x_] = max_filter_1d_0(tmp_parameter,n_x,A_x_,n_w); tmp_t = toc(tmp_t);
disp(sprintf(' %% %% single step: %0.4fs',tmp_t));
tmp_t = tic(); 
F_x_ = A_x_;
for nw=0:n_w-1;
[~,F_x_] = max_filter_1d_0(tmp_parameter,n_x,F_x_);
end;%for nw=0:n_w-1;
tmp_t = toc(tmp_t);
disp(sprintf(' %% %% iterative: %0.4fs',tmp_t));
tmp_parameter = struct('type','parameter','str_strategy','imdilate');
tmp_t = tic(); [~,G_x_] = max_filter_1d_0(tmp_parameter,n_x,A_x_,n_w); tmp_t = toc(tmp_t); disp(sprintf(' %% %% imdilate: %0.4fs',tmp_t));
disp(sprintf(' %% %% E_x_ vs F_x_: %0.16f',fnorm(E_x_-F_x_)/fnorm(E_x_)));
disp(sprintf(' %% %% E_x_ vs G_x_: %0.16f',fnorm(E_x_-G_x_)/fnorm(E_x_)));
end;%for str_strategy_={'loop+max0','loop+max1'};
%%%%;
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); A_x_=[]; end; na=na+1;
if (nargin<1+na); n_w=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_periodic'); parameter.flag_periodic=0; end;
flag_periodic=parameter.flag_periodic;
if ~isfield(parameter,'str_strategy'); parameter.str_strategy='imdilate'; end;
str_strategy=parameter.str_strategy;

if isempty(n_w); n_w = 1; end;

if flag_periodic==0;
A_min = min(A_x_,[],'all');
B_x_ = [A_min*ones(n_w,1) ; reshape(A_x_,[n_x,1]) ; A_min*ones(n_w,1)];
parameter_periodic = parameter; parameter_periodic.flag_periodic = 1;
[parameter_periodic,D_x_] = max_filter_1d_0(parameter_periodic,2*n_w+n_x,B_x_,n_w);
E_x_ = reshape(D_x_(n_w+1:n_w+n_x),size(A_x_));
parameter = parameter_periodic; parameter.flag_periodic = 0;
end;%if flag_periodic==0;

if flag_periodic==1;
if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
A_min = min(A_x_,[],'all');
%%%%%%%%;
if strcmp(str_strategy,'imdilate');
E_x_ = reshape(imdilate(reshape(A_x_,[n_x,1]),ones(1+2*n_w,1)),size(A_x_));
end;%if strcmp(str_strategy,'imdilate');
%%%%%%%%;
if strcmp(str_strategy,'loop+max0');
B_x_ = [reshape(A_x_,[1,n_x])];
n_s = 1+2*n_w;
C_wx__ = zeros(n_s,n_x);
for ns=0:n_s-1;
nw = ns-n_w;
C_wx__(1+ns,:) = circshift(B_x_,nw);
end;%for ns=0:n_s-1;
D_x_ = max(C_wx__,[],1);
E_x_ = reshape(D_x_,size(A_x_));
end;%if strcmp(str_strategy,'loop+max0');
%%%%%%%%;
if strcmp(str_strategy,'loop+max1');
B_x_ = [reshape(A_x_,[n_x,1])];
n_s = 1+2*n_w;
C_xw__ = zeros(n_x,n_s);
for ns=0:n_s-1;
nw = ns-n_w;
C_xw__(:,1+ns) = circshift(B_x_,nw);
end;%for ns=0:n_s-1;
D_x_ = max(C_xw__,[],2);
E_x_ = reshape(D_x_,size(A_x_));
end;%if strcmp(str_strategy,'loop+max1');
%%%%%%%%;
end;%if flag_periodic==1;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
