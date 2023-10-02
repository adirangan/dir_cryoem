function ...
[ ...
 parameter ...
,E_xy__ ...
] = ...
max_filter_2d_0( ...
 parameter ...
,n_x ...
,n_y ...
,A_xy__ ...
,n_w_x ... 
,n_w_y ... 
);
%%%%%%%%;
% Just for testing purposes. ;
% Use imdilate instead. ;
%%%%%%%%;
str_thisfunction = 'max_filter_2d_0';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
rng(0);
n_x = 16; n_y = 12; A_xy__ = round(128*randn(n_x,n_y));
%%%%;
for flag_periodic_x=[0,1];
for flag_periodic_y=[0,1];
for str_strategy_={'loop+max1'};
str_strategy = str_strategy_{1};
%disp(sprintf(' %% flag_periodic %d,%d; str_strategy %s',flag_periodic_x,flag_periodic_y,str_strategy));
tmp_parameter = struct('type','parameter','flag_periodic_x',flag_periodic_x,'flag_periodic_y',flag_periodic_y,'str_strategy',str_strategy);
tmp_parameter_bf = tmp_parameter; tmp_parameter_bf.str_strategy='brute_force';
n_w_x = 0; n_w_y = 0; 
%%;
for ntest=0:5-1;
if ntest==0; n_w_x = 1; n_w_y = 0; end;
if ntest==1; n_w_x = 0; n_w_y = 1; end;
if ntest==2; n_w_x = 1; n_w_y = 1; end;
if ntest==3; n_w_x = 2; n_w_y = 3; end;
if ntest==4; n_w_x = 3; n_w_y = 2; end;
[~,E_xy__] = max_filter_2d_0(tmp_parameter,n_x,n_y,A_xy__,n_w_x,n_w_y);
[~,F_xy__] = max_filter_2d_0(tmp_parameter_bf,n_x,n_y,A_xy__,n_w_x,n_w_y);
if (fnorm(F_xy__-E_xy__)<=1e-12);
disp(sprintf(' %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %% '));
disp(sprintf(' %% flag_periodic %d,%d; str_strategy %s: <-- passed',flag_periodic_x,flag_periodic_y,str_strategy));
end;%if (fnorm(F_xy__-E_xy__)<=1e-12);
if (fnorm(F_xy__-E_xy__)> 1e-12);
disp(sprintf(' %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %% '));
disp(sprintf(' %% flag_periodic %d,%d; str_strategy %s',flag_periodic_x,flag_periodic_y,str_strategy));
disp(sprintf(' %% %% n_w %d,%d: A_xy__',n_w_x,n_w_y));
disp(num2str(A_xy__,' %+.3d'));
disp(sprintf(' %% %% n_w %d,%d: F_xy__',n_w_x,n_w_y));
disp(num2str(F_xy__,' %+.3d'));
disp(sprintf(' %% %% n_w %d,%d: E_xy__',n_w_x,n_w_y));
disp(num2str(E_xy__,' %+.3d'));
disp(sprintf(' %% %% n_w %d,%d: F_xy__ - E_xy__',n_w_x,n_w_y));
disp(num2str(F_xy__ - E_xy__,' %+.3d'));
disp('test failed: returning'); return;
end;%if (fnorm(F_xy__-E_xy__)> 1e-12);
end;%for ntest=0:2-1;
disp(sprintf(' %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %% '));
%%;
end;%for str_strategy_={'loop+max1'};
end;%for flag_periodic_y=[0,1];
end;%for flag_periodic_x=[0,1];
%%%%;
for str_strategy_={'loop+max1'};
str_strategy = str_strategy_{1};
disp(sprintf(' %% str_strategy: %s',str_strategy));
tmp_parameter = struct('type','parameter','str_strategy',str_strategy);
n_x = 1024*5; n_y = 1024*5; A_xy__ = round(128*randn(n_x,n_y));
n_w_x = 4; n_w_y = n_w_x;
n_GB_1 = n_x*n_y*max(n_w_x*n_w_y)*8/1e9;
disp(sprintf(' %% n_x %d n_y %d n_w_x %d n_w_y %d <-- single step memory %.2f GB',n_x,n_y,n_w_x,n_w_y,n_GB_1));
if n_GB_1< 4;
tmp_t = tic(); [~,E_xy__] = max_filter_2d_0(tmp_parameter,n_x,n_y,A_xy__,n_w_x,n_w_y); tmp_t = toc(tmp_t);
disp(sprintf(' %% %% single step: %0.4fs',tmp_t));
end;%if n_GB_1< 4;
n_GB_0 = n_x*n_y*max(1*1)*8/1e9;
disp(sprintf(' %% n_x %d n_y %d n_w_x %d n_w_y %d <-- iterative memory %.2f GB',n_x,n_y,n_w_x,n_w_y,n_GB_0));
if n_GB_0< 4;
tmp_t = tic(); 
F_xy__ = A_xy__;
for nw_x=0:n_w_x-1;
[~,F_xy__] = max_filter_2d_0(tmp_parameter,n_x,n_y,F_xy__);
end;%for nw_x=0:n_w_x-1;
tmp_t = toc(tmp_t);
disp(sprintf(' %% %% iterative: %0.4fs',tmp_t));
end;%if n_GB_0< 4;
if (n_GB_1< 4 & n_GB_0< 4); disp(sprintf(' %% %% E_xy__ vs F_xy__: %0.16f',fnorm(E_xy__-F_xy__)/fnorm(E_xy__))); end;
end;%for str_strategy_={'loop+max1'};
%%%%;
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); n_y=[]; end; na=na+1;
if (nargin<1+na); A_xy__=[]; end; na=na+1;
if (nargin<1+na); n_w_x=[]; end; na=na+1;
if (nargin<1+na); n_w_y=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_periodic_x'); parameter.flag_periodic_x=0; end;
flag_periodic_x=parameter.flag_periodic_x;
if ~isfield(parameter,'flag_periodic_y'); parameter.flag_periodic_y=0; end;
flag_periodic_y=parameter.flag_periodic_y;
if ~isfield(parameter,'str_strategy'); parameter.str_strategy='loop+max1'; end;
str_strategy=parameter.str_strategy;

if isempty(n_w_x); n_w_x = 1; end;
if isempty(n_w_y); n_w_y = 1; end;

%%%%%%%%%%%%%%%%;
if  strcmp(str_strategy,'brute_force');
A_min = min(A_xy__,[],'all');
E_xy__ = zeros(size(A_xy__));
for ny=0:n_y-1;
for nx=0:n_x-1;
tmp_A = A_min;
for nw_y=-n_w_y:+n_w_y;
ny2 = ny+nw_y; 
if flag_periodic_y==1; ny2 = periodize(ny2,0,n_y); end;
if flag_periodic_y==0; ny2 = max(0,min(n_y-1,ny2)); end;
for nw_x=-n_w_x:+n_w_x;
nx2 = nx+nw_x; 
if flag_periodic_x==1; nx2 = periodize(nx2,0,n_x); end;
if flag_periodic_x==0; nx2 = max(0,min(n_x-1,nx2)); end;
tmp_A = max(tmp_A,A_xy__(1+nx2,1+ny2));
end;%for nw_x=-n_w_x:+n_w_x;
end;%for nw_y=-n_w_y:+n_w_y;
E_xy__(1+nx,1+ny) = tmp_A;
end;%for nx=0:n_x-1;
end;%for ny=0:n_y-1;
end;%if  strcmp(str_strategy,'brute_force');
%%%%%%%%%%%%%%%%;
if ~strcmp(str_strategy,'brute_force');
if (flag_periodic_x==0) | (flag_periodic_y==0);
n_p_x = (flag_periodic_x==0)*n_w_x;
n_p_y = (flag_periodic_y==0)*n_w_y;
A_min = min(A_xy__,[],'all');
B_xy__ = A_min*ones(2*n_p_x+n_x,2*n_p_y+n_y);
B_xy__(n_p_x+1:n_p_x+n_x,n_p_y+1:n_p_y+n_y) = A_xy__;
parameter_periodic = parameter; parameter_periodic.flag_periodic_x = 1; parameter_periodic.flag_periodic_y = 1;
[parameter_periodic,D_xy__] = max_filter_2d_0(parameter_periodic,2*n_p_x+n_x,2*n_p_y+n_y,B_xy__,n_w_x,n_w_y);
E_xy__ = reshape(D_xy__(n_p_x+1:n_p_x+n_x,n_p_y+1:n_p_y+n_y),size(A_xy__));
parameter = parameter_periodic; parameter.flag_periodic_x = flag_periodic_x; parameter.flag_periodic_y = flag_periodic_y;
end;%if (flag_periodic_x==0) | (flag_periodic_y==0);
%%%%%%%%;
if (flag_periodic_x==1) & (flag_periodic_y==1);
if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
A_min = min(A_xy__,[],'all');
n_w_x = min(n_x,n_w_x);
n_w_y = min(n_y,n_w_y);
n_B_x = 2*n_w_x+n_x; n_B = n_B_x*n_y;
ij_x_ = n_w_x+1:n_w_x+n_x;
B_xy__ = zeros(n_B_x,n_y);
B_xy__(ij_x_,:) = A_xy__(1:n_x,:); %<-- center <-- center. ;
B_xy__(1:n_w_x,:) = A_xy__(n_x-n_w_x+1:n_x,:); %<-- x periodic. ;
B_xy__(n_w_x+n_x+1:n_w_x+n_x+n_w_x,:) = A_xy__(1:n_w_x,:); %<-- x periodic. ;
[~,C_xy__] = max_filter_1d_0(parameter,n_B,B_xy__,n_w_x); C_xy__ = C_xy__(n_w_x+1:n_w_x+n_x,:);
C_yx__ = transpose(C_xy__);
n_B_y = 2*n_w_y+n_y; n_B = n_B_y*n_x;
ij_y_ = n_w_y+1:n_w_y+n_y;
B_yx__ = zeros(n_B_y,n_x);
B_yx__(ij_y_,:) = C_yx__(1:n_y,:); %<-- center <-- center. ;
B_yx__(1:n_w_y,:) = C_yx__(n_y-n_w_y+1:n_y,:); %<-- y periodic. ;
B_yx__(n_w_y+n_y+1:n_w_y+n_y+n_w_y,:) = C_yx__(1:n_w_y,:); %<-- y periodic. ;
[~,C_yx__] = max_filter_1d_0(parameter,n_B,B_yx__,n_w_y); C_yx__ = C_yx__(n_w_y+1:n_w_y+n_y,:);
C_xy__ = transpose(C_yx__);
E_xy__ = C_xy__;
end;%if (flag_periodic_x==1) & (flag_periodic_y==1);
%%%%%%%%;
end;%if ~strcmp(str_strategy,'brute_force');
%%%%%%%%%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
