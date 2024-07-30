function [o0,o1,o2,o3,o4,o5] = tmp17(i0,i1,i2,i3,i4,i5);

if (nargin<1);
%%%%;
disp(sprintf(' %% external: calling ''tmp17(1);'''));
tmp17(1);
%%%%;
disp(sprintf(' %% external: calling ''[a1,a2,a3] = tmp17(1);'''));
[a1,a2,a3] = tmp17(1);
%%%%;
disp(sprintf(' %% external: calling ''[a1,a2,a3,~,~,~] = tmp17(1);'''));
[a1,a2,a3,~,~,~] = tmp17(1);
%%%%;
disp('returning');return;
end;%if (nargin<1);

disp(sprintf(' %% internal: tmp17 called with nargin %d nargout %d',nargin,nargout));
o0=[];
o1=[];
o2=[];
o3=[];
o4=[];
o5=[];
