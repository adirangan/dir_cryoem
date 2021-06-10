function output = decenter2(input);

if (nargin<1);
n_x = 8; n_y = 10;
input = zeros(n_x,n_y);
for ny=1:n_y; for nx=1:n_x;
input(nx,ny) = (nx-1) + i*(ny-1);
end;end;%for ny=1:n_y; for nx=1:n_x;
disp('original');
disp(num2str(input));
output = recenter2(input);
disp('recentered');
disp(num2str(output));
output2 = decenter2(output);
disp('decentered');
disp(num2str(output2));
disp('returning'); return;
end;%if (nargin<1);

[nrows,ncols] = size(input);
new_rows = [ceil(nrows/2)+1:nrows , 1:ceil(nrows/2)];
new_cols = [ceil(ncols/2)+1:ncols , 1:ceil(ncols/2)];
output = input(new_rows,new_cols);
