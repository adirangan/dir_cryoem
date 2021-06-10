function output = decenter3(input);

if (nargin<1);
n_x = 8; n_y = 10; n_z = 12;
n_x = 9; n_y = 11; n_z = 13;
input = zeros(n_x,n_y,n_z);
for nz=1:n_z; for ny=1:n_y; for nx=1:n_x;
input(nx,ny,nz) = (nx-1) + i*(ny-1) + sqrt(2)*(nz-1);
end;end;end;%for nz=1:n_z; for ny=1:n_y; for nx=1:n_x;
disp('original');
disp(num2str(input));
output = recenter3(input);
disp('recentered');
disp(num2str(output));
output2 = decenter3(output);
disp('decentered');
disp(num2str(output2));
disp('returning'); return;
end;%if (nargin<1);

[nrows,ncols,nlyrs] = size(input);
new_rows = [ceil(nrows/2)+1:nrows , 1:ceil(nrows/2)];
new_cols = [ceil(ncols/2)+1:ncols , 1:ceil(ncols/2)];
new_lyrs = [ceil(nlyrs/2)+1:nlyrs , 1:ceil(nlyrs/2)];
output = input(new_rows,new_cols,new_lyrs);

