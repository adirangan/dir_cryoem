function output = recenter2(input);
% test with: ;
%{

  n_x = 7; n_y = 9;
  input = zeros(n_x,n_y);
  for ny=1:n_y; for nx=1:n_x;
  input(nx,ny) = (nx-1) + i*(ny-1);
  end;end;%for ny=1:n_y; for nx=1:n_x;
  disp('original');
  disp(num2str(input));
  output = recenter2(input);
  disp('recentered');
  disp(num2str(output));
  output2 = recenter2(output);
  disp('recentered again');
  disp(num2str(output2));
    
  %}
[nrows,ncols] = size(input);
new_rows = [floor(nrows/2)+1:nrows , 1:floor(nrows/2)];
new_cols = [floor(ncols/2)+1:ncols , 1:floor(ncols/2)];
output = input(new_rows,new_cols);
