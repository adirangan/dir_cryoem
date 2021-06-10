function output = recenter(input);
[nrows,ncols] = size(input);
new_rows = [floor(nrows/2)+1:nrows , 1:floor(nrows/2)];
new_cols = [floor(ncols/2)+1:ncols , 1:floor(ncols/2)];
output = input(new_rows,new_cols);
