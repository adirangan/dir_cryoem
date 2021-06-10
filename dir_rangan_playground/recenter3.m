function output = recenter3(input);
[nrows,ncols,nlyrs] = size(input);
new_rows = [floor(nrows/2)+1:nrows , 1:floor(nrows/2)];
new_cols = [floor(ncols/2)+1:ncols , 1:floor(ncols/2)];
new_lyrs = [floor(nlyrs/2)+1:nlyrs , 1:floor(nlyrs/2)];
output = input(new_rows,new_cols,new_lyrs);
