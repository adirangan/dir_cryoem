function output = decenter1(input);
[nrows] = length(input);
new_rows = [ceil(nrows/2)+1:nrows , 1:ceil(nrows/2)];
output = input(new_rows);
