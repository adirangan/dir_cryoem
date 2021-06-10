function output = recenter1(input);
[nrows] = length(input);
new_rows = [floor(nrows/2)+1:nrows , 1:floor(nrows/2)];
output = input(new_rows);
