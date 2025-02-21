function ax = subplot_t(p_row,p_col,ij);
index = ij-1;
nr = mod(index,p_row);
index = index - nr;
index = index/p_row;
nc = mod(index,p_col);

ax = subplot(p_row,p_col,1 + nc + nr*p_col);



