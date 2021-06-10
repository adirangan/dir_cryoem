function S_x_c_ = get_F2_x_c_(n_x,grid_x_c_,max_x_c,n_y,grid_y_c_,max_y_c,param_1,param_2);
for ny=0:n_y-1;
if (ny<floor(n_y/2)); Y_x_c = grid_y_c_(1+ny); else Y_x_c = grid_y_c_(1+ny) - max_y_c; end;
for nx=0:n_x-1;
if (nx<floor(n_x/2)); X_x_c = grid_x_c_(1+nx); else X_x_c = grid_x_c_(1+nx) - max_x_c; end;
F_x_c = get_F2_x_c(max_x_c,X_x_c,Y_x_c,param_1,param_2);
S_x_c_(1+nx+ny*n_x) = F_x_c;
end;%for nx=0:n_x-1;
end;%for ny=0:n_y-1;
