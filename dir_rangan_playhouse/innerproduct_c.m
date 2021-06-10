function C_c = innerproduct_c(n_x,grid_x_c_,n_y,grid_y_c_,T_c_,M_c_);
% assumes that M_c_ is the same size and dimensions as T_c_ ;
% assumes uniform grid and periodic boundary conditions ;
dx = mean(diff(grid_x_c_));
dy = mean(diff(grid_y_c_));
C_c = sum(sum(conj(T_c_).*(M_c_)))*dx*dy;
