function M_c_ = transl_c_to_c(n_x,max_x_c,n_y,max_y_c,S_c_,delta_x,delta_y);
% Assumes that M_c_ is the same size and dimensions as S_c_;
T_c_ = recenter2(S_c_);
N_c_ = zeros(n_x,n_y);
for ny=0:n_y-1;
for nx=0:n_x-1;
X_c = 0.0 + nx*max_x_c/n_x - delta_x;
Y_c = 0.0 + ny*max_y_c/n_y - delta_y;
C_c = interp2(n_x,0.0,max_x_c,n_y,0.0,max_y_c,T_c_,X_c,Y_c);
N_c_(1+nx+ny*n_x) = C_c;
end;%for;
end;%for;
M_c_ = decenter2(N_c_);

