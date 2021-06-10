function S_p_ = interp_c_to_p(n_x,max_x_c,n_y,max_y_c,S_c_,n_r,grid_p_,n_w_,n_A);

if (nargin<9);
test_transforms_dr();
disp('returning');return;
end;%if (nargin<9);

S_p_ = zeros(n_A,1);
T_c_ = recenter2(S_c_);
max_r_c = sqrt(max_x_c.^2 + max_y_c.^2);
ic=0;
for nr=0:n_r-1;
R_c = grid_p_(1+nr);
for nw=0:n_w_(1+nr)-1;
W_c = 0.0 + nw*(2*pi)/(n_w_(1+nr));
X_c = R_c*cos(W_c);
Y_c = R_c*sin(W_c);
C_c = interp2(n_x,-max_x_c/2.0,max_x_c/2.0,n_y,-max_y_c/2.0,max_y_c/2.0,T_c_,X_c,Y_c);
S_p_(1+ic) = C_c;
ic = ic + 1;
end;%for nw=0:n_w_(1+nr)-1;
end;%for nr=0:n_r-1;

