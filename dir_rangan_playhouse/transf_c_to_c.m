function M_c_ = transf_c_to_c(n_x,max_x_c,n_y,max_y_c,S_c_,delta_x,delta_y);

if (nargin<7);
test_transforms_dr.m();
disp('returning'); return;
end;%if (nargin<7);

T_c_ = recenter2(S_c_);
N_c_ = zeros(n_x,n_y);
for ny=0:n_y-1;
for nx=0:n_x-1;
X_c = 0.0 + nx*max_x_c/n_x - max_x_c/2.0;
Y_c = 0.0 + ny*max_y_c/n_y - max_y_c/2.0;
L_c = (X_c * delta_x) + (Y_c * delta_y);
C_c = exp(-i*2*pi*L_c);
N_c_(1+nx,1+ny) = C_c*T_c_(1+nx,1+ny);
end;%for nx=0:n_x-1;
end;%for ny=0:n_y-1;
M_c_ = decenter2(N_c_);
