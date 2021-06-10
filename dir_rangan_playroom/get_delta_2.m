function [n_delta_v_out,delta_x_,delta_y_] = get_delta_2(delta_r_max,n_delta_v_0in);
% This function outputs an array of displacements delta_x_, delta_y_ ;
% that are supported on a disc of radius delta_r_max. ;
% The number of requested displacents n_delta_v_0in is a lower-bound on the actual number produced. ;

if n_delta_v_0in<=1;
n_delta_v_out = 1; delta_x_ = 0; delta_y_ = 0;
end;%if n_delta_v_0in<=1;

if n_delta_v_0in> 1;
n_x = 1 + floor(sqrt(n_delta_v_0in)); continue_flag=1;
while continue_flag;
x_ = linspace(-delta_r_max,+delta_r_max,n_x);
y_ = linspace(-delta_r_max,+delta_r_max,n_x);
[X_,Y_] = ndgrid(x_,y_);
R_ = sqrt(X_.^2 + Y_.^2);
tmp_ij_ = find(R_(:)<=delta_r_max);
if numel(tmp_ij_)>=n_delta_v_0in; continue_flag=0; end;
if numel(tmp_ij_)< n_delta_v_0in; n_x = n_x + 1; continue_flag=1; end;
end;%while;
n_delta_v_out = numel(tmp_ij_);
delta_x_ = X_(tmp_ij_);
delta_y_ = Y_(tmp_ij_);
if (mod(n_x,2)==0); n_delta_v_out = n_delta_v_out + 1; end;
if (mod(n_x,2)==0); delta_x_ = [0;delta_x_]; end;
if (mod(n_x,2)==0); delta_y_ = [0;delta_y_]; end;
end;%if n_delta_v_0in> 1;
