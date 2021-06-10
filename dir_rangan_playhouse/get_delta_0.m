function [delta_x_,delta_y_] = get_delta_0(N_pixels_in,n_r,half_diameter_x_c,n_delta_x,n_delta_y);
for ndx=0:n_delta_x-1;
if (n_delta_x>1);
delta_x = (-N_pixels_in + ndx*2*N_pixels_in/(n_delta_x-1))/n_r*half_diameter_x_c;
 else;
delta_x = 0.0;
end;%if;
delta_x_(1+ndx) = delta_x;
for ndy=0:n_delta_y-1;
if (n_delta_y>1);
delta_y = (-N_pixels_in + ndy*2*N_pixels_in/(n_delta_y-1))/n_r*half_diameter_x_c;
 else;
delta_y = 0.0;
end;%if;
delta_y_(1+ndy) = delta_y;
end;%for;
end;%for;

