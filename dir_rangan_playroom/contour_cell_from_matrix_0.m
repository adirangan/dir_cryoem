function [n_contour,contour_cell___,contour_level_] = contour_cell_from_matrix_0(contour_matrix__);

n_col = size(contour_matrix__,2);
ncontour=0; na=0;
while (na<n_col-1);
n_point = contour_matrix__(1+1,1+na);
ncontour = ncontour+1;
na = na+1+n_point;
end;%while (na<n_col-1);

n_contour = ncontour;
contour_level_ = zeros(n_contour,1);
contour_cell___ = cell(n_contour,1);

ncontour=0; na=0;
while (na<n_col-1);
n_point = contour_matrix__(1+1,1+na);
contour_level = contour_matrix__(1+0,1+na);
contour_x_ = contour_matrix__(1+0,1+na + [1:n_point]);
contour_y_ = contour_matrix__(1+1,1+na + [1:n_point]);
contour_level_(1+ncontour) = contour_level;
contour_cell___{1+ncontour} = [transpose(contour_x_) , transpose(contour_y_)];
ncontour = ncontour+1;
na = na+1+n_point;
end;%while (na<n_col-1);

