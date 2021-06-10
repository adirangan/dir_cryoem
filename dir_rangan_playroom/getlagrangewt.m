function val = getlagrangewt(delta_polar_a,delta_azimu_b,n_order,index_polar_a,index_azimu_b,spacing_polar_a,spacing_azimu_b);
% compute interpolation_weight at (delta_polar_a,delta_azimu_b) from node (index_polar_a,index_azimu_b). ;
% Input: ;
% delta_polar_a = real offest from base point (assumed to be the origin). ;
% delta_azimu_b = real offest from base point (assumed to be the origin). ;
% n_order = integer interpolation order. ;
% index_polar_a = integer grid point index. lagrange polynomial will be evaluated at (delta_polar_a,delta_azimu_b). ;
% index_azimu_b = integer grid point index. lagrange polynomial will be evaluated at (delta_polar_a,delta_azimu_b). ;
% spacing_polar_a = real grid spacing in polar_a. ;
% spacing_azimu_b = real grid spacing in azimu_b. ;
% ;
% Output: ;
% val = value of lagrange interpolant associated with (index_polar_a,index_azimu_b). ;

delta_polar_a_location = delta_polar_a/spacing_polar_a;
delta_azimu_b_location = delta_azimu_b/spacing_azimu_b;

product_polar_a = 1.0;
for norder=-floor(n_order/2):floor(n_order/2);
if (index_polar_a~=norder);
product_polar_a = product_polar_a*(delta_polar_a - norder*spacing_polar_a)/(index_polar_a*spacing_polar_a - norder*spacing_polar_a); 
end;%if (index_polar_a~=norder);
end;%for norder=-floor(n_order/2):floor(n_order/2);

product_azimu_b = 1.0;
for norder=-floor(n_order/2):floor(n_order/2);
if (index_azimu_b~=norder);
product_azimu_b = product_azimu_b*(delta_azimu_b - norder*spacing_azimu_b)/(index_azimu_b*spacing_azimu_b - norder*spacing_azimu_b); 
end;%if (index_azimu_b~=norder);
end;%for norder=-floor(n_order/2):floor(n_order/2);

val = product_polar_a*product_azimu_b;


