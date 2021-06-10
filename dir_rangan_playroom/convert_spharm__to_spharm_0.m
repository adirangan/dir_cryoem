function a_k_Y_ = convert_spharm__to_spharm_0(l_max,a_k_Y__);
% converts a stacked array of spherical-harmonic-expansion into a linearly-ordered array. ;
% More specifically, if a_k_Y__ stores the components of a so that l_val corresponds to the rows, ;
% and m_val corresponds to the columns, ;
% i.e., a_k_Y(l_val;m_val) = a_k_Y__(1+l_val,1+l_max+m_val); %<-- note shift by l_max. ;
% then a_k_Y_ stores the components of a such that m_val varies quickly and l_val varies slowly. ;
% (i.e., a_k_Y(l_val;m_val) = a_k_Y_(1+l_val^2+l_val+m_val) ). ;

n_l_val = 1 + 1*l_max;
n_lm_max = (1+l_max)^2;
n_m_val = 1 + 2*l_max;
a_k_Y_ = zeros(n_lm_max,1);
for l_val=0:l_max;
a_k_Y_(1+l_val^2+l_val+[-l_val:+l_val]) = a_k_Y__(1+l_val,1+l_max+[-l_val:+l_val]);
end;%for l_val=0:l_max;
