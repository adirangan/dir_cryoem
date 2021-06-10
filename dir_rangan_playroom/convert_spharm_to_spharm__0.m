function a_k_Y__ = convert_spharm_to_spharm__0(l_max,a_k_Y_);
% converts a linearly-ordered spherical-harmonic-expansion into a stacked array. ;
% More specifically, if a_k_Y_ stores the components of a such that m_val varies ;
% quickly and l_val varies slowly (i.e., a_k_Y(l_val;m_val) = a_k_Y_(1+l_val^2+l_val+m_val) ) ;
% then a_k_Y__ stores the components so that l_val corresponds to the rows, ;
% and m_val corresponds to the columns. ;
% Specifically: a_k_Y(l_val;m_val) = a_k_Y__(1+l_val,1+l_max+m_val); %<-- note shift by l_max. ;

n_l_val = 1 + 1*l_max;
n_lm_max = (1+l_max)^2;
n_m_val = 1 + 2*l_max;

a_k_Y__ = zeros(n_l_val,n_m_val);
nlm_max=0;
for l_val=0:l_max;
for m_val=-l_val:+l_val;
a_k_Y__(1+l_val,1+l_max+m_val) = a_k_Y_(1+nlm_max); %<-- note shift by l_max. ;
nlm_max=nlm_max+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
assert(nlm_max==n_lm_max);
