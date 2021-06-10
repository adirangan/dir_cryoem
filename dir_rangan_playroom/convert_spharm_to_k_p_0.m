function ...
[ ...
 n_all ...
,n_sub_ ...
,k_p_r_all_ ...
,azimu_b_all_ ...
,polar_a_all_ ...
,weight_all_ ...
,a_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
] = ...
convert_spharm_to_k_p_0( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,l_val_ ...
,a_ ...
,eq_d ...
);
% evaluates spherical-harmonic-expansion a_ on a collection of points on spherical shells determined by k_p_r_. ;
% These points are themselves represented using polar (spherical) coordinates (i.e., azimu_b_ and polar_a_). ;
% The quadrature weights associated with these points are also calculated. ; 
% ;
% inputs: ;
% ;
% verbose = integer verbosity_level ;
% n_k_p_r = integer maximum k ;
% k_p_r_ = real array of length n_k_p_r; k_p_r_(nk_p_r) = k_p_r_value for shell nk_p_r ;
% l_val_ = integer array of length n_k_p_r; l_val_(nk_p_r) = spherical harmonic order on shell nk_p_r; l_val_(nk_p_r) corresponds to n_lm_(nk_p_r) = (l_val_(nk_p_r)+1)^2 coefficients ;
% a_ = complex array of length \sum_{nk_p_r} (n_lm_(nk_p_r)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% eq_d = distance between points on the equator (passed to sample_shell_0) ;
% ;
% outputs: ;
% ;
% n_all = integer total number of points ;
% n_sub_ = integer array of starting indices associated with each k-value ;
% k_p_r_all_ = real array of k-values for each point ;
% azimu_b_all_ = real array of azimu_b-values for each point ;
% polar_a_all_ = real array of polar_a-values for each point ;
% weight_all_ = real array of quadrature weights for each point ;
% a_all_ = complex array of a-values for each point ;
% k_c_0_all_ = real array of k_c_0 values (cartesian) for each point ;
% k_c_1_all_ = real array of k_c_1 values (cartesian) for each point ;
% k_c_2_all_ = real array of k_c_2 values (cartesian) for each point ;

n_lm_ = (l_val_+1).^2;
k_p_r_max = k_p_r_(end);
l_val_max = l_val_(end);
m_max_ = -l_val_max : +l_val_max;
n_m_max = length(m_max_);

n_all_ = zeros(n_k_p_r,1);
for nk_p_r=1:n_k_p_r;
k_p_r = k_p_r_(nk_p_r); 
n_all_(nk_p_r) =  sample_shell_0(k_p_r,eq_d,'L');
end%for nk_p_r=1:n_k_p_r;
n_all = sum(n_all_);
k_p_r_all_ = zeros(n_all,1); azimu_b_all_ = zeros(n_all,1); polar_a_all_ = zeros(n_all,1); weight_all_ = zeros(n_all,1);
a_all_ = zeros(n_all,1);

n_sub_ = zeros(n_k_p_r,1);
n_sub = 0;
for nk_p_r=1:n_k_p_r;
k_p_r = k_p_r_(nk_p_r); 
n_sub_(nk_p_r) = n_sub;
[length_sub,azimu_b_sub_,polar_a_sub_,weight_sub_] = sample_shell_0(k_p_r,eq_d,'L');
l_val_max = l_val_(nk_p_r); n_lm = n_lm_(nk_p_r); ix_base = sum(n_lm_(1:nk_p_r-1));
a_k_ = a_(ix_base + (1:n_lm)); 
l_ = []; m_ = []; for nl=0:l_val_max; l_ = [l_ , nl*ones(1,2*nl+1) ]; m_ = [m_ , [-nl:+nl] ]; end;%for nl=0:l_val_max;
if (verbose>1); disp(sprintf(' %% nk_p_r %d k_p_r %d/%d: l_val_max %d n_lm %d ix_base %d',nk_p_r,k_p_r,k_p_r_max,l_val_max,n_lm,ix_base)); end;
for nl=0:l_val_max;
l_val = nl;
if (verbose>2); disp(sprintf(' %% nk_p_r %d k %d/%d: nl %d l_val %d',nk_p_r,k_p_r,k_p_r_max,nl,l_val)); end;
Llm__=legendre(l_val,cos(polar_a_sub_),'unnorm');
A_ = zeros(1,length_sub);
for m_val = -l_val:+l_val;
ix = 1+l_val*(l_val+1)+m_val;
m_abs = abs(m_val);
if (length_sub>0); if (l_val>0); Llm_ = squeeze(Llm__(1+m_abs,:,:)); end; if (l_val==0); Llm_ = Llm__(:,:); end; end;
a1=((2*l_val+1)/(4*pi));
a2=exp(lfactorial(l_val-m_abs) - lfactorial(l_val+m_abs));
c=sqrt(a1*a2);
s=(-1)^((m_val<0)*m_val); % needed to preserve condon-shortley phase. ;
%s=1; % original phase ;
if (length_sub>0); 
Ylm_ = s*c*Llm_.*exp(+i*m_val*transpose(azimu_b_sub_)); 
A_ = A_ + s*a_k_(ix)*Ylm_;
end;%if (length_sub>0); 
end;%for m_val = -l_val:+l_val;
if (length_sub>0); 
a_all_(1 + n_sub + (0:length_sub-1)) = a_all_(1 + n_sub + (0:length_sub-1)) + transpose(A_);
end;%if (length_sub>0); 
end;%for nl=0:l_val_max;
k_p_r_all_(1 + n_sub + (0:length_sub-1)) = k_p_r*ones(length_sub,1);
azimu_b_all_(1 + n_sub + (0:length_sub-1)) = azimu_b_sub_;
polar_a_all_(1 + n_sub + (0:length_sub-1)) = polar_a_sub_;
weight_all_(1 + n_sub + (0:length_sub-1)) = weight_sub_;
n_sub = n_sub + length_sub;
end;%for nk_p_r=1:n_k_p_r;
k_c_0_all_ = k_p_r_all_ .* cos(azimu_b_all_) .* sin(polar_a_all_);
k_c_1_all_ = k_p_r_all_ .* sin(azimu_b_all_) .* sin(polar_a_all_);
k_c_2_all_ = k_p_r_all_ .* cos(polar_a_all_);



