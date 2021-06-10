function [b_all_] = transf_k_p_to_k_p_0(verbose,n_all,k_p_r_all_,azimu_b_all_,polar_a_all_,a_all_,delta_);
% Applies translation in fourier space to list of points. ;
% ;
% inputs: ;
% ;
% verbose = integer verbosity_level ;
% n_all = integer total number of points ;
% k_p_r_all_ = real array of k-values for each point ;
% azimu_b_all_ = real array of azimu_b-values for each point ;
% polar_a_all_ = real array of polar_a-values for each point ;
% a_all_ = complex array of a-values for each point ;
% delta_ = real array of real-space displacements; [delta_x,delta_y,delta_z] ;
% ;
% outputs: ;
% ;
% b_all_ = complex array of b-values for each point ;

if (verbose>1); disp(sprintf(' %% translating by [%0.2f %0.2f %0.2f]',delta_)); end;
b_all_ = a_all_ ;
if (norm(delta_)>0);
kx_all_ = k_p_r_all_ .* cos(azimu_b_all_) .* sin(polar_a_all_);
ky_all_ = k_p_r_all_ .* sin(azimu_b_all_) .* sin(polar_a_all_);
kz_all_ = k_p_r_all_ .* cos(polar_a_all_);
for nall=1:n_all;
kd = kx_all_(nall)*delta_(1) + ky_all_(nall)*delta_(2) + kz_all_(nall)*delta_(3);
b_all_(nall) = a_all_(nall) * exp(i*kd);
end;%for nall=1:n_all;
end;%if (norm(delta_)>0);
