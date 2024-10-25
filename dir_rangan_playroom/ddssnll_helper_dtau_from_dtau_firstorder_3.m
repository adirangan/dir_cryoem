%%%%%%%%;
% Now if dtau is empty we fill dtau with dtau_firstorder. ;
%%%%%%%%;
if isempty(dtau_euler_polar_a_M_); dtau_euler_polar_a_M_ = real(dtau_firstorder_M3__(:,1+0)); end;
if isempty(dtau_euler_azimu_b_M_); dtau_euler_azimu_b_M_ = real(dtau_firstorder_M3__(:,1+1)); end;
if isempty(dtau_euler_gamma_z_M_); dtau_euler_gamma_z_M_ = real(dtau_firstorder_M3__(:,1+2)); end;
