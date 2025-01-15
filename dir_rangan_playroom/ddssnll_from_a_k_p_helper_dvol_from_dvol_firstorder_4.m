%%%%%%%%;
% Now if dvol_a is empty we fill dtau with dtau_firstorder. ;
%%%%%%%%;
if isempty(dvol_a_k_p_qk_); dvol_a_k_p_qk_ = dvol_a_firstorder_k_p_qk_; dvol_a_R_k_p_qk_ = []; end;
if isempty(dvol_a_k_p_qk__); dvol_a_k_p_qk__ = dvol_a_firstorder_k_p_qk__; dvol_a_R_k_p_qk__ = []; end;

