function [polar_a_periodize_,azimu_b_periodize_] = periodize_polar_a_azimu_b_0(polar_a_,azimu_b_);

%%%%%%%%;
% periodize azimu_b_ (standard). ;
%%%%%%%%;
azimu_b_periodize_ = periodize(azimu_b_,0,2*pi);
%%%%%%%%;
% periodize polar_a_ (note behavior at poles). ;
%%%%%%%%;
polar_a_periodize_ = periodize(polar_a_,-1*pi/2,+3*pi/2);
tmp_index_ = efind(polar_a_periodize_< 0*pi);
polar_a_periodize_(1+tmp_index_) = 0*pi - (polar_a_periodize_(1+tmp_index_) - 0*pi);
azimu_b_periodize_(1+tmp_index_) = +azimu_b_periodize_(1+tmp_index_) + 1*pi;
tmp_index_ = efind(polar_a_periodize_> 1*pi);
polar_a_periodize_(1+tmp_index_) = 1*pi - (polar_a_periodize_(1+tmp_index_) - 1*pi);
azimu_b_periodize_(1+tmp_index_) = +azimu_b_periodize_(1+tmp_index_) - 1*pi;
%%%%%%%%;
% periodize azimu_b_ (standard). ;
%%%%%%%%;
azimu_b_periodize_ = periodize(azimu_b_periodize_,0,2*pi);
