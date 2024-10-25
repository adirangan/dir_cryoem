if  isempty(dtau_ssnll_q2d_3M__) & ~isempty(dtau_ssnll_q2d_M3__);
dtau_ssnll_q2d_3M__ = permute(dtau_ssnll_q2d_M3__,[2,1]);
end;%if  isempty(dtau_ssnll_q2d_3M__) & ~isempty(dtau_ssnll_q2d_M3__);
if  isempty(dtau_dtau_ssnll_q2d_33M___) & ~isempty(dtau_dtau_ssnll_q2d_M33___);
dtau_dtau_ssnll_q2d_33M___ = permute(dtau_dtau_ssnll_q2d_M33___,[2,3,1]);
end;%if  isempty(dtau_dtau_ssnll_q2d_33M___) & ~isempty(dtau_dtau_ssnll_q2d_M33___);
if  isempty(dtau_dvol_ssnll_q2d_3M__) & ~isempty(dtau_dvol_ssnll_q2d_M3__);
dtau_dvol_ssnll_q2d_3M__ = permute(dtau_dvol_ssnll_q2d_M3__,[2,1]);
end;%if  isempty(dtau_dvol_ssnll_q2d_3M__) & ~isempty(dtau_dvol_ssnll_q2d_M3__);

