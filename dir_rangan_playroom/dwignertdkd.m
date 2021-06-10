function dWtdkd = dwignertdkd(m_val,l_val_a,l_val_b);
% calculates the derivative of wignert(kd;m_val,l_val_a,l_val_b) with respect to kd. ;
dWtdkd=0;
if (abs(l_val_a-l_val_b)==1);
l_pos = max(l_val_a,l_val_b); l_neg = min(l_val_a,l_val_b); m_abs = abs(m_val);
dWtdkd = 2*pi*i*sqrt(exp(lfactorial(l_pos-m_abs)+lfactorial(l_pos+m_abs)-lfactorial(l_neg-m_abs)-lfactorial(l_neg+m_abs))/(2*l_pos+1)/(2*l_neg+1));
end;%if (abs(l_val_a-l_val_b)==1);
