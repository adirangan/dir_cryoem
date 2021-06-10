function c = legendre_normalization_factor_0(l_val,m_val);
a1=((2*l_val+1)/(4*pi));
m_abs = abs(m_val);
a2=exp(lfactorial(l_val-m_abs) - lfactorial(l_val+m_abs));
c=sqrt(a1.*a2);
