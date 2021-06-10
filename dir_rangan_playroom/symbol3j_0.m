function output = symbol3j_0(l_A,l_B,l_C,m_A,m_B,m_C);
% integral associated with: ;
% \int_{sphere} conjg(Y_{l_C}^{m_C}) .* Y_{l_A}^{m_A} .* Y_{l_B}^{m_B}. ;
tmp_wx = wigner3j_1(l_C,l_A,l_B,-m_C,+m_A,+m_B);
tmp_w0 = 0.0d0; if (mod(l_C+l_A+l_B,2)==0); tmp_w0 = wigner3j_1(l_C,l_A,l_B,0,0,0); end;
tmp_zz = sqrt((2*l_A+1)*(2*l_B+1)*(2*l_C+1))/sqrt(4*pi);
%tmp_ss = (-1).^(l_A+l_B+l_C) * (-1).^((m_C>0)*m_C);
%tmp_ss = (-1).^(l_A+l_B+l_C);
%tmp_ss = 1;
tmp_ss = (-1).^(l_A+l_B+l_C) * (-1).^( (m_A*m_C<0)*(m_A+m_C) + (m_A*m_C>=0)*(max(abs(m_A),abs(m_C))) );
output = tmp_wx*tmp_w0*tmp_zz*tmp_ss;

