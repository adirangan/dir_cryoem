      subroutine get_Jtaylor(l_max,n_J,f_,b_,p_)
c$$$      This calculates the terms in the taylor-expansion of 
c$$$      exp(-i*z*cos(phi)).
c$$$      In general, there will be l+1 terms of order l.
c$$$      The j^th term of order l (with j in 0,..,l) is of the form:
c$$$      f*b*exp(i*p*phi)*z^l , where
c$$$      f = (-i/2)^l / factorial(l)
c$$$      b = nchoosek(l,j)
c$$$      p = -l + 2*j
c$$$      
c$$$      Inputs:
c$$$      l_max = order of expansion 
c$$$      (e.g., l_max = 0 --> constant, 1 term)
c$$$      (e.g., l_max = 1 --> linear, 3 terms)
c$$$      (e.g., l_max = 2 --> quadratic, 6 terms)
c$$$
c$$$      Outputs:
c$$$      n_J = total number of terms (i.e., (l_max+1)*(l_max+2)/2)
c$$$      f_ = complex array of factors f
c$$$      b_ = integer array of binomial-coefficients b
c$$$      p_ = integer array of powers p
c$$$      
c$$$      Note: output arrays f_, b_ and p_ should be preallocated
c$$$      Note: index j on level l has index l*(l+1)/2+j, 
c$$$            assuming that the indexing runs from 0,..,n_J-1
      integer l_max,n_J,b_(*),p_(*)
      complex *16 f_(*)
c$$$      integer l_factorial,l_start
      integer jc,l,j,j_factorial,lmj_factorial
      complex *16 ci,f
      ci = ( 0.0 , -0.5 )
      f = ( 1.0 , 0.0 )
      n_J = ((l_max+1) * (l_max+2)) / 2
c$$$      l_factorial = 1
      jc = 1
      do l=0,l_max
         if (l.gt.0) then
            f = f*ci/max(1,l)
         end if
c$$$         l_start = (l * (l + 1)) / 2
c$$$         l_factorial = l_factorial * max(1,l)
         j_factorial = 1
         lmj_factorial = 1
         do j=0,l
            j_factorial = j_factorial * max(1,j)
c$$$            f_(1 + l_start + j) = f
c$$$            p_(1 + l_start + j) = -l + 2*j
c$$$            b_(1 + l_start + j) = lmj_factorial / j_factorial
            f_(jc) = f
            p_(jc) = -l + 2*j
            b_(jc) = lmj_factorial / j_factorial
            lmj_factorial = lmj_factorial * max(1,l-j)
            jc = jc + 1
         enddo
      enddo
      end
