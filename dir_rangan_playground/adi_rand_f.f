!> Doxygen comment: ;\n
!> create a pseudo-random real *8 from rseed: ;\n
!> this is useful when you want your ;\n
!> random numbers to be fixed from one build to another. ;\n
      real *8 function adi_rand_f(rseed)
      implicit none
      integer *4 rseed
      include 'adi_rand_f_excerpt_define_constant.f'
      rseed = mod(rseed * POW2RPOWPLUSRADD,POW22RPOWMINUSONE)
      adi_rand_f = max(0.0d0,min(1.0d0,max(1,rseed)
     $     *POW22RPOWMINUSONE_INVERSE))
c$$$      adi_rand_f = max(0.0d0,min(1.0d0,max(1,rseed)/(1.0d0
c$$$     $     *POW22RPOWMINUSONE)))
      end
