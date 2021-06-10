!> Doxygen comment: ;\n
!> This header file defines some of the constants used in adi_rand_f.f. ;
      integer *4 POW2RPOWPLUSRADD
      parameter (POW2RPOWPLUSRADD=35)
      integer *4 POW22RPOWMINUSONE
c$$$      parameter (POW22RPOWMINUSONE=2147483647)
      parameter (POW22RPOWMINUSONE=8388607)
      real *8 POW22RPOWMINUSONE_INVERSE
      parameter (POW22RPOWMINUSONE_INVERSE=1.192093037616376592679d-7)
