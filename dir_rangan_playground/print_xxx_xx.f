!> Doxygen comment: ;\n
!> several printing functions for various data types. ;\n
      subroutine print_all_i4(n_a,a_,prefix_string)
      implicit none
      integer n_a
      integer *4 a_(0:n_a-1)
      character(len=*) prefix_string
      character(len=1024) format_string
      integer na
      write(format_string,'(A,I0,A)') '(A, ' , n_a , '(I0,1X))'
      write(6,format_string) prefix_string , (a_(na),na=0,n_a-1)
      end

      subroutine print_all_i4__(n_r,n_c,a_,prefix_string)
      implicit none
      integer n_r,n_c
      integer *4 a_(0:n_r*n_c-1)
      character(len=*) prefix_string
      character(len=1024) format_string
      integer nr,nc
      write(format_string,'(A,I0,A)') '(A, ' , n_c , '(I0,1X))'
      do nr=0,n_r-1
         write(6,format_string) prefix_string , (a_(nc),nc=nr+0*n_r,nr
     $        +(n_c-1)*n_r,n_r)
      enddo ! do nr=0,n_r-1
      end

      subroutine print_sub_i4(n_a,a_,prefix_string)
      implicit none
      integer n_a
      integer *4 a_(0:n_a-1)
      character(len=*) prefix_string
      character(len=1024) format_string
      integer na
      if (n_a.lt.4) then
         call print_all_i4(n_a,a_,prefix_string)
      else
         write(format_string,'(A)') '(A,I0,1X,I0,1X,A,I0,A,I0,1X,I0) '
         write(6,format_string) prefix_string , a_(0) , a_(1) , ' ..['
     $        ,n_a,'].. ' , a_(n_a-2) , a_(n_a-1)
      end if ! if (n_a.lt.4) then
      end

      subroutine print_all_l2(n_a,a_,prefix_string)
      implicit none
      integer n_a
      logical a_(0:n_a-1)
      character(len=*) prefix_string
      character(len=1024) format_string
      integer na
      write(format_string,'(A,I0,A)') '(A, ' , n_a , '(L1,1X))'
      write(6,format_string) prefix_string , (a_(na),na=0,n_a-1)
      end

      subroutine print_all_l2__(n_r,n_c,a_,prefix_string)
      implicit none
      integer n_r,n_c
      logical a_(0:n_r*n_c-1)
      character(len=*) prefix_string
      character(len=1024) format_string
      integer nr,nc
      write(format_string,'(A,I0,A)') '(A, ' , n_c , '(L1,1X))'
      do nr=0,n_r-1
         write(6,format_string) prefix_string , (a_(nc),nc=nr+0*n_r,nr
     $        +(n_c-1)*n_r,n_r)
      enddo ! do nr=0,n_r-1
      end

      subroutine print_sub_l2(n_a,a_,prefix_string)
      implicit none
      integer n_a
      logical a_(0:n_a-1)
      character(len=*) prefix_string
      character(len=1024) format_string
      integer na
      if (n_a.lt.4) then
         call print_all_l2(n_a,a_,prefix_string)
      else
         write(format_string,'(A)') '(A,L1,1X,L1,1X,A,I0,A,L1,1X,L1) '
         write(6,format_string) prefix_string , a_(0) , a_(1) , ' ..['
     $        ,n_a,'].. ' , a_(n_a-2) , a_(n_a-1)
      end if ! if (n_a.lt.4) then
      end

      subroutine print_all_r8(n_a,a_,prefix_string)
      implicit none
      integer n_a
      real *8 a_(0:n_a-1)
      character(len=*) prefix_string
      character(len=1024) format_string
      integer na
      write(format_string,'(A,I0,A)') '(A, ' , n_a , '(F21.12,1X))'
      write(6,format_string) prefix_string , (a_(na),na=0,n_a-1)
      end

      subroutine print_sub_r8(n_a,a_,prefix_string)
      implicit none
      integer n_a
      real *8 a_(0:n_a-1)
      character(len=*) prefix_string
      character(len=1024) format_string
      integer na
      if (n_a.lt.4) then
         call print_all_r8(n_a,a_,prefix_string)
      else
         write(format_string,'(A)')
     $        '(A,F21.12,1X,F21.12,1X,A,I0,A,F21.12,1X,F21.12) '
         write(6,format_string) prefix_string , a_(0) , a_(1) , ' ..['
     $        ,n_a,'].. ' , a_(n_a-2) , a_(n_a-1)
      end if ! if (n_a.lt.4) then
      end

      subroutine print_mid_r8(n_a,a_,prefix_string)
      implicit none
      integer n_a
      real *8 a_(0:n_a-1)
      character(len=*) prefix_string
      character(len=1024) format_string
      integer na2,na
      if (n_a.lt.8) then
         call print_all_r8(n_a,a_,prefix_string)
      else
         na2 = n_a/2
         write(format_string,'(A,A)')
     $        '(A,F21.12,1X,F21.12,A,F21.12,1X,F21.12,'
     $        ,'F21.12,A,F21.12,1X,F21.12) '
         write(6,format_string) prefix_string , a_(0) , a_(1) , ' ... '
     $        ,a_(na2-1),a_(na2-0),a_(na2+1), ' ... ' , a_(n_a-2) ,
     $        a_(n_a-1)
      end if ! if (n_a.lt.4) then
      end

      subroutine print_all_c16(n_a,a_,prefix_string)
      implicit none
      integer n_a
      complex *16 a_(0:n_a-1)
      character(len=*) prefix_string
      character(len=1024) format_string
      integer na
      write(format_string,'(A,I0,A)') '(A, ' , n_a ,
     $     '(F21.12,1X,F21.12,4X))'
      write(6,format_string) prefix_string , (a_(na),na=0,n_a-1)
      end

      subroutine print_sub_c16(n_a,a_,prefix_string)
      implicit none
      integer n_a
      complex *16 a_(0:n_a-1)
      character(len=*) prefix_string
      character(len=1024) format_string
      integer na
      if (n_a.lt.4) then
         call print_all_c16(n_a,a_,prefix_string)
      else
         write(format_string,'(A,A,A,A,A,A,A)') '(A,'
     $        ,'A,F21.12,1X,F21.12,A,' ,'A,F21.12,1X,F21.12,A,' ,
     $        'A,I0,A,' , 'A,F21.12,1X,F21.12,A,' ,'A,F21.12,1X,F21.12
     $        ,A' , ') '
         write(6,format_string) prefix_string , ' (' , a_(0) , ') ' ,
     $        ' (' , a_(1) , ') ' , '..[',n_a,']..' , ' (' , a_(n_a-2) ,
     $        ') ' , ' (' , a_(n_a-1) , ') ' 
      end if ! if (n_a.lt.4) then
      end

      subroutine print_mid_c16(n_a,a_,prefix_string)
      implicit none
      integer n_a
      complex *16 a_(0:n_a-1)
      character(len=*) prefix_string
      character(len=1024) format_string
      integer na2,na
      if (n_a.lt.4) then
         call print_all_c16(n_a,a_,prefix_string)
      else
         na2 = n_a/2
         write(format_string,'(9A)') 
     $        '(A,'
     $        ,'A,F21.12,1X,F21.12,A,' 
     $        ,'A,' 
     $        ,'A,F21.12,1X,F21.12,A,' 
     $        ,'A,F21.12,1X,F21.12,A,' 
     $        ,'A,F21.12,1X,F21.12,A,' 
     $        ,'A,' 
     $        ,'A,F21.12,1X,F21.12,A' 
     $        , ') '
         write(6,format_string) 
     $        prefix_string 
     $        , ' (' 
     $        , a_(0) 
     $        , ') ' 
     $        , '...' 
     $        , ' (' 
     $        , a_(na2-1) 
     $        , ') ' 
     $        , ' (' 
     $        , a_(na2) 
     $        , ') ' 
     $        , ' (' 
     $        , a_(na2+1) 
     $        , ') ' 
     $        ,'...' 
     $        , ' (' 
     $        , a_(n_a-1) 
     $        , ') ' 
      end if ! if (n_a.lt.4) then
      end

