      subroutine write_all_i4(n_a,a_,n_p,prefix_string)
      implicit none
      integer n_a,n_p
      integer *4 a_(0:n_a-1)
      character(len=n_p) prefix_string
      character(len=1024) format_string
      integer na
      write(format_string,'(A,I0,A)') '(A, ' , n_a , '(I0,1X))'
      write(6,format_string) prefix_string , (a_(na),na=0,n_a-1)
      end

      subroutine write_all_i4__(n_r,n_c,a_,n_p,prefix_string)
      implicit none
      integer n_r,n_c,n_p
      integer *4 a_(0:n_r*n_c-1)
      character(len=n_p) prefix_string
      character(len=1024) format_string
      integer nr,nc
      write(format_string,'(A,I0,A)') '(A, ' , n_c , '(I0,1X))'
      do nr=0,n_r-1
         write(6,format_string) prefix_string , (a_(nc),nc=nr+0*n_r,nr
     $        +(n_c-1)*n_r,n_r)
      enddo ! do nr=0,n_r-1
      end

      subroutine write_sub_i4(n_a,a_,n_p,prefix_string)
      implicit none
      integer n_a,n_p
      integer *4 a_(0:n_a-1)
      character(len=n_p) prefix_string
      character(len=1024) format_string
      integer na
      if (n_a.lt.4) then
         call write_all_i4(n_a,a_,n_p,prefix_string)
      else
         write(format_string,'(A)') '(A,I0,1X,I0,1X,A,I0,1X,I0) '
         write(6,format_string) prefix_string , a_(0) , a_(1) ,
     $        ' ... ' , a_(n_a-2) , a_(n_a-1)
      end if ! if (n_a.lt.4) then
      end

      subroutine write_all_l2(n_a,a_,n_p,prefix_string)
      implicit none
      integer n_a,n_p
      logical a_(0:n_a-1)
      character(len=n_p) prefix_string
      character(len=1024) format_string
      integer na
      write(format_string,'(A,I0,A)') '(A, ' , n_a , '(L1,1X))'
      write(6,format_string) prefix_string , (a_(na),na=0,n_a-1)
      end

      subroutine write_all_l2__(n_r,n_c,a_,n_p,prefix_string)
      implicit none
      integer n_r,n_c,n_p
      logical a_(0:n_r*n_c-1)
      character(len=n_p) prefix_string
      character(len=1024) format_string
      integer nr,nc
      write(format_string,'(A,I0,A)') '(A, ' , n_c , '(L1,1X))'
      do nr=0,n_r-1
         write(6,format_string) prefix_string , (a_(nc),nc=nr+0*n_r,nr
     $        +(n_c-1)*n_r,n_r)
      enddo ! do nr=0,n_r-1
      end

      subroutine write_sub_l2(n_a,a_,n_p,prefix_string)
      implicit none
      integer n_a,n_p
      logical a_(0:n_a-1)
      character(len=n_p) prefix_string
      character(len=1024) format_string
      integer na
      if (n_a.lt.4) then
         call write_all_l2(n_a,a_,n_p,prefix_string)
      else
         write(format_string,'(A)') '(A,L1,1X,L1,1X,A,L1,1X,L1) '
         write(6,format_string) prefix_string , a_(0) , a_(1) ,
     $        ' ... ' , a_(n_a-2) , a_(n_a-1)
      end if ! if (n_a.lt.4) then
      end

      subroutine write_all_r8(n_a,a_,n_p,prefix_string)
      implicit none
      integer n_a,n_p
      real *8 a_(0:n_a-1)
      character(len=n_p) prefix_string
      character(len=1024) format_string
      integer na
      write(format_string,'(A,I0,A)') '(A, ' , n_a , '(F10.5,1X))'
      write(6,format_string) prefix_string , (a_(na),na=0,n_a-1)
      end

      subroutine write_sub_r8(n_a,a_,n_p,prefix_string)
      implicit none
      integer n_a,n_p
      real *8 a_(0:n_a-1)
      character(len=n_p) prefix_string
      character(len=1024) format_string
      integer na
      if (n_a.lt.4) then
         call write_all_r8(n_a,a_,n_p,prefix_string)
      else
         write(format_string,'(A)')
     $        '(A,F16.8,1X,F16.8,1X,A,F16.8,1X,F16.8) '
         write(6,format_string) prefix_string , a_(0) , a_(1) ,
     $        ' ... ' , a_(n_a-2) , a_(n_a-1)
      end if ! if (n_a.lt.4) then
      end


      subroutine write_all_c16(n_a,a_,n_p,prefix_string)
      implicit none
      integer n_a,n_p
      complex *16 a_(0:n_a-1)
      character(len=n_p) prefix_string
      character(len=1024) format_string
      integer na
      write(format_string,'(A,I0,A)') '(A, ' , n_a ,
     $     '(F16.8,1X,F16.8,4X))'
      write(6,format_string) prefix_string , (a_(na),na=0,n_a-1)
      end

      subroutine write_sub_c16(n_a,a_,n_p,prefix_string)
      implicit none
      integer n_a,n_p
      complex *16 a_(0:n_a-1)
      character(len=n_p) prefix_string
      character(len=1024) format_string
      integer na
      if (n_a.lt.4) then
         call write_all_c16(n_a,a_,n_p,prefix_string)
      else
         write(format_string,'(A,A,A,A,A,A,A)') '(A,'
     $        ,'A,F16.8,1X,F16.8,A,' ,'A,F16.8,1X,F16.8,A,' , 'A,' , 
     $        'A,F16.8,1X,F16.8,A,' ,'A,F16.8,1X,F16.8,A' , ') '
         write(6,format_string) prefix_string , ' (' , a_(0) , ') ' 
     $        , ' (' , a_(1) , ') ' , '...' , 
     $        ' (' , a_(n_a-2) , ') ' , ' (' , a_(n_a-1) , ') ' 
      end if ! if (n_a.lt.4) then
      end

c$$$

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
         write(format_string,'(A)') '(A,I0,1X,I0,1X,A,I0,1X,I0) '
         write(6,format_string) prefix_string , a_(0) , a_(1) ,
     $        ' ... ' , a_(n_a-2) , a_(n_a-1)
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
         write(format_string,'(A)') '(A,L1,1X,L1,1X,A,L1,1X,L1) '
         write(6,format_string) prefix_string , a_(0) , a_(1) ,
     $        ' ... ' , a_(n_a-2) , a_(n_a-1)
      end if ! if (n_a.lt.4) then
      end

      subroutine print_all_r8(n_a,a_,prefix_string)
      implicit none
      integer n_a
      real *8 a_(0:n_a-1)
      character(len=*) prefix_string
      character(len=1024) format_string
      integer na
      write(format_string,'(A,I0,A)') '(A, ' , n_a , '(F10.5,1X))'
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
     $        '(A,F16.8,1X,F16.8,1X,A,F16.8,1X,F16.8) '
         write(6,format_string) prefix_string , a_(0) , a_(1) ,
     $        ' ... ' , a_(n_a-2) , a_(n_a-1)
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
     $     '(F16.8,1X,F16.8,4X))'
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
     $        ,'A,F16.8,1X,F16.8,A,' ,'A,F16.8,1X,F16.8,A,' , 'A,' , 
     $        'A,F16.8,1X,F16.8,A,' ,'A,F16.8,1X,F16.8,A' , ') '
         write(6,format_string) prefix_string , ' (' , a_(0) , ') ' 
     $        , ' (' , a_(1) , ') ' , '...' , 
     $        ' (' , a_(n_a-2) , ') ' , ' (' , a_(n_a-1) , ') ' 
      end if ! if (n_a.lt.4) then
      end

