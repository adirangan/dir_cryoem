      subroutine print_alpha_SM_0(n_M,n_SM_max,n_SM_,alpha_SM__
     $     ,prefix_string)
      implicit none
      integer *4 n_M
      integer *4 n_SM_max
      integer *4 n_SM_(0:0)
      real *8 alpha_SM__(0:0)
      character(len=*) prefix_string
      include 'excerpt_define_nalpha.f'
      character(len=1024) prefix_2
      character(len=1024) format_string
      integer nM,nS,nT
      write(prefix_2,'(A,A)') prefix_string , ' length: '      
      call print_all_i4(n_M,n_SM_,trim(prefix_2))
c$$$  %%%%%%%%
      do nM=0,n_M-1
         if (n_SM_(nM).gt.0) then
         do nS=0,n_SM_(nM)-1
            write(prefix_2,'(A,A,I0,A)') prefix_string , ' entry: ' , nS
     $           , ': '
            nT = n_alpha*(nS + n_SM_max*nM)
            write(6,'(A,I0)') ' nT: ' , nT
            call print_all_r8(n_alpha,alpha_SM__(nT),trim(prefix_2))
         enddo !do nS=0,n_SM_(nM)-1
         end if !if (n_SM_(nM).gt.0) then
      enddo !do nM=0,n_M-1
c$$$      %%%%%%%%
      write(prefix_2,'(A,A)') prefix_string , ' M_index: '
      do nM=0,n_M-1
         if (n_SM_(nM).gt.0) then
         nT = nalpha_M_index + n_alpha*n_SM_max*nM
         write(format_string,'(A,I0,A)') '(A, ' , n_SM_(nM) ,
     $        '(F24.16,1X))'
         write(6,format_string) trim(prefix_2) , (alpha_SM__(nS),nS
     $        =nT,nT+n_alpha*(n_SM_(nM)-1),n_alpha)
         end if !if (n_SM_(nM).gt.0) then
      enddo !do nM=0,n_M-1
c$$$      %%%%%%%%
      write(prefix_2,'(A,A)') prefix_string , ' S_index: '
      do nM=0,n_M-1
         if (n_SM_(nM).gt.0) then
         nT = nalpha_S_index + n_alpha*n_SM_max*nM
         write(format_string,'(A,I0,A)') '(A, ' , n_SM_(nM) ,
     $        '(F24.16,1X))'
         write(6,format_string) trim(prefix_2) , (alpha_SM__(nS),nS
     $        =nT,nT+n_alpha*(n_SM_(nM)-1),n_alpha)
         end if !if (n_SM_(nM).gt.0) then
      enddo !do nM=0,n_M-1
c$$$      %%%%%%%%
      write(prefix_2,'(A,A)') prefix_string , ' CTF_R_S: '
      do nM=0,n_M-1
         if (n_SM_(nM).gt.0) then
         nT = nalpha_CTF_R_S + n_alpha*n_SM_max*nM
         write(format_string,'(A,I0,A)') '(A, ' , n_SM_(nM) ,
     $        '(F24.16,1X))'
         write(6,format_string) trim(prefix_2) , (alpha_SM__(nS),nS
     $        =nT,nT+n_alpha*(n_SM_(nM)-1),n_alpha)
         end if !if (n_SM_(nM).gt.0) then
      enddo !do nM=0,n_M-1
c$$$      %%%%%%%%
      write(prefix_2,'(A,A)') prefix_string , ' C_Z_opt: '
      do nM=0,n_M-1
         if (n_SM_(nM).gt.0) then
         nT = nalpha_C_Z_opt + n_alpha*n_SM_max*nM
         write(format_string,'(A,I0,A)') '(A, ' , n_SM_(nM) ,
     $        '(F24.16,1X))'
         write(6,format_string) trim(prefix_2) , (alpha_SM__(nS),nS
     $        =nT,nT+n_alpha*(n_SM_(nM)-1),n_alpha)
         end if !if (n_SM_(nM).gt.0) then
      enddo !do nM=0,n_M-1
c$$$      %%%%%%%%

      end !subroutine
