      subroutine test_innerproduct_3d_stage_0(verbose,n_beta,beta_,n_k
     $     ,k_,n_l_,a_,b_,X_ ,fftw_plan_frwd,fftw_plan_back
     $     ,fftw_in1_,fftw_out_)
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      integer *4 n_beta,n_k,n_l_(0:n_k-1)
      real *8 beta_(0:n_beta-1),k_(0:n_k-1)
      complex *16 a_(0:0),b_(0:0),X_(0:0)
c$$$      fftw plans and workspace should be passed in as input
c$$$      integer *8, allocatable :: fftw_plan_frwd_(:)
c$$$      integer *8, allocatable :: fftw_plan_back_(:)
c$$$      complex *16, allocatable :: fftw_in1__(:)
c$$$      complex *16, allocatable :: fftw_out__(:)
      integer *8 fftw_plan_frwd(0:0)
      integer *8 fftw_plan_back(0:0)
      complex *16 fftw_in1_(0:0)
      complex *16 fftw_out_(0:0)
c$$$      indices
      integer *4 max_i4_f
      integer *4 n_l_max,n_m_max,n_A,nl,nm
      integer *4 nbeta
      integer *4 nmn,mn,nmp,mp,nk,n_l,n_lm
      integer *4 a_base,b_base,a_mn,b_mp
      logical mn_flag,mp_flag
      real *8 beta,k
      complex *16 a,b,C
c$$$      wignerd matrix
      integer *4, allocatable :: n_wignerd_(:)
      integer *4 n_wignerd,nwignerd
      real *8, allocatable :: wignerd_(:)
      real *8 wignerd
      character(len=64) format_string

      if (verbose.gt.1) then
         write (6,'(A)') ' [entering test_innerproduct_3d_stage_0]'
      end if

      n_l_max = max_i4_f(n_k,n_l_)
      n_m_max = 1 + 2*n_l_max
      n_A = 0
      do nl=0,n_l_max
         nm = 1 + 2*nl
         n_A = n_A + nm*nm
      enddo !do nl=0,n_l_max

      if (verbose.gt.1) then
         write (6,'(A,I0)') ' n_A: ' , n_A
      end if

      allocate(n_wignerd_(0:n_l_max))
      allocate(wignerd_(0:n_A-1))
      do nbeta = 0,n_beta-1
         beta = beta_(nbeta)
         call wignerd_b(n_l_max,beta,wignerd_,n_wignerd_)
         if (verbose.gt.2) then
            write(format_string,'(A,I0,A)') '(',n_l_max+1,'(I0,1X))'
            write(6,format_string) (n_wignerd_(nwignerd),nwignerd=0
     $           ,n_l_max)
            do nl=0,n_l_max
               nm = 1+2*nl
               n_wignerd = nm*nm
               write(6,'(A,I0,A,I0,A,I0)') ' nl ' , nl , ' nm ' , nm ,
     $              ' n_wignerd ' , n_wignerd
               write(6,'(A,I0,A,I0)') ' n_wignerd_(nl) ' ,
     $              n_wignerd_(nl) ,' n_wignerd_(nl)+n_wignerd-1 ' ,
     $              n_wignerd_(nl)+n_wignerd-1
               write(format_string,'(A,I0,A)') '(',nm,'(F10.6,1X))'
               write(6,format_string) (wignerd_(nwignerd),nwignerd
     $              =n_wignerd_(nl),n_wignerd_(nl)+n_wignerd-1)
            enddo !do nl=0,n_l_max
         end if !if (verbose.gt.2) then
         do nmn = 0,n_m_max-1
            mn = -n_l_max + nmn
            if (verbose.gt.2) then
               write(6,'(A,I0,A,I0)') ' nmn ' , nmn , ' mn ' , mn
            end if !if (verbose.gt.2) then
            do nmp = 0,n_m_max-1
               mp = -n_l_max + nmp
               if (verbose.gt.2) then
                  write(6,'(A,I0,A,I0)') ' nmp ' , nmp , ' mp ' , mp
               end if !if (verbose.gt.2) then
               C = cmplx(0.0d0,0.0d0)
               a_base = 0
               b_base = 0
               do nk = 0,n_k-1
                  k = k_(nk)
                  n_l = n_l_(nk)
                  n_lm = (1+n_l)*(1+n_l)
                  if (verbose.gt.2) then
                     write(6,'(A,I0,A,F8.4)') ' nk ' , nk , ' k ' , k
                     write(6,'(A,I0,A,I0)') ' n_l ' , n_l , ' n_lm ' ,
     $                    n_lm
                     write(6,'(A,I0,A,I0)') ' a_base ' , a_base ,
     $                    ' b_base ' , b_base
                  end if !if (verbose.gt.2) then
                  if ((abs(mn).le.n_l) .and. (abs(mp).le.n_l)) then
                  do nl = 0,n_l
                     n_wignerd = n_wignerd_(nl)
                     mn_flag = .false.
                     if (abs(mn).le.nl) then
                        mn_flag = .true.
                        a_mn = nl*(1+nl) + mn
                     end if !if (abs(mn).le.nl) then
                     mp_flag = .false.
                     if (abs(mp).le.nl) then
                        mp_flag = .true.
                        b_mp = nl*(1+nl) + mp
                     end if !if (abs(mp).le.nl) then
                     if (verbose.gt.2) then
                        write(6,'(A,I0,A,I0,A,I0)') ' n_wignerd ' ,
     $                       n_wignerd , ' a_mn ' , a_mn , ' b_mp ' ,
     $                       b_mp
                        write(6,'(A,L2,A,L2)') ' mn_flag ' , mn_flag ,
     $                       ' mp_flag ' , mp_flag
                     end if     !if (verbose.gt.2) then
                     if (mn_flag .and. mp_flag) then
                        if (verbose.gt.3) then
                           write(6,'(A,I0,A,I0)') ' nmn ' , nmn , ' mn '
     $                          , mn
                        end if  !if (verbose.gt.3) then
                        if (verbose.gt.3) then
                           write(6,'(A,I0,A,I0)') ' nmp ' , nmp , ' mp '
     $                          , mp
                        end if  !if (verbose.gt.3) then
                        if (verbose.gt.3) then
                           write(6,'(A,I0,A,F8.4)') ' nk ' , nk , ' k '
     $                          , k
                           write(6,'(A,I0,A,I0)') ' n_l ' , n_l ,
     $                          ' n_lm ' ,n_lm
                           write(6,'(A,I0,A,I0)') ' a_base ' , a_base ,
     $                          ' b_base ' , b_base
                        end if  !if (verbose.gt.3) then
                        if (verbose.gt.3) then
                           write(6,'(A,I0,A,I0,A,I0)') ' n_wignerd ' ,
     $                          n_wignerd , ' a_mn ' , a_mn , ' b_mp ' ,
     $                          b_mp
                        end if  !if (verbose.gt.3) then
                        nwignerd = n_wignerd + (nl+mn) + (nl+mp)*(1+2
     $                       *nl)
                        wignerd = wignerd_(nwignerd)
                        a = a_(a_base+a_mn)
                        b = b_(b_base+b_mp)
                        if (verbose.gt.3) then
                           write(6,'(A,I0,A,F8.4)') ' nwignerd ' ,
     $                          nwignerd , ' wignerd ' ,
     $                          wignerd
                           write(6,'(A,2(F8.4,1X),A,2(F8.4,1X))') ' a '
     $                          , a , ' b ' , b
                        end if  !if (verbose.gt.3) then                        
                        C = C + k*k*dconjg(a)*wignerd*b
                        if (verbose.gt.3) then
                           write(6,'(A,2(F8.4,1X))') 'C: ' , C
                        end if !if (verbose.gt.3) then
                     end if !if (mn_flag .and. mp_flag) then
                  enddo !do nl = 0,n_l
                  end if !if ((abs(mn).le.n_l) .and. (abs(mp).le.n_l)) then
                  a_base = a_base + n_lm
                  b_base = b_base + n_lm
               enddo !do nk = 0,n_k-1
               fftw_in1_(nmn + nmp*(n_m_max)) = C
            enddo !do nmp = 0,n_m_max-1
         enddo ! do nmn = 0,n_m_max-1
         call recenter2_c16(n_m_max,n_m_max,fftw_in1_,fftw_in1_)
         call dfftw_execute_(fftw_plan_frwd)
         call recenter2_c16(n_m_max,n_m_max,fftw_out_,fftw_out_)
         call cp1_c16(n_m_max*n_m_max,fftw_out_,X_(nbeta*n_m_max
     $        *n_m_max))
      enddo !do nbeta = 0,n_beta-1

      deallocate(wignerd_)
      deallocate(n_wignerd_)

      if (verbose.gt.1) then
         write (6,'(A)') ' [finished test_innerproduct_3d_stage_0]'
      end if

      end
