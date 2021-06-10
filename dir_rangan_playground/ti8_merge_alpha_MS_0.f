      subroutine ti8_merge_alpha_MS_0(
     $     verbose
     $     ,n_S_tot
     $     ,n_S
     $     ,nS_add
     $     ,n_MS_max
     $     ,n_omp_sub_0in
     $     ,n_MS_omp__
     $     ,n_MS_
     $     ,alpha_MS_omp__
     $     ,alpha_MS__
     $     )
      implicit none
      integer verbose
      integer *4 n_S_tot
      integer *4 n_S
      integer *4 nS_add
      integer *4 n_MS_max
      integer *4 n_omp_sub_0in
      integer *4 n_MS_omp__(0:0)
      integer *4 n_MS_(0:0)
      real *8 alpha_MS_omp__(0:0)
      real *8 alpha_MS__(0:0)
      integer *4 nS,nS_tab
      integer *4 nM,n_MS
      integer *4 nomp,tab_pre,tab_cur,tab_pos
      include 'excerpt_define_nalpha.f'
      if (verbose.gt.1) then
         write(6,'(A,5(A,I0))') ' [entering ti8_merge_alpha_MS_0] '
     $        ,' n_S_tot: ' , n_S_tot
     $        ,' n_S: ' , n_S
     $        ,' nS_add: ' , nS_add
     $        ,' n_MS_max: ' , n_MS_max
     $        ,' n_omp_sub_0in: ' , n_omp_sub_0in
      end if !if (verbose.gt.1) then
      do nS=0,n_S-1
         nS_tab = nS + nS_add
         if (verbose.gt.1) then
            write(6,'(2(A,I0))')
     $           ' nS: ' , nS
     $           ,' nS_tab: ' , nS_tab
         end if !if (verbose.gt.1) then
         do nomp=0,n_omp_sub_0in-1
            tab_pre = (nS_tab + n_S_tot*nomp)*n_alpha*n_MS_max
            if (verbose.gt.1) then
               write(6,'(4(A,I0))')
     $              ' nS: ' , nS
     $              ,' nS_tab: ' , nS_tab
     $              ,' nomp: ' , nomp
     $              ,' tab_pre: ' , tab_pre
            end if !if (verbose.gt.1) then
            n_MS = n_MS_omp__(nS_tab + n_S_tot*nomp)
            tab_pos = tab_pre + n_MS*n_alpha-1
            if (verbose.gt.1) then
               write(6,'(6(A,I0))')
     $              ' nS: ' , nS
     $              ,' nS_tab: ' , nS_tab
     $              ,' nomp: ' , nomp
     $              ,' tab_pre: ' , tab_pre
     $              ,' n_MS: ' , n_MS
     $              ,' tab_pos: ' , tab_pos
            end if !if (verbose.gt.1) then
            do nM=0,n_MS-1
               tab_cur = tab_pre + nM*n_alpha
               if (verbose.gt.1) then
               write(6,'(2(A,I0))') ' nM: ' , nM , ' tab_cur: ' ,
     $              tab_cur
               call print_all_r8(n_alpha,alpha_MS_omp__(tab_cur)
     $              ,' alpha_MS_: ')
               end if !if (verbose.gt.1) then
               call alpha_SM_update_1(
     $              n_MS_max
     $              ,n_MS_(nS_tab)
     $              ,alpha_MS__(nS_tab*n_alpha*n_MS_max)
     $              ,alpha_MS_omp__(tab_cur)
     $              )
            enddo !do nM=0,n_MS-1
         enddo !do nomp=0,n_omp_sub_0in-1
      enddo !do nS=0,n_S-1
      if (verbose.gt.1) then
         write(6,'(A)') ' [finished ti8_merge_alpha_MS_0] '
      end if !if (verbose.gt.1) then
      end !subroutine
