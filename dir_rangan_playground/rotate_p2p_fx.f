!> Doxygen comment: ;\n
!> Applies rotation in ?_p (i.e., polar) coordinates. ;\n
!> This uses fourier interpolation. ;\n
!> Assumes that M_p_ is the same size and dimensions as S_p_. ;\n
      subroutine rotate_p2p_fx(n_r,fftw_plan_frwd_,fftw_plan_back_,n_w_
     $     ,n_A,fftw_0in_,fftw_out_,S_p_,gamma,M_p_)
c$$$      Assumes that M_p_ is the same size and dimensions as S_p_ ;
c$$$      Also assumes that fftw plans are passed in. ;
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 gamma
      integer *8 fftw_plan_frwd_(0:n_r-1)
      integer *8 fftw_plan_back_(0:n_r-1)
      complex *16 fftw_0in_(0:n_A-1),fftw_out_(0:n_A-1)
      complex *16 S_p_(0:n_A-1),M_p_(0:n_A-1)
      complex *16, allocatable :: C_(:)
      complex *16 C
      real *8 pi;
      integer nr,ic,icstart,nw,n_w_max,nq,q;
      real *8 al2_c16_f
      real *8 l2_pre,l2_pos

      ic = 0
      do nr=0,n_r-1
         ic = ic + n_w_(nr)
      enddo !do nr=0,n_r-1
      if (ic.ne.n_A) then
         write(6,'(2(A,I0))') ' Warning, n_A: ' , n_A , ' ic: ' , ic
      end if !if (ic.ne.n_A) then

      n_w_max = n_w_(n_r-1)
      allocate(C_(0:n_w_max-1));
      do nw=0,n_w_max-1
         q = nw - n_w_max/2
         C = dcmplx( +dcos(q*gamma) , -dsin(q*gamma) )
         C_(nw) = C
      enddo !do nw=0,n_w_max-1

      ic=0
      do nr=0,n_r-1
         if (verbose.gt.0) then
            write(6,'(2(A,I0))') ' nr: ' , nr , ' n_w_(nr): ' , n_w_(nr)
         end if !if (verbose.gt.0) then
         if (n_w_(nr).gt.0) then            
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling cp1_c16: '
            end if !if (verbose.gt.2) then
            call cp1_c16(n_w_(nr),S_p_(ic),fftw_0in_(ic))
            if (verbose.gt.0) then
               write(6,'(A,F16.8)') ' fftw_0in_ l2: ' ,
     $              al2_c16_f(n_w_(nr),fftw_0in_(ic))
            end if !if (verbose.gt.0) then
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; executing fftw_plan: '
            end if !if (verbose.gt.2) then
            call dfftw_execute_(fftw_plan_frwd_(nr))
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling affine_c16: '
            end if !if (verbose.gt.2) then
            call af1_c16(n_w_(nr),1.0d0*dcmplx(1.0d0/dsqrt(1.0d0
     $           *n_w_(nr))),1.0d0*dcmplx(0.0d0),fftw_out_(ic)
     $           ,fftw_out_(ic))
            if (verbose.gt.0) then
               write(6,'(A,F16.8)') ' fftw_out_ l2: ' ,
     $              al2_c16_f(n_w_(nr),fftw_out_(ic))
            end if !if (verbose.gt.0) then
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; applying rotation: '
            end if !if (verbose.gt.2) then
            do nq=0,n_w_(nr)-1
               q = nq
               if (q.gt.n_w_(nr)/2-1) then 
                  q = q - n_w_(nr)
               end if !if (q.ge.n_w_(nr)/2-1) then 
               C = C_(n_w_max/2 + q)
               fftw_out_(ic + nq) = fftw_out_(ic + nq) * C
            enddo !do nq=0,n_w_(nr)-1
            if (verbose.gt.0) then
               write(6,'(A,F16.8)') ' fftw_out_ l2: ' ,
     $              al2_c16_f(n_w_(nr),fftw_out_(ic))
            end if !if (verbose.gt.0) then
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; executing fftw_plan: '
            end if !if (verbose.gt.2) then
            call dfftw_execute_(fftw_plan_back_(nr))
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling affine_c16: '
            end if !if (verbose.gt.2) then
            call af1_c16(n_w_(nr),1.0d0*dcmplx(1.0d0/dsqrt(1.0d0
     $           *n_w_(nr))),1.0d0*dcmplx(0.0d0),fftw_0in_(ic)
     $           ,fftw_0in_(ic))
            if (verbose.gt.0) then
               write(6,'(A,F16.8)') ' fftw_0in_ l2: ' ,
     $              al2_c16_f(n_w_(nr),fftw_0in_(ic))
            end if !if (verbose.gt.0) then
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling cp1_c16: '
            end if !if (verbose.gt.2) then
            call cp1_c16(n_w_(nr),fftw_0in_(ic),M_p_(ic))
            if (verbose.gt.2) then
               write(6,*) '% % setting ic ',ic,' to ic ',ic + n_w_(nr)
            end if !if (verbose.gt.2) then
            ic = ic + n_w_(nr)
         end if !if (n_w_(nr).gt.0) then
      enddo !do nr=0,n_r-1
      
      deallocate(C_)

      end

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

      subroutine rotate_p2p_fz(n_r,n_w_,n_A,S_p_,gamma,M_p_)
c$$$      Assumes that M_p_ is the same size and dimensions as S_p_ ;
c$$$      Generates fftw plans internally. ;
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      data verbose / 0 /
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 gamma
      integer *8, allocatable :: fftw_plan_frwd_(:)
      integer *8, allocatable :: fftw_plan_back_(:)
      complex *16, allocatable :: fftw_0in_(:)
      complex *16, allocatable :: fftw_out_(:)
      complex *16 S_p_(0:n_A-1),M_p_(0:n_A-1)
      complex *16, allocatable :: C_(:)
      complex *16 C
      real *8 pi;
      integer na,nr,ic,icstart,nw,n_w_max,nq,q;
      real *8 al2_c16_f
      real *8 l2_pre,l2_pos

      if (verbose.gt.0) then
         write(6,'(A)') ' [entering rotate_p2p_fz] '
      end if !if (verbose.gt.0) then

      na = 0
      do nr=0,n_r-1
         na = na + n_w_(nr)
      enddo !do nr=0,n_r-1
      if (na.ne.n_A) then
         write(6,'(2(A,I0))') ' Warning, n_A: ' , n_A , ' na: ' , na
      end if !if (na.ne.n_A) then

      if (verbose.gt.1) then
         write(6,'(A)') ' Generating fftw_plans for local use'
      end if !if (verbose.gt.1) then
      allocate(fftw_plan_frwd_(0:n_r-1))
      allocate(fftw_plan_back_(0:n_r-1))
      allocate(fftw_0in_(0:n_A-1))
      allocate(fftw_out_(0:n_A-1))
      na = 0
      do nr=0,n_r-1
         call dfftw_plan_dft_1d_(fftw_plan_frwd_(nr),n_w_(nr)
     $        ,fftw_0in_(na),fftw_out_(na),FFTW_FORWARD,FFTW_MEASURE) 
         call dfftw_plan_dft_1d_(fftw_plan_back_(nr),n_w_(nr)
     $        ,fftw_out_(na),fftw_0in_(na),FFTW_BACKWARD,FFTW_MEASURE) 
         na = na + n_w_(nr)
      enddo !do nr=0,n_r-1

      if (verbose.gt.1) then
         write(6,'(A)') ' generating complex exponentials '
      end if !if (verbose.gt.1) then
      n_w_max = n_w_(n_r-1)
      allocate(C_(0:n_w_max-1));
      do nw=0,n_w_max-1
         q = nw - n_w_max/2
         C = dcmplx( +dcos(q*gamma) , -dsin(q*gamma) )
         C_(nw) = C
      enddo !do nw=0,n_w_max-1

      ic=0
      do nr=0,n_r-1
         if (verbose.gt.0) then
            write(6,'(2(A,I0))') ' nr: ' , nr , ' n_w_(nr): ' , n_w_(nr)
         end if !if (verbose.gt.0) then
         if (n_w_(nr).gt.0) then            
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling cp1_c16: '
            end if !if (verbose.gt.2) then
            call cp1_c16(n_w_(nr),S_p_(ic),fftw_0in_(ic))
            if (verbose.gt.0) then
               write(6,'(A,F16.8)') ' fftw_0in_ l2: ' ,
     $              al2_c16_f(n_w_(nr),fftw_0in_(ic))
            end if !if (verbose.gt.0) then
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; executing fftw_plan: '
            end if !if (verbose.gt.2) then
            call dfftw_execute_(fftw_plan_frwd_(nr))
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling affine_c16: '
            end if !if (verbose.gt.2) then
            call af1_c16(n_w_(nr),1.0d0*dcmplx(1.0d0/dsqrt(1.0d0
     $           *n_w_(nr))),1.0d0*dcmplx(0.0d0),fftw_out_(ic)
     $           ,fftw_out_(ic))
            if (verbose.gt.0) then
               write(6,'(A,F16.8)') ' fftw_out_ l2: ' ,
     $              al2_c16_f(n_w_(nr),fftw_out_(ic))
            end if !if (verbose.gt.0) then
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; applying rotation: '
            end if !if (verbose.gt.2) then
            do nq=0,n_w_(nr)-1
               q = nq
               if (q.gt.n_w_(nr)/2-1) then 
                  q = q - n_w_(nr)
               end if !if (q.ge.n_w_(nr)/2-1) then 
               C = C_(n_w_max/2 + q)
               fftw_out_(ic + nq) = fftw_out_(ic + nq) * C
            enddo !do nq=0,n_w_(nr)-1
            if (verbose.gt.0) then
               write(6,'(A,F16.8)') ' fftw_out_ l2: ' ,
     $              al2_c16_f(n_w_(nr),fftw_out_(ic))
            end if !if (verbose.gt.0) then
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; executing fftw_plan: '
            end if !if (verbose.gt.2) then
            call dfftw_execute_(fftw_plan_back_(nr))
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling affine_c16: '
            end if !if (verbose.gt.2) then
            call af1_c16(n_w_(nr),1.0d0*dcmplx(1.0d0/dsqrt(1.0d0
     $           *n_w_(nr))),1.0d0*dcmplx(0.0d0),fftw_0in_(ic)
     $           ,fftw_0in_(ic))
            if (verbose.gt.0) then
               write(6,'(A,F16.8)') ' fftw_0in_ l2: ' ,
     $              al2_c16_f(n_w_(nr),fftw_0in_(ic))
            end if !if (verbose.gt.0) then
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling cp1_c16: '
            end if !if (verbose.gt.2) then
            call cp1_c16(n_w_(nr),fftw_0in_(ic),M_p_(ic))
            if (verbose.gt.2) then
               write(6,*) '% % setting ic ',ic,' to ic ',ic + n_w_(nr)
            end if !if (verbose.gt.2) then
            ic = ic + n_w_(nr)
         end if !if (n_w_(nr).gt.0) then
      enddo !do nr=0,n_r-1
      
      deallocate(C_)

      if (verbose.gt.1) then
         write(6,'(A)') ' Destroying fftw_plans for local use.'
      end if !if (verbose.gt.1) then
      do nr=0,n_r-1
         call dfftw_destroy_plan(fftw_plan_frwd_(nr))
         call dfftw_destroy_plan(fftw_plan_back_(nr))
      enddo !do nr=0,n_r-1
      deallocate(fftw_plan_frwd_)
      deallocate(fftw_plan_back_)
      deallocate(fftw_0in_)
      deallocate(fftw_out_)

      if (verbose.gt.0) then
         write(6,'(A)') ' [finished rotate_p2p_fz] '
      end if !if (verbose.gt.0) then

      end

c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$% Notes: test with the following in ti8_dr.f: ;
c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$
c$$$      ns=0
c$$$      nS_sample = I_S_sample_(ns)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      call cp1_c16(n_A,S_k_p__(nS_sample*n_A),T0_k_p_)
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_A
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_T0_.mda'
c$$$      call MDA_write_c16(MDA_n_d,MDA_d_,T0_k_p_,MDA_string)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      gamma_z = pi/3;
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      call rotate_p_to_p(n_r,n_w_,n_A,T0_k_p_,+gamma_z,T1_k_p_)
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_A
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_T1_.mda'
c$$$      call MDA_write_c16(MDA_n_d,MDA_d_,T1_k_p_,MDA_string)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      call rotate_p2p_l6(n_r,n_w_,n_A,T0_k_p_,+gamma_z,T2_k_p_)
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_A
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_T2_.mda'
c$$$      call MDA_write_c16(MDA_n_d,MDA_d_,T2_k_p_,MDA_string)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      call rotate_p2p_fx(n_r,fftw_plan_frwd_,fftw_plan_back_,n_w_
c$$$     $     ,n_A,fftw_0in_,fftw_out_,T0_k_p_,+gamma_z,T3_k_p_)
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_A
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_T3_.mda'
c$$$      call MDA_write_c16(MDA_n_d,MDA_d_,T3_k_p_,MDA_string)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      call rotate_p_to_p(n_r,n_w_,n_A,T1_k_p_,-gamma_z,T0_k_p_)
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_A
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_U1_.mda'
c$$$      call MDA_write_c16(MDA_n_d,MDA_d_,T0_k_p_,MDA_string)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      call rotate_p2p_l6(n_r,n_w_,n_A,T2_k_p_,-gamma_z,T0_k_p_)
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_A
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_U2_.mda'
c$$$      call MDA_write_c16(MDA_n_d,MDA_d_,T0_k_p_,MDA_string)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      call rotate_p2p_fx(n_r,fftw_plan_frwd_,fftw_plan_back_,n_w_
c$$$     $     ,n_A,fftw_0in_,fftw_out_,T3_k_p_,-gamma_z,T0_k_p_)
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_A
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_U3_.mda'
c$$$      call MDA_write_c16(MDA_n_d,MDA_d_,T0_k_p_,MDA_string)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_r
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_n_w_.mda'
c$$$      call MDA_write_i4(MDA_n_d,MDA_d_,n_w_,MDA_string)
c$$$      goto 10

c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$% Notes: along with the following matlab code.; 
c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$
c$$$n_w_ = MDA_read_i4('dir_mda/tmp_n_w_.mda'); 
c$$$n_r = length(n_w_); grid_p_ = 1:n_r; n_A = sum(n_w_);
c$$$T0_ = MDA_read_c16('dir_mda/tmp_T0_.mda');
c$$$T1_ = MDA_read_c16('dir_mda/tmp_T1_.mda');
c$$$T2_ = MDA_read_c16('dir_mda/tmp_T2_.mda');
c$$$T3_ = MDA_read_c16('dir_mda/tmp_T3_.mda');
c$$$U1_ = MDA_read_c16('dir_mda/tmp_U1_.mda');
c$$$U2_ = MDA_read_c16('dir_mda/tmp_U2_.mda');
c$$$U3_ = MDA_read_c16('dir_mda/tmp_U3_.mda');
c$$$
c$$$disp_flag=0;
c$$$if disp_flag;
c$$$%%%%%%%%%%%%%%%%;
c$$$figure(1);
c$$$pcols = 4;
c$$$np=1;
c$$$%%%%%%%%%%%%%%%%;
c$$$subplot(2,pcols,np+0*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,real(T0_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('real(T0_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$subplot(2,pcols,np+1*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,imag(T0_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('imag(T0_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$np=np+1;
c$$$%%%%%%%%%%%%%%%%;
c$$$subplot(2,pcols,np+0*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,real(T1_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('real(T1_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$subplot(2,pcols,np+1*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,imag(T1_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('imag(T1_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$np=np+1;
c$$$%%%%%%%%%%%%%%%%;
c$$$subplot(2,pcols,np+0*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,real(T2_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('real(T2_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$subplot(2,pcols,np+1*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,imag(T2_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('imag(T2_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$np=np+1;
c$$$%%%%%%%%%%%%%%%%;
c$$$subplot(2,pcols,np+0*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,real(T3_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('real(T3_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$subplot(2,pcols,np+1*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,imag(T3_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('imag(T3_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$np=np+1;
c$$$%%%%%%%%%%%%%%%%;
c$$$end;%if disp_flag;
c$$$
c$$$E1_ = T0_ - U1_;
c$$$E2_ = T0_ - U2_;
c$$$E3_ = T0_ - U3_;
c$$$
c$$$disp_flag=1;
c$$$if disp_flag;
c$$$%%%%%%%%%%%%%%%%;
c$$$figure(1);
c$$$pcols = 3;
c$$$np=1;
c$$$%%%%%%%%%%%%%%%%;
c$$$subplot(2,pcols,np+0*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,real(E1_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('real(E1_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$subplot(2,pcols,np+1*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,imag(E1_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('imag(E1_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$np=np+1;
c$$$%%%%%%%%%%%%%%%%;
c$$$subplot(2,pcols,np+0*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,real(E2_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('real(E2_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$subplot(2,pcols,np+1*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,imag(E2_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('imag(E2_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$np=np+1;
c$$$%%%%%%%%%%%%%%%%;
c$$$subplot(2,pcols,np+0*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,real(E3_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('real(E3_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$subplot(2,pcols,np+1*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,imag(E3_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('imag(E3_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$np=np+1;
c$$$%%%%%%%%%%%%%%%%;
c$$$end;%if disp_flag;
