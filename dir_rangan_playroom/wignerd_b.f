      subroutine wignerd_b(n_l,beta,W_,n_W_)
      implicit none
      integer verbose
      data verbose / 0 /
      integer *4 n_l
      real *8 beta
      real *8 W_(0:0)
      integer *4 n_W_(0:0)
      integer *4 nl,nlp,nlc,nmn,nmp,nw
      integer *4 w_base,w_ij,v_base,v_ij
      real *8 cb,sb
      real *8 W1,W2,W3,A1,A2,A3,B1,B2,B3,C1,C2,C3
      real *8 W_tmp
      real *8 smn,smp
      character(len=64) str_tmp
      cb = cos(beta/2)
      sb = sin(beta/2)
      n_W_(0) = 0
      W_(0) = 1.0d0
      do nl=1,n_l
         if (nl.eq.1) then
            n_W_(nl) = 1
         end if !if (nl.eq.1) then
         nlp = nl-1
         nlc = nl
         v_base = n_W_(nlp)
         w_base = n_W_(nlc)
         if (verbose.gt.0) then
            write(6,'(A)') ' '
            write(6,'(A,I0,A,I0,A,I0)') ' nl ' , nl , ' v_base ' ,
     $           v_base ,' w_base ' , w_base
            write(6,'(A)') ' '
         end if !verbose
         do nmn = -nlc,+nlc
            do nmp = -nlc,+nlc
               W_tmp = 0.0d0
c$$$               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$               Recurrence A
c$$$               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               if ((nlc.ne.-nmp) .and. (nlc.ne.1-nmp)) then
                  write(str_tmp,'(A)') 'A'
                  W1 = 0;
                  if ((abs(nmp-1).le.nlp) .and. (abs(nmn-1).le.nlp))
     $                 then
                     v_ij = nlp+(nmp-1) + (nlp+(nmn-1))*(1+2*nlp)
                     W1 = W_(v_base + v_ij)
                  end if !if ((abs(nmp-1).le.nlp) .and. (abs(nmn-1).le.nlp)) then
                  W2 = 0;
                  if ((abs(nmp-1).le.nlp) .and. (abs(nmn-0).le.nlp))
     $                 then
                     v_ij = nlp+(nmp-1) + (nlp+(nmn-0))*(1+2*nlp)
                     W2 = W_(v_base + v_ij)
                  end if !if ((abs(nmp-1).le.nlp) .and. (abs(nmn-0).le.nlp)) then
                  W3 = 0;
                  if ((abs(nmp-1).le.nlp) .and. (abs(nmn+1).le.nlp))
     $                 then
                     v_ij = nlp+(nmp-1) + (nlp+(nmn+1))*(1+2*nlp)
                     W3 = W_(v_base + v_ij)
                  end if !if ((abs(nmp-1).le.nlp) .and. (abs(nmn+1).le.nlp)) then
                  call wignerd_A1(nl,nmn,nmp,A1)
                  call wignerd_A2(nl,nmn,nmp,A2)
                  call wignerd_A3(nl,nmn,nmp,A3)
                  W_tmp = cb*cb*A1*W1 - 2*cb*sb*A2*W2 + sb*sb*A3*W3;
               end if !if ((nlc.ne.-nmp) .and. (nlc.ne.1-nmp)) then
c$$$               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$               Recurrence B
c$$$               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               if ((nlc.ne.+nmp) .and. (nlc.ne.1+nmp)) then
                  write(str_tmp,'(A)') 'B'
                  W1 = 0;
                  if ((abs(nmp+1).le.nlp) .and. (abs(nmn-1).le.nlp))
     $                 then
                     v_ij = nlp+(nmp+1) + (nlp+(nmn-1))*(1+2*nlp)
                     W1 = W_(v_base + v_ij)
                  end if !if ((abs(nmp+1).le.nlp) .and. (abs(nmn-1).le.nlp)) then
                  W2 = 0;
                  if ((abs(nmp+1).le.nlp) .and. (abs(nmn-0).le.nlp))
     $                 then
                     v_ij = nlp+(nmp+1) + (nlp+(nmn-0))*(1+2*nlp)
                     W2 = W_(v_base + v_ij)
                  end if !if ((abs(nmp+1).le.nlp) .and. (abs(nmn-0).le.nlp)) then
                  W3 = 0;
                  if ((abs(nmp+1).le.nlp) .and. (abs(nmn+1).le.nlp))
     $                 then
                     v_ij = nlp+(nmp+1) + (nlp+(nmn+1))*(1+2*nlp)
                     W3 = W_(v_base + v_ij)
                  end if !if ((abs(nmp+1).le.nlp) .and. (abs(nmn+1).le.nlp)) then
                  call wignerd_B1(nl,nmn,nmp,B1)
                  call wignerd_B2(nl,nmn,nmp,B2)
                  call wignerd_B3(nl,nmn,nmp,B3)
                  W_tmp = sb*sb*B1*W1 + 2*sb*cb*B2*W2 + cb*cb*B3*W3;
               end if !if ((nlc.ne.+nmp) .and. (nlc.ne.1+nmp)) then
c$$$               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$               Recurrence C
c$$$               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               if ((nlc.ne.-nmp) .and. (nlc.ne.+nmp)) then
                  write(str_tmp,'(A)') 'C'
                  W1 = 0;
                  if ((abs(nmp-0).le.nlp) .and. (abs(nmn-1).le.nlp))
     $                 then
                     v_ij = nlp+(nmp-0) + (nlp+(nmn-1))*(1+2*nlp)
                     W1 = W_(v_base + v_ij)
                  end if !if ((abs(nmp-0).le.nlp) .and. (abs(nmn-1).le.nlp)) then
                  W2 = 0;
                  if ((abs(nmp-0).le.nlp) .and. (abs(nmn-0).le.nlp))
     $                 then
                     v_ij = nlp+(nmp-0) + (nlp+(nmn-0))*(1+2*nlp)
                     W2 = W_(v_base + v_ij)
                  end if !if ((abs(nmp-0).le.nlp) .and. (abs(nmn-0).le.nlp)) then
                  W3 = 0;
                  if ((abs(nmp-0).le.nlp) .and. (abs(nmn+1).le.nlp))
     $                 then
                     v_ij = nlp+(nmp-0) + (nlp+(nmn+1))*(1+2*nlp)
                     W3 = W_(v_base + v_ij)
                  end if !if ((abs(nmp-0).le.nlp) .and. (abs(nmn+1).le.nlp)) then
                  call wignerd_C1(nl,nmn,nmp,C1)
                  call wignerd_C2(nl,nmn,nmp,C2)
                  call wignerd_C3(nl,nmn,nmp,C3)
                  W_tmp = sb*cb*C1*W1 + (cb*cb-sb*sb)*C2*W2 - sb*cb*C3
     $                 *W3;
               end if !if ((nlc.ne.-nmp) .and. (nlc.ne.+nmp)) then
               w_ij = (nlc+nmp) + (nlc+nmn)*(1+2*nlc)
               if (verbose.gt.0) then
                  write(6 ,'(A , I0 , A , I0 , A , F10.6 , A , A , A , F
     $10.6 , 1X , F10.6 , 1X ,F10.6)') ' setting W_(' , nlc+nmp , ',' ,
     $                 nlc + nmn , ') <-- ' , W_tmp ,' using ' ,
     $                 trim(str_tmp) , ' ' , W1 , W2 , W3
               end if !verbose 
               W_(w_base + w_ij) = W_tmp
            enddo !do nmp = -nlc,+nlc
         enddo !do nmn = -nlc,+nlc
         if (nl.lt.n_l) then
            n_W_(nlc+1) = n_W_(nlc) + (1+2*nlc)*(1+2*nlc)
         end if !if (nl.lt.n_l) then
      enddo !do nl=1,n_l

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      Fix phase
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      do nl=1,n_l
         nlc = nl
         w_base = n_W_(nlc)
         do nmn = -nlc,+nlc
            smn = 1.0d0
            if ((nmn.lt.0) .and. (abs(mod(nmn,2)).eq.1)) then
               smn = -1.0d0
            end if !if ((nmn.lt.0) .and. (abs(mod(nmn,2)).eq.1)) then
            do nmp = -nlc,+nlc
               smp = 1.0d0
               if ((nmp.lt.0) .and. (abs(mod(nmp,2)).eq.1)) then
                  smp = -1.0d0
               end if !if ((nmp.lt.0) .and. (abs(mod(nmp,2)).eq.1)) then
               w_ij = (nlc+nmp) + (nlc+nmn)*(1+2*nlc)
               W_(w_base + w_ij) = W_(w_base + w_ij)*smn*smp
            enddo !do nmp = -nlc,+nlc
         enddo !do nmn = -nlc,+nlc
      enddo !do nl=1,n_l
      
      end

      subroutine wignerd_A1(nl,nmn,nmp,output)
      implicit none
      integer *4 nl,nmn,nmp
      real *8 output
      real *8 numerator,denominator
      numerator   = (nl+nmn)*(nl+nmn-1)
      denominator = (nl+nmp)*(nl+nmp-1)
      output = dsqrt(numerator/denominator)
      end

      subroutine wignerd_A2(nl,nmn,nmp,output)
      implicit none
      integer *4 nl,nmn,nmp
      real *8 output
      real *8 numerator,denominator
      numerator   = (nl+nmn)*(nl-nmn)
      denominator = (nl+nmp)*(nl+nmp-1)
      output = dsqrt(numerator/denominator)
      end

      subroutine wignerd_A3(nl,nmn,nmp,output)
      implicit none
      integer *4 nl,nmn,nmp
      real *8 output
      real *8 numerator,denominator
      numerator   = (nl-nmn)*(nl-nmn-1)
      denominator = (nl+nmp)*(nl+nmp-1)
      output = dsqrt(numerator/denominator)
      end

      subroutine wignerd_B1(nl,nmn,nmp,output)
      implicit none
      integer *4 nl,nmn,nmp
      real *8 output
      real *8 numerator,denominator
      numerator   = (nl+nmn)*(nl+nmn-1)
      denominator = (nl-nmp)*(nl-nmp-1)
      output = dsqrt(numerator/denominator)
      end

      subroutine wignerd_B2(nl,nmn,nmp,output)
      implicit none
      integer *4 nl,nmn,nmp
      real *8 output
      real *8 numerator,denominator
      numerator   = (nl+nmn)*(nl-nmn)
      denominator = (nl-nmp)*(nl-nmp-1)
      output = dsqrt(numerator/denominator)
      end

      subroutine wignerd_B3(nl,nmn,nmp,output)
      implicit none
      integer *4 nl,nmn,nmp
      real *8 output
      real *8 numerator,denominator
      numerator   = (nl-nmn)*(nl-nmn-1)
      denominator = (nl-nmp)*(nl-nmp-1)
      output = dsqrt(numerator/denominator)
      end

      subroutine wignerd_C1(nl,nmn,nmp,output)
      implicit none
      integer *4 nl,nmn,nmp
      real *8 output
      real *8 numerator,denominator
      numerator   = (nl+nmn)*(nl+nmn-1)
      denominator = (nl-nmp)*(nl+nmp)
      output = dsqrt(numerator/denominator)
      end

      subroutine wignerd_C2(nl,nmn,nmp,output)
      implicit none
      integer *4 nl,nmn,nmp
      real *8 output
      real *8 numerator,denominator
      numerator   = (nl+nmn)*(nl-nmn)
      denominator = (nl-nmp)*(nl+nmp)
      output = dsqrt(numerator/denominator)
      end

      subroutine wignerd_C3(nl,nmn,nmp,output)
      implicit none
      integer *4 nl,nmn,nmp
      real *8 output
      real *8 numerator,denominator
c$$$      numerator   = (nl-nmn)*(nl-nmn+1)
      numerator   = (nl-nmn)*(nl-nmn-1)
      denominator = (nl-nmp)*(nl+nmp)
      output = dsqrt(numerator/denominator)
      end
      




