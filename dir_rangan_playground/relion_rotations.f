
      subroutine relion_euler_rot_mat(phi,theta,psi,RMat)
      
      implicit none
      
c     RELION rotation angles
c     rotation around z axis in lab frame by angle phi
c     0 < phi < 2.0d0*pi
c     the rotation results in new axis x',y',z'
      real *8 phi
c     rotation around the new y' axis by angle theta
c     0 < theta < pi
c     the rotation results in new axis x'',y'',z''
      real *8 theta
c     rotation around the new z'' axis by angle psi
c     0 < psi < 2.0d0*pi
      real *8 psi
c     output is the rotation matrix
      real *8 RMat(3,3)


      real *8 RotMat(3,3)
      real *8 TRotMat(3,3)

      integer i,j

      RotMat(1,1) = dcos(psi)*dcos(theta)*dcos(phi)-dsin(psi)*dsin(phi)
      RotMat(2,1) = -dsin(psi)*dcos(theta)*dcos(phi)-dcos(psi)*dsin(phi)
      RotMat(3,1) = dsin(theta)*dcos(phi)
      RotMat(1,2) = dcos(psi)*dcos(theta)*dsin(phi)+dsin(psi)*dcos(phi)
      RotMat(2,2) = -dsin(psi)*dcos(theta)*dsin(phi)+dcos(psi)*dcos(phi)
      RotMat(3,2) = dsin(theta)*dsin(phi)
      RotMat(1,3) = -dcos(psi)*dsin(theta)
      RotMat(2,3) = dsin(psi)*dsin(theta)
      RotMat(3,3) = dcos(theta)

c     transpose
      do i=1,3
         do j = 1,3
            TRotMat(i,j) = RotMat(j,i)
         enddo
      enddo

      do i=1,3
         do j = 1,3
c            RMat(i,j) = RotMat(i,j)
            RMat(i,j) = TRotMat(i,j)
         enddo
      enddo

      

      return
      end

      subroutine relion_euler_rot_mat_verbose(phi,theta,psi)
      
      implicit none
      
c     RELION rotation angles
c     rotation around z axis in lab frame by angle phi
c     0 < phi < 2.0d0*pi
c     the rotation results in new axis x',y',z'
      real *8 phi
c     rotation around the new y' axis by angle theta
c     0 < theta < pi
c     the rotation results in new axis x'',y'',z''
      real *8 theta
c     rotation around the new z'' axis by angle psi
c     0 < psi < 2.0d0*pi
      real *8 psi
      
      real *8 RotMat(3,3)
      real *8 RotMat1(3,3)
      real *8 RotMat2(3,3)
      real *8 RotMat3(3,3)
      real *8 Rpsi(3,3)
      real *8 Rtheta(3,3)
      real *8 Rphi(3,3)

      real *8 TRotMat(3,3)
      real *8 TRotMat1(3,3)
      real *8 TRotMat2(3,3)
      real *8 TRotMat3(3,3)
      real *8 TRpsi(3,3)
      real *8 TRtheta(3,3)
      real *8 TRphi(3,3)

      integer i,j


      Rphi(1,1) = dcos(phi)
      Rphi(2,1) = -dsin(phi)
      Rphi(3,1) = 0.0d0
      Rphi(1,2) = dsin(phi)
      Rphi(2,2) = dcos(phi)
      Rphi(3,2) = 0.0d0
      Rphi(1,3) = 0.0d0
      Rphi(2,3) = 0.0d0
      Rphi(3,3) = 1.0d0

      Rtheta(1,1) = dcos(theta)
      Rtheta(2,1) = 0.0d0
      Rtheta(3,1) = dsin(theta)
      Rtheta(1,2) = 0.0d0
      Rtheta(2,2) = 1.0d0
      Rtheta(3,2) = 0.0d0
      Rtheta(1,3) = -dsin(theta)
      Rtheta(2,3) = 0.0d0
      Rtheta(3,3) = dcos(theta)

      Rpsi(1,1) = dcos(psi)
      Rpsi(2,1) = -dsin(psi)
      Rpsi(3,1) = 0.0d0
      Rpsi(1,2) = dsin(psi)
      Rpsi(2,2) = dcos(psi)
      Rpsi(3,2) = 0.0d0
      Rpsi(1,3) = 0.0d0
      Rpsi(2,3) = 0.0d0
      Rpsi(3,3) = 1.0d0

      RotMat1(1,1) = Rtheta(1,1)*Rphi(1,1)+Rtheta(1,2)*Rphi(2,1)
     1     +Rtheta(1,3)*Rphi(3,1)
      RotMat1(2,1) = Rtheta(2,1)*Rphi(1,1)+Rtheta(2,2)*Rphi(2,1)
     1     +Rtheta(2,3)*Rphi(3,1)
      RotMat1(3,1) = Rtheta(3,1)*Rphi(1,1)+Rtheta(3,2)*Rphi(2,1)
     1     +Rtheta(3,3)*Rphi(3,1)
      RotMat1(1,2) = Rtheta(1,1)*Rphi(1,2)+Rtheta(1,2)*Rphi(2,2)
     1     +Rtheta(1,3)*Rphi(3,2)
      RotMat1(2,2) = Rtheta(2,1)*Rphi(1,2)+Rtheta(2,2)*Rphi(2,2)
     1     +Rtheta(2,3)*Rphi(3,2)
      RotMat1(3,2) = Rtheta(3,1)*Rphi(1,2)+Rtheta(3,2)*Rphi(2,2)
     1     +Rtheta(3,3)*Rphi(3,2)
      RotMat1(1,3) = Rtheta(1,1)*Rphi(1,3)+Rtheta(1,2)*Rphi(2,3)
     1     +Rtheta(1,3)*Rphi(3,3)
      RotMat1(2,3) = Rtheta(2,1)*Rphi(1,3)+Rtheta(2,2)*Rphi(2,3)
     1     +Rtheta(2,3)*Rphi(3,3)
      RotMat1(3,3) = Rtheta(3,1)*Rphi(1,3)+Rtheta(3,2)*Rphi(2,3)
     1     +Rtheta(3,3)*Rphi(3,3)

      RotMat3(1,1) = Rpsi(1,1)*RotMat1(1,1)+Rpsi(1,2)*RotMat1(2,1)
     1     +Rpsi(1,3)*RotMat1(3,1)
      RotMat3(2,1) = Rpsi(2,1)*RotMat1(1,1)+Rpsi(2,2)*RotMat1(2,1)
     1     +Rpsi(2,3)*RotMat1(3,1)
      RotMat3(3,1) = Rpsi(3,1)*RotMat1(1,1)+Rpsi(3,2)*RotMat1(2,1)
     1     +Rpsi(3,3)*RotMat1(3,1)
      RotMat3(1,2) = Rpsi(1,1)*RotMat1(1,2)+Rpsi(1,2)*RotMat1(2,2)
     1     +Rpsi(1,3)*RotMat1(3,2)
      RotMat3(2,2) = Rpsi(2,1)*RotMat1(1,2)+Rpsi(2,2)*RotMat1(2,2)
     1     +Rpsi(2,3)*RotMat1(3,2)
      RotMat3(3,2) = Rpsi(3,1)*RotMat1(1,2)+Rpsi(3,2)*RotMat1(2,2)
     1     +Rpsi(3,3)*RotMat1(3,2)
      RotMat3(1,3) = Rpsi(1,1)*RotMat1(1,3)+Rpsi(1,2)*RotMat1(2,3)
     1     +Rpsi(1,3)*RotMat1(3,3)
      RotMat3(2,3) = Rpsi(2,1)*RotMat1(1,3)+Rpsi(2,2)*RotMat1(2,3)
     1     +Rpsi(2,3)*RotMat1(3,3)
      RotMat3(3,3) = Rpsi(3,1)*RotMat1(1,3)+Rpsi(3,2)*RotMat1(2,3)
     1     +Rpsi(3,3)*RotMat1(3,3)

      RotMat(1,1) = dcos(psi)*dcos(theta)*dcos(phi)-dsin(psi)*dsin(phi)
      RotMat(2,1) = -dsin(psi)*dcos(theta)*dcos(phi)-dcos(psi)*dsin(phi)
      RotMat(3,1) = dsin(theta)*dcos(phi)
      RotMat(1,2) = dcos(psi)*dcos(theta)*dsin(phi)+dsin(psi)*dcos(phi)
      RotMat(2,2) = -dsin(psi)*dcos(theta)*dsin(phi)+dcos(psi)*dcos(phi)
      RotMat(3,2) = dsin(theta)*dsin(phi)
      RotMat(1,3) = -dcos(psi)*dsin(theta)
      RotMat(2,3) = dsin(psi)*dsin(theta)
      RotMat(3,3) = dcos(theta)

c      write(6,*) RotMat(1,1), RotMat3(1,1)
c      write(6,*) RotMat(2,1), RotMat3(2,1)
c      write(6,*) RotMat(3,1), RotMat3(3,1)
c      write(6,*) RotMat(1,2), RotMat3(1,2)
c      write(6,*) RotMat(2,2), RotMat3(2,2)
c      write(6,*) RotMat(3,2), RotMat3(3,2)
c      write(6,*) RotMat(1,3), RotMat3(1,3)
c      write(6,*) RotMat(2,3), RotMat3(2,3)
c      write(6,*) RotMat(3,3), RotMat3(3,3)

c     transpose
      do i=1,3
         do j = 1,3
            TRotMat(i,j) = RotMat(j,i)
         enddo
      enddo
c
      do i=1,3
         do j = 1,3
            TRphi(i,j) = Rphi(j,i)
         enddo
      enddo
c
      do i=1,3
         do j = 1,3
            TRtheta(i,j) = Rtheta(j,i)
         enddo
      enddo
c
      do i=1,3
         do j = 1,3
            TRpsi(i,j) = Rpsi(j,i)
         enddo
      enddo
      
      TRotMat1(1,1) = TRtheta(1,1)*TRpsi(1,1)+TRtheta(1,2)*TRpsi(2,1)
     1     +TRtheta(1,3)*TRpsi(3,1)
      TRotMat1(2,1) = TRtheta(2,1)*TRpsi(1,1)+TRtheta(2,2)*TRpsi(2,1)
     1     +TRtheta(2,3)*TRpsi(3,1)
      TRotMat1(3,1) = TRtheta(3,1)*TRpsi(1,1)+TRtheta(3,2)*TRpsi(2,1)
     1     +TRtheta(3,3)*TRpsi(3,1)
      TRotMat1(1,2) = TRtheta(1,1)*TRpsi(1,2)+TRtheta(1,2)*TRpsi(2,2)
     1     +TRtheta(1,3)*TRpsi(3,2)
      TRotMat1(2,2) = TRtheta(2,1)*TRpsi(1,2)+TRtheta(2,2)*TRpsi(2,2)
     1     +TRtheta(2,3)*TRpsi(3,2)
      TRotMat1(3,2) = TRtheta(3,1)*TRpsi(1,2)+TRtheta(3,2)*TRpsi(2,2)
     1     +TRtheta(3,3)*TRpsi(3,2)
      TRotMat1(1,3) = TRtheta(1,1)*TRpsi(1,3)+TRtheta(1,2)*TRpsi(2,3)
     1     +TRtheta(1,3)*TRpsi(3,3)
      TRotMat1(2,3) = TRtheta(2,1)*TRpsi(1,3)+TRtheta(2,2)*TRpsi(2,3)
     1     +TRtheta(2,3)*TRpsi(3,3)
      TRotMat1(3,3) = TRtheta(3,1)*TRpsi(1,3)+TRtheta(3,2)*TRpsi(2,3)
     1     +TRtheta(3,3)*TRpsi(3,3)

      TRotMat3(1,1) = TRphi(1,1)*TRotMat1(1,1)+TRphi(1,2)*TRotMat1(2,1)
     1     +TRphi(1,3)*TRotMat1(3,1)
      TRotMat3(2,1) = TRphi(2,1)*TRotMat1(1,1)+TRphi(2,2)*TRotMat1(2,1)
     1     +Rphi(2,3)*TRotMat1(3,1)
      TRotMat3(3,1) = TRphi(3,1)*TRotMat1(1,1)+TRphi(3,2)*TRotMat1(2,1)
     1     +TRphi(3,3)*TRotMat1(3,1)
      TRotMat3(1,2) = TRphi(1,1)*TRotMat1(1,2)+TRphi(1,2)*TRotMat1(2,2)
     1     +TRphi(1,3)*TRotMat1(3,2)
      TRotMat3(2,2) = TRphi(2,1)*TRotMat1(1,2)+TRphi(2,2)*TRotMat1(2,2)
     1     +TRphi(2,3)*TRotMat1(3,2)
      TRotMat3(3,2) = TRphi(3,1)*TRotMat1(1,2)+TRphi(3,2)*TRotMat1(2,2)
     1     +Rphi(3,3)*TRotMat1(3,2)
      TRotMat3(1,3) = TRphi(1,1)*TRotMat1(1,3)+TRphi(1,2)*TRotMat1(2,3)
     1     +TRphi(1,3)*TRotMat1(3,3)
      TRotMat3(2,3) = TRphi(2,1)*RotMat1(1,3)+TRphi(2,2)*TRotMat1(2,3)
     1     +TRphi(2,3)*TRotMat1(3,3)
      TRotMat3(3,3) = TRphi(3,1)*TRotMat1(1,3)+TRphi(3,2)*TRotMat1(2,3)
     1     +TRphi(3,3)*TRotMat1(3,3)

      do i = 1,3
         write(6,*) TRotMat(i,1),TRotMat(i,2),TRotMat(i,2)
      enddo


c      write(6,*) TRotMat(1,1), TRotMat3(1,1)
c      write(6,*) TRotMat(2,1), TRotMat3(2,1)
c      write(6,*) TRotMat(3,1), TRotMat3(3,1)
c      write(6,*) TRotMat(1,2), TRotMat3(1,2)
c      write(6,*) TRotMat(2,2), TRotMat3(2,2)
c      write(6,*) TRotMat(3,2), TRotMat3(3,2)
c      write(6,*) TRotMat(1,3), TRotMat3(1,3)
c      write(6,*) TRotMat(2,3), TRotMat3(2,3)
c      write(6,*) TRotMat(3,3), TRotMat3(3,3)



      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c conversions to our rotations
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine rot_mat_to_euler(RotMat,angles)
      implicit none
      
      real *8 RotMat(3,3)
      real *8 angles(3)
      real *8 pi
      
      pi=4.0d0*datan(1.0d0)
      
      angles(1) = RotMat(3,3)
      angles(2) = datan2(RotMat(2,3),RotMat(1,3))
      if(angles(2) < 0) then
         angles(2) = angles(2)+2.0d0*pi
      endif
      angles(3) = datan2(RotMat(3,2),-RotMat(3,1))
      if(angles(3) < 0) then
         angles(3) = angles(3)+2.0d0*pi
      endif

      return
      end


      subroutine euler_to_rot_mat(angles,RotMat)
      implicit none

      real *8 angles(3)
      real *8 RotMat(3,3)

      real *8 R1(3,3),R2(3,3),R3(3,3)
      real *8 RotMat1(3,3)
      real *8 calphai,salphai
      real *8 cbetaj,sbetaj
      real *8 cgamma,sgamma

      calphai=angles(1)
      salphai=dsqrt(1-calphai**2)
      cbetaj=dcos(angles(2))
      sbetaj=dsin(angles(2))
      cgamma=dcos(angles(3))
      sgamma=dsin(angles(3))

      R3(1,1) = cgamma;
      R3(1,2) = -sgamma;
      R3(1,3) = 0.0d0;
      R3(2,1) = sgamma;
      R3(2,2) = cgamma;
      R3(2,3) = 0.0d0;
      R3(3,1) = 0.0d0;
      R3(3,2) = 0.0d0;
      R3(3,3) = 1;

      R1(1,1) = calphai
      R1(1,2) = 0.0d0
      R1(1,3) = salphai
      R1(2,1) = 0.0d0;
      R1(2,2) = 1;
      R1(2,3) = 0.0d0;
      R1(3,1) = -salphai
      R1(3,2) = 0.0d0
      R1(3,3) = calphai         
            
      R2(1,1) = cbetaj;
      R2(1,2) = -sbetaj;
      R2(1,3) = 0.0d0;
      R2(2,1) = sbetaj;
      R2(2,2) = cbetaj;
      R2(2,3) = 0.0d0;
      R2(3,1) = 0.0d0;
      R2(3,2) = 0.0d0;
      R2(3,3) = 1;

      RotMat1(1,1) = calphai*cbetaj
      RotMat1(1,2) = -sbetaj
      RotMat1(1,3) = salphai*cbetaj
      RotMat1(2,1) = calphai*sbetaj
      RotMat1(2,2) = cbetaj
      RotMat1(2,3) = salphai*sbetaj
      RotMat1(3,1) = -salphai
      RotMat1(3,2) = 0.0d0
      RotMat1(3,3) = calphai

      RotMat(1,1) = calphai*cbetaj*cgamma-sbetaj*sgamma
      RotMat(1,2) = -calphai*cbetaj*sgamma-sbetaj*cgamma
      RotMat(1,3) = salphai*cbetaj
      RotMat(2,1) = calphai*sbetaj*cgamma+cbetaj*sgamma
      RotMat(2,2) = -calphai*sbetaj*sgamma+cbetaj*cgamma
      RotMat(2,3) = salphai*sbetaj
      RotMat(3,1) = -salphai*cgamma
      RotMat(3,2) = salphai*sgamma
      RotMat(3,3) = calphai
      
      return
      end
