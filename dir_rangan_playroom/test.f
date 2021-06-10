      implicit none
      integer k
      k = 3;
      write(6,'(A,I0)') ' -k/2 ' , -k/2
      write(6,'(A,I0)') ' -(k/2) ' , -(k/2)
      write(6,'(A,I0)') ' (-k)/2 ' , (-k)/2
      k = 4;
      write(6,'(A,I0)') ' -k/2 ' , -k/2
      write(6,'(A,I0)') ' -(k/2) ' , -(k/2)
      write(6,'(A,I0)') ' (-k)/2 ' , (-k)/2
      stop
      end
