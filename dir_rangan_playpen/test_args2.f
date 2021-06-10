      program test_args2
c     gfortran -w -o test_args2.out test_args2.f ; ./test_args2.out 523 1 45;
      character(len=100) :: argstring
      integer i1,i2,i3
      call get_command_argument(1,argstring)
      read (argstring, '(I10)') i1
      call get_command_argument(2,argstring)
      read (argstring, '(I10)') i2
      call get_command_argument(3,argstring)
      read (argstring, '(I10)') i3
      write(6,*) i1,i2,i3
      stop
      end
