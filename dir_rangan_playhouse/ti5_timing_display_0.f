      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 
