program example3
use iotk_module
implicit none

real(KIND=4) :: dat1(10)
real(KIND=8) :: dat2(10) 

dat1 = real( (/0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0/) , kind=kind(dat1) )

! WRITING FILE: (a binary one)

call iotk_open_write(10,"example3.dat",binary=.true.)
  call iotk_write_dat(10, "dat1", dat1)
call iotk_close_write(10)

! READING FILE: 
! there is no need to say the library that the file is binary
call iotk_open_read(10,"example3.dat")
  ! kind conversion is transparent
  call iotk_scan_dat (10,"dat1", dat2)
  write(0,*) dat2
call iotk_close_read(10)

write(0,*) kind(dat1),kind(dat2)

end program example3
