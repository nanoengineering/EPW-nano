program test5
use iotk_module
use iotk_xtox_interf
implicit none
integer :: i,j
character(2) :: c
character(30) :: str

goto 2

open(10,file="pippo.txt")
write(10,"(a)") "<Root>"
write(10,"(a)") "<prova i='3'>"
write(10,"(a)") "</prova i='4'>"
write(10,"(a)") "</Root>"
close(10)
goto 1

call iotk_open_write(10,"pippo1.txt",root="ee~~")
call iotk_close_write(10)

call iotk_open_read(10,"pippo.txt")
write(0,*) "-"
call iotk_scan_begin(10,"prova")
write(0,*) "-"
call iotk_scan_end(10,"prova")
write(0,*) "-"
call iotk_close_read(10)
write(0,*) "-"

1 continue
do j = 1 , 10
write(0,*) j
call iotk_open_write(11,file="pippo.txt")
write(0,*) j,"+"
do i = 1 ,100
call iotk_write_dat(11,"pippo","aa")
end do
write(0,*) j,"+"
call iotk_close_write(11)

write(0,*) j
call iotk_open_read(11,file="pippo.txt")
do i = 1 ,100
call iotk_scan_dat(11,"pippo",c)
end do
call iotk_close_read(11)
end do

2 continue

do i = 1 , 1000000
  str = iotk_itoa(i)
end do
write(0,*) str

write(0,*) "A"//iotk_index(-123)//"B"
write(0,*) "A"//iotk_index((/1,-2,-4524254/))//"B"


end program test5
