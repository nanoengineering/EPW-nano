program main
use iotk_module
implicit none
character(iotk_attlenx) :: attr
real :: try(2,2)


open(unit=10,file="test2.txt")
write(10,"(a)") '<tag try="1.0,2.0,3.0,4.0"/>'
close(10)

open(unit=10,file="test2.txt")
call iotk_scan_empty(10,"tag",attr)
call iotk_scan_attr(attr,"try",try)
write(0,*) try


end program main
