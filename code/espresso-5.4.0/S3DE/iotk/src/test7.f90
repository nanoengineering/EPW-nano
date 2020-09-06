PROGRAM test
USE iotk_module
USE iotk_misc_interf
USE iotk_str_interf
IMPLICIT NONE
character(iotk_attlenx) :: attr
integer :: ierr

call iotk_write_attr(attr,"ciccio",1,first=.true.)
call iotk_write_attr(attr,"pippo","sei")
call iotk_write_attr(attr,"aa","nove")
call iotk_write_attr(attr,"WWWW","nove")

attr='cc=  "s"   ee  = "qqqq"  arg ="ee"   '//iotk_eos
write(*,"(a)") "%"//attr(1:iotk_strlen(attr))//"%"

call iotk_delete_attr(attr,"cc",ierr)
if(ierr/=0) write (*,*) "ARG"
write(*,"(a)") "%"//attr(1:iotk_strlen(attr))//"%"

!
! check newline for attributes
!
call iotk_write_attr(attr,"ciccio",1,first=.true.)
call iotk_write_attr(attr,"pippo","sei",newline=.true.)
call iotk_write_attr(attr,"aa","nove",newline=.true.)
call iotk_write_attr(attr,"WWWW","nove")
write(*,"(a)") "%"//attr(1:iotk_strlen(attr))//"%"

END PROGRAM


