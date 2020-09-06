program test6
use iotk_module
implicit none
character(32000) :: line
integer :: ierr,l
integer :: index,n(3)
logical :: ll(3)
real    :: rr(3)
complex :: cc(3)
character(20) :: words(3),words2(3),stringa,stringa2
character(iotk_attlenx) :: attr

call iotk_copyfile(source="uno",dest="due")

stop

words(1) = "Primo"
words(2) = "Secondo"
words(3) = "Contorno"
stringa  = "    WER  R  "

call iotk_open_write(10,"test6b.txt")
call iotk_write_attr(attr,"prova",stringa,first=.true.)
call iotk_write_empty(10,"empty",attr)
call iotk_write_dat(10,"menu",words)
call iotk_close_write(10)

call iotk_open_read(10,"test6.txt")
call iotk_scan_dat(10,"menu",words2)
call iotk_scan_empty(10,"empty",attr)
call iotk_scan_attr(attr,"prova",stringa2)
call iotk_close_read(10)

write(0,"(a)") words(1)//"%"
write(0,"(a)") words(2)//"%"
write(0,"(a)") words(3)//"%"
write(0,"(a)") "#"//stringa//"#"
write(0,"(a)") words2(1)//"%"
write(0,"(a)") words2(2)//"%"
write(0,"(a)") words2(3)//"%"
write(0,"(a)") "#"//stringa2//"#"

stop

n=255
index = 0
call iotk_read(n,"10000  1, ,   ",index,ierr)
write(0,*) n
write(0,*) index,ierr

ll = .false.
index = 0
call iotk_read(ll,"             .true. false ",index,ierr)
write(0,*) ll
write(0,*) index,ierr

rr = 0.0
index = 0
call iotk_read(rr,"    ,,,,        -3.01   4.2222E3, 1.2",index,ierr)
write(0,*) rr
write(0,*) index,ierr
call iotk_error_print(ierr,0)

cc = 0.0
index = 0
call iotk_read(cc,"1.0 2.0 3.0 4.0 5.0 ",index,ierr)
write(0,*) cc
write(0,*) index,ierr
call iotk_error_print(ierr,0)

stop

open(10,file="test6.in")
open(11,file="test6.out")

#if 0
do
  call iotk_getline(10,line,length=l)
  write(11,"(i5,a)") l,line(1:l)
end do
#else

l=1
call iotk_getline(10,line,length=l)
write(11,"(a)") line(1:l)

call iotk_getline(10,line,length=l)
write(11,"(a)") line(1:l)

backspace(10)
backspace(10)

call iotk_getline(10,line,length=l)
write(11,"(a)") line(1:l)

call iotk_getline(10,line,length=l)
write(11,"(a)") line(1:l)
#endif

close(10)
close(11)

end program test6
