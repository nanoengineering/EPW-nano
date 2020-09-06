program test
use iotk_module
implicit none

integer :: i
call iotk_set(error_warn_overflow=.true.)

do i = 1 , 1

write(0,*) "Test for IOTK library"

call test1(.false.,"test1.xml")
write(0,*) "Textual OK"
call test1(.true., "test1.dat")
write(0,*) "Binary OK"

call test1(.false.,"test1l.xml",link=.true.)
write(0,*) "Textual with links OK"
call test1(.true., "test1l.dat",link=.true.)
write(0,*) "Binary with links OK"


call test1(.false.,"test1dl.xml",link=.true.,deferred=.true.)
write(0,*) "Textual with deferred links OK"
call test1(.true., "test1dl.dat",link=.true.,deferred=.true.)
write(0,*) "Binary with deferred links OK"

call test1(.false.,"test1lr.xml",link=.true.,raw=.true.)
write(0,*) "Textual with raw links OK"
call test1(.true., "test1lr.dat",link=.true.,raw=.true.)
write(0,*) "Binary with raw links OK"

call test1(.false.,"test1dlr.xml",link=.true.,raw=.true.,deferred=.true.)
write(0,*) "Textual with raw deferred links OK"
call test1(.true., "test1dlr.dat",link=.true.,raw=.true.,deferred=.true.)
write(0,*) "Binary with raw deferred links OK"


end do


contains

subroutine test1(binary,file,link,deferred,raw)
  logical, intent(in) :: binary
  character(*), intent(in) :: file
  logical, optional, intent(in) :: link
  logical, optional, intent(in) :: deferred
  logical, optional, intent(in) :: raw
  logical :: lraw,ldeferred,llink
!  integer, parameter :: dp=selected_real_kind(10),i8=selected_int_kind(10)
  integer, parameter :: dp=selected_real_kind(5),i8=selected_int_kind(5)
  character(30) :: char_array(3)
  character(25) :: char_array_read(3)
  logical       :: logical_array(30),logical_array_read(30)
  integer       :: integer_array(30)
  integer(i8)   :: integer_array_read(30)
  real(dp)      :: real_array(30)
  real          :: real_array_read(30)
  complex       :: complex_array(10,2)
  complex(dp)   :: complex_array_read(10,2)
  character(iotk_attlenx) :: attr
  character(100):: escapes,escapes_read
  character(100):: escapes1,escapes1_read
  character(100):: escapes2,escapes2_read
  character(100):: value

  logical :: logical_att(3),logical_att_read(3)
  integer :: integer_att(3),integer_att_read(3)
  real    :: real_att(2,2),real_att_read(2,2)
  complex :: complex_att(3),complex_att_read(3)

  logical :: stream

  llink=.false.
  ldeferred=.false.
  lraw=.false.
  if(present(raw)) lraw=raw
  if(present(deferred)) ldeferred=deferred
  if(present(link)) llink=link
  
  char_array(1) = "First && <<>line"
  char_array(2) = "Second line'''" // '"' // " ss"
  char_array(3) = "Third line"
  logical_array = .true.
  integer_array = 7
  real_array    = real( 3.14 , kind=kind(real_array) )
  complex_array = (2.0,1.0)
  escapes = "try&try<try>try'"//'"'
  escapes1 = "'ert'"
  escapes2 = '"ert"'
  logical_att=(/.true.,.false.,.true./)
  integer_att=(/1,2,-443/)
  real_att=reshape((/1.0,2.0,5.0,7.0/),(/2,2/))
  complex_att=(-1.0,-2.0)
  

  call iotk_open_write(unit=10,file=file,binary=binary)
  call iotk_write_begin(10,"Outer")
  if(ldeferred) then
    call iotk_link(10,name="Inner1",file="Inner1.xml",binary=.false.)
  else
    call iotk_write_attr(attr,"escapes",trim(escapes),first=.true.)
    call iotk_write_attr(attr,"escapes1",trim(escapes1))
    call iotk_write_attr(attr,"escapes2",trim(escapes2))
    if(llink) call iotk_link(10,"Inner1",file="Inner1.xml",binary=.false.,create=.true.)
    call iotk_write_begin(10,"Inner1",attr)
    if(llink) call iotk_link(10,"Dat1",file="Dat1.xml",create=.true.,raw=lraw)
    call iotk_write_dat  (10,"Dat1",char_array)
    if(llink) call iotk_link(10,"Dat2",file="Dat2.xml",create=.true.,raw=lraw)
    call iotk_write_dat  (10,"Dat2",logical_array,columns=4,sep="^^")
    if(llink) call iotk_link(10,"Dat3",file="Dat3.xml",create=.true.,raw=lraw)
    call iotk_write_dat  (10,"Dat3",integer_array)
    if(llink) call iotk_link(10,"Dat4",file="Dat4.xml",create=.true.,raw=lraw)
    call iotk_write_dat  (10,"Dat4",real_array)
    if(llink) call iotk_link(10,"Dat5",file="Dat5.xml",create=.true.,raw=lraw)
    call iotk_write_dat  (10,"Dat5",complex_array,columns=3)
    call iotk_write_end  (10,"Inner1")
  end if
  attr=" "
  call iotk_write_attr(attr,"log",logical_att)
  call iotk_write_attr(attr,"int",integer_att)
  call iotk_write_attr(attr,"rea",real_att)
  call iotk_write_attr(attr,"com",complex_att)
  call iotk_write_begin(10,"Inner2",attr)
  if(llink) call iotk_link(10,"Dat5",file="Dat5.dat",binary=.true.,create=.true.)
  call iotk_write_dat  (10,"Dat5",char_array)
  if(llink) call iotk_link(10,"Dat4",file="Dat4.dat",binary=.true.,create=.true.)
  call iotk_write_dat  (10,"Dat4",logical_array,fmt="(10l1)")
  if(ldeferred) then
    call iotk_link(10,name="Dat3",file="Dat3.dat",binary=.true.,raw=lraw)
  else
    if(llink) call iotk_link(10,"Dat3",file="Dat3.dat",binary=.true.,create=.true.,raw=lraw)
    call iotk_write_dat  (10,"Dat3",integer_array,fmt="(3i8)")
  end if
  if(llink) call iotk_link(10,"Dat2",file="Dat2.dat",binary=.true.,create=.true.)
  call iotk_write_dat  (10,"Dat2",real_array,fmt="(f25.14)")
  if(llink) call iotk_link(10,"Dat1",file="Dat1.dat",binary=.true.,create=.true.)
  call iotk_write_attr (attr,"attribute","value",first=.true.)
  call iotk_write_dat  (10,"Dat1",complex_array,fmt="(f25.14)",attr=attr)
  call iotk_write_end  (10,"Inner2")
  call iotk_write_end  (10,"Outer")
  call iotk_close_write(10)


  stream=.false.
#ifdef __IOTK_STREAM
  stream=.true.
#endif

  call iotk_open_read(unit=10,file=file,stream=stream)
  call iotk_scan_begin(10,"Outer")
  call iotk_scan_begin(10,"Inner2",attr)
  logical_att_read = .false.
  integer_att_read = 0
  real_att_read = 0.0
  complex_att_read = 0.0
  call iotk_scan_attr(attr,"log",logical_att_read)
  call iotk_scan_attr(attr,"int",integer_att_read)
  call iotk_scan_attr(attr,"rea",real_att_read)
  call iotk_scan_attr(attr,"com",complex_att_read)
  if(any(logical_att_read.neqv.logical_att)) stop "logical_att"
  if(any(integer_att_read/=integer_att)) stop "integer_att"
  if(any(abs(real_att_read-real_att)>0.00001)) stop "real_att"
  if(any(complex_att_read/=complex_att)) stop "complex_att"
  char_array_read = " "
  logical_array_read = .false.
  integer_array_read = -1
  real_array_read = 0.0
  complex_array_read = 0.0
  attr = " "
  call iotk_scan_dat  (10,"Dat1",complex_array_read,attr=attr)
  call iotk_scan_attr (attr,"attribute",value)
  if(value/="value") stop "value"
  call iotk_scan_dat  (10,"Dat2",real_array_read)
  call iotk_scan_dat  (10,"Dat3",integer_array_read)
  call iotk_scan_dat  (10,"Dat4",logical_array_read)
  call iotk_scan_dat  (10,"Dat5",char_array_read)
  if(any(char_array_read/=char_array)) stop "char"
  if(any(logical_array_read.neqv.logical_array)) stop "logical"
  if(any(integer_array_read/=integer_array)) stop "integer"
  if(any(abs(real_array_read-real_array)>0.00001)) stop "real"
  if(any(complex_array_read/=complex_array)) stop "complex"
  call iotk_scan_end  (10,"Inner2")
  attr=" "
  escapes_read=""
  escapes1_read=""
  escapes2_read=""
  call iotk_scan_begin(10,"Inner1",attr)
  call iotk_scan_attr(attr,"escapes",escapes_read)
  call iotk_scan_attr(attr,"escapes1",escapes1_read)
  call iotk_scan_attr(attr,"escapes2",escapes2_read)
  if(escapes/=escapes_read) stop "escape"
  if(escapes1/=escapes1_read) stop "escape1"
  if(escapes2/=escapes2_read) stop "escape2"
  char_array_read = " "
  logical_array_read = .false.
  integer_array_read = -1
  real_array_read = 0.0
  complex_array_read = 0.0
  call iotk_scan_dat  (10,"Dat5",complex_array_read)
  call iotk_scan_dat  (10,"Dat4",real_array_read)
  call iotk_scan_dat  (10,"Dat3",integer_array_read)
  call iotk_scan_dat  (10,"Dat2",logical_array_read)
  call iotk_scan_dat  (10,"Dat1",char_array_read)
  if(any(char_array_read/=char_array)) stop "char1"
  if(any(logical_array_read.neqv.logical_array)) stop "logical1"
  if(any(integer_array_read/=integer_array)) stop "integer1"
  if(any(abs(real_array_read-real_array)>0.00001)) stop "real1"
  if(any(complex_array_read/=complex_array)) stop "complex1"
  call iotk_scan_end  (10,"Inner1")
  call iotk_scan_end  (10,"Outer")
  call iotk_close_read(10)
end subroutine test1

end program test


