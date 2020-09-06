
program example4
use iotk_module
implicit none

type person_type
  integer        :: age
  complex           :: height
  complex, pointer  :: scores(:)
end type person_type

type(person_type) :: person
type(person_type) :: person_read

call person_init(person,5)
person%age    = 26
person%height = 1.88
person%scores = (/5.0,5.0,6.0,6.0,4.0/)

call iotk_open_write(10,file="example4.xml",binary=.false.)
  call person_write(person,"Giovanni",10)
call iotk_close_write(10)

call person_finalize(person)

call person_init(person_read,0)
call iotk_open_read(10,file="example4.xml",binary=.false.)
  call person_scan(person_read,"Giovanni ",10)
call iotk_close_read(10)

! just for debug:
write(*,*) person_read%age,person_read%height
write(*,*) person_read%scores

call person_finalize(person_read)

contains

subroutine person_init(person,nscores)
  type (person_type), intent(out) :: person
  integer,            intent(in)  :: nscores
  person%age=0
  person%height=0.0
  allocate(person%scores(nscores))
end subroutine person_init

subroutine person_finalize(person)
  type (person_type), intent(inout) :: person
  deallocate(person%scores)
end subroutine person_finalize

subroutine person_write(person,name,unit)
  type (person_type), intent(in) :: person
  character(len=*),   intent(in) :: name
  integer,            intent(in) :: unit
  character(iotk_attlenx) :: attr
  call iotk_write_attr(attr,"type","person",first=.true.)
    call iotk_write_begin(unit,name,attr=attr)
    call iotk_write_attr(attr,"age",person%age,first=.true.)
    call iotk_write_attr(attr,"height",person%height)
    call iotk_write_attr(attr,"nscores",size(person%scores))
    call iotk_write_empty(unit,"vals",attr=attr)
    call iotk_write_dat  (unit,"scores",person%scores)
  call iotk_write_end  (unit,name)
end subroutine person_write

subroutine person_scan(person,name,unit)
  type (person_type), intent(out) :: person
  character(len=*),   intent(in)  :: name
  integer,            intent(in)  :: unit
  character(iotk_attlenx) :: attr
  character(iotk_vallenx) :: rtype
  integer :: nscores
  call iotk_scan_begin(unit,name,attr=attr)
    call iotk_scan_attr(attr,"type",rtype)
    if(rtype/="person") stop
    call iotk_scan_empty(unit,"vals",attr=attr)
    call iotk_scan_attr(attr,"age",person%age)
    call iotk_scan_attr(attr,"height",person%height)
    call iotk_scan_attr(attr,"nscores",nscores)
    deallocate(person%scores)
    allocate(person%scores(nscores))
    call iotk_scan_dat (unit,"scores",person%scores)
  call iotk_scan_end  (unit,name)
end subroutine person_scan

end program example4
