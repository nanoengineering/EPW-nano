program example1
use iotk_module
implicit none

character(iotk_attlenx) :: attr
character(100) :: type
integer :: i,age,isize
logical :: found

! WRITING FILE:

call iotk_open_write(10,"example2.xml")
  call iotk_write_attr (attr,"size",2,first=.true.)
  call iotk_write_begin(10,"Animals",attr=attr)
! attributes are written on string attr BEFORE writing the tag
    call iotk_write_attr (attr,"name","Anna",first=.true.)
    call iotk_write_attr (attr,"age",1)
    call iotk_write_attr (attr,"type","cat")
    call iotk_write_empty(10,"Animal"//iotk_index(1),attr)
    call iotk_write_attr (attr,"name","Peggy",first=.true.)
    call iotk_write_attr (attr,"age",6)
    call iotk_write_attr (attr,"type","dog")
    call iotk_write_empty(10,"Animal"//iotk_index(2),attr)
  call iotk_write_end  (10,"Animals")
call iotk_close_write(10)

! READING FILE:

call iotk_open_read(10,"example2.xml")
  call iotk_scan_begin(10,"Animals",attr=attr)
  ! attributes are read from string attr AFTER reading the tag
  call iotk_scan_attr (attr,"size",isize)
    do i = 1,isize
      call iotk_scan_empty(10,"Animal"//iotk_index(i),attr=attr)
      call iotk_scan_attr (attr,"type",type)
      write(0,*) trim(type) ! this is for debug
      call iotk_scan_attr (attr,"age",age,found=found)
           ! the 'found' keyword can be used for optional attributes
      if(.not.found) write (0,*) "ERRATO!"
      write(0,*) age !this is for debug
      call iotk_scan_attr (attr,"notes",age,found=found)
      if(found) write (0,*) "ERRATO!"
    end do
  call iotk_scan_end(10,"Animals")
call iotk_close_read(10)

end program example1
