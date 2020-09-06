program example1
use iotk_module
implicit none

logical :: tmp_logical
real    :: tmp_real4(4)

! WRITING A FILE:

call iotk_open_write(unit=10,file="example1.xml")
  call iotk_write_comment(10,"This is an example file for iotk library")
  call iotk_write_begin(10,"First_set")
    write(iotk_phys_unit(10),*) 0.0,1.0,2.0,3.0
! In this example Fortran IO is used for data.
! However, keep in mind that it is preferred to use the library
! for that, see later examples
! Note also that you should use the iotk_phys_unit function to have
! the physical unit for Fortran IO.
  call iotk_write_end  (10,"First_set")
  call iotk_write_begin(10,"Second_set_more_complex")
    write(iotk_phys_unit(10),*) 1.0,2.0,3.0,4.0
    call iotk_write_begin(10,"First_set")
      write(iotk_phys_unit(10),*) .true.
    call iotk_write_end  (10,"First_set")
  call iotk_write_end  (10,"Second_set_more_complex")
call iotk_close_write(10)

! READING FILE BACK:

! Note that order is changed
! Note also that the library works right even if tags in
! different points of the hierarchy have the same name

call iotk_open_read(10,"example1.xml")
  call iotk_scan_begin(10,"Second_set_more_complex")
    call iotk_scan_begin(10,"First_set")
      read(iotk_phys_unit(10),*) tmp_logical
    call iotk_scan_end  (10,"First_set")
  call iotk_scan_end  (10,"Second_set_more_complex")
  call iotk_scan_begin(10,"First_set")
    read(iotk_phys_unit(10),*) tmp_real4
  call iotk_scan_end  (10,"First_set")
call iotk_close_read(10)

! Just for debug:
write(0,*) tmp_logical
write(0,*) tmp_real4

end program example1
