program test4
use iotk_base
use iotk_error_interf


integer :: ierr,ierr1,ierr2,ierr3,ierr4
ierr = 0
ierr1= 0
ierr2= 0
ierr3= 0
ierr3= 4


write(0,*) "#########################################"
call iotk_error_issue(ierr,"main",__FILE__,__LINE__)
call iotk_error_msg  (ierr,"Error in IOTK")
call iotk_error_issue(ierr1,"main",__FILE__,__LINE__)
call iotk_error_msg  (ierr1,"IO error")
call iotk_error_write(ierr1,"iostat",30)
call iotk_error_issue(ierr2,"main",__FILE__,__LINE__)
call iotk_error_msg  (ierr2,"Other error")
call iotk_error_write(ierr2,"iostat",32)
call iotk_error_issue(ierr3,"main",__FILE__,__LINE__)
call iotk_error_msg  (ierr3,"Ultimo")
call iotk_error_write(ierr3,"iostat",34)
call iotk_error_issue(ierr4,"main",__FILE__,__LINE__)
call iotk_error_msg  (ierr4,"Ultimissimo")

write(0,*) "===",ierr,iotk_error_pool_order(ierr)
call iotk_error_print(ierr,0)
write(0,*) "===",ierr1,iotk_error_pool_order(ierr1)
call iotk_error_print(ierr1,0)
write(0,*) "===",ierr2,iotk_error_pool_order(ierr2)
call iotk_error_print(ierr2,0)
write(0,*) "===",ierr3,iotk_error_pool_order(ierr3)
call iotk_error_print(ierr3,0)
write(0,*) "===",ierr4,iotk_error_pool_order(ierr4)
call iotk_error_print(ierr4,0)

write(0,*) iotk_error_pool_pending(),ierr,ierr2

call iotk_error_clear(ierr)

write(0,*) iotk_error_pool_pending()

call iotk_error_print(16,0)
call iotk_error_print(17,0)
write(0,*) iotk_error_check(0)
write(0,*) iotk_error_check(1)
write(0,*) iotk_error_check(2)
write(0,*) iotk_error_check(3)
write(0,*) iotk_error_check(16)
write(0,*) iotk_error_check(17)

call iotk_error_handler(ierr2)

end program test4
