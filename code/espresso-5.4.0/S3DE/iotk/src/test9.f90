program test9
use iotk_unit_list_module
use iotk_module
implicit none

type(iotk_unit_list) :: unit_list
type(iotk_unit), pointer :: ptr

call iotk_unit_list_init(unit_list)

call iotk_unit_list_search(unit_list,ptr,unit=10)
write(0,*) "should be false",associated(ptr)

call iotk_unit_list_add(unit_list,ptr)
ptr%unit=10
ptr%root="aa"

call iotk_unit_list_add(unit_list,ptr)
ptr%unit=11
ptr%root="aa"

nullify(ptr)

call iotk_unit_list_search(unit_list,ptr,unit=10)

write(0,*) ptr%unit
write(0,*) ptr%root

call iotk_unit_list_del(unit_list,ptr)

call iotk_unit_list_destroy(unit_list)





end program test9
