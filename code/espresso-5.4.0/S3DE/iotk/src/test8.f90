program test8
use iotk_module
use iotk_stream_interf
implicit none
integer, allocatable :: val(:,:)
integer, allocatable :: val1(:,:)
integer :: i,rate,ttt1,ttt2,j
logical :: stream

integer :: rsize
integer :: rnum

rsize = 2500000
rnum  = 20

allocate(val(rsize,rnum))
allocate(val1(rsize,rnum))

call system_clock(count_rate=rate)

call system_clock(count=ttt1)

do i=1,rnum
  val(:,i)=(/ (i,j,j=1,rsize) /)
end do
val=77

!call write()

call system_clock(count=ttt2)
write(0,*) "WRITE:",(ttt2-ttt1)/real(rate)
ttt1=ttt2

!call read1()

call system_clock(count=ttt2)
write(0,*) "READ1:",(ttt2-ttt1)/real(rate)
ttt1=ttt2

call read2()

call system_clock(count=ttt2)
write(0,*) "READ2:",(ttt2-ttt1)/real(rate)
ttt1=ttt2

contains

subroutine write()

call iotk_open_write(10,file="aaa.dat",binary=.true.)
do i=1,rnum
  call iotk_write_dat(10,"vv"//iotk_index(i),val(:,i))
end do
call iotk_close_write(10)

end subroutine write

subroutine read1()

val1=0.0
call iotk_open_read(10,file="aaa.dat",binary=.true.)
do i=1,rnum
  call iotk_scan_dat(10,"vv"//iotk_index(i),val1(:,i))
end do
call iotk_close_read(10)
if(any(val1/=val)) stop "EE"

end subroutine read1

subroutine read2()

val1=0.0
stream=.false.
#ifdef __IOTK_STREAMS
stream=.true.
#endif

call iotk_open_read(10,file="aaa.dat",binary=.true.,stream=.true.)
do i=rnum,1,-1
  call iotk_scan_dat(10,"vv"//iotk_index(i),val1(:,i))
end do
call iotk_close_read(10)
if(any(val1/=val)) stop "AA"

end subroutine read2



end program test8
