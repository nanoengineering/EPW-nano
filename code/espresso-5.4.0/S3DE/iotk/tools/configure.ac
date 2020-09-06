AC_INIT(IOTK, 0.1)

# LIST OF USED ENVIRONMENT VARIABLES:
# IOTK_INTEGER_KINDS
# IOTK_LOGICAL_KINDS
# IOTK_REAL_KINDS
# colon separated lists of kinds
# a kind ALL expands to hopefully all kinds

# Default integer and logical are added by default
IOTK_INTEGER_KINDS="kind(1):$IOTK_INTEGER_KINDS"
IOTK_LOGICAL_KINDS="kind(.true.):$IOTK_LOGICAL_KINDS"

# Save the default IFS
IFS_BACKUP="$IFS"

# Define the newline character
NEWLINE="
"

# Define the backslash character

AC_PROG_FC(ifort g95 f95 pgf90 xlf90)

# Selecting the language
AC_LANG(Fortran)

# Selecting the source extension
AC_FC_SRCEXT(f90)

AC_MSG_NOTICE([checking intrinsic FORTRAN kinds])

# Loop over types
for TYPE in INTEGER LOGICAL REAL
do
  eval LIST=\"\$IOTK_${TYPE}_KINDS\"

# Expansion of the 'ALL' keyword
  NLIST=""
  IFS=":"
  for KIND in $LIST
  do
    case "$KIND" in
    ALL) NLIST="$NLIST:1:2:3:4:5:6:7:8:10:12:16:24:32" ;;
    *)   NLIST="$NLIST:$KIND" ;;
    esac
  done
  IFS="$IFS_BACKUP"
  LIST="$NLIST"
  
  NLIST=""
  IFS=":"
  for KIND in $LIST
  do
  IFS="$IFS_BACKUP"
  
  if test -n "$KIND"
  then
    AC_MSG_CHECKING([intrinsic type $TYPE(kind=$KIND)])
    AC_RUN_IFELSE(
    [[
    program pippo
      implicit none
      $TYPE (kind=$KIND) :: try
      open(10,file="test.txt")
      if(kind(try)==$KIND) write(10,"(a,i8)") "KIND=",kind(try)
      close(10)
    end program pippo
    ]],
    DONE=yes,
    DONE=""
    )

    if test -n "$DONE"
    then
      if grep -q KIND= < test.txt
      then
        AC_MSG_RESULT([yes])
        NLIST="$NLIST:"`sed -n 's/KIND= *//p' < test.txt`
      else
        AC_MSG_RESULT([no])
      fi
    else
      AC_MSG_RESULT([no])
    fi
    if test -f test.txt ; then rm test.txt ; fi
  fi

  done
  IFS="$IFS_BACKUP"
  LIST="$NLIST"

  AC_MSG_NOTICE([Eliminating duplicates])

  NLIST=""
  IFS=":"
  for KIND in $LIST
  do
    if test -n "$KIND"
    then
      WAS_FOUND=""
      for KIND1 in $NLIST
      do
        if test -n "$KIND1"
        then
          test "$KIND1" -eq "$KIND" && WAS_FOUND=yes
        fi
      done
      test -n "$WAS_FOUND" || NLIST="$NLIST:$KIND"
    fi
  done
  IFS="$IFS_BACKUP"
  LIST="$NLIST"

  IFS=":"
  COUNT=0
  for KIND in $LIST
  do
    if test -n "$KIND"
    then
      COUNT=`expr $COUNT + 1`
      eval "IOTK_$TYPE$COUNT='#define __IOTK_$TYPE$COUNT $KIND'"
    fi
  done
  IFS="$IFS_BACKUP"
done

##############
# KNOWN BUGS #
##############

IOTK_BUGS=""
AC_MSG_NOTICE([checking for known bugs])

# BUG 1
AC_MSG_CHECKING([bug in non-advancing read])
AC_RUN_IFELSE(
[[
program bug1
implicit none
character(10) :: line
open(file="bug1.txt",unit=10)
open(file="result.txt",unit=11)
write(10,"(a)") "A"
write(10,"(a)") "B"
write(10,"(a)") "C"
close(10)
open(file="bug1.txt",unit=10)
read(10,*)
read(10,"(a)",advance="no",eor=1) line
write(11,*) "FAILED" ! This line should not be reached
stop
1 continue
if(line /= "B") then
  write(11,*) "FAILED" ! This line should not be reached
  stop
end if
backspace(10)
read(10,"(a)",advance="no",eor=2) line
write(11,*) "FAILED" ! This line should not be reached
stop
2 continue
if(line /= "B") then
  write(11,*) "FAILED" ! This line should not be reached
  stop
end if
close(10)
write(11,*) "PASSED"
close(11)
end program bug1
]],DONE=yes,DONE="")

if test -n "$DONE"
then
  if grep -q PASSED < result.txt
  then
    AC_MSG_RESULT([no])
  else
    AC_MSG_RESULT([yes])
    IOTK_WORKAROUND1='#define __IOTK_WORKAROUND1'
  fi
  rm bug1.txt result.txt
else
  AC_MSG_RESULT([yes])
  IOTK_WORKAROUND1='#define __IOTK_WORKAROUND1'
fi

# BUG 2
AC_MSG_CHECKING([bug in converting between logical(2) and logical(8) (if they exist)])

AC_COMPILE_IFELSE(
[[
program main
implicit none
logical(2) :: l2
logical(8) :: l8
l8 = .true.
l2 = l8
write(*,*) l2
end program
]],DONE="yes",DONE="")

if test -n "$DONE"
then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([yes])
  IOTK_WORKAROUND2='#define __IOTK_WORKAROUND2'
fi

# BUG 3
# to do

AC_MSG_CHECKING([bug in pack])
AC_RUN_IFELSE(
[[
program bug3
implicit none
integer, parameter :: size1 = 10000
integer, parameter :: size2 = 1000
integer :: out1(size1*size2)
integer :: out2(size1*size2)
integer :: in(size1,size2)
in=1
call mypack(out1,in,size(in))
out2=pack(in,mask=.true.)
open(10,file="test.txt")
if(all(out1==out2)) write(10,*) "PASSED"
close(10)
end program bug3

subroutine mypack(out,in,n)
implicit none
integer, intent(in)  :: n
integer, intent(in)  :: in(n)
integer, intent(out) :: out(n)
  out = in
end subroutine mypack
]],DONE="yes",DONE="")

if test -n "$DONE"
then
  if grep -q PASSED test.txt
  then
    AC_MSG_RESULT([no])
  else
    AC_MSG_RESULT([yes])
    IOTK_WORKAROUND3='#define __IOTK_WORKAROUND3'
  fi
else
  AC_MSG_RESULT([yes])
  IOTK_WORKAROUND3='#define __IOTK_WORKAROUND3'
fi
if test -f test.txt ; then rm test.txt ; fi


# BUG 4
# to do

AC_MSG_CHECKING([other bug in pack])
AC_COMPILE_IFELSE(
[[
program bug4
implicit none
call sub((/"a","b","c"/))
write(*,*) "PASSED"
contains

subroutine sub(str)
character(len=*), intent(in) :: str(:)
write(*,*) pack(str,mask=.true.)
end subroutine sub

end program bug4
]],DONE="yes",DONE="")

if test -n "$DONE"
then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([yes])
  IOTK_WORKAROUND4='#define __IOTK_WORKAROUND4'
fi


# BUG 5
AC_MSG_CHECKING([huge(1_1) does not compile as a parameter])

AC_COMPILE_IFELSE(
[[
program bug5
integer, parameter :: i=huge(1_1)
write(*,*) "PASSED",i
end program bug5
]],DONE="yes",DONE="")

if test -n "$DONE"
then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([yes])
  IOTK_WORKAROUND5='#define __IOTK_WORKAROUND5'
fi

# BUG 6
AC_MSG_CHECKING([bug in intent])
AC_COMPILE_IFELSE(
[[
subroutine pippo(arg)
implicit none
character(len=*), intent(out) :: arg
character(len=len(arg)), allocatable :: cc(:)
arg="ss"
end subroutine pippo
]],DONE="yes",DONE="")
if test -n "$DONE"
then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([yes])
  IOTK_WORKAROUND6='#define __IOTK_WORKAROUND6'
fi

AC_SUBST(IOTK_INTEGER)
AC_SUBST(IOTK_LOGICAL)
AC_SUBST(IOTK_REAL)
AC_SUBST(IOTK_BUGS)

AC_SUBST(IOTK_INTEGER1)
AC_SUBST(IOTK_INTEGER2)
AC_SUBST(IOTK_INTEGER3)
AC_SUBST(IOTK_INTEGER4)
AC_SUBST(IOTK_LOGICAL1)
AC_SUBST(IOTK_LOGICAL2)
AC_SUBST(IOTK_LOGICAL3)
AC_SUBST(IOTK_LOGICAL4)
AC_SUBST(IOTK_REAL1)
AC_SUBST(IOTK_REAL2)
AC_SUBST(IOTK_REAL3)
AC_SUBST(IOTK_REAL4)
AC_SUBST(IOTK_WORKAROUND1)
AC_SUBST(IOTK_WORKAROUND2)
AC_SUBST(IOTK_WORKAROUND3)
AC_SUBST(IOTK_WORKAROUND4)
AC_SUBST(IOTK_WORKAROUND5)
AC_SUBST(IOTK_WORKAROUND6)
AC_SUBST(IOTK_WORKAROUND7)
AC_SUBST(IOTK_WORKAROUND8)
AC_SUBST(IOTK_WORKAROUND9)

for outfile in `find . -name iotk_config.h.in`
do
  file=`echo $outfile | sed 's/.in$//'`
  AC_CONFIG_FILES($file)
done


AC_OUTPUT

