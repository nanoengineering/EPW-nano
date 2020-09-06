#! /bin/bash
# Generate info target
# Usage: makeinfo.sh
# makefile is read on stdin
# info target is written on stdout

shopt -s extglob

#
# this check has been removed for portability issues (AIX bash)
#
#[ $# -ne 0 ] && { echo "Incorrect usage" ; exit 1 ; }

INFO=
while read LINE ; do
  if [[ $LINE == \#*([[:space:]])INFO* ]] ; then
    INFO="${LINE##\#*([[:space:]])INFO*([[:space:]])}"
  elif [[ $LINE == +([^=#]):* && "$INFO" ]] ; then
    TARGET="${LINE%%:*}"
    echo -n "	@echo \"  $TARGET "
    count=$(( 30 - ${#TARGET} ))
    for((i=0;i<count;i++)) ; do echo -n " " ; done
    echo "$INFO\""
    INFO=
  else
    INFO=
  fi
done
echo -n "	@echo " 


