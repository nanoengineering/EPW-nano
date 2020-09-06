#! /bin/bash
# Generates outer makefile form inner makefile
# Usage: makeouter.sh PATH
# where PATH is the relative path of the inner makefile
# with respect to the outer makefile
# inner makefile is read on stdin
# outer makefile is written on stdout

shopt -s extglob

[ $# -ne 1 ] && { echo "Incorrect usage" ; exit 1 ; }

while read LINE ; do
  if [[ $LINE == +([^=#\"]):* ]] ; then
    TARGET="${LINE%%*([[:space:]]):*}"
    [[ $TARGET == *.o ]] && continue
    echo "$TARGET:"
    echo "	cd $1 ; make $TARGET"
    echo
  fi
done


