#! /bin/bash

source IDE/IDE.conf

shopt -s extglob

directory=$1

tarlist="$(find $directory -type f -o -type l)"
newtarlist=
#cleanlist=$(cd $directory/src ; ls $SUFFIXES_CLEAN 2>/dev/null)
for file in $tarlist ; do
  found=
  for expression in $SUFFIXES_CLEAN ; do
    [[ $file == $directory/src/$expression ]] && found=yes && break
  done
  [ $found ] || { [[ $file == $directory/tmp/[^.]* ]] && found=yes ; }
  [ $found ] || { [[ $file == $directory/bin/*.x ]] && found=yes ; }
  [ $found ] || newtarlist="$newtarlist $file"
done
tarlist="$newtarlist"

echo "$tarlist"

