#! /bin/bash 

MANUAL="\
 USAGE
    packls.sh [-h]
 returns the list of available packages l
"
shopt -s extglob

while getopts :h OPT
do
   case $OPT in
   (h)  echo "$MANUAL" ; exit 0 ;;
   (:)  echo "error: $OPTARG requires an argument" ; exit 1 ;;
   (?)  echo "error: unkwown option $OPTARG" ; exit 1 ;;
   esac
done

#
# load DEFAULTS
if [ ! -e ./IDE/IDE.conf ] ; then
   echo "ERROR: Unable to find ./IDE/IDE.conf" ; exit 1
fi
. ./IDE/IDE.conf


#
# load OPTIONS
options=
test -r OPTIONS && options="$(cat OPTIONS)"

#
# load dir list and skip selected dir
skip_dirs="$(echo "$options" | IDE/getoption -c SKIP )
IDE
CVS
CVSROOT
CONFIG"
for lib in $NAMERULE
do
  found=
  for lib2 in $skip_dirs ; do
    [[ $lib = $lib2 ]] && { found=yes ; break ; }
  done
  test -z "$found" &&
  test -d $lib &&
  [[ $lib != @(IDE|CONFIG|tmp) ]] &&
  IDE_LIBS="$IDE_LIBS$lib
"
done

echo "$IDE_LIBS"
exit 0

