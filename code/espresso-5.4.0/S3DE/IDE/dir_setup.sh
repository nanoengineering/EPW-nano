#! /bin/bash 
#
# dir_setup.sh
#
MANUAL="\
USAGE
  dir_setup.sh [-h] [ [-f] <dir>] 

  -h    :   print this manual
  -f    :   force the creation of all defined directories 
  <dir> :   the directory to be set

  The script sets up the required directories and files in <dir>
  according to the \${IDEHOME}/IDE/IDE.conf file."

#
# input check
shopt -s extglob
if [ "$1" = "-h" ] ; then echo "$MANUAL" ; exit 0 ; fi
force=""
if [ "$1" = "-f" ] ; then force=yes ; shift 1 ; fi
if [ "$#" != "1" ] ; then 
    echo " Invalid number of arguments, type dir_setup.sh -h for help" 
    exit 1
fi

directory=$1
if [ ! -d $directory ] ; then 
    echo " Input argument is not a directory, type dir_setup.sh -h for help" 
    exit 1
fi

#
# load DEFAULTS
if [ ! -e ./IDE/IDE.conf ] ; then
   echo "ERROR: Unable to find ./IDE/IDE.conf" ; exit 1
fi
. ./IDE/IDE.conf

list="\
$STDDIR_SOURCE 
$STDDIR_INCLUDE 
$STDDIR_BINARY 
$STDDIR_TMP 
$STDDIR_DOC
"

#
# check if <dir> exists and it is not a file 
for dir in $list
do
   if [ -e $directory/$dir -a ! -d  $directory/$dir ] ; then
      echo "ERROR: $dir exists and it is NOT a directory"
      exit 1
   fi
done

#
# create dirs
for dir in $list
do
    test -d $directory/$dir || mkdir $directory/$dir
    test -e $directory/$dir/.touch || touch $directory/$dir/.touch
done




