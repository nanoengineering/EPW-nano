#! /bin/sh 
# moduledep.sh -- script that computes dependencies on Fortran 90 modules
# 13.08.04     -- also dipendencies on #include files are serched

# files whose dependencies must be computed
touch __tmp__.f90
sources=`ls *.f90 `
if [ "$sources" = ""  ] ; then
   exit 0
fi

# files that may contain modules
# extra directories can be specified on the command line
touch __tmp__.h
sources_all="$sources `ls *.h`"
for dir in $*
do
  touch $dir/__tmp__.f90
  touch $dir/__tmp__.h
  sources_all="$sources_all `ls $dir/*.f90 $dir/*.h `"
done

# remove previous contents
rm -f  moduledep.tmp1 moduledep.tmp2 


# create list of module dependencies
# each line is of the form:
# file_name.o : @module_name@
# cast all module names to lowercase because Fortran is case insensitive
egrep -i "^ *use " $sources |             # look for "USE name"
sed 's/f90:/o /
     s/,/ /' |                            # replace extension, insert space
#                                         #   and remove trailing comma
awk '{print $1 " : @" tolower($3) "@"}' | # create dependency entry
sort | uniq > moduledep.tmp1              # remove duplicates


# create alist of dependencies on .h files
# the form is as before :
# filename.o : @include_file.h@
# case is mantained for include file names
egrep -i "^# *include" $sources |         # look for #include statement
egrep -v "<.*>"  |                        # avoid lines like "#include <pippo.h>"
sed 's/f90:/o /
     s/,/ /
     s/\"/ /g
     s/#/ /
     s/include/ /' | tr -d "\'"   |
#                                         # replace extension, insert space
#                                         #   remove trailing comma
#                                         #   remove double and single quotes 
#                                         #   remove '# include' statements
#                                         
awk '{print $1 " : @" $2 "@"}' |          # create dependency entry
sort | uniq >> moduledep.tmp1             # remove duplicates


# create list of available modules
# for each module, create a line of the form:
# s/@module_name@/file_name/
egrep -i "^ *module " $sources_all |           # look for "MODULE name"
sed 's/f90:/o /
     s/\//\\\//g' |                            # replace extension, insert
#                                              #   space and escape slashes
awk '{print "s/@" tolower($3) "@/" $1 "/" }' | # create substitution line
sort | uniq > moduledep.tmp2                   # remove duplicates

# The line is the following for include files:
# s/@include_file.h@/complete_file_name/

rm -rf moduledep.pathnames moduledep.names
LIST=`
echo $sources_all          |                       
tr " " "\n"                |                   # everithing on a different line
grep -v "\.f90"            |                   # delete all f90 files
grep -v "__tmp__.h"`                           # delete __tmp__.h from the LIST
echo "$LIST"               |                   # insert space and escape slashes
sed 's/\//\\\//g'  >  moduledep.pathnames

echo "$LIST" | 
tr "\/" "\n" |                                 # all the / become new-line
grep "\.h"  > moduledep.names                  # get the simple names of the include files

paste moduledep.names moduledep.pathnames |    # merging the names and the names+path
awk '{print "s/@" $1 "@/" $2 "/" }'       |        
sort | uniq >> moduledep.tmp2

rm -f moduledep.names moduledep.pathnames

# replace module names with file names
# by applying the file of substitution patterns just created
sed -f moduledep.tmp2 moduledep.tmp1 |
awk '{if ($1 != $3) print}' |          # remove self dependencies
sort | uniq                            # remove duplicates


# remove false module names arosen from 
# MODULE PROCEDURE interfaces

egrep -v "@procedure@" moduledep.tmp2 > __tmp__
mv __tmp__ moduledep.tmp2


# removing __tmp__.f90 __tmp__.h from all dirs
# including the current dir
test -e __tmp__.f90 && rm __tmp__.f90
test -e __tmp__.h   && rm __tmp__.h
for dir in $*
do
  test -e $dir/__tmp__.f90 && rm $dir/__tmp__.f90
  test -e $dir/__tmp__.h   && rm $dir/__tmp__.h
done


rm -f moduledep.tmp1 moduledep.tmp2 # remove temporary files

