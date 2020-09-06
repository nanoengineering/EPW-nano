#! /bin/bash 

MANUAL="\
 USAGE
    cvssnap.sh [-h] [ -n <file> ] [ -d dir...dirN]
      the script produces a second script file able to reproduce the CVS dir 
      structure and and automatically updating all the files with the current 
      revision is produced. 
      If no dir is given ./ is assumed.  
"
shopt -s extglob
IDEHOME=$(pwd)

filename=
list=
while getopts :hn:d: OPT
do
   case $OPT in
   (h)  echo "$MANUAL" ; exit 0 ;;
   (n)  filename=$OPTARG ;;
   (d)  list=$OPTARG ;;
   (:)  echo "error: $OPTARG requires an argument" ; exit 1 ;;
   (?)  echo "error: unkwown option $OPTARG" ; exit 1 ;;
   esac
done

if [ -z "$filename" ] ; then
    echo "ERROR: -n <file> should be specified"
    exit 1
fi

if [ -n "$filename" ] ; then 
   test -e $filename && rm $filename
   touch $filename
   chmod ug+x $filename
fi

# by default, void list is set to the current dir
[ -z "$list" ] && list="./" 

#
# if the script is generated here write the header 
#
HEADER="\
#!/bin/sh 
#
# Today is $(date)
# The script has been automatically generated.
#
MANUAL=\"\\
  USAGE:
    $filename [-h] [-u cvsuser]  

  The present script extract a snapshot of a project or part 
  of it from the CVS repository. The -u option overwrite the name
  of the cvs user set in the snapshot. \"

while getopts :hu: OPT
do
   case \$OPT in
   (h)  echo \"\$MANUAL\" ; exit 0 ;;
   (u)  newuser=\$OPTARG ;;
   (:)  echo \"error: \$OPTARG requires an argument\" ; exit 1 ;;
   (?)  echo \"error: unkwown option \$OPTARG\" ; exit 1 ;;
   esac
done

#
# here start the generation of the CVS tree
CVSROOT_LOCAL=\$CVSROOT
"

echo "$HEADER" > $filename

#
# storing the CVS directories
#
for dir in $list
do
    cvslist=$(find $dir -name CVS)   
    for cvsdir in $cvslist
    do
        dirname=${cvsdir#"./"}
        dirname=${dirname%CVS}
        #
        # the list of file in the format
        #   revision  filename
        #
        filelist="$filelist $( awk -v DIR=$dirname  -F /   \
                '{ if ( $1 != "D" )  print DIR"@"$2"@"$3 }' $cvsdir/Entries )"
        #
        # reproducing the CVS dirs
        #
        $IDEHOME/IDE/arch.sh --script $filename $cvsdir/*

        #
        # changing the CVS/Root files if needed
        #
        TMP="\

# Change the CVSROOT environment in CVS/Root file if needed
if [ -n \"\$newuser\" ] ; then
      rootstr=\`cat $cvsdir/Root \`
      new_rootstr=\`echo \$rootstr | sed "s/:[a-z]*@/:\$newuser@/" \`
      rm $cvsdir/Root 
      echo \${new_rootstr} > $cvsdir/Root
      CVSROOT_LOCAL=\$new_rootstr
fi
        
"
        echo "$TMP" >> $filename

        #
        # adding the sticky tags to the Entries file
        #
        TMP="\

# Add the sticky tags to the CVS/Entries file
    tmp=\"\$( awk -F / ' { if ( \$1 != \"D\" ) { sticky=\"T\"\$3 } else { sticky=\" \" } }
                   { print \$0sticky }' $cvsdir/Entries )\"
    rm $cvsdir/Entries
    echo \"\$tmp\" > $cvsdir/Entries
    
"
        echo "$TMP" >> $filename

    done
done


#
# updating the files
#
TMP="\

#Retrieving files with the correct revisions
#
# NOTA BENE: when the snapshotted dirs have different root it probably does not work...
#
cvs -d \$CVSROOT_LOCAL login

for dir in $list
do
    cd \$dir
    cvs update
    cd -
done

cvs -d \$CVSROOT_LOCAL logout 
exit 0

"
echo "$TMP" >> $filename

exit 0

