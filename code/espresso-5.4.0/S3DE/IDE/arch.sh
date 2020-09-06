#! /bin/bash 

NEXT=
SCRIPT=
FILES=
AFTER=
for OPT
do
  case "$NEXT" in
  (--after-extraction|-a) AFTER="$OPT" ; NEXT= ;;
  (--script|-s) SCRIPT="$OPT" ; NEXT= ;;
  (*)
    case "$OPT" in
    (--help) ;;
    (--after-extraction|-a) NEXT="$OPT" ;;
    (--script|-s) NEXT="$OPT" ;;
    (-*) echo "$0: Unknown option $OPT" ; exit 1 ;;
    (*) FILES="$FILES"$'\n'"$OPT" ;;
    esac
  ;;
  esac
done

test -z "$SCRIPT" && { echo "$0: specify the script file with -s" ; exit 1 ; }

#test -e "$SCRIPT" && rm "$SCRIPT"
#touch "$SCRIPT"
test -e "$SCRIPT" || touch "$SCRIPT"

IFS=$'\n'
for FILE in $FILES ; do
  test -r "$FILE" || { echo "$0: file not found: $FILE" ; continue ; }
  IFS=/
  ACCDIR=.
  for DIR in ${FILE%\/*} ; do
    ACCDIR="$ACCDIR/$DIR"
    echo "test -e \"$ACCDIR\" || mkdir \"$ACCDIR\"" >> "$SCRIPT"
  done
  IFS=" "$'\n'
  EOF="$RANDOM$RANDOM$RANDOM$RANDOM$RANDOM$RANDOM"
  echo "cat << \"$EOF\" > \"$FILE\"" >> "$SCRIPT"
  cat "$FILE" >> "$SCRIPT"
  echo "$EOF" >> "$SCRIPT"
done

if [ -n "$AFTER" ] ; then
if [ -r "$AFTER" ] ; then
  cat "$AFTER" >> "$SCRIPT"
else
  echo "$0: file not found $AFTER"
fi
fi

# for simple direct use
chmod u+x $SCRIPT

