#! /bin/bash

source IDE/IDE.conf

shopt -s extglob

directory=$1
DIRECTORY=$(echo $directory | tr "[a-z]" "[A-Z]")
version="$(IDE/getoption -c VERSION < $directory/README)"
if [ -r $directory/include/${directory}_version.sh ] ; then
  . $directory/include/${directory}_version.sh
  eval written_version=\$__${DIRECTORY}_VERSION
  test $version = $written_version && exit
fi
test -z "$version" && exit
[[ "$version" != $VERSRULE ]] && exit

version_major="${version%%.*}"
version_minor="${version#+([[:alnum:]]).}"
version_minor="${version_minor%%.*}"
version_patch="${version#+([[:alnum:]]).+([[:alnum:]]).}"
version_patch="${version_patch%%[[:alpha:]]*}"
version_extra="${version#+([[:alnum:]]).+([[:alnum:]]).+([[:alnum:]])}"

echo > $directory/include/${directory}_version.h "\
#define __${DIRECTORY}_VERSION \"$version\"
#define __${DIRECTORY}_VERSION_MAJOR $version_major
#define __${DIRECTORY}_VERSION_MINOR $version_minor
#define __${DIRECTORY}_VERSION_PATCH $version_patch
#define __${DIRECTORY}_VERSION_EXTRA \"$version_extra\""

echo > $directory/include/${directory}_version.sh "\
__${DIRECTORY}_VERSION=$version
__${DIRECTORY}_VERSION_MAJOR=$version_major
__${DIRECTORY}_VERSION_MINOR=$version_minor
__${DIRECTORY}_VERSION_PATCH=$version_patch
__${DIRECTORY}_VERSION_EXTRA=$version_extra"


