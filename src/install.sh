#!/bin/bash

USAGE=""

if [ $# == 0 ]; then
  echo "Must specify argument"
  echo "Type ./install.sh help for more information"
  exit 1
fi

TARGET=$1
EXECDIR=./build
MAKEDIR=./makefiles

if [ $TARGET == "help" ]
then
  printf '\nTo compile type: ./install.sh TARGET\n\nWhere TARGET indicates that Makefile.TARGET that will be used \n'
  printf 'Available machines with a TARGET (in the makefiles/ directory) are:\n    '
  for filename in makefiles/Makefile.*; do
    printf ${filename#"makefiles/Makefile."},
    printf ' '
  done
  printf '\n\nOther usage:\n\n ./install.sh help: prints this help \n ./install.sh clean: deletes code, binary, and object files\n ./install.sh realclean: clean + deletes Makefiles\n'
  exit 0
fi
if [ $TARGET == "clean" ]
then
	rm -f $EXECDIR/*.o $EXECDIR/gomc $EXECDIR/*.cpp $EXECDIR/*.h
elif [ $TARGET == "realclean" ]
then
	rm -rf $EXECDIR
else
  mkdir -p $EXECDIR
	cp -p transport/*.cpp transport/*.h $EXECDIR
	cp -p hydro/*.cpp hydro/*.h $EXECDIR
  cp -p hydro/*/*.cpp hydro/*/*.h $EXECDIR
	cp -p main/*.cpp main/*.h $EXECDIR
  cp -p utils/*.cpp utils/*.h $EXECDIR
	cp -p opacity/*.cpp opacity/*.h $EXECDIR
	cp -p grid/*.cpp grid/*.h $EXECDIR
  cp -p $MAKEDIR/make.exec $EXECDIR/make.exec
  cp -p $MAKEDIR/Makefile.$TARGET $EXECDIR/Makefile
  cd $EXECDIR;

  default:
    if make  ; then \
      printf '\n .-._.-._.-._.-._.-. \n'
      printf '(                   ) \n'
      printf ' )    IT WORKED    (  \n'
      printf '(                   ) \n'
      printf ' ._.-._.-._.-._.-._. \n\n'
      printf 'The executable ./sedona6.ex should now be in your src/ directory\n\n' ; \
  else printf "\n *** Compilation Failed *** \n\n" ; fi

  cd ..
  cp $EXECDIR/gomc sedona6.ex



fi
