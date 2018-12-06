#!/bin/bash

USAGE="usage: ./install.sh TARGET, where TARGET indicates that Makefile.TARGET will be executed
                    help: prints usage
                    clean: deletes code, binary, and object files
                    realclean: clean + deletes Makefiles"

if [ $# == 0 ]; then
  echo "Must specify argument"
  echo "$USAGE"
  exit 1
fi

TARGET=$1
EXECDIR=./build
MAKEDIR=./makefiles

if [ $TARGET == "help" ]
then
  echo $USAGE
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
  cd $EXECDIR; make
  cd ..
  cp $EXECDIR/gomc sedona6.ex
fi
