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
EXECDIR=./EXEC

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
	rm -f $EXECDIR/*.o $EXECDIR/*.cpp $EXECDIR/*.h $EXECDIR/gomc $EXECDIR/Makefile $EXECDIR/make.inc $EXECDIR/make.exec
else
	cp -p sedona.cpp  $EXECDIR
	cp -p transport/*.cpp transport/*.h $EXECDIR	
	cp -p hydro/*.cpp hydro/*.h $EXECDIR	
	cp -p helper/*.cpp helper/*.h $EXECDIR		
	cp -p opacity/*.cpp opacity/*.h $EXECDIR	
	cp -p grid/*.cpp grid/*.h $EXECDIR	
  cp -p make.exec $EXECDIR/make.exec
  cp -p Makefile.$TARGET $EXECDIR/Makefile
	cd $EXECDIR; make
  cd ..
  cp $EXECDIR/gomc sedona6
fi
