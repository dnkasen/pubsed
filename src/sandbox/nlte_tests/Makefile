appname := a.out

CXX := mpicxx
#CXXFLAGS := -std=c++11
CXXFLAGS = -O3 -DNDEBUG -DMPICH_IGNORE_CXX_SEEK -std=c++0x

GSL_DIR = /usr/local/Cellar/gsl/2.6
# location of gsl
GSL_INC=$(GSL_DIR)/include
GSL_LIB=$(GSL_DIR)/lib

LUA_DIR = /usr/local/Cellar/lua/5.3.5_1

# location of lua
LUA_INC=$(LUA_DIR)/include
LUA_LIB=$(LUA_DIR)/lib

# location of hdf5

HDF5_DIR = /usr/local/Cellar/hdf5/1.12.0

HDF_INC=$(HDF5_DIR)/include
HDF_LIB=$(HDF5_DIR)/lib

EIGEN_DIR=/Users/dkhatami/Documents/codes/sedona_nlte/sedona6/src/external/eigen-3.3.7/Eigen/
EIGEN_INC=$(EIGEN_DIR)

CCOPT = -I$(GSL_INC) -I$(LUA_INC) -I$(HDF_INC) -I$(EIGEN_INC)
CLOPT = $(CCOPT) -L$(GSL_LIB) -L$(LUA_LIB) -L$(HDF_LIB) -llua -lgsl -lgslcblas -lhdf5 -lhdf5_hl -ldl -I$(GSL_INC) -I$(LUA_INC) -I$(HDF_INC)

srcfiles := $(shell find . -name "*.cpp")
objects  := $(patsubst %.C, %.o, $(srcfiles))

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(appname) $(objects) $(LDLIBS) $(CCOPT) $(CLOPT)

depend: .depend

.depend: $(srcfiles)
	rm -f ./.depend a.out
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	rm -f *.o

dist-clean: clean
	rm -f *~ .depend

include .depend
