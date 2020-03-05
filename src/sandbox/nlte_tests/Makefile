appname := a.out

CXX := mpicxx
CXXFLAGS := -std=c++11


# location of gsl
GSL_INC=$(GSL_DIR)/include
GSL_LIB=$(GSL_DIR)/lib

# location of lua
LUA_INC=$(LUA_DIR)/src
LUA_LIB=$(LUA_DIR)/src

# location of hdf5
HDF_INC=$(HDF5_DIR)/include
HDF_LIB=$(HDF5_DIR)/lib

CCOPT = -I$(GSL_INC) -I$(LUA_INC) -I$(HDF_INC)
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
