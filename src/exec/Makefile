#####
# SEDONA makefile
#######

# put your make.inc here
include make.inc.tagore

TARGET  =  gomc

SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

CCOPT = -I$(GSL_INC) -I$(LUA_INC) -I$(HDF_INC)
CLOPT = -L$(GSL_LIB) -L$(LUA_LIB) -L$(HDF_LIB) -llua -lgsl -lgslcblas -lhdf5 -lhdf5_hl


$(TARGET): $(OBJECTS)
	$(CXX) -o $(TARGET) $(OBJECTS) $(CLOPT)

.cpp.o: 
	$(CXX) $(CCOPT) -I ./exec -c -o $@ $<

clean:
	rm -f $(OBJECTS)

realclean: clean
	rm -f $(TARGET)
