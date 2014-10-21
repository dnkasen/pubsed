#####
# SEDONA makefile
#######

# put your make.inc here
MAKE_INC = make.inc.tagore

EXECDIR = ./exec

all:
	cp -p transport/*.cpp $(EXECDIR)	
	cp -p hydro/*.cpp $(EXECDIR)	
	cp -p helper/*.cpp $(EXECDIR)		
	cp -p opacity/*.cpp $(EXECDIR)	
	cp -p grid/*.cpp $(EXECDIR)	
	cp -p $(MAKE_INC) $(EXECDIR)/make.inc
	cd $(EXECDIR); make

clean:
	rm -f $(EXECDIR)/*.o gomc

