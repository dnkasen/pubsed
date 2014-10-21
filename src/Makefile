#####
# SEDONA makefile
#######

# put your make.inc here
MAKE_INC = make.inc.tagore

EXECDIR = ./exec

all:
	cp -p main.cpp $(EXECDIR)
	cp -p transport/*.cpp transport/*.h $(EXECDIR)	
	cp -p hydro/*.cpp hydro/*.h $(EXECDIR)	
	cp -p helper/*.cpp helper/*.h $(EXECDIR)		
	cp -p opacity/*.cpp opacity/*.h $(EXECDIR)	
	cp -p grid/*.cpp grid/*.h $(EXECDIR)	
	cp -p $(MAKE_INC) $(EXECDIR)/make.inc
	cp make.exec $(EXECDIR)/Makefile
	cd $(EXECDIR); make

clean:
	rm -f $(EXECDIR)/*.o gomc

