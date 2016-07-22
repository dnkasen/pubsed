#####
# SEDONA makefile
#######

# put your make.inc here
MAKE_INC = make.inc.tagore

# executable directory where compiliation takes place
EXECDIR = ./EXEC

all:
	cp -p sedona.cpp $(EXECDIR)
	cp -p transport/*.cpp transport/*.h $(EXECDIR)	
	cp -p hydro/*.cpp hydro/*.h $(EXECDIR)	
	cp -p helper/*.cpp helper/*.h $(EXECDIR)		
	cp -p opacity/*.cpp opacity/*.h $(EXECDIR)	
	cp -p grid/*.cpp grid/*.h $(EXECDIR)	
	cp -p $(MAKE_INC) $(EXECDIR)/make.inc
	cp make.exec $(EXECDIR)/Makefile
	cd $(EXECDIR); make

clean:
	rm -f $(EXECDIR)/*.o $(EXECDIR)/gomc

realclean:
	rm -f $(EXECDIR)/*.o $(EXECDIR)/*.cpp $(EXECDIR)/*.h $(EXECDIR)/gomc $(EXECDIR)/Makefile $(EXECDIR)/make.inc
