CFLAGS = $(shell root-config --cflags)
INCLUDE = -Iinclude
LIBS = $(shell root-config --libs)
ROOLIBS = -lFoam -lRooFitCore -lRooFit -lRooStats -lMinuit -lMathMore

DSRhoFit : src/DSRhoFit.cpp
	g++ $(LIBS) $(ROOLIBS) -o DSRhoFit src/DSRhoFit.cpp $(INCLUDE) $(CFLAGS)


clean : 
	rm -f DSRhoFit

