CFLAGS = $(shell root-config --cflags)
INCLUDE = -Iinclude
LIBS = $(shell root-config --libs)
ROOLIBS = -lFoam -lRooFitCore -lRooFit -lRooStats -lMinuit# -lMathMore

DSRhoFit : src/DSRhoGraphs.cpp
	g++ -o bin/Release/DSRhoGraphs src/DSRhoGraphs.cpp src/ObservablesCollectionCartesian.cpp $(INCLUDE) $(CFLAGS) $(LIBS) $(ROOLIBS)

clean : 
	rm -f bin/Release/DSRhoGraphs

