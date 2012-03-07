CFLAGS = $(shell root-config --cflags)
INCLUDE = -Iinclude
LIBS = $(shell root-config --libs)
ROOLIBS = -lFoam -lRooFitCore -lRooFit -lRooStats -lMinuit# -lMathMore

DSRhoFit : src/DSRhoDataMining.cpp src/Particle.cpp
	g++ $(LIBS) $(ROOLIBS) -o bin/Release/DSRhoDataMining src/DSRhoDataMining.cpp src/Particle.cpp $(INCLUDE) $(CFLAGS)

clean : 
	rm -f bin/Release/DSRhoDataMining

