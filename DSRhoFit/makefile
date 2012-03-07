CFLAGS = $(shell root-config --cflags)
INCLUDE = -Iinclude
LIBS = $(shell root-config --libs)
ROOLIBS = -lFoam -lRooFitCore -lRooFit -lRooStats -lMinuit# -lMathMore

DSRhoFit : src/DSRhoFit.cpp src/FitterTrans.cpp src/DSRhoPDF.cpp
	g++ $(LIBS) $(ROOLIBS) -o bin/Release/DSRhoFit src/DSRhoFit.cpp src/FitterTrans.cpp src/DSRhoPDF.cpp $(INCLUDE) $(CFLAGS)

clean : 
	rm -f bin/Release/DSRhoFit

