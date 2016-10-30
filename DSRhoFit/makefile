CFLAGS = $(shell root-config --cflags)
INCLUDE = -Iinclude
LIBS = $(shell root-config --libs)
ROOLIBS = -lFoam -lRooFitCore -lRooFit -lRooStats -lMinuit

DSRhoFit : src/DSRhoFit.cpp src/FitterTrans.cpp src/DSRhoPDF.cpp
	g++ -o bin/Release/DSRhoFit src/DSRhoFit.cpp src/FitterTrans.cpp src/DSRhoPDF.cpp src/DSRhoPDFTIndep.cpp src/FitterTransTIndep.cpp $(INCLUDE) $(CFLAGS) $(LIBS) $(ROOLIBS)

clean : 
	rm -f bin/Release/DSRhoFit

