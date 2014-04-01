CFLAGS = $(shell root-config --cflags)
LIBS = $(shell root-config --libs)
ROOLIBS = -lFoam -lRooFitCore -lRooFit -lRooStats -lMinuit -lMinuit2

all : overlay sum

overlay : overlay.cpp
	g++ -o overlay overlay.cpp $(CFLAGS) $(LIBS) $(ROOLIBS)

sum : sum.cpp
	g++ -o sum sum.cpp $(CFLAGS) $(LIBS) $(ROOLIBS)


