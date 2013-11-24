# This makefile is tailored to the ipnp30 custom ROOT setup. For other computers use 'make -f makefile.default'.

CFLAGS = -pthread -m64 -I/home/cervenkov/root/include
INCLUDE = -Iinclude
LIBS = -L/home/cervenkov/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic
ROOLIBS = -lFoam -lRooFitCore -lRooFit -lRooStats -lMinuit -lMinuit2# -lMathMore

DSRhoFit : src/DSRhoFit.cpp src/FitterTrans.cpp src/DSRhoPDF.cpp
	g++ -o bin/Release/DSRhoFitCartesian src/DSRhoFit.cpp src/FitterTrans.cpp src/DSRhoPDF.cpp src/DSRhoPDFTIndep.cpp src/FitterTransTIndep.cpp -Wl,-rpath -Wl,$(HOME)/root/lib $(INCLUDE) $(CFLAGS) $(LIBS) $(ROOLIBS)

clean : 
	rm -f bin/Release/DSRhoFitCartesian

