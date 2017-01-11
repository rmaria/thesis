ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -Wall -fPIC -Wno-deprecated

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit -lRooFit -lRooFitCore -lFoam -lTMVA

CXXFLAGS      += $(ROOTCFLAGS)
CXX           += -I./
LIBS           = $(ROOTLIBS) 

GLIBS          = $(filter-out -lNew, $(NGLIBS))

CXX	      += -I./lib/
OUTLIB	      = ./lib/
.SUFFIXES: .cc,.C
.PREFIXES: ./lib/

#----------------------------------------------------#

main: lib/myDecay.o  

lib/myDecay.o:
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)myDecay.o $<
lib/myDecay.o: src/myDecay.cc

clean:
	rm -f src/*~
	rm -f *~
	rm -f $(OUTLIB)*.o

#----------------------------------------------------#
