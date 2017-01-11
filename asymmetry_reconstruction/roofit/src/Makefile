ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -Wall -fPIC -Wno-deprecated

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit -lRooFit -lRooFitCore -lFoam


CXXFLAGS      += $(ROOTCFLAGS)
CXX           += -I./
LIBS           = $(ROOTLIBS) 

GLIBS          = $(filter-out -lNew, $(NGLIBS))

CXX	      += -I./lib/
OUTLIB	      = ./lib/
.SUFFIXES: .cc,.C
.PREFIXES: ./lib/

#----------------------------------------------------#

main: macro.cc lib/RooDDecay.o
	$(CXX) $(CXXFLAGS) -o AppExample $(OUTLIB)/*.o  $(GLIBS) $<

#Main code:
lib/RooDDecay.o:
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)RooDDecay.o $<
lib/RooDDecay.o: src/RooDDecay.cc

clean:
	rm -f src/*~
	rm -f *~
	rm -f $(OUTLIB)*.o
	rm -f AppExample

#----------------------------------------------------#
