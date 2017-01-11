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

MYFUNC = /Users/rmaria/work/programs/roofit/learning_again/test22.06.15/mypdf/lib/myDecay.o

#----------------------------------------------------#

main: macro.cc lib/Example.o
	$(CXX) $(CXXFLAGS) -o AppExample $(OUTLIB)/*.o  $(MYFUNC) $(GLIBS) $<
#	$(CXX) $(CXXFLAGS) -o AppExample $(OUTLIB)/*.o  $(GLIBS) $<

#Main code:
lib/Example.o:
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)Example.o $<
lib/Example.o: src/Example.cc

clean:
	rm -f src/*~
	rm -f *~
	rm -f $(OUTLIB)*.o
	rm -f AppExample

#----------------------------------------------------#
