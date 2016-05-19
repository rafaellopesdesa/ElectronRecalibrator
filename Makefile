#
# stuff to make
#
LD=g++
SOURCES=$(wildcard *.cxx)
OBJECTS=$(SOURCES:.cxx=.o)
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
LIB=recalibrator.so
SOFLAGS=-O -fPIC -shared
COREINC=-I./

$(LIB): $(OBJECTS) 
	$(LD) $(SOFLAGS) $(OBJECTS) $(ROOTLIBS) -lRooFit -lRooFitCore -lTMVA -lEG -lGenVector -o $@

%.o:	%.cxx
	$(CXX) $(COREINC) $(ROOTCFLAGS) -O2 -Wall -fPIC -c $< -o $@

#
# target to build
# likelihood id library
#

all: $(LIB) 
clean:
	rm -f *.o \
	rm -f *.d \
	rm -f *.so \
