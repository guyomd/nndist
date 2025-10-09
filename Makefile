CXX = g++
LIBS = strutils.o arrayutils.o fileutils.o nnanalysis.o
INCDIR = ./include
SRCDIR = ./src
BINDIR = ./bin

OPTS = -O3 -fopenmp -I$(INCDIR)

all: nndistance nndecluster clean
	@echo "\nDone. NNDISTANCE and NNDECLUSTER binaries stored in subdirectory $(BINDIR)\n"

nndistance: $(LIBS) nndistance.o
	$(CXX) $(OPTS) $^ -o $(BINDIR)/$@ 

nndecluster: $(LIBS) nndeclust.o
	$(CXX) $(OPTS) $^ -o $(BINDIR)/$@ 

strutils.o: $(SRCDIR)/strutils.cc
	$(CXX) $(OPTS) -c $(SRCDIR)/strutils.cc

arrayutils.o: $(SRCDIR)/arrayutils.cc 
	$(CXX) $(OPTS) -c $(SRCDIR)/arrayutils.cc

fileutils.o: $(SRCDIR)/fileutils.cc 
	$(CXX) $(OPTS) -c $(SRCDIR)/fileutils.cc

nnanalysis.o: $(SRCDIR)/nnanalysis.cc 
	$(CXX) $(OPTS) -c $(SRCDIR)/nnanalysis.cc 

nndistance.o: $(SRCDIR)/nndistance.cc $(SRCDIR)/strutils.cc $(SRCDIR)/arrayutils.cc $(SRCDIR)/fileutils.cc $(SRCDIR)/nnanalysis.cc
	$(CXX) $(OPTS) -c $(SRCDIR)/nndistance.cc $(SRCDIR)/strutils.cc $(SRCDIR)/arrayutils.cc $(SRCDIR)/fileutils.cc $(SRCDIR)/nnanalysis.cc

nndeclust.o: $(SRCDIR)/nndeclust.cc $(SRCDIR)/strutils.cc $(SRCDIR)/arrayutils.cc $(SRCDIR)/fileutils.cc $(SRCDIR)/nnanalysis.cc 
	$(CXX) $(OPTS) -c $(SRCDIR)/nndeclust.cc $(SRCDIR)/strutils.cc $(SRCDIR)/arrayutils.cc $(SRCDIR)/fileutils.cc $(SRCDIR)/nnanalysis.cc


clean:
	rm strutils.o arrayutils.o fileutils.o nnanalysis.o nndistance.o nndeclust.o
