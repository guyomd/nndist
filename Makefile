CXX = g++
LIBS = strutils.o arrayutils.o fileutils.o nnanalysis.o nnstats.o
INCDIR = ./include
BOOSTDIR = /home/b94678/Code/Cpp/boost_1_82_0 
SRCDIR = ./src
BINDIR = ./bin

OPTS = -O3 -fopenmp -I$(INCDIR) -I$(BOOSTDIR)

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

nnstats.o: $(SRCDIR)/nnstats.cc 
	$(CXX) $(OPTS) -c $(SRCDIR)/nnstats.cc 

nndistance.o: $(SRCDIR)/nndistance.cc $(SRCDIR)/strutils.cc $(SRCDIR)/arrayutils.cc $(SRCDIR)/fileutils.cc $(SRCDIR)/nnanalysis.cc $(SRCDIR)/nnstats.cc
	$(CXX) $(OPTS) -c $(SRCDIR)/nndistance.cc $(SRCDIR)/strutils.cc $(SRCDIR)/arrayutils.cc $(SRCDIR)/fileutils.cc $(SRCDIR)/nnanalysis.cc $(SRCDIR)/nnstats.cc

nndeclust.o: $(SRCDIR)/nndeclust.cc $(SRCDIR)/strutils.cc $(SRCDIR)/arrayutils.cc $(SRCDIR)/fileutils.cc $(SRCDIR)/nnanalysis.cc $(SRCDIR)/nnstats.cc
	$(CXX) $(OPTS) -c $(SRCDIR)/nndeclust.cc $(SRCDIR)/strutils.cc $(SRCDIR)/arrayutils.cc $(SRCDIR)/fileutils.cc $(SRCDIR)/nnanalysis.cc $(SRCDIR)/nnstats.cc


clean:
	rm strutils.o arrayutils.o fileutils.o nnanalysis.o nnstats.o nndistance.o nndeclust.o
