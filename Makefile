#Makefile

MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build

.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR); fi;
	touch build/.base

vpath %.cpp source:test
vpath %.o build
vpath .base build

# C compiler:
#CC       = gcc -std=c++17
CC	 = mpic++ -std=c++17

# Additional libraries
LBFLAG = -lstdc++ -lm -lfftw3 -lmpi_cxx -lboost_mpi -lboost_serialization -Ofast #-lmpi

# Optimization flags
OPTIM_FLAG = #-O3

# where to find include files *.h
INCLUDES = -I../include 
HEADERFILES = $(wildcard ./include/*.h)

%.o:  %.cpp .base $(HEADERFILES)
	cd $(WRKDIR);$(CC) $(LBFLAG) $(INCLUDES) -c ../$< -o $*.o


SOURCE = field.o scalar.o vector.o lie.o group.o lattice.o timeslice.o evolver.o
SOURCE_FIELD = field.o scalar.o vector.o lie.o group.o

TEST_FIELD = test_field.o
TEST_TIMESLICE = test_timeslice.o
TEST_EVOLVER = test_evolver.o

MAIN = ReNLiD.o

all: ReNLiD

ReNLiD: $(MAIN) $(SOURCE)
	$(CC) -o ReNLiD $(addprefix build/,$(notdir $^)) $(LBFLAG)

test_field: $(SOURCE_FIELD) $(TEST_FIELD)
	$(CC) -o $@ $(addprefix build/,$(notdir $^)) $(LBFLAG)

test_timeslice: $(SOURCE) $(TEST_TIMESLICE)
	$(CC) -o $@ $(addprefix build/,$(notdir $^)) $(LBFLAG)

test_evolver: $(SOURCE) $(TEST_EVOLVER)
	$(CC) -o $@ $(addprefix build/,$(notdir $^)) $(LBFLAG) $(OPTIM_FLAG)

clean: .base
	rm -rf $(WRKDIR);
