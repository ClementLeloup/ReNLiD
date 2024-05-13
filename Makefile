#Makefile

MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build

.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR); fi;
	touch build/.base

vpath %.cpp source
vpath %.o build
vpath .base build

# C compiler:
#CC       = gcc -std=c++17
CC	 = mpic++ -std=c++17

# Additional libraries
LBFLAG = -lstdc++ --lfftw3 -lm

# where to find include files *.h
INCLUDES = -I../include
HEADERFILES = $(wildcard ./include/*.h)

%.o:  %.cpp .base $(HEADERFILES)
	cd $(WRKDIR);$(CC) $(LBFLAG) $(INCLUDES) -c ../$< -o $*.o


SOURCE = field.o scalar.o vector.o lie.o group.o lattice.o

MAIN = ReNLiD.o

all: ReNLiD

ReNLiD: $(MAIN) $(SOURCE)
	$(CC) -o ReNLiD $(addprefix build/,$(notdir $^)) $(LBFLAG)

clean: .base
	rm -rf $(WRKDIR);
