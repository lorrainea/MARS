MF=     Makefile
 
CC=     g++
 
CFLAGS= -g -D_USE_OMP -O3 -fomit-frame-pointer -funroll-loops -pthread

# Detect OS for platform-specific settings
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    # macOS
    RPATH_FLAG := -Wl,-rpath,$(PWD)/libsdsl/lib
    # Dynamically get Homebrew prefix (libomp path)
    BREW_PREFIX := $(shell brew --prefix libomp)
    # OpenMP flags for macOS
    OPENMP_FLAGS := -Xpreprocessor -fopenmp -I$(BREW_PREFIX)/include
    OPENMP_LIBS := -L$(BREW_PREFIX)/lib -lomp
    CFLAGS += $(OPENMP_FLAGS)
    LFLAGS += $(OPENMP_LIBS)
else
    # Linux and others
    RPATH_FLAG := -Wl,-rpath=$(PWD)/libsdsl/lib
    CFLAGS += -fopenmp
    LFLAGS += -fopenmp
endif

 
LFLAGS= -std=c++11 -I ./ -I ./libsdsl/include/ -I ./seqan/include/ -L ./libsdsl/lib/ -lsdsl -ldivsufsort -ldivsufsort64 $(RPATH_FLAG) $(OPENMP_LIBS)
 
EXE=    mars
 
SRC=    mars.cc matrices.cc utils.cc sacsc.cc ced.cc nj.cc progAlignment.cc cyclic.cc RestrictedLevenshtein.cc bb.cc heap.cc edlib.cc
 
HD=     EBLOSUM62.h EDNAFULL.h mars.h sacsc.h ced.h nj.h RestrictedLevenshtein.h heap.h Makefile

# 
# No need to edit below this line 
# 
 
.SUFFIXES: 
.SUFFIXES: .cc .o 
 
OBJ=    $(SRC:.cc=.o) 
 
.cc.o: 
	$(CC) $(CFLAGS) -c $(LFLAGS) $< 
 
all:    $(EXE) 
 
$(EXE): $(OBJ) 
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS) 
 
$(OBJ): $(MF) $(HD) 
 
clean: 
	rm -f $(OBJ) $(EXE) *~

clean-all: 
	rm -f $(OBJ) $(EXE) *~
	rm -r libsdsl
	rm -r sdsl-lite
