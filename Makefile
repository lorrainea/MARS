MF=     Makefile
 
CC=     g++
 
CFLAGS= -g -fopenmp -D_USE_OMP -msse4.2 -O3 -fomit-frame-pointer -funroll-loops  
 
LFLAGS= -std=c++11 -I ./ -I ./libsdsl/include/ -L ./libsdsl/lib/ -lsdsl -ldivsufsort -ldivsufsort64 -Wl,-rpath=$(PWD)/libsdsl/lib
 
EXE=    mars
 
SRC=    mars.cc utils.cc sacsc.cc ced.cc nj.cc progAlignment.cc
 
HD=     mars.h sacsc.h ced.h nj.h Makefile
 
# 
# No need to edit below this line 
# 
 
.SUFFIXES: 
.SUFFIXES: .cc .o 
 
OBJ=    $(SRC:.cc=.o) 
 
.cc.o: 
	$(CC) $(CFLAGS)-c $(LFLAGS) $< 
 
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
