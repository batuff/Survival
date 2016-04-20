# Makefile: use this to create the executable (Linux and MacOS X)

CC = g++
CCFLAGS = -O3 -Wall -W -fopenmp -std=c++11
INCLUDE = ./include
SRC = ./src
EXT_LIB = ./usr/local/lib/
BIN = ./
INCLUDE_PATHS = -I$(INCLUDE)
#debug: -g
LDFLAGS = -L$(EXT_LIB) -lgsl -lgslcblas -lm


ifeq ($(OS),Windows_NT)
    CCFLAGS += -D WIN32
    ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
        CCFLAGS += -D AMD64
    endif
    ifeq ($(PROCESSOR_ARCHITECTURE),x86)
        CCFLAGS += -D IA32
    endif
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        CCFLAGS += -D LINUX
    endif
    ifeq ($(UNAME_S),Darwin)
        CCFLAGS += -D OSX
    endif
    UNAME_P := $(shell uname -p)
    ifeq ($(UNAME_P),x86_64)
        CCFLAGS += -D AMD64
    endif
    ifneq ($(filter %86,$(UNAME_P)),)
        CCFLAGS += -D IA32
    endif
    ifneq ($(filter arm%,$(UNAME_P)),)
        CCFLAGS += -D ARM
    endif
endif


.PHONY : all clean

all : $(BIN)/survival

clean :
	@rm -f $(SRC)/*.o *~ $(BIN)/survival

ifeq (,$(findstring $(MAKECMDGOALS),clean))
-include dep.mk
endif

dep.mk : $(SRC)/*.cpp $(INCLUDE)/*.h
	g++ -MM $(INCLUDE_PATHS) $(SRC)/*.cpp > dep.mk
	
$(BIN)/survival : $(SRC)/main.o $(SRC)/Particles.o $(SRC)/Tracks.o $(SRC)/Track_Scholz2000.o $(SRC)/Track_Elsasser2007.o $(SRC)/Track_Elsasser2008.o $(SRC)/Track_KieferChatterjee.o $(SRC)/CellLine.o $(SRC)/Nucleus_Pixel.o $(SRC)/Nucleus_Integral.o $(SRC)/Nucleus_MonteCarlo.o $(SRC)/Nucleus_MKM.o $(SRC)/Calculus.o $(SRC)/Nucleus_tMKM.o $(SRC)/Nucleus_Integral_t.o $(SRC)/usefulFunctions.o
	$(CC) $(CCFLAGS) $(INCLUDE_PATHS) $^ $(LDFLAGS)  -o $@

%.o : %.cpp
	$(CC) $(CCFLAGS) $(INCLUDE_PATHS) -c $< -o $@

