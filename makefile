
LT = /Users/ycwu/Library/Mathematica/Applications/LoopTools/x86_64-Darwin/lib
#LT = /Users/mac/work/LoopTools/x86_64-Darwin/lib

MAINDIR := $(shell pwd)
SRCDIR := src
INCDIR := include
OBJDIR := obj
# CXX = $(LT)/../bin/f++
CXX = clang++
FLAGS = -I$(LT)/../include -I$(INCDIR)
LIBS = -L$(LT) -looptools -Wl,-no_compact_unwind -L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin16/6.3.0 -L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin16/6.3.0/../../.. -lgfortran -lSystem -lgcc_ext.10.5 -lgcc -lquadmath -lm -lgcc_ext.10.5 -lgcc -lSystem -L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin16/6.3.0 -L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin16/6.3.0/../../.. -lgfortran -lSystem -lgcc_ext.10.5 -lgcc -lquadmath -lm -lgcc_ext.10.5 -lgcc -lSystem -m64

gsllibs = $(shell gsl-config --libs)
gslcflags = $(shell gsl-config --cflags)

muSTULibPath := $(MAINDIR)/../HiggsSignalStrength_STU
muSTUFLAGS := -I$(muSTULibPath)/
muSTULIBS := -L$(muSTULibPath)/ -lHiggsSignalStrengthSTU -Wl,-rpath,$(muSTULibPath)

FLAGS += $(gslcflags) $(muSTUFLAGS)
LIBS += $(gsllibs) $(muSTULIBS)

SRC = $(wildcard $(SRCDIR)/*.cpp)
OBJ = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRC))

all: SMCSHiggsCoupCalc.x

.PHONY: clean

.SECONDARY: $(OBJ)

%.x:%.cpp $(OBJ)
	$(CXX) $(FLAGS) $(LIBS) -o $@ $< $(OBJ) $(LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(FLAGS) -c $< -o $@


clean:
	rm -f *.x
	rm -f $(OBJDIR)/*.o

