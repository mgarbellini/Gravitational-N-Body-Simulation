# F. Borando - francesco.borando@studenti.unimi.it
# M. Garbellini - matteo.garbellini@studenti.unimi.it
# Università degli Studi di Milano
# Dept. of Physics
# January 2019
# GRAVITATIONAL N-BODY SIMULATION
#
# Description: Rules to build the code for Runge-Kutta simulations.
# 

CXX := g++

EXECUTABLES := tree_code_RK_energy tree_code_RK

CXXFLAGS := $(CXXFLAGS) $(shell root-config --cflags) 

ROOT_LIBS := $(shell root-config --libs)

all: $(EXECUTABLES)
tree_code_%: tree_code_%.o \
	$(CXX) -o $(@F) $(@F).o

tools_%: tools_%.o ; \
	$(CXX) -o $(@F) $(@F).o $(ROOT_LIBS)

%: %.cpp ; \
        $(LINK.cpp) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean: ; /bin/rm -f *.o

#like the zen master taught
mrproper: clean ; /bin/rm -f $(EXECUTABLES)
