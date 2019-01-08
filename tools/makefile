# F. Borando - francesco.borando@studenti.unimi.it
# M. Garbellini - matteo.garbellini@studenti.unimi.it
# Università degli Studi di Milano
# Dept. of Physics
# January 2019
# GRAVITATIONAL N-BODY SIMULATION
#
# Description: Rules to build the code for analysis tools (includes ROOT).
# 

CXX := g++

EXECUTABLES := tools_plot_errors tools_plot_energy tools_energy_data_fixer tools_dispersion

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
