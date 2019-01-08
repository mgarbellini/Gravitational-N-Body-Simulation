# F. Borando - francesco.borando@studenti.unimi.it
# M. Garbellini - matteo.garbellini@studenti.unimi.it
# Università degli Studi di Milano
# Dept. of Physics
# January 2019
# GRAVITATIONAL N-BODY SIMULATION
#
# Description: Rules to build the code of graphic simulation.
# 

CXX := g++

EXECUTABLES := glut_Graphics

CXXFLAGS := $(CXXFLAGS) $(shell root-config --cflags) 

ROOT_LIBS := $(shell root-config --libs)

GL_LIBS := $(-framework OpenGL -framework GLUT -lglew -lglut)

all: $(EXECUTABLES)
tree_code_%: tree_code_%.o \
	$(CXX) -o $(@F) $(@F).o

tools_%: tools_%.o ; \
	$(CXX) -o $(@F) $(@F).o $(ROOT_LIBS)

glut_%: glut_%.o ; \
	$(CXX) -o $(@F) $(@F).o $(ROOT_LIBS)

%: %.cpp ; \
        $(LINK.cpp) $^ $(LOADLIBES) $(LDLIBS) -o $@m

clean: ; /bin/rm -f *.o

#like the zen master taught
mrproper: clean ; /bin/rm -f $(EXECUTABLES)