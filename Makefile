INC = -I./include
LIBDIR = lib
SRCDIR = src
BUILDDIR = build
LDFLAGS = -lxdrfile -larmadillo -lboost_program_options
STAT_LIB = $(LIBDIR)/libz_vec.a $(LIBDIR)/libz_gromacs.a $(LIBDIR)/libz_molecule.a $(LIBDIR)/libz_atom.a $(LIBDIR)/libz_atom_group.a $(LIBDIR)/libz_string.a $(LIBDIR)/libz_sim_params.a $(LIBDIR)/libz_tcf.a
CFLAGS += -g -Wall -std=c++11
LIBFLAGS = -c -static -fpic
TARGET = bin/z_inst_surf

all:	build vec gromacs molecule atom atom_group string sim_params tcf inst_surf histogram

build:
	mkdir -p $(BUILDDIR)	

vec:
	$(CXX) $(CFLAGS) $(LIBFLAGS) $(INC) -o $(BUILDDIR)/z_vec.o $(SRCDIR)/z_vec.cpp $(MODE)
	ar rvs $(LIBDIR)/libz_vec.a $(BUILDDIR)/z_vec.o

gromacs:
	$(CXX) $(CFLAGS) $(LIBFLAGS) $(INC) -o $(BUILDDIR)/z_gromacs.o $(SRCDIR)/z_gromacs.cpp  $(MODE)
	ar rvs $(LIBDIR)/libz_gromacs.a $(BUILDDIR)/z_gromacs.o

molecule:
	$(CXX) $(CFLAGS) $(LIBFLAGS) $(INC) -o $(BUILDDIR)/z_molecule.o $(SRCDIR)/z_molecule.cpp $(MODE)
	ar rvs $(LIBDIR)/libz_molecule.a $(BUILDDIR)/z_molecule.o

atom:
	$(CXX) $(CFLAGS) $(LIBFLAGS) $(INC) -o $(BUILDDIR)/z_atom.o $(SRCDIR)/z_atom.cpp $(MODE)
	ar rvs $(LIBDIR)/libz_atom.a $(BUILDDIR)/z_atom.o

atom_group:
	$(CXX) $(CFLAGS) $(LIBFLAGS) $(INC) -o $(BUILDDIR)/z_atom_group.o $(SRCDIR)/z_atom_group.cpp $(MODE)
	ar rvs $(LIBDIR)/libz_atom_group.a $(BUILDDIR)/z_atom_group.o

string:
	$(CXX) $(CFLAGS) $(LIBFLAGS) $(INC) -o $(BUILDDIR)/z_string.o $(SRCDIR)/z_string.cpp $(MODE)
	ar rvs $(LIBDIR)/libz_string.a $(BUILDDIR)/z_string.o

sim_params:
	$(CXX) $(CFLAGS) $(LIBFLAGS) $(INC) -o $(BUILDDIR)/z_sim_params.o $(SRCDIR)/z_sim_params.cpp $(MODE)
	ar rvs $(LIBDIR)/libz_sim_params.a $(BUILDDIR)/z_sim_params.o

tcf:
	$(CXX) $(CFLAGS) $(LIBFLAGS) $(INC) -o $(BUILDDIR)/z_tcf.o $(SRCDIR)/z_tcf.cpp $(MODE)
	ar rvs $(LIBDIR)/libz_tcf.a $(BUILDDIR)/z_tcf.o

histogram:
	$(CXX) $(CFLAGS) $(LIBFLAGS) $(INC) -o $(BUILDDIR)/z_histogram.o $(SRCDIR)/z_histogram.cpp $(MODE)
	ar rvs $(LIBDIR)/libz_histogram.a $(BUILDDIR)/z_histogram.o

inst_surf:
	$(CXX) $(CFLAGS) $(INC) $(SRCDIR)/main.cpp $(STAT_LIB) -o $(TARGET) $(LDFLAGS) $(MODE)

clean:
	rm lib/*.a build/*.o
