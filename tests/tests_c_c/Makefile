ELM_BASE = ../../src
KERNEL_FOLDER = cc
ELM_ROOT = $(ELM_BASE)/$(KERNEL_FOLDER)

#
# COMPILERS
# =================================
CXX	      = mpic++
CXXFLAGS      += -Wall -Wshadow -O3 -std=c++14

include $(ELM_BASE)/config/Makefile.config
include $(ELM_BASE)/config/Makefile.rules

TESTS = test_CanopyHydrology_module

TEST1 = mpi_CanopyHydrology_module

EXEC_TESTS = CanopyHydrology_module

VPATH = ../utils

.PHONY: library test compute

default: main

include $(ELM_BASE)/config/Makefile.logo

main: logo $(TESTS)

test: main $(EXEC_TESTS)
	python ../compare_to_gold.py $(TESTS)

CanopyHydrology_module: test_CanopyHydrology_module
	mpiexec -n 6 ./test_CanopyHydrology_module 

test_%: %.cc.o utils.hh array.hh readers.hh readers.cc.o
	$(CXX) -o $@ readers.cc.o $< $(NETCDF_FLAGS)  

clean:
	$(RM) test_* mpi_* *.o 

allclean: clean

