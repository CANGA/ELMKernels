ELM_BASE = ../../src
KERNEL_FOLDER = cc
ELM_ROOT = $(ELM_BASE)/$(KERNEL_FOLDER)

#
# COMPILERS
# =================================
HAVE_MPI = 1
HAVE_PNETCDF = 1

include $(ELM_BASE)/config/Makefile.config
include $(ELM_BASE)/config/Makefile.rules

TESTS = test_CanopyHydrology_module
EXEC_TESTS = CanopyHydrology_module
CXXFLAGS += -Wall -Wshadow -O3 -std=c++14

VPATH = ../utils

.PHONY: default library test compute

default: all

all: logo report library $(TESTS)

include $(ELM_BASE)/config/Makefile.logo

test: all $(EXEC_TESTS)
	python ../compare_to_gold.py $(TESTS)

report:
	@echo "Building with:"
	@echo "---------------"
	@echo "CXX = $(CXX)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "INC_FLAGS = $(INC_FLAGS)"
	@echo ""

CanopyHydrology_module: test_CanopyHydrology_module
	mpiexec -n 6 $(MPI_EXEC_GLOBAL_ARGS) ./test_CanopyHydrology_module 
	#jsrun -n1 -c20 -a20 ./test_CanopyHydrology_module 

test_%: %.cc.o utils.hh utils.cc.o array.hh read_pnetcdf.hh read_input.hh read_input.cc.o
	$(CXX) -o $@ utils.cc.o read_input.cc.o $< $(CXX_LD_FLAGS)

clean:
	$(RM) test_* mpi_* *.o 

allclean: clean

