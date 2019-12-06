ELM_BASE = ../../src
KERNEL_FOLDER = cc
ELM_ROOT = $(ELM_BASE)/$(KERNEL_FOLDER)

include $(ELM_BASE)/config/Makefile.config
include $(ELM_BASE)/config/Makefile.rules

TESTS = test_CanopyHydrology_module

TEST1 = mpi_CanopyHydrology_module

EXEC_TESTS = CanopyHydrology_module

VPATH = ../utils

.PHONY: links library test compute

default: test

include $(ELM_BASE)/config/Makefile.logo


test: debug logo links library $(TESTS) $(EXEC_TESTS)
	python ../compare_to_gold.py $(TESTS)

debug:
	$(eval CXXFLAGS += -DUNIT_TEST)

CanopyHydrology_module: test_CanopyHydrology_module
	mpiexec --oversubscribe -n 6 ./test_CanopyHydrology_module > test_CanopyHydrology_module.stdout

test_%: %.cc.o utils.hh array.hh readers.hh readers.cc.o
	$(CXX) -o $@ readers.cc.o $< $(NETCDF_FLAGS) 

compute: logo links library $(TEST1)

clean:
	$(RM) test_* mpi_* *.o 

allclean: clean

links:
	@echo "making in links"
	$(MAKE) -C ../links links
