ELM_BASE = ../../src
KERNEL_FOLDER = cc
ELM_ROOT = $(ELM_BASE)/$(KERNEL_FOLDER)

include $(ELM_BASE)/config/Makefile.config
include $(ELM_BASE)/config/Makefile.rules

TESTS = test_CanopyHydrology_kern1_single \
        test_CanopyHydrology_kern1_multiple \
        test_CanopyHydrology_module

TEST1 = mpi_CanopyHydrology_kern1_single \
        mpi_CanopyHydrology_kern1_multiple \
        mpi_CanopyHydrology_module

EXEC_TESTS = CanopyHydrology_kern1_single \
             CanopyHydrology_kern1_multiple \
             CanopyHydrology_module

.PHONY: links library test compute

default: test

include $(ELM_BASE)/config/Makefile.logo


test: debug logo links library $(TESTS) $(EXEC_TESTS)
	python ../compare_to_gold.py $(TESTS)

debug:
	$(eval CXXFLAGS += -DDEBUG)

CanopyHydrology_kern1_single: test_CanopyHydrology_kern1_single
	./test_CanopyHydrology_kern1_single > test_CanopyHydrology_kern1_single.stdout

CanopyHydrology_kern1_multiple: test_CanopyHydrology_kern1_multiple
	./test_CanopyHydrology_kern1_multiple > test_CanopyHydrology_kern1_multiple.stdout

CanopyHydrology_module: test_CanopyHydrology_module
	./test_CanopyHydrology_module > test_CanopyHydrology_module.stdout

test_%: %.cpp.o readers.hh utils.hh 
	$(CXX) -o $@ $< $(NETCDF_FLAGS)

compute: logo links library $(TEST1)

mpi_%: %.cpp.o readers.hh utils.hh
	$(MPICXX) -DMPICOMP -o $@ $< $(NETCDF_FLAGS) $(MPI_INC_FLAGS) -L$(MPI_LIB_ROOT)/lib -lmpi

clean:
	$(RM) test_* mpi_* *.o 

allclean: clean

links:
	@echo "making in links"
	$(MAKE) -C ../links links
