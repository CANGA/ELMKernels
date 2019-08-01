ELM_BASE = ../../src
KERNEL_FOLDER = fortran
ELM_ROOT = $(ELM_BASE)/$(KERNEL_FOLDER)

FC = gfortran

include $(ELM_BASE)/config/Makefile.config
include $(ELM_BASE)/config/Makefile.rules

FC_FLAGS += -I$(ELM_ROOT) -I$(ELM_ROOT)/$(ELM_UTILS_DIR) -I$(NETCDF_ROOT)/include

TESTS = test_CanopyHydrology_kern1_single \
        test_CanopyHydrology_kern1_multiple \
	test_CanopyHydrology_module

EXEC_TESTS = CanopyHydrology_kern1_single \
             CanopyHydrology_kern1_multiple \
             CanopyHydrology_module

.PHONY: links library test

default: all

all: library $(TESTS)

test: $(EXEC_TESTS)
	python ../compare_to_gold.py $(TESTS)


CanopyHydrology_kern1_single: test_CanopyHydrology_kern1_single
	./test_CanopyHydrology_kern1_single > test_CanopyHydrology_kern1_single.stdout

CanopyHydrology_kern1_multiple: test_CanopyHydrology_kern1_multiple
	./test_CanopyHydrology_kern1_multiple > test_CanopyHydrology_kern1_multiple.stdout

CanopyHydrology_module: test_CanopyHydrology_module
	./test_CanopyHydrology_module > test_CanopyHydrology_module.stdout


test_%: %.F90.o library
	$(FC) -o $@ $< $(FC_LD_FLAGS)


clean:
	@$(ELM_CLEAN)
	$(RM) test_*

allclean:
	@$(ELM_CLEAN)
	$(RM) test_*
	$(MAKE) -C $(ELM_ROOT) allclean

links:
	@echo "making in links"
	$(MAKE) -C ../links links

library:
	$(MAKE) -C $(ELM_ROOT)
