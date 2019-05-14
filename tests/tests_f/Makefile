SRCDIR	      =	../../src/fortran/
OBJECT        = ../../src/

include $(OBJECT)config/Makefile.config

FC_FLAGS += -I$(SRCDIR) -I$(SRCDIR)$(ELM_UTILS_DIR) -I$(NETCDF_ROOT)/include


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
	./test_CanopyHydrology_kern1_single &> test_CanopyHydrology_kern1_single.stdout

CanopyHydrology_kern1_multiple: test_CanopyHydrology_kern1_multiple
	./test_CanopyHydrology_kern1_multiple &> test_CanopyHydrology_kern1_multiple.stdout

CanopyHydrology_module: test_CanopyHydrology_module
	./test_CanopyHydrology_module &> test_CanopyHydrology_module.stdout


test_%: %.F90.o library
	$(FC) -o $@ $< $(FC_LD_FLAGS)


clean:
	@$(ELM_CLEAN)
	$(RM) test_*

allclean:
	@$(ELM_CLEAN)
	$(RM) test_*
	$(MAKE) -C $(OBJECT) allclean

links:
	@echo "making in links"
	$(MAKE) -C ../links links

library:
	$(MAKE) -C $(OBJECT) fortran
