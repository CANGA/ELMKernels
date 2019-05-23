SRCDIR	      =	../../src/fortran/
CSRCDIR	      =	../../src/cc_fortran_wrappers/
OBJECT        = ../../src/

include $(OBJECT)config/Makefile.config


INC_FLAGS = -I$(NETCDF_ROOT)/include -I$(CSRCDIR) -I../tests_c


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


test_%: %.cc.o readers.hh utils.hh library
	$(CXX) -o $@ $< $(CC_LD_FLAGS)


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
	$(MAKE) -C $(OBJECT) cc_fortran_wrappers 
