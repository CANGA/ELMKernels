OBJECT        = ../../src/
KERNEL_LANG	  = cuda
SRCDIR	      =	$(OBJECT)$(KERNEL_LANG)

include $(OBJECT)config/Makefile.config

INC_FLAGS = -I$(AMANZI_TPLS_DIR)/include -I../../src/cuda/ -I../tests_cuda

TESTS = test_CanopyHydrology_kern1_multiple

EXEC_TESTS = CanopyHydrology_kern1_multiple


.PHONY: links library test

default: all

all: links library $(TESTS)

test: $(EXEC_TESTS)
	python ../compare_to_gold.py $(TESTS)


CanopyHydrology_kern1_single: test_CanopyHydrology_kern1_single
	./test_CanopyHydrology_kern1_multiple&> test_CanopyHydrology_kern1_single.stdout


test_%: %.cu.o readers.hh utils.hh library
	$(NVCC) -o $@ $< $(CC_LD_FLAGS)


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
	$(MAKE) -C $(OBJECT) all