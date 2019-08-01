OBJECT        = ../../src/
KERNEL_LANG	  = cc
SRCDIR	      =	$(OBJECT)$(KERNEL_LANG)
CUDA_FLAGS += -DNATURE=__global__
NVCC = nvcc
include $(OBJECT)config/Makefile.config
include $(OBJECT)config/Makefile.rules

INC_FLAGS ?= -I$(NETCDF_ROOT)/include -I$(SRCDIR) -I../tests_cuda

TESTS = test_CanopyHydrology_kern1_single \
        test_CanopyHydrology_kern1_multiple \
        test_CanopyHydrology_module

EXEC_TESTS = CanopyHydrology_kern1_single \
             CanopyHydrology_kern1_multiple \
             CanopyHydrology_module


.PHONY: links library test

default: all

all: links library $(TESTS)

test: $(EXEC_TESTS)
	python ../compare_to_gold.py $(TESTS)


CanopyHydrology_kern1_single: test_CanopyHydrology_kern1_single
	./test_CanopyHydrology_kern1_single > test_CanopyHydrology_kern1_single.stdout

CanopyHydrology_kern1_multiple: test_CanopyHydrology_kern1_multiple
	./test_CanopyHydrology_kern1_multiple > test_CanopyHydrology_kern1_multiple.stdout

CanopyHydrology_module: test_CanopyHydrology_module
	./test_CanopyHydrology_module > test_CanopyHydrology_module.stdout

test_%: %.cu.o readers.hh utils.hh
	$(NVCC) -o $@ $< $(CUDA_LD_FLAGS)


clean:
	@$(ELM_CLEAN)
	$(RM) test_*

allclean:
	@$(ELM_CLEAN)
	$(RM) test_*
	$(MAKE) -C $(SRCDIR) allclean

links:
	@echo "making in links"
	$(MAKE) -C ../links links