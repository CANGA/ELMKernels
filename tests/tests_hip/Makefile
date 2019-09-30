OBJECT        = ../../src/
KERNEL_LANG	  = cc
SRCDIR	      =	$(OBJECT)$(KERNEL_LANG)
HIP_FLAGS += -DNATURE=__global__
HIPCC =/home/7hp/Downloads/HIP/build/bin/hipcc
HIPFY=/home/7hp/Downloads/HIP/build/bin/hipify-perl

include $(OBJECT)config/Makefile.config
include $(OBJECT)config/Makefile.rules

INC_FLAGS ?= -I$(NETCDF_ROOT)/include -I$(SRCDIR) 

TESTS = test_CanopyHydrology_kern1_single \
        test_CanopyHydrology_kern1_multiple \
        test_CanopyHydrology_module

EXEC_TESTS = CanopyHydrology_kern1_single \
             CanopyHydrology_kern1_multiple \
             CanopyHydrology_module


.PHONY: links library test

default: all

all: links library $(TESTS)

hipfy: 
	$(HIPFY) CanopyHydrology_kern1_single.cu > CanopyHydrology_kern1_single.cpp
	$(HIPFY) CanopyHydrology_kern1_multiple.cu > CanopyHydrology_kern1_multiple.cpp
	$(HIPFY) CanopyHydrology_module.cu > CanopyHydrology_module.cpp

test: readers.hh utils.hh #$(EXEC_TESTS)
# 	python ../compare_to_gold.py $(TESTS)

# CanopyHydrology_kern1_single: test_CanopyHydrology_kern1_single
# 	./test_CanopyHydrology_kern1_single > test_CanopyHydrology_kern1_single.stdout

# CanopyHydrology_kern1_multiple: test_CanopyHydrology_kern1_multiple
# 	./test_CanopyHydrology_kern1_multiple > test_CanopyHydrology_kern1_multiple.stdout

# CanopyHydrology_module: test_CanopyHydrology_module
# 	./test_CanopyHydrology_module > test_CanopyHydrology_module.stdout

# test_%: readers.hh utils.hh
	#$(HIPCC) CanopyHydrology_kern1_single.cpp $(INC_FLAGS) $(HIP_LD_FLAGS) -o CanopyHydrology_kern1_single.out
	$(HIPCC) $(HIP_FLAGS) CanopyHydrology_kern1_multiple.cpp $(INC_FLAGS) $(HIP_LD_FLAGS) -o CanopyHydrology_kern_multiple.out
	#$(HIPCC) CanopyHydrology_module.cpp $(INC_FLAGS) $(HIP_LD_FLAGS) -o CanopyHydrology_module.out


clean:
	@$(ELM_CLEAN)
	$(RM) test_* *.cpp

allclean:
	@$(ELM_CLEAN)
	$(RM) test_*
	$(MAKE) -C $(SRCDIR) allclean

links:
	@echo "making in links"
	$(MAKE) -C ../links links