ELM_BASE = ../../src
KERNEL_FOLDER = cc
ELM_ROOT = $(ELM_BASE)/$(KERNEL_FOLDER)

# Flags for directing the runtime makefile what to include
DEBUG           ?= 1		# Include debugging symbols
OUTPUT_LEVEL    ?= LEVEL_DEBUG	# Compile time logging level
SHARED_LOWLEVEL ?= 0		# Use shared-memory runtime (not recommended)
USE_CUDA        ?= 0		# Include CUDA support (requires CUDA)
USE_GASNET      ?= 0		# Include GASNet support (requires GASNet)
USE_HDF         ?= 0		# Include HDF5 support (requires HDF5)
ALT_MAPPERS     ?= 0		# Include alternative mappers (not recommended)

# include the user's config
include $(ELM_BASE)/config/Makefile.config

# legion variables
CC_FLAGS= -DLEGION_SPY
NVCC_FLAGS	?=
GASNET_FLAGS	?=

# note: INC_FLAGS is defined by legion
CXXFLAGS += $(NETCDF_INC_FLAGS) $(INC_FLAGS)

CXX_LD_FLAGS += $(NETCDF_FLAGS)

#
# Load Legion logic, but not rules
#
NO_BUILD_RULES = 1
include $(LG_RT_DIR)/runtime.mk

# now include our rules
include $(ELM_BASE)/config/Makefile.rules

TESTS = test_CanopyHydrology_kern1_multiple \
        test_CanopyHydrology_module


EXEC_TESTS = CanopyHydrology_kern1_multiple \
             CanopyHydrology_module

#
# Legion rules
#
.PHONY: links test

default: all


all: links $(TESTS)

test: $(EXEC_TESTS)
	python ../compare_to_gold.py $(TESTS)


CanopyHydrology_kern1_multiple: test_CanopyHydrology_kern1_multiple
	./test_CanopyHydrology_kern1_multiple > test_CanopyHydrology_kern1_multiple.stdout
# 	./test_CanopyHydrology_kern1_multiple -lg:spy -logfile spy_multiple%.log
# 	$(LG_RT_DIR)/../tools/legion_spy.py -dez spy_multiple*.log

CanopyHydrology_module: test_CanopyHydrology_module
	./test_CanopyHydrology_module > test_CanopyHydrology_module.stdout
# 	./test_CanopyHydrology_module -lg:spy -logfile spy_module%.log
# 	$(LG_RT_DIR)/../tools/legion_spy.py -dez spy_module*.log


sandbox: test_sandbox_domain_template_magic
	./test_sandbox_domain_template_magic

# If we're using CUDA we have to link with nvcc
test_%: %.cpp.o tasks.cpp.o $(SLIB_LEGION) $(SLIB_REALM)
	@echo "---> Linking objects into one binary: $()"
	$(CXX) -o $@ $< tasks.cpp.o $(CXX_LD_FLAGS) $(LEGION_LIBS) $(LEGION_LD_FLAGS) $(GASNET_FLAGS)

$(SLIB_LEGION) : $(LEGION_OBJS) $(MAPPER_OBJS)
	rm -f $@
	$(AR) rc $@ $^

$(SLIB_REALM) : $(REALM_OBJS)
	rm -f $@
	$(AR) rc $@ $^

$(GEN_OBJS) : %.cpp.o : %.cpp

$(ASM_OBJS) : %.S.o : %.S

$(REALM_OBJS) : %.cc.o : %.cc

$(LEGION_OBJS) : %.cc.o : %.cc

$(MAPPER_OBJS) : %.cc.o : %.cc

$(GEN_GPU_OBJS) : %.cu.o : %.cu

$(GPU_RUNTIME_OBJS): %.cu.o : %.cu

clean :
	$(ELM_CLEAN)
	rm -f test_*

allclean : clean
	$(MAKE) -C $(OBJECT) allclean
	$(RM) -f $(SLIB_LEGION) $(SLIB_REALM) $(GEN_OBJS) $(GEN_GPU_OBJS) $(REALM_OBJS) $(LEGION_OBJS) $(GPU_RUNTIME_OBJS) $(MAPPER_OBJS) $(ASM_OBJS)

ifeq ($(strip $(USE_LLVM)),1)
llvmjit_internal.cpp.o : CC_FLAGS += $(LLVM_CXXFLAGS)
%/llvmjit_internal.cpp.o : CC_FLAGS += $(LLVM_CXXFLAGS)
endif



links:
	@echo "making in links"
	$(MAKE) -C ../links links


library:
	$(MAKE) -C $(OBJECT) all

