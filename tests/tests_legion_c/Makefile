ifndef LG_RT_DIR
$(error LG_RT_DIR variable is not defined, aborting build)
endif

# from Makefile.config
OBJECT        = ../../src/
KERNEL_LANG	  = cc
SRCDIR	      =	$(OBJECT)$(KERNEL_LANG)
# Flags for directing the runtime makefile what to include
DEBUG           ?= 1		# Include debugging symbols
OUTPUT_LEVEL    ?= LEVEL_DEBUG	# Compile time logging level
SHARED_LOWLEVEL ?= 0		# Use shared-memory runtime (not recommended)
USE_CUDA        ?= 0		# Include CUDA support (requires CUDA)
USE_GASNET      ?= 0		# Include GASNet support (requires GASNet)
USE_HDF         ?= 0		# Include HDF5 support (requires HDF5)
ALT_MAPPERS     ?= 0		# Include alternative mappers (not recommended)

# Put the binary file name here

# You can modify these variables, some will be appended to by the runtime makefile
GEN_SRC		?= tasks.cc
INC_FLAGS	?= -I$(AMANZI_TPLS_DIR)/include -I../../src/cc -I../tests_c_legion
NVCC_FLAGS	?=
GASNET_FLAGS	?=

include $(OBJECT)config/Makefile.config

TESTS = test_CanopyHydrology_kern1_multiple \
        test_CanopyHydrology_module


EXEC_TESTS = CanopyHydrology_kern1_multiple \
             CanopyHydrology_module



#
# Load Legion logic, but not rules
#
NO_BUILD_RULES = 1
include $(LG_RT_DIR)/runtime.mk


#
# Legion rules
#
.PHONY: links library test

default: all

all: links library $(TESTS)

test: $(EXEC_TESTS)
	python ../compare_to_gold.py $(TESTS)


CanopyHydrology_kern1_multiple: test_CanopyHydrology_kern1_multiple
	./test_CanopyHydrology_kern1_multiple &> test_CanopyHydrology_kern1_multiple.stdout

CanopyHydrology_module: test_CanopyHydrology_module
	./test_CanopyHydrology_module &> test_CanopyHydrology_module.stdout


sandbox: test_sandbox_domain_template_magic
	./test_sandbox_domain_template_magic

# If we're using CUDA we have to link with nvcc
test_%: %.cc.o $(GEN_OBJS) $(GEN_GPU_OBJS) $(SLIB_LEGION) $(SLIB_REALM)
	@echo "---> Linking objects into one binary: $()"
	$(CXX) -o $@ $< $(GEN_OBJS) $(GEN_GPU_OBJS) $(CC_LD_FLAGS) $(LEGION_LIBS) $(LEGION_LD_FLAGS) $(GASNET_FLAGS)

$(SLIB_LEGION) : $(LEGION_OBJS) $(MAPPER_OBJS)
	rm -f $@
	$(AR) rc $@ $^

$(SLIB_REALM) : $(REALM_OBJS)
	rm -f $@
	$(AR) rc $@ $^

$(GEN_OBJS) : %.cc.o : %.cc
	$(CXX) -o $@ -c $< $(CC_FLAGS) $(CPP_FLAGS) $(INC_FLAGS)

$(ASM_OBJS) : %.S.o : %.S
	$(CXX) -o $@ -c $< $(CC_FLAGS) $(INC_FLAGS)

$(REALM_OBJS) : %.cc.o : %.cc
	$(CXX) -o $@ -c $< $(CC_FLAGS) $(INC_FLAGS)

$(LEGION_OBJS) : %.cc.o : %.cc
	$(CXX) -o $@ -c $< $(CC_FLAGS) $(INC_FLAGS)

$(MAPPER_OBJS) : %.cc.o : %.cc
	$(CXX) -o $@ -c $< $(CC_FLAGS) $(INC_FLAGS)

$(GEN_GPU_OBJS) : %.cu.o : %.cu
	$(NVCC) -o $@ -c $< $(NVCC_FLAGS) $(INC_FLAGS)

$(GPU_RUNTIME_OBJS): %.cu.o : %.cu
	$(NVCC) -o $@ -c $< $(NVCC_FLAGS) $(INC_FLAGS)

# disable gmake's default rule for building % from %.o
#% : %.o

%.cc.o : %.cc domains.hh readers.hh utils.hh
	$(CXX) $(CC_FLAGS) $(CPP_FLAGS) $(INC_FLAGS) -c $< -o $@

clean :
	$(ELM_CLEAN)

allclean : clean
	$(MAKE) -C $(OBJECT) allclean
	$(RM) -f $(SLIB_LEGION) $(SLIB_REALM) $(GEN_OBJS) $(GEN_GPU_OBJS) $(REALM_OBJS) $(LEGION_OBJS) $(GPU_RUNTIME_OBJS) $(MAPPER_OBJS) $(ASM_OBJS)

ifeq ($(strip $(USE_LLVM)),1)
llvmjit_internal.cc.o : CC_FLAGS += $(LLVM_CXXFLAGS)
%/llvmjit_internal.cc.o : CC_FLAGS += $(LLVM_CXXFLAGS)
endif



links:
	@echo "making in links"
	$(MAKE) -C ../links links


library:
	$(MAKE) -C $(OBJECT) all

