ifndef LG_RT_DIR
$(error LG_RT_DIR variable is not defined, aborting build)
endif

# from Makefile.config
ELM_BASE        = ../../src
KERNEL_FOLDER = cc
ELM_ROOT	= $(ELM_BASE)/$(KERNEL_FOLDER)
NATURE		=

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
OUTFILE = CanopyHydrology_module
GEN_SRC	= tasks.cc readers_seq.cc CanopyHydrology_module.cc

CXXFLAGS        += -DLEGION_SPY

include $(ELM_BASE)/config/Makefile.config
include $(ELM_BASE)/config/Makefile.rules

INC_FLAGS += -I ../utils $(NETCDF_INC_FLAGS)
CC_FLAGS += $(CXXFLAGS)
LD_FLAGS += $(CXX_LD_FLAGS)

TESTS = test_CanopyHydrology_module
EXEC_TESTS = CanopyHydrology_module

VPATH = ../utils

#
# Load Legion logic, but not rules
#

.PHONY: default
default: all
	@echo "GEN OBJS = $(GEN_OBJS)"

all: $(TESTS)

test: $(EXEC_TESTS)
	../compare_to_gold.py $(TESTS)


CanopyHydrology_module: test_CanopyHydrology_module
	./test_CanopyHydrology_module > test_CanopyHydrology_module.stdout
# 	./test_CanopyHydrology_module -lg:spy -logfile spy_module%.log
# 	$(LG_RT_DIR)/../tools/legion_spy.py -dez spy_module*.log

include $(LG_RT_DIR)/runtime.mk

# test: $(EXEC_TESTS)
# 	python ../compare_to_gold.py $(TESTS)

# CanopyHydrology_module: test_CanopyHydrology_module
# 	./test_CanopyHydrology_module > test_CanopyHydrology_module.stdout


allclean: clean
	$(MAKE) -C $(OBJECT) allclean
	$(RM) -f $(SLIB_LEGION) $(SLIB_REALM) $(GEN_OBJS) $(GEN_GPU_OBJS) $(REALM_OBJS) $(LEGION_OBJS) $(GPU_RUNTIME_OBJS) $(MAPPER_OBJS) $(ASM_OBJS)


