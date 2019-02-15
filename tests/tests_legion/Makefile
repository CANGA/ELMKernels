ifndef LG_RT_DIR
$(error LG_RT_DIR variable is not defined, aborting build)
endif

# from Makefile.config
OBJECT        = ../../src/
NETCDF_ROOT   = /Users/uec/codes/ats/amanzi-tpls/install-0.88
STD_LIB_ROOT  = /usr



# Flags for directing the runtime makefile what to include
DEBUG           ?= 1		# Include debugging symbols
OUTPUT_LEVEL    ?= LEVEL_DEBUG	# Compile time logging level
SHARED_LOWLEVEL ?= 0		# Use shared-memory runtime (not recommended)
USE_CUDA        ?= 0		# Include CUDA support (requires CUDA)
USE_GASNET      ?= 0		# Include GASNet support (requires GASNet)
USE_HDF         ?= 0		# Include HDF5 support (requires HDF5)
ALT_MAPPERS     ?= 0		# Include alternative mappers (not recommended)

# Put the binary file name here
OUTFILE		?= CanopyHydrology_kern1_multiple_read_phenology
# List all the application source files here
GEN_SRC		?= CanopyHydrology_kern1_multiple_read_phenology.cc	# .cc files

# You can modify these variables, some will be appended to by the runtime makefile
INC_FLAGS	?= -I$(AMANZI_TPLS_DIR)/include -I../../src/cpp -I../tests_c
CC_FLAGS	?= -std=c++11
NVCC_FLAGS	?=
GASNET_FLAGS	?=
LD_FLAGS	?= -L$(OBJECT)fortran -lelm -L$(NETCDF_ROOT)/lib -lnetcdf -L$(STD_LIB_ROOT)/lib -lc++ 

###########################################################################
#
#   Don't change anything below here
#   
###########################################################################

include $(LG_RT_DIR)/runtime.mk

