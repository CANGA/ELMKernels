ELM_BASE = ../../src
KERNEL_FOLDER = cc
ELM_ROOT = $(ELM_BASE)/$(KERNEL_FOLDER)

KOKKOS_DEBUG = "yes"
EXE_NAME = test_CanopyHydrology_module

PARALLEL_NETCDF = TRUE
include $(ELM_BASE)/config/Makefile.config

#
# Set up Kokkos flags custom by device
#
ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))

CXXFLAGS += -Wall -Wshadow -ccbin ${NVCC_WRAPPER_DEFAULT_COMPILER}
EXE = ${EXE_NAME}.cuda

else

CXX	      = mpic++
CXXFLAGS      += -Wall -Wshadow -std=c++11
EXE = ${EXE_NAME}.host

endif

CXXFLAGS += -DNATURE=KOKKOS_INLINE_FUNCTION 

# build the tests
TEST = test_CanopyHydrology_module
SRC = CanopyHydrology_module.cc     
OBJ =  $(SRC:.cc=.cc.o)	



default: main

# include Kokkos
include $(KOKKOS_ROOT)/Makefile.kokkos
include $(ELM_BASE)/config/Makefile.rules


LINK = ${CXX}
VPATH = ../utils

include $(ELM_BASE)/config/Makefile.logo

.PHONY: logo all

main: logo report $(EXE)

report: 
	@echo "Kokkos setup"
	@echo "----------------"
	@echo "LINK = $(LINK)"
	@echo "CXX = $(CXX)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "CXX_LD_FLAGS = $(CXX_LD_FLAGS)"
	@echo "KOKKOS_CXXFLAGS = $(KOKKOS_CXXFLAGS)"
	@echo "KOKKOS_LDFLAGS = $(KOKKOS_LDFLAGS)"
	@echo "KOKKOS_LIBS = $(KOKKOS_LIBS)"

build_msg:
	@echo "Building"
	@echo "----------------"

build: build_msg $(OBJ) $(KOKKOS_LINK_DEPENDS) utils.hh kokkos_utils.hh array.hh readers.hh readers.cc.o

$(EXE): report build 
	@echo "Linking"
	@echo "----------------"
	$(LINK) readers.cc.o $(OBJ) $(CXX_LD_FLAGS) $(KOKKOS_LDFLAGS) $(KOKKOS_LIBS) -o $(EXE)

clean:
	@$(ELM_CLEAN)
	rm -f test_* *.o *.cuda *.host KokkosCore_config.h KokosCore_config.tmp

test : main
	${MPIEXEC} ${MPI_NP}6 ./$(EXE) > test_CanopyHydrology_module_multiple.stdout
	python ../compare_to_gold.py $(TESTS)

allclean: clean
	$(MAKE) -C $(SRCDIR) clean
