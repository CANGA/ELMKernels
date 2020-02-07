ELM_BASE = ../../src
KERNEL_FOLDER = cc
ELM_ROOT = $(ELM_BASE)/$(KERNEL_FOLDER)

KOKKOS_DEBUG = "yes"
EXE_NAME = test_CanopyHydrology_module

include $(ELM_BASE)/config/Makefile.config

# Get the Kokkos Device
ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
EXE = ${EXE_NAME}.cuda
else
EXE = ${EXE_NAME}.host
endif

CXXFLAGS += -DNATURE=KOKKOS_INLINE_FUNCTION 

default: all

# include Kokkos
include $(KOKKOS_ROOT)/Makefile.kokkos
CXXFLAGS += ${KOKKOS_CXXFLAGS}
include $(ELM_BASE)/config/Makefile.rules


# build the tests
TEST = test_CanopyHydrology_module
SRC = CanopyHydrology_module.cc     
OBJ =  $(SRC:.cc=.cc.o)	

LINK = ${CXX}
VPATH = ../utils

include $(ELM_BASE)/config/Makefile.logo

.PHONY: logo all

all: logo build test

report: 
	@echo "LINKING"
	@echo "----------------"
	@echo "LINK = $(LINK)"
	@echo "KOKKOS_LDFLAGS = $(KOKKOS_LDFLAGS)"
	@echo "OBJ = $(OBJ)"
	@echo "NETCDF_FLAGS = $(NETCDF_FLAGS)"
	@echo "KOKKOS_LIBS = $(KOKKOS_LIBS)"
	@echo "EXE = $(EXE)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "CXX_LD_FLAGS = $(CXX_LD_FLAGS)"
	@echo "KOKKOS_CXXFLAGS = $(KOKKOS_CXXFLAGS)"
	@echo "KOKKOS_CPPFLAGS = $(KOKKOS_CPPFLAGS)"

build: $(OBJ) $(KOKKOS_LINK_DEPENDS) utils.hh kokkos_utils.hh array.hh readers.hh readers.cc.o $(EXE)

$(EXE): build report
	$(LINK) readers.cc.o $(OBJ) $(NETCDF_FLAGS) $(KOKKOS_LDFLAGS) $(KOKKOS_LIBS) $(CXX_LD_FLAGS) -o $(EXE)

clean:
	@$(ELM_CLEAN)
	rm -f test_* *.o *.cuda *.host

test : $(SRC)
	./$(EXE) > test_CanopyHydrology_module_multiple.stdout
	python ../compare_to_gold.py $(TESTS)

allclean: clean
	$(MAKE) -C $(SRCDIR) clean
