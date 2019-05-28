KOKKOS_DEVICES = "Serial"
EXE_NAME1 = "test_CanopyHydrology_kern1_multiple"
EXE_NAME1 = "test_CanopyHydrology_module"
OBJECT  = ../../../src/
KERNEL_LANG	  = cc_serial
SRCDIR	      =	$(OBJECT)$(KERNEL_LANG)
include $(OBJECT)config/Makefile.config

SRC1 = CanopyHydrology_kern1_multiple.cpp
SRC2 = CanopyHydrology_module.cpp

default: build1 build2
	echo "Start Build"


ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
#EXE = ${EXE_NAME}.cuda
KOKKOS_ARCH = "BSW,Pascal60"
KOKKOS_CUDA_OPTIONS = "enable_lambda"
else
CXX = g++
#EXE = ${EXE_NAME}.host
KOKKOS_ARCH = "BSW"
endif

CXXFLAGS = -g -O0
LINK = ${CXX}
LINKFLAGS = -lnetcdf
EXTRA_PATH = -I/usr/local/include

DEPFLAGS = -M

OBJ1 =  $(SRC1:.cpp=.o)
OBJ2 =  $(SRC2:.cpp=.o)
LIB = -I$(NETCDF_ROOT)/include  -I$(SRCDIR)

include $(KOKKOS_PATH)/Makefile.kokkos

.PHONY: links library test

default: all

all: links library $(EXE_NAME1) $(EXE_NAME2)

build1: $(SRC1)
	python ../../compare_to_gold.py $(EXE_NAME1)

CanopyHydrology_kern1_multiple: test_CanopyHydrology_kern1_multiple
	./test_CanopyHydrology_kern1_multiple.host > test_CanopyHydrology_kern1_multiple.stdout

build2: $(SRC2)
	python ../../compare_to_gold.py $(EXE_NAME2)

CanopyHydrology_module: test_CanopyHydrology_module
	./test_CanopyHydrology_module.host > test_CanopyHydrology_module.stdout


sandbox: test_sandbox_domain_template_magic
	./test_sandbox_domain_template_magic

$(EXE_NAME1): $(OBJ1) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(CC_LD_FLAGS) $(OBJ1) $(KOKKOS_LIBS) $(LIB) -o $(EXE_NAME1) $(LINKFLAGS) $(EXTRA_PATH)
$(EXE_NAME2): $(OBJ2) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(CC_LD_FLAGS) $(OBJ2) $(KOKKOS_LIBS) $(LIB) -o $(EXE_NAME2) $(LINKFLAGS) $(EXTRA_PATH)

clean: kokkos-clean
	rm -f *.o *.cuda *.host test_*

links:
	@echo "making in links"
	$(MAKE) -C ../../links links

library:
	$(MAKE) -C $(OBJECT) cc_serial

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) -I$(SRCDIR) $(EXTRA_PATH) $(EXTRA_INC) -c $<

test: $(EXE_NAME1)
	./$(EXE_NAME1)
