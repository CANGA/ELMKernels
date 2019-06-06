KOKKOS_PATH = ${HOME}/Downloads/kokkos
KOKKOS_DEVICES = "Cuda"
KOKKOS_DEBUG = "yes"
EXE_NAME1 = "test_CanopyHydrology_kern1_multiple"
EXE_NAME2 = "test_CanopyHydrology_module"
OBJECT        = ../../src/
KERNEL_FOLDER = cc_serial
SRCDIR	      =	$(OBJECT)$(KERNEL_FOLDER)

include $(OBJECT)config/Makefile.config

TESTS = test_CanopyHydrology_kern1_multiple \
        test_CanopyHydrology_module
SRC1 = CanopyHydrology_kern1_multiple.cpp
SRC2 = CanopyHydrology_module.cpp      

default: build1 build2 test
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE1 = ${EXE_NAME1}.cuda
EXE2 = ${EXE_NAME2}.cuda
KOKKOS_ARCH = "BSW,Pascal60"
KOKKOS_CUDA_OPTIONS = "enable_lambda,force_uvm"
else
CXX = g++
EXE1 = ${EXE_NAME1}.host
EXE2 = ${EXE_NAME2}.host
KOKKOS_ARCH = "BSW"
endif

CXXFLAGS = -g -O0
LINK = ${CXX}
LINKFLAGS = -lnetcdf
EXTRA_PATH = -I$(NETCDF_ROOT)/include

DEPFLAGS = -M
OBJ1 =  $(SRC1:.cpp=.o)
OBJ2 =  $(SRC2:.cpp=.o)
EXTRA_INC = -I$(NETCDF_ROOT)/include
LIB = -I$(NETCDF_ROOT)/include

include $(KOKKOS_PATH)/Makefile.kokkos

.PHONY: links

build1: links $(EXE1)

$(EXE1): $(OBJ1) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ1) $(KOKKOS_LIBS) $(LIB) $(CXX_LD_FLAGS) -o $(EXE1)

build2: links $(EXE2)

$(EXE2): $(OBJ2) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ2) $(KOKKOS_LIBS) $(LIB) $(CXX_LD_FLAGS) -o $(EXE2)
clean:
	@$(ELM_CLEAN)
	rm -f test_* *.o *.cuda *.host

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

test : $(SRC1) $(SRC2)
	./$(EXE1) > test_CanopyHydrology_kern1_multiple.stdout
	./$(EXE2) > test_CanopyHydrology_module.stdout
	python ../compare_to_gold.py $(TESTS)

links:
	@echo "making in links"
	$(MAKE) -C ../links links