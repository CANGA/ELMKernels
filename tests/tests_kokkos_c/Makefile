KOKKOS_PATH = ${HOME}/Downloads/kokkos
KOKKOS_DEBUG = "yes"
EXE_NAME1 = "test_CanopyHydrology_kern1_single"
EXE_NAME2 = "test_CanopyHydrology_kern1_multiple"
EXE_NAME3 = "test_CanopyHydrology_module"
OBJECT        = ../../src/
KERNEL_FOLDER = cc
SRCDIR	      =	$(OBJECT)$(KERNEL_FOLDER)
KOKKOS_DEVICES = $1
myDEVICE     = Cuda


# ifeq ($(KOKKOS_DEVICES),$(myDEVICE))
#     NATURE = KOKKOS_INLINE_FUNCTION
# else
#     NATURE = 
# endif


include $(OBJECT)config/Makefile.config

TESTS = test_CanopyHydrology_kern1_single \
        test_CanopyHydrology_kern1_multiple \
        test_CanopyHydrology_module

SRC1 = CanopyHydrology_kern1_single.cpp        
SRC2 = CanopyHydrology_kern1_multiple.cpp
SRC3 = CanopyHydrology_module.cpp      

serial: build1 build2 build3 test
	KOKKOS_DEVICES = "Serial"
	echo "Start Build"
openmp: build1 build2 build3 test
	KOKKOS_DEVICES = "OpenMP"
	echo "Start Build"
cuda: build1 build2 build3 test
	KOKKOS_DEVICES = "Cuda"
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE1 = ${EXE_NAME1}.cuda
EXE2 = ${EXE_NAME2}.cuda
EXE3 = ${EXE_NAME3}.cuda
KOKKOS_ARCH = "BSW,Pascal60"
KOKKOS_CUDA_OPTIONS = "enable_lambda"
else
CXX = g++
EXE1 = ${EXE_NAME1}.host
EXE2 = ${EXE_NAME2}.host
EXE3 = ${EXE_NAME3}.host
KOKKOS_ARCH = "BSW"
endif

CXXFLAGS = -g -O0
LINK = ${CXX}
LINKFLAGS = -lnetcdf
EXTRA_PATH = -I$(NETCDF_ROOT)/include

DEPFLAGS = -M
OBJ1 =  $(SRC1:.cpp=.o)
OBJ2 =  $(SRC2:.cpp=.o)
OBJ3 =  $(SRC3:.cpp=.o)	
EXTRA_INC = -I$(NETCDF_ROOT)/include  -I$(SRCDIR)
LIB = -I$(NETCDF_ROOT)/include  -I$(SRCDIR)

include $(KOKKOS_PATH)/Makefile.kokkos

.PHONY: links library 

build1: links library $(EXE1)

$(EXE1): $(OBJ1) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ1) $(KOKKOS_LIBS) $(LIB) $(CXX_LD_FLAGS) -o $(EXE1)

build2: links library $(EXE2)

$(EXE2): $(OBJ2) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ2) $(KOKKOS_LIBS) $(LIB) $(CXX_LD_FLAGS) -o $(EXE2)

build3: links library $(EXE3)

$(EXE3): $(OBJ3) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ3) $(KOKKOS_LIBS) $(LIB) $(CXX_LD_FLAGS) -o $(EXE3)

clean:
	@$(ELM_CLEAN)
	rm -f test_* *.o *.cuda *.host

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

test : $(SRC1) $(SRC2) $(SRC3)
	./$(EXE1) > test_CanopyHydrology_kern1_single.stdout
	./$(EXE2) > test_CanopyHydrology_kern1_multiple.stdout
	./$(EXE3) > test_CanopyHydrology_module.stdout
	python ../compare_to_gold.py $(TESTS)

links:
	@echo "making in links"
	$(MAKE) -C ../links links

library:
	$(MAKE) -C $(OBJECT) cc