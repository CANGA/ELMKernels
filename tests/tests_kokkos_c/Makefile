KOKKOS_PATH = ${HOME}/Downloads/kokkos
KOKKOS_DEVICES = "Cuda"
EXE_NAME = "CanopyHydrology_kern1_multiple"
OBJECT        = ../../src/
include $(OBJECT)config/Makefile.config

SRC = CanopyHydrology_kern1_multiple.cpp

default: build
	echo "Start Build"


ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE = ${EXE_NAME}.cuda
KOKKOS_ARCH = "HSW,Pascal60"
KOKKOS_CUDA_OPTIONS = "enable_lambda,force_uvm"
else
CXX = g++
EXE = ${EXE_NAME}.host
KOKKOS_ARCH = "HSW"
endif

CXXFLAGS = -g -O0
LINK = ${CXX}
LINKFLAGS = -lnetcdf  #-I../test_c_test
EXTRA_PATH = -I/usr/local/include

DEPFLAGS = -M

OBJ =  $(SRC:.cpp=.o)
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) ../../src/cpp/CanopyHydrology_Interception.cpp.o $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE) $(LINKFLAGS) $(EXTRA_PATH)

clean: kokkos-clean
	rm -f *.o *.cuda *.host

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) -I../../src/cpp $(EXTRA_PATH) $(EXTRA_INC) -c $<

test: $(EXE)
	./$(EXE)