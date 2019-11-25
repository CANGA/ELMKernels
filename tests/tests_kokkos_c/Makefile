ELM_BASE = ../../src
KERNEL_FOLDER = cc
ELM_ROOT = $(ELM_BASE)/$(KERNEL_FOLDER)

KOKKOS_DEBUG = "yes"
EXE_NAME1 = "test_CanopyHydrology_kern1_single"
EXE_NAME2 = "test_CanopyHydrology_kern1_multiple"
EXE_NAME3 = "test_CanopyHydrology_module"
EXE_NAME4 = "test_CanopyHydrology_module_multiple"
KOKKOS_DEVICES = ${myDEVICE}

# Get the Kokkos Device
ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))

CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE1 = ${EXE_NAME1}.cuda
EXE2 = ${EXE_NAME2}.cuda
EXE3 = ${EXE_NAME3}.cuda
EXE4 = ${EXE_NAME4}.cuda
KOKKOS_ARCH = "BSW,Pascal60"
KOKKOS_CUDA_OPTIONS = "enable_lambda"

else

EXE1 = ${EXE_NAME1}.host
EXE2 = ${EXE_NAME2}.host
EXE3 = ${EXE_NAME3}.host
EXE4 = ${EXE_NAME4}.host
KOKKOS_ARCH = "BSW"

endif

include $(ELM_BASE)/config/Makefile.config
include $(ELM_BASE)/config/Makefile.rules

# FC = mpifort
# CXX = mpicxx
# CC = mpicc
CXXFLAGS += -DNATURE=KOKKOS_INLINE_FUNCTION 

default: Serial
# include Kokkos
include $(KOKKOS_PATH)/Makefile.kokkos

# build the tests
TESTS = test_CanopyHydrology_kern1_single \
        test_CanopyHydrology_kern1_multiple \
        test_CanopyHydrology_module

TESTS1 = test_CanopyHydrology_kern1_multiple \
         test_CanopyHydrology_module

SRC1 = CanopyHydrology_kern1_single.cpp        
SRC2 = CanopyHydrology_kern1_multiple.cpp
SRC3 = CanopyHydrology_module.cpp     
SRC4 = CanopyHydrology_module_multiple.cpp   

LINK = ${CXX}

OBJ1 =  $(SRC1:.cpp=.cpp.o)
OBJ2 =  $(SRC2:.cpp=.cpp.o)
OBJ3 =  $(SRC3:.cpp=.cpp.o)	
OBJ4 =  $(SRC4:.cpp=.cpp.o)

include $(ELM_BASE)/config/Makefile.logo

.PHONY: logo links Serial OpenMP Cuda

Serial: logo build1 build2 build3 build4 test
	echo "Serial"
	echo "Start Build"
OpenMP: logo build1 build2 build3 build4 test
	echo "OpenMP"
	echo "Start Build"
Cuda: logo build2 build3 build4 test1
	echo "Cuda"
	echo "Start Build"

Report: 
	@echo "LINKING"
	@echo "----------------"
	@echo "LINK = $(LINK)"
	@echo "KOKKOS_LDFLAGS = $(KOKKOS_LDFLAGS)"
	@echo "OBJ1 = $(OBJ1)"
	@echo "NETCDF_FLAGS = $(NETCDF_FLAGS)"
	@echo "KOKKOS_LIBS = $(KOKKOS_LIBS)"
	@echo "EXE = $(EXE1)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "CXX_LD_FLAGS = $(CXX_LD_FLAGS)"
	@echo "KOKKOS_CXXFLAGS = $(KOKKOS_CXXFLAGS)"
	@echo "KOKKOS_CPPFLAGS = $(KOKKOS_CPPFLAGS)"

build1: links $(EXE1)

$(EXE1): $(KOKKOS_CPP_DEPENDS) $(OBJ1) $(KOKKOS_LINK_DEPENDS)
	# @echo "LINKING"
	# @echo "----------------"
	# @echo "LINK = $(LINK)"
	# @echo "KOKKOS_LDFLAGS = $(KOKKOS_LDFLAGS)"
	# @echo "OBJ1 = $(OBJ1)"
	# @echo "NETCDF_FLAGS = $(NETCDF_FLAGS)"
	# @echo "KOKKOS_LIBS = $(KOKKOS_LIBS)"
	# @echo "EXE = $(EXE1)"
	# @echo "CXXFLAGS = $(CXXFLAGS)"
	# @echo "KOKKOS_CXXFLAGS = $(KOKKOS_CXXFLAGS)"
	# @echo "KOKKOS_CPPFLAGS = $(KOKKOS_CPPFLAGS)"
	$(LINK) $(OBJ1) $(NETCDF_FLAGS) $(KOKKOS_LDFLAGS) $(KOKKOS_LIBS) $(CXX_LD_FLAGS) -o $(EXE1)

build2: links Report $(EXE2)

$(EXE2): $(OBJ2) $(KOKKOS_LINK_DEPENDS)
	@echo "LINKING"
	@echo "----------------"
	@echo "LINK = $(LINK)"
	@echo "KOKKOS_LDFLAGS = $(KOKKOS_LDFLAGS)"
	@echo "OBJ1 = $(OBJ1)"
	@echo "NETCDF_FLAGS = $(NETCDF_FLAGS)"
	@echo "KOKKOS_LIBS = $(KOKKOS_LIBS)"
	@echo "EXE = $(EXE1)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "KOKKOS_CXXFLAGS = $(KOKKOS_CXXFLAGS)"
	@echo "KOKKOS_CPPFLAGS = $(KOKKOS_CPPFLAGS)"
	$(LINK) $(OBJ2) $(NETCDF_FLAGS) $(KOKKOS_LDFLAGS) $(KOKKOS_LIBS) $(CXX_LD_FLAGS) -o $(EXE2)

build3: links $(EXE3)

$(EXE3): $(OBJ3) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(OBJ3) $(NETCDF_FLAGS) $(KOKKOS_LDFLAGS) $(KOKKOS_LIBS) $(CXX_LD_FLAGS) -o $(EXE3)

build4: links $(EXE4)

$(EXE4): $(OBJ4) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(OBJ4) $(NETCDF_FLAGS) $(KOKKOS_LDFLAGS) $(KOKKOS_LIBS) $(CXX_LD_FLAGS) -o $(EXE4)

clean:
	@$(ELM_CLEAN)
	rm -f test_* *.o *.cuda *.host

test : $(SRC1) $(SRC2) $(SRC3) $(SRC4)
	./$(EXE_NAME1).host > test_CanopyHydrology_kern1_single.stdout
	./$(EXE_NAME2).host > test_CanopyHydrology_kern1_multiple.stdout
	./$(EXE_NAME3).host > test_CanopyHydrology_module.stdout
	./$(EXE_NAME4).host > test_CanopyHydrology_module_multiple.stdout
	python ../compare_to_gold.py $(TESTS)

test1 : $(SRC2) $(SRC3) $(SRC4)
	./${EXE_NAME2}.cuda > test_CanopyHydrology_kern1_multiple.stdout
	./${EXE_NAME3}.cuda > test_CanopyHydrology_module.stdout
	./${EXE_NAME4}.cuda > test_CanopyHydrology_module_multiple.stdout
	python ../compare_to_gold.py $(TESTS1)

links:
	@echo "making in links"
	$(MAKE) -C ../links links


allclean: clean
	$(MAKE) -C $(SRCDIR) clean
