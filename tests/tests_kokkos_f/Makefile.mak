CSRCDIR	      =	../../src/cpp/
KOKKOS_PATH = ${HOME}/Downloads/kokkos
KOKKOS_DEVICES = "Cuda,OpenMP"
OBJECT        = ../../src/
include $(OBJECT)config/Makefile.config

TESTS = test_CanopyHydrology_kern1_single \
        test_CanopyHydrology_kern1_multiple 

EXEC_TESTS = CanopyHydrology_kern1_single \
             CanopyHydrology_kern1_multiple

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE = ${TESTS }.cuda
KOKKOS_ARCH = "HSW,Pascal60"
KOKKOS_CUDA_OPTIONS = "enable_lambda,force_uvm"
else
CXX = g++
EXE = ${TESTS }.host
KOKKOS_ARCH = "HSW"
endif

CXXFLAGS = -g -O0
LINK = ${CXX}
LINKFLAGS = -lnetcdf -I../../src/cpp -I../tests_c
EXTRA_PATH = -I/usr/local/include

DEPFLAGS = -M

OBJ = $(SRC:.cpp=.o)
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos




.PHONY: links library test

default: all

all: links library $(TESTS)

test: $(EXEC_TESTS)
	python ../compare_to_gold.py $(TESTS)


CanopyHydrology_kern1_single: test_CanopyHydrology_kern1_single
	./test_CanopyHydrology_kern1_single &> test_CanopyHydrology_kern1_single.stdout

CanopyHydrology_kern1_multiple: test_CanopyHydrology_kern1_multiple
	./test_CanopyHydrology_kern1_multiple &> test_CanopyHydrology_kern1_multiple.stdout

test_%: $(OBJ) $(KOKKOS_LINK_DEPENDS) readers.hh utils.hh library
	$(LINK) $(KOKKOS_LDFLAGS) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE) $(LINKFLAGS) $(EXTRA_PATH)

%.o : %.cpp $(KOKKOS_CPP_DEPENDS) domains.hh readers.hh utils.hh
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(EXTRA_INC) -c $<

clean:
	@$(ELM_CLEAN)
	$(RM) test_*

allclean:
	@$(ELM_CLEAN)
	$(RM) test_*
	$(MAKE) -C $(OBJECT) allclean

links:
	@echo "making in links"
	$(MAKE) -C ../links links

library:
	$(MAKE) -C $(OBJECT) all
