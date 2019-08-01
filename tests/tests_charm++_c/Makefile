-include /home/7hp/Downloads/charm/examples/common.mk
CHARMC=/home/7hp/Downloads/charm/bin/charmc $(OPTS)
ELM_BASE = ../../src
KERNEL_FOLDER = cc
ELM_ROOT = $(ELM_BASE)/$(KERNEL_FOLDER)

include $(ELM_BASE)/config/Makefile.config
include $(ELM_BASE)/config/Makefile.rules

TESTS = test_CanopyHydrology_kern1_single \
         test_CanopyHydrology_kern1_multiple \
         test_CanopyHydrology_module

EXEC_TESTS = CanopyHydrology_kern1_single \
             CanopyHydrology_kern1_multiple \
             CanopyHydrology_module


.PHONY: links library test

default: all

all: links library $(TESTS)

test: $(EXEC_TESTS)
	python ../compare_to_gold.py $(TESTS)


CanopyHydrology_kern1_single: test_CanopyHydrology_kern1_single
	./test_CanopyHydrology_kern1_single > test_CanopyHydrology_kern1_single.stdout

CanopyHydrology_kern1_multiple: test_CanopyHydrology_kern1_multiple
	./test_CanopyHydrology_kern1_multiple > test_CanopyHydrology_kern1_multiple.stdout

CanopyHydrology_module: test_CanopyHydrology_module
	./test_CanopyHydrology_module > test_CanopyHydrology_module.stdout

OBJ1 = CanopyHydrology_kern1_single.o 
OBJ2 = CanopyHydrology_kern1_multiple.o
OBJ3 = CanopyHydrology_module.o

all: links library test_CanopyHydrology_kern1_single test_CanopyHydrology_kern1_multiple test_CanopyHydrology_module

test_CanopyHydrology_kern1_single: $(OBJ1)
	$(CHARMC) $(NETCDF_FLAGS) -language charm++ -o test_CanopyHydrology_kern1_single $(OBJ1) 

CanopyHydrology_kern1_single.decl.h: CanopyHydrology_kern1_single.ci
	$(CHARMC)  CanopyHydrology_kern1_single.ci 

CanopyHydrology_kern1_single.o: CanopyHydrology_kern1_single.C CanopyHydrology_kern1_single.decl.h readers.hh utils.hh
	$(CHARMC) -c CanopyHydrology_kern1_single.C $@ $< $(NETCDF_FLAGS) -I$(ELM_ROOT) -lelm -L$(NETCDF_ROOT) -lnetcdf -L$(STD_LIB_ROOT)/lib -lstdc++ -L$(FC_LIB_ROOT) -lgfortran

test_CanopyHydrology_kern1_multiple: $(OBJ2)
	$(CHARMC) $(NETCDF_FLAGS) -language charm++ -o test_CanopyHydrology_kern1_multiple $(OBJ2) 

CanopyHydrology_kern1_multiple.decl.h: CanopyHydrology_kern1_multiple.ci
	$(CHARMC)  CanopyHydrology_kern1_multiple.ci 

CanopyHydrology_kern1_multiple.o: CanopyHydrology_kern1_multiple.C CanopyHydrology_kern1_multiple.decl.h readers.hh utils.hh
	$(CHARMC) -c CanopyHydrology_kern1_multiple.C $@ $< $(NETCDF_FLAGS) -I$(ELM_ROOT) -lelm -L$(NETCDF_ROOT) -lnetcdf -L$(STD_LIB_ROOT)/lib -lstdc++ -L$(FC_LIB_ROOT) -lgfortran

test_CanopyHydrology_module: $(OBJ3)
	$(CHARMC) $(NETCDF_FLAGS) -language charm++ -o test_CanopyHydrology_module $(OBJ3) 

CanopyHydrology_module.decl.h: CanopyHydrology_module.ci
	$(CHARMC)  CanopyHydrology_module.ci 

CanopyHydrology_module.o: CanopyHydrology_module.C CanopyHydrology_module.decl.h readers.hh utils.hh
	$(CHARMC) -c CanopyHydrology_module.C $@ $< $(NETCDF_FLAGS) -I$(ELM_ROOT) -lelm -L$(NETCDF_ROOT) -lnetcdf -L$(STD_LIB_ROOT)/lib -lstdc++ -L$(FC_LIB_ROOT) -lgfortran


clean:
	$(RM) test_* *.o 
	rm -f conv-host pgm *.def.h *.decl.h *.o *~ charmrun charmrun.exe pgm.exe pgm.pdb pgm.ilk

allclean: clean

links:
	@echo "making in links"
	$(MAKE) -C ../links links

