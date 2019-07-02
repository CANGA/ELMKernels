-include /home/7hp/Downloads/charm/common.mk
CHARMC=/home/7hp/Downloads/charm/bin/charmc $(OPTS)
ELM_BASE = ../../src
KERNEL_FOLDER = cc
ELM_ROOT = $(ELM_BASE)/$(KERNEL_FOLDER)

include $(ELM_BASE)/config/Makefile.config

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


test_%: %.o
	$(CHARMC) %.o -o $@ $< $(NETCDF_FLAGS) -language charm++

%.o : %.C readers.hh utils.hh %.ci
	$(CHARMC) -c $@ $< $(NETCDF_FLAGS)


clean:
	$(RM) test_* *.o 
	rm -f conv-host pgm *.def.h *.decl.h *.o *~ charmrun charmrun.exe pgm.exe pgm.pdb pgm.ilk

allclean: clean

links:
	@echo "making in links"
	$(MAKE) -C ../links links

