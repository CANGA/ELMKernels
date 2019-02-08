.PHONY: links library CanopyHydrology_kern1_single

SRCDIR	      =	../src/fortran/
CSRCDIR	      =	../src/cpp/
OBJECT        = ../src/

default: all

include $(OBJECT)config/Makefile.config

CPPFLAGS += -I$(CSRCDIR) -I$(NETCDF_ROOT)/include

all: links library CanopyHydrology_kern1_single

CanopyHydrology_kern1_single: test_CanopyHydrology_kern1_single

test_%: %.o
	$(LINKER) $(FCFLAGS) -L$(OBJECT)fortran -lelm -L$(NETCDF_ROOT)/lib -lnetcdf -L$(STD_LIB_ROOT)/lib -lc++ $< -o $@

test: test_CanopyHydrology_kern1_single
	./test_CanopyHydrology_kern1_single


clean:
	@$(ELM_CLEAN)
	$(RM) test_*

allclean:
	@$(ELM_CLEAN)
	$(RM) test_*
	$(MAKE) -C $(OBJECT) allclean

links:
	@echo "making in links"
	$(MAKE) -C links links

library:
	$(MAKE) -C $(OBJECT) all
