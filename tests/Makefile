.PHONY: links library

SRCDIR	      =	../src/
OBJECT        = ../src/

default: all

include $(OBJECT)config/Makefile.config
FCFLAGS += -I$(SRCDIR) -I$(SRCDIR)$(ELM_UTILS_DIR) -I$(NETCDF_ROOT)/include

all: links library CanopyHydrology_kern1_single

CanopyHydrology_kern1_single: test_CanopyHydrology_kern1_single
#	./test_$@
	echo 'fail'

test_%: %.o
	$(LINKER) $(FCFLAGS) -L$(OBJECT) -lelm -L$(NETCDF_ROOT)/lib -lnetcdff $< -o $@

clean:
	@$(ELM_CLEAN)
	$(RM) test_*

links:
	@echo "making in links"
	$(MAKE) -C links links

library:
	$(MAKE) -C ../src all
