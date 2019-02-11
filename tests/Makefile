SRCDIR	      =	../src/fortran/
OBJECT        = ../src/

include $(OBJECT)config/Makefile.config
FCFLAGS += -I$(SRCDIR) -I$(SRCDIR)$(ELM_UTILS_DIR) -I$(NETCDF_ROOT)/include


TESTS = test_CanopyHydrology_kern1_single \
        test_CanopyHydrology_kern1_multiple

EXEC_TESTS = CanopyHydrology_kern1_single \
             CanopyHydrology_kern1_multiple

.PHONY: links library test

default: all

all: links library $(TESTS)

test: $(EXEC_TESTS)
	echo "done"


CanopyHydrology_kern1_single: test_CanopyHydrology_kern1_single
	./test_CanopyHydrology_kern1_single

CanopyHydrology_kern1_multiple: test_CanopyHydrology_kern1_multiple
	./test_CanopyHydrology_kern1_multiple



test_%: %.o
	$(LINKER) $(FCFLAGS) -L$(OBJECT)fortran -lelm -L$(NETCDF_ROOT)/lib -lnetcdff $< -o $@


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
