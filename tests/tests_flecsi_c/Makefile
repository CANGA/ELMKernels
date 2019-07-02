# from Makefile.config
OBJECT        = ../../src/
KERNEL_LANG	  = cc_kokkos
SRCDIR	      =	$(OBJECT)$(KERNEL_LANG)

include $(OBJECT)config/Makefile.config

TESTS = test_CanopyHydrology_kern1_multiple \
        test_CanopyHydrology_module


EXEC_TESTS = CanopyHydrology_kern1_multiple \
             CanopyHydrology_module

#------------------------------------------------------------------------------#

all:
	$(eval test := $(shell pwd | sed 's,^.*\/,,g;s,[0-9][0-9]-,,g'))
	@echo "/*------------------------------------------------------------* "
	@echo " * FleCSI Tutorial"
	@echo " *------------------------------------------------------------*/"
	@echo ""
	@echo "To build the tutorial example, execute:"
	@echo "   $$ flecsit compile $(test).cc"
	@echo ""
	@echo "You may need to load the flecsi-tutorial environment module:"
	@echo "   $$ module load flecsi-tutorial"
	@echo ""
	@echo "OR source the bash or csh configuration file:"
	@echo "   $$ source INSTALL_PATH_TO_FLECSI/bin/flecsi-tutorial.{sh,csh}"
	@echo ""
	@echo "For convenience, this Makefile can be used to generate pdf"
	@echo "files from the tutorial output by running:"
	@echo "   $$ make pdfs"
	@echo ""
	@echo "This Makefile can also be used to clean-up after you have"
	@echo "built an example by running:"
	@echo "   $$ make clean"
	@echo ""

texfiles = $(wildcard *.tex)
pdffiles = $(patsubst %.tex, %.pdf, $(texfiles))

%.pdf: %.tex
	@pdflatex $< 2>&1 > /dev/null

pdfs: $(pdffiles)

clean:
	$(eval test := $(shell pwd | sed 's,^.*\/,,g;s,[0-9][0-9]-,,g'))
	rm -rf *.tex *.aux *.log *.pdf $(test) tmp*


.PHONY: links library test

default: all
	echo "GEN OBJS = $(GEN_OBJS)"


all: links library $(TESTS)

test: $(EXEC_TESTS)
	python ../compare_to_gold.py $(TESTS)


CanopyHydrology_kern1_multiple: test_CanopyHydrology_kern1_multiple
	./test_CanopyHydrology_kern1_multiple > test_CanopyHydrology_kern1_multiple.stdout
# 	./test_CanopyHydrology_kern1_multiple -lg:spy -logfile spy_multiple%.log
# 	$(LG_RT_DIR)/../tools/legion_spy.py -dez spy_multiple*.log

CanopyHydrology_module: test_CanopyHydrology_module
	./test_CanopyHydrology_module > test_CanopyHydrology_module.stdout
	./test_CanopyHydrology_module -lg:spy -logfile spy_module%.log
	$(LG_RT_DIR)/../tools/legion_spy.py -dez spy_module*.log


sandbox: test_sandbox_domain_template_magic
	./test_sandbox_domain_template_magic

# If we're using CUDA we have to link with nvcc
test_%: %.cc.o $(GEN_OBJS) $(GEN_GPU_OBJS) $(SLIB_LEGION) $(SLIB_REALM)
	@echo "---> Linking objects into one binary: $()"
	$(CXX) -o $@ $< $(GEN_OBJS) $(GEN_GPU_OBJS) $(CC_LD_FLAGS) $(LEGION_LIBS) $(LEGION_LD_FLAGS) $(GASNET_FLAGS)

%.cc.o : %.cc domains.hh readers.hh utils.hh
	$(CXX) $(CC_FLAGS) $(CPP_FLAGS) $(INC_FLAGS) -c $< -o $@

clean :
	$(ELM_CLEAN)
	rm -f test_*

allclean : clean
	$(MAKE) -C $(OBJECT) allclean
	$(RM) -f $(SLIB_LEGION) $(SLIB_REALM) $(GEN_OBJS) $(GEN_GPU_OBJS) $(REALM_OBJS) $(LEGION_OBJS) $(GPU_RUNTIME_OBJS) $(MAPPER_OBJS) $(ASM_OBJS)

links:
	@echo "making in links"
	$(MAKE) -C ../links links


library:
	$(MAKE) -C $(OBJECT) all

