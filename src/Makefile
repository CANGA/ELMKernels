SRCDIR	      =	.
OBJECT        = .

.PHONY: cc_fortran_wrappers cc fortran

default: all

all: cc_fortran_wrappers fortran cc

cc_fortran_wrappers: fortran
	$(MAKE) -C cc_fortran_wrappers all

cc:
	$(MAKE) -C cc all

fortran:
	$(MAKE) -C fortran all

clean:
	$(MAKE) -C cc clean
	$(MAKE) -C fortran clean

allclean:
	$(MAKE) -C cc clean
	$(MAKE) -C fortran clean



