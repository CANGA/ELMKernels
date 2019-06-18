SRCDIR	      =	.
OBJECT        = .

.PHONY: cc_fortran_wrappers cc fortran

default: all

all: cc_fortran_wrappers fortran cc

cc_fortran_wrappers: fortran
	$(MAKE) -C cc_fortran_wrappers all

cc:
	echo "header only"

fortran:
	$(MAKE) -C fortran all

clean:
	$(MAKE) -C fortran clean

allclean:
	$(MAKE) -C fortran clean



