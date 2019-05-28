SRCDIR	      =	.
OBJECT        = .

.PHONY: cc_fortran_wrappers cc_serial cc_kokkos fortran

default: all

all: cc_fortran_wrappers fortran cc_serial cc_kokkos

cc_fortran_wrappers: fortran
	$(MAKE) -C cc_fortran_wrappers all

cc_serial:
	$(MAKE) -C cc_serial all

cc_kokkos:
	$(MAKE) -C cc_kokkos all

fortran:
	$(MAKE) -C fortran all


clean:
	$(MAKE) -C cc clean
	$(MAKE) -C fortran clean

allclean:
	$(MAKE) -C cc_serial clean
	$(MAKE) -C cc_kokkos clean
	$(MAKE) -C fortran clean



