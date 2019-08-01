SRCDIR	      =	.
OBJECT        = .

.PHONY: tests_c_c tests_legion_c tests_kokkos_c tests_c_f tests_f tests_legion_f tests_cuda tests_charm++_c tests_raja_c

default: all

all: tests_c_c tests_legion_c tests_kokkos_c tests_c_f tests_f tests_legion_f tests_cuda tests_charm++_c tests_raja_c

tests_c_c:
	$(MAKE) -C tests_c_c all

tests_cuda:
	$(MAKE) -C tests_cuda all

tests_legion_c:
	$(MAKE) -C tests_legion_c all

tests_kokkos_c:
	$(MAKE) -C tests_kokkos_c all

tests_c_f:
	$(MAKE) -C tests_c_f all

tests_f:
	$(MAKE) -C tests_f all

tests_legion_f:
	$(MAKE) -C tests_legion_f all

tests_charm++_c:
	$(MAKE) -C tests_charm++_c all

tests_raja_c:
	$(MAKE) -C tests_raja_c all

clean:
	$(MAKE) -C tests_c_c clean
	$(MAKE) -C tests_cuda clean
	$(MAKE) -C tests_legion_c clean
	$(MAKE) -C tests_kokkos_c clean
	$(MAKE) -C tests_c_f clean
	$(MAKE) -C tests_f clean
	$(MAKE) -C tests_legion_f clean
	$(MAKE) -C tests_charm++_c clean
	$(MAKE) -C tests_raja_c clean

allclean:
	$(MAKE) -C tests_c_c clean
	$(MAKE) -C tests_cuda clean
	$(MAKE) -C tests_legion_c clean
	$(MAKE) -C tests_kokkos_c clean
	$(MAKE) -C tests_c_f clean
	$(MAKE) -C tests_f clean
	$(MAKE) -C tests_legion_f clean
	$(MAKE) -C tests_charm++_c clean
	$(MAKE) -C tests_raja_c clean



