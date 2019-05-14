SRCDIR	      =	.
OBJECT        = .

.PHONY: cc fortran

default: all

all: cc fortran

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



