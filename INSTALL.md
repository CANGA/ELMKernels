Installing ELM Kernels
===========================

First, this is an extremely hacky, but relatively simple Make-based build system.

All builds occur in-source.

First, edit the file: ```src/config/Makefile.config```

Then, making the library should simply require a:

```
cd src
make
```

Or alternatively build a C++ or Fortran test via:

```
cd tests_c
make all test
```

or

```
cd tests
make all test
```

respectively.

