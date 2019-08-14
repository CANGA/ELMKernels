Installing ELM Kernels
===========================

First, this is an extremely hacky, but relatively simple Make-based build system.

All builds occur in-source.

First, edit the file: ```src/config/Makefile.config```

Then, various program model mini-apps can be built, e.g.:

```
cd tests_c_c
make all test
```


