ELM Kernels
================


Kernels lifted from ELM, enucleated of their data structures, and then
re-wrapped in a variety of drivers.

To build, see INSTALL.md in the top level directory.


ELM modules covered:

CanopyHydrologyMod
--------------------

Kernels include:

* Interception: (rename me) partitions incoming precip into throughfall/precip



Work Plan
=================

Single kernel and single grid cell tests/questions:
-----------------------------------------------------------------------------

- how do write kernels such that they can be called by:         (Ethan)
  * fortran
  * c
  * Legion
  * Kokkos

- how do we wrap a kernel for Legion                    (Himanshu)
  * sequence of blocked tasks
  * sequence of non-blocked tasks (based on futures)
        



Single kernel and multiple grid cell tests/questions:   (Himanshu / Ethan)
-----------------------------------------------------------------------------

- work out how to use IndexSpace/FieldSpace             



Multiple kernel and multiple grid cell tests (module test): (Ethan / Himanshu)
------------------------------------------------------------------------------

- IndexSpace/FieldSpace + actual task graph that isn't trivial


Module level part 2:            (Dali / Ethan / Himanshu)
-----------------------------------------------------------------------------

- Implement driver with ACTUAL ELM data layout (with filter / lake PFTs / ???) (masking??)
- compare alternative data layouts
  * ragged layouts (no filter)
  * fixed layout with no filter (maybe better for threading models e.g. GPUs?)


automation...
-----------------------------------------------------------------------------

- automated kernel generation
- automated C kernel wrapper
- automated Legion wrappers?


and beyond...
-----------------------------------------------------------------------------

- integrate with ATS
