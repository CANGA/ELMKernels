
if (ENABLE_CC)
  #add_subdirectory(cc)
endif()

if (ENABLE_KOKKOS)
  add_subdirectory(kokkos)
endif()

