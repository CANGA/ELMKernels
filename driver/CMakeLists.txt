
if (ENABLE_CC)
  add_subdirectory(cc)
endif()

if (ENABLE_Kokkos)
  add_subdirectory(kokkos)
endif()

# include_directories(${NetCDF_INCLUDE_DIR})
# add_executable (ELM_driver driver.cc)
# target_link_libraries (ELM_driver LINK_PUBLIC physics ${NetCDF_LIBRARIES})
# install(TARGETS ELM_driver)

# add_test(NAME ELM_driver COMMAND ELM_driver)
