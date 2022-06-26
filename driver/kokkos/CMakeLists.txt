
add_executable (ELM_kokkos new_kokkos_driver.cc kokkos_elm_initialize.cc)

target_include_directories (ELM_kokkos PUBLIC ${ELM_PHYSICS_SOURCE_DIR} ${ELM_UTILS_SOURCE_DIR} ${KOKKOS_INCLUDE_DIR} ${NetCDF_INCLUDE_DIR} )
target_link_libraries (ELM_kokkos LINK_PUBLIC elm_physics elm_utils Kokkos::kokkos ${NetCDF_LIBRARIES})

install(TARGETS ELM_kokkos)


