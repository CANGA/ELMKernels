include_directories(${ELM_UTILS_SOURCE_DIR})
include_directories(${KOKKOS_INCLUDE_DIR})

add_executable (ELM_kokkos driver.cc)

target_link_libraries (ELM_kokkos LINK_PUBLIC elm_physics Kokkos::kokkos )
install(TARGETS ELM_kokkos)

add_test(NAME ELM_kokkos COMMAND ELM_kokkos)
