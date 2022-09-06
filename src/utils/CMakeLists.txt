project(ELM_UTILS)

set(utils_src_files
read_input.cc
read_test_input.cc
utils.cc
)

set(utils_inc_files
array.hh
compile_options.hh
data_types.hh
date_time.hh
helper_functions.hh
invoke_kernel.hh
kokkos_utils.hh
legion_utils.hh
min_max_sum.hh
mpi_types.hh
read_input.hh
read_netcdf.hh
read_pnetcdf.hh
read_test_input.hh
utils.hh
)

add_library (elm_utils
${utils_src_files}
${utils_inc_files}
)

if (ENABLE_KOKKOS)
target_include_directories (elm_utils PUBLIC ${ELM_UTILS_SOURCE_DIR} ${NetCDF_INCLUDE_DIR} ${KOKKOS_INCLUDE_DIR})
target_link_libraries (elm_utils LINK_PUBLIC ${NetCDF_C_LIBRARIES} Kokkos::kokkos)
else()
target_include_directories (elm_utils PUBLIC ${ELM_UTILS_SOURCE_DIR} ${NetCDF_INCLUDE_DIR})
target_link_libraries (elm_utils LINK_PUBLIC ${NetCDF_C_LIBRARIES})
endif()

install(TARGETS elm_utils)
