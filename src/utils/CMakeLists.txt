project(ELM_UTILS)

add_library (elm_utils read_test_input.cc read_input.cc utils.cc)
target_include_directories (elm_utils PUBLIC ${ELM_UTILS_SOURCE_DIR}  ${NetCDF_INCLUDE_DIR})

target_link_libraries (elm_utils LINK_PUBLIC ${NetCDF_C_LIBRARIES})

install(TARGETS elm_utils)
