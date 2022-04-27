
add_executable (whole_system whole_system.cc)

target_include_directories (whole_system PUBLIC ${ELM_PHYSICS_SOURCE_DIR} ${ELM_UTILS_SOURCE_DIR} ${NetCDF_INCLUDE_DIR})
target_link_libraries (whole_system LINK_PUBLIC elm_physics elm_utils ${NetCDF_LIBRARIES})

install(TARGETS whole_system)
