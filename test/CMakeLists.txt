include_directories (${ELM_PHYSICS_SOURCE_DIR} ${ELM_UTILS_SOURCE_DIR})

add_executable (test_CanHydro test_CanHydro.cc)
target_link_libraries (test_CanHydro LINK_PUBLIC elm_physics elm_utils)
install(TARGETS test_CanHydro)
