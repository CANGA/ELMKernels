include_directories (${ELM_PHYSICS_SOURCE_DIR} ${ELM_UTILS_SOURCE_DIR})

add_executable (test_CanHydro test_CanHydro.cc)
target_link_libraries (test_CanHydro LINK_PUBLIC elm_physics elm_utils)
install(TARGETS test_CanHydro)

add_executable (test_CanSunShade test_CanSunShade.cc)
target_link_libraries (test_CanSunShade LINK_PUBLIC elm_physics elm_utils)
install(TARGETS test_CanSunShade)

add_executable (test_SurfRad test_SurfRad.cc)
target_link_libraries (test_SurfRad LINK_PUBLIC elm_physics elm_utils)
install(TARGETS test_SurfRad)
