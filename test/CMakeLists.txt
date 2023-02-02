include_directories (${ELM_PHYSICS_SOURCE_DIR} ${ELM_DATA_SOURCE_DIR} ${ELM_UTILS_SOURCE_DIR})

add_compile_definitions(TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/")

add_executable (test_CanHydro test_CanHydro.cc)
target_link_libraries (test_CanHydro LINK_PUBLIC elm_physics elm_data elm_utils)
install(TARGETS test_CanHydro)

add_executable (test_CanSunShade test_CanSunShade.cc)
target_link_libraries (test_CanSunShade LINK_PUBLIC elm_physics elm_data elm_utils)
install(TARGETS test_CanSunShade)

add_executable (test_SurfRad test_SurfRad.cc)
target_link_libraries (test_SurfRad LINK_PUBLIC elm_physics elm_data elm_utils)
install(TARGETS test_SurfRad)

add_executable (test_CanTemp test_CanTemp.cc)
target_link_libraries (test_CanTemp LINK_PUBLIC elm_physics elm_data elm_utils)
install(TARGETS test_CanTemp)

add_executable (test_BGFlux test_BGFlux.cc)
target_link_libraries (test_BGFlux LINK_PUBLIC elm_physics elm_data elm_utils)
install(TARGETS test_BGFlux)

add_executable (test_CanFlux test_CanFlux.cc)
target_link_libraries (test_CanFlux LINK_PUBLIC elm_physics elm_data elm_utils)
install(TARGETS test_CanFlux)

add_executable (test_SurfAlb test_SurfAlb.cc)
target_link_libraries (test_SurfAlb LINK_PUBLIC elm_physics elm_data elm_utils)
install(TARGETS test_SurfAlb)

#add_executable (test_SurfAlb_input test_SurfAlb_input.cc)
#target_link_libraries (test_SurfAlb_input LINK_PUBLIC elm_physics elm_data elm_utils)
#install(TARGETS test_SurfAlb_input)

#add_executable (test_CanHydro_newphysics test_CanHydro_newphysics.cc)
#target_link_libraries (test_CanHydro_newphysics LINK_PUBLIC elm_physics elm_data elm_utils)
#install(TARGETS test_CanHydro_newphysics)



