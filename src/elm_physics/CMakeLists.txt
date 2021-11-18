project(ELM_PHYSICS)

add_library (elm_physics FrictionVelocity.cc QSat.cc BareGroundFluxes.cc CanopyHydrology.cc CanopyTemperature.cc InitTopography.cc Photosynthesis.cc
SurfaceAlbedo.cc CanopyHydrology_newphysics.cc CanHydro_newphys_internal_logic.cc CanHydro_newphys_external_logic.cc)
target_include_directories (elm_physics PUBLIC ${ELM_PHYSICS_SOURCE_DIR} ${ELM_UTILS_SOURCE_DIR})

target_link_libraries (elm_physics LINK_PUBLIC elm_utils)
install(TARGETS elm_physics)
