project(ELM_PHYSICS)

add_library (elm_physics FrictionVelocity.cc QSat.cc BareGroundFluxes.cc)
target_include_directories (elm_physics PUBLIC ${ELM_PHYSICS_SOURCE_DIR})
install(TARGETS elm_physics)
