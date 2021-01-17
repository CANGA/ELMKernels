add_library (physics physics_dummy.cc CanopyHydrology.cc coldstart_topography.cc coldstart_snow.cc SurfaceRadiation.cc SurfaceResistance.cc QSat.cc CanopyTemperature.cc FrictionVelocity.cc BareGroundFluxes.cc CanopyFluxes.cc Photosynthesis.cc)
target_include_directories (physics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
