add_library (physics physics_dummy.cc CanopyHydrology.cc coldstart_topography.cc coldstart_snow.cc SurfaceRadiation.cc SurfaceResistanceMod_functions.cc QSatMod_functions.cc CanopyTemperature.cc FrictionVelocityMod_functions.cc)
target_include_directories (physics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
