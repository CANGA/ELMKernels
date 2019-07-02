
find_package(RAJA)
cmake -DRAJA_DIR=/home/7hp/Downloads/RAJA/share/raja/cmake

#raja_add_executable(
#  NAME test_CanopyHydrology_kern1_multiple
#  SOURCES CanopyHydrology_kern1_multiple.cpp)

#raja_add_executable(
#  NAME test_CanopyHydrology_module
#  SOURCES CanopyHydrology_module.cpp)
add_executable(test_CanopyHydrology_kern1_single.exe
               CanopyHydrology_kern1_single.cpp)

target_link_libraries(test_CanopyHydrology_kern1_single.exe RAJA)

target_include_directories(test_CanopyHydrology_kern1_single.exe ${RAJA_INCLUDE_DIR}

if(ENABLE_TARGET_OPENMP)
raja_add_executable(
  NAME test_CanopyHydrology_kern1_single
  SOURCES CanopyHydrology_kern1_single.cpp)
endif()