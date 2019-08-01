cmake_minimum_required(VERSION 3.10)

set(RAJA_DIR "/home/7hp/Downloads/RAJA/build/share/raja/cmake")
set(RAJA_INCLUDE_DIR "/home/7hp/Downloads/RAJA/build/include")
set(RAJA_LIB_DIR "/home/7hp/Downloads/RAJA/build/lib/")
set(ELM_ROOT "/home/7hp/Downloads/elm_kernels/src/cc/")
set (NETCDF "/usr/local")
set(CMAKE_MODULE_PATH $RAJA_DIR)

install(FILES
    ${RAJA_INCLUDE_DIR}/RAJA/config.hpp
    ${RAJA_INCLUDE_DIR}/RAJA/module.modulemap
    ${RAJA_INCLUDE_DIR}/RAJA/module.private.modulemap
    DESTINATION "include/RAJA/")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(RAJA)
include_directories(
    /home/7hp/Downloads/elm_kernels/src/cc/
    /usr/local/include
    )
add_executable(
  test_CanopyHydrology_kern1_multiple
  CanopyHydrology_kern1_multiple.cpp)
target_link_libraries(test_CanopyHydrology_kern1_multiple -L/usr/local/lib -lnetcdf)
target_include_directories(test_CanopyHydrology_kern1_multiple PUBLIC ${RAJA_INCLUDE_DIR} $(ELM_ROOT))

add_executable(
  test_CanopyHydrology_kern1_single
  CanopyHydrology_kern1_single.cpp)
target_link_libraries(test_CanopyHydrology_kern1_single -L/usr/local/lib -lnetcdf)
target_include_directories(test_CanopyHydrology_kern1_single PUBLIC ${RAJA_INCLUDE_DIR})

#add_executable(
#  test_CanopyHydrology_module
#  CanopyHydrology_module.cpp)
#target_link_libraries(test_CanopyHydrology_module -L/usr/local/lib -lnetcdf)
#target_include_directories(test_CanopyHydrology_module PUBLIC ${RAJA_INCLUDE_DIR} $(ELM_ROOT))

execute_process(
  COMMAND ./test_CanopyHydrology_kern1_single
  OUTPUT_FILE "test_CanopyHydrology_kern1_single.stdout"
)

execute_process(
  COMMAND ./test_CanopyHydrology_kern1_multiple
  OUTPUT_FILE "test_CanopyHydrology_kern1_multiple.stdout"
)

enable_testing()
add_test(NAME RAJAtest
  COMMAND python ../compare_to_gold.py test_CanopyHydrology_kern1_single test_CanopyHydrology_kern1_multiple
)

add_custom_target(cleanall
  COMMAND rm -r CMakeFiles/ Testing/ 
  COMMAND rm test_* CMakeCache.txt cmake_* CTestTestfile.cmake Makefile
)