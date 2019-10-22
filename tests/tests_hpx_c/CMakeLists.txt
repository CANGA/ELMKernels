cmake_minimum_required(VERSION 3.10)

set(HPX_DIR "/home/7hp/Downloads/hpx/build")
set(HPX_INCLUDE_DIR "/home/7hp/Downloads/hpx/build/include")
set(HPX_LIB_DIR "/home/7hp/Downloads/hpx/build/lib/")
set(ELM_ROOT "/home/7hp/Downloads/elm_kernels/src/cc/")
set(NETCDF "/usr/local")
set(HPX_IGNORE_COMPILER_COMPATIBILITY On)
set(HPX_CXX_COMPILER "/usr/bin/c++") 
set(CMAKE_CXX_COMPILER "/usr/bin/c++")
set(HPX_CXX_COMPILER_ID "GNU")
set(CMAKE_CXX_COMPILER_ID "GNU")
set(CMAKE_MODULE_PATH $HPX_DIR)

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
endif()

message(" ")
message(" ${Blue}EeeeeeeeeeeeeeE  LLlL            Mm                        Mm       Kk           kk                                                          LllL${ColourReset}")
message(" ${Blue}EeeeeeeeeeeeeeE  LLlL            Mm Mm                  Mm Mm       Kk         kk                                                            LllL${ColourReset}")
message(" ${Blue}Eee              LLlL            Mm   Mm              Mm   Mm       Kk       kk                                                              LllL${ColourReset}")
message(" ${Blue}Eee              LllL            Mm     Mm          Mm     Mm       Kk     kk                                                                LllL${ColourReset}")
message(" ${Blue}Eee              LllL            Mm       Mm      Mm       Mm       Kk    kk                                                                 LllL${ColourReset}")
message(" ${Blue}Eee              LllL            Mm         Mm  Mm         Mm       Kk  kk                 eeeeeeeee  RrRrRrRrRr  Nnn         Nn  eeeeeeeee  LllL${ColourReset}")
message(" ${Blue}EeeeeeeeE        LllL            Mm           Mm           Mm       Kk kk kk               eeeeeeeee  Rr      Rr  Nn nn       Nn  eeeeeeeee  LllL${ColourReset}")
message(" ${Blue}EeeeeeeeE        LllL            Mm                        Mm       Kkk     kk             eeee       Rr      Rr  Nn  nn      Nn  eeee       LllL${ColourReset}")
message(" ${Blue}Eee              LllL            Mm                        Mm       Kk        kk           eeeeeeee   RrRrRrRrRr  Nn   nn     Nn  eeeeeeee   LllL${ColourReset}")
message(" ${Blue}Eee              LllL            Mm                        Mm       Kk          kk         eeeeeeee   RrR         Nn    nn    Nn  eeeeeeee   LllL${ColourReset}")
message(" ${Blue}Eee              LllL            Mm                        Mm       Kk            kk       eeee       Rr Rr       Nn     nn   Nn  eeee       LllL${ColourReset}")
message(" ${Blue}Eee              LllL            Mm                        Mm       Kk              kk     eeee       Rr   Rr     Nn      nn  Nn  eeee       LllL${ColourReset}")
message(" ${Blue}EeeeeeeeeeeeeeE  LllllllllllLllL Mm                        Mm       Kk                kk   eeeeeeeee  Rr     Rr   Nn       nn Nn  eeeeeeeee  LlllllllLllL${ColourReset}")
message(" ${Blue}EeeeeeeeeeeeeeE  LllllllllllLllL Mm                        Mm       Kk                  kk eeeeeeeee  Rr       Rr Nn        nnNn  eeeeeeeee  LlllllllLllL${ColourReset}")
message(" ")
message(" ")
message("     ${Green} Copyright 2019, UT Battelle / Oak Ridge National Laboratory ${ColourReset}")
message(" ")
message(" ${Magenta} Collaboration ${ColourReset}")
message("     ${Cyan}Oak Ridge Leadership Computing Facility ${ColourReset}")
message("     ${Cyan}United States Department of Energy ${ColourReset}")
message("     ${Cyan}The Energy Exascale Earth System Model (E3SM) Project ${ColourReset}")
message("     ${Cyan}Coupling Approaches for Next-Generation Architectures (CANGA) - Scientific Discovery through Advanced Computing (SciDAC) ${ColourReset}")
message(" ")
message("${Yellow} Contact ${ColourReset} ${Red}coonet@ornl.gov , pillaihk@ornl.gov ${ColourReset}")

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

project(my_hpx_project CXX)
find_package(HPX REQUIRED)
include_directories(
    /home/7hp/Downloads/elm_kernels/src/cc/
    /usr/local/include
    )
link_directories(${HPX_LIBRARY_DIR})

add_compile_options(-std=c++11)

add_hpx_executable(test_CanopyHydrology_kern1_multiple
    ESSENTIAL
    SOURCES CanopyHydrology_kern1_multiple.cpp
    LINK_FLAGS "-I/usr/local/include -I/home/7hp/Downloads/elm_kernels/src/cc/"
    COMPONENT_DEPENDENCIES iostreams init
    DEPENDENCIES -g -L/usr/local/lib -lnetcdf)

add_hpx_executable(test_CanopyHydrology_kern1_single
    ESSENTIAL
    SOURCES CanopyHydrology_kern1_single.cpp
    LINK_FLAGS "-I/usr/local/include -I/home/7hp/Downloads/elm_kernels/src/cc/"
    COMPONENT_DEPENDENCIES iostreams init 
    DEPENDENCIES -L/usr/local/lib -lnetcdf)


add_hpx_executable(test_CanopyHydrology_module
  ESSENTIAL
  SOURCES CanopyHydrology_module.cpp
  LINK_FLAGS "-I/usr/local/include -I/home/7hp/Downloads/elm_kernels/src/cc/"
  COMPONENT_DEPENDENCIES iostreams init 
  DEPENDENCIES -L/usr/local/lib -lnetcdf)


#execute_process(
#  COMMAND ./test_CanopyHydrology_kern1_single
#  OUTPUT_FILE "test_CanopyHydrology_kern1_single.stdout"
#)

execute_process(
  COMMAND ./test_CanopyHydrology_kern1_multiple
  OUTPUT_FILE "test_CanopyHydrology_kern1_multiple.stdout"
)

execute_process(
  COMMAND ./test_CanopyHydrology_module
  OUTPUT_FILE "test_CanopyHydrology_module.stdout"
)

enable_testing()
add_test(NAME HPXtest
  COMMAND python ../compare_to_gold.py test_CanopyHydrology_kern1_multiple test_CanopyHydrology_module
)

add_custom_target(cleanall
  COMMAND rm -r CMakeFiles/ Testing/ 
  COMMAND rm test_* CMakeCache.txt cmake_* CTestTestfile.cmake Makefile *.soln *.stdout
)