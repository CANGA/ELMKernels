cmake_minimum_required(VERSION 3.10)

project(ELM_physics)

add_subdirectory (src)
add_subdirectory (driver)

install(TARGETS ELM_driver DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)
install(TARGETS physics DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib) 