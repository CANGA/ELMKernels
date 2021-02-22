include_directories(${ELM_UTILS_SOURCE_DIR})
add_executable (ELM_cc driver.cc)

target_link_libraries (ELM_cc LINK_PUBLIC elm_physics elm_utils)
install(TARGETS ELM_cc)

add_test(NAME ELM_cc COMMAND ELM_cc)
