include_directories(${ELM_UTILS_SOURCE_DIR})

add_executable (date_time_test date_time_test.cc)
target_link_libraries (date_time_test LINK_PUBLIC elm_physics elm_utils)

add_executable (ELM_cc driver.cc)
target_link_libraries (ELM_cc LINK_PUBLIC elm_physics elm_utils)

add_executable (whole_system whole_system.cc)
target_link_libraries (whole_system LINK_PUBLIC elm_physics elm_utils)


install(TARGETS ELM_cc date_time_test whole_system)
