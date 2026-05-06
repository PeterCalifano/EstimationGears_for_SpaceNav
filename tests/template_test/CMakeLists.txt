# Exclude EXCLUDED_LIST from the list of tests
set(EXCLUDED_LIST "test_to_exclude")  # Specify files to exclude WE
if(NOT SPDLOG_ENABLED)
    list(APPEND EXCLUDED_LIST "testLogging" "testLogging.cpp")
endif()

add_tests(${project_name} EXCLUDED_LIST TESTS_LIST ${CUDA_COMPILE_TARGET} CATCH2_TEST_PROPERTIES Catch2::Catch2WithMain)

set(TESTS_LIST ${TESTS_LIST} PARENT_SCOPE) # Export list of tests to parent scope (applied ONLY here)
