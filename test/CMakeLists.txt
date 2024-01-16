if (NOT TARGET "test-drive::test-drive")
  find_package("test-drive" REQUIRED)
endif()

macro(ADDTEST name)
    add_executable(test_${name} test_${name}.f90 ${test-srcs})
    target_link_libraries(test_${name} "${PROJECT_NAME}" "test-drive::test-drive"
                         )
    target_compile_definitions(
      "test_${name}"
      PRIVATE
      "-D_DP=1"
      "-D_METIS=${WITH_METIS}"
      "-D_PARDISO=1"
      "-D_SPAINV=1"
      "-D_VERBOSE=0"
    )
    add_test(NAME ${name}
             COMMAND $<TARGET_FILE:test_${name}> ${CMAKE_CURRENT_BINARY_DIR}
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro(ADDTEST)

list(
  APPEND fyppFlags
  "-I${PROJECT_SOURCE_DIR}/src"
)

set(
    test-srcs
    modtest_common.f90
    modtest_coo.f90
    modtest_crs64.f90
    modtest_crs.f90
    modtest_random.f90
)

ADDTEST(sparse)
