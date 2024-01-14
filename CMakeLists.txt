cmake_minimum_required(VERSION 3.20)

# Include overwrites before setting up the project
set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_CURRENT_SOURCE_DIR}/config/DefaultFlags.cmake)

get_directory_property(is-subproject PARENT_DIRECTORY)

project(
  "libsparse"
  LANGUAGES "Fortran"
  VERSION "0.1.0"
  DESCRIPTION "Sparse matrices"
)

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# General configuration information
add_subdirectory("config")

if(NOT TARGET "OpenMP::OpenMP_Fortran")
  find_package("OpenMP" REQUIRED)
endif()

#MKL
message(STATUS "ENV{MKLROOT}: $ENV{MKLROOT}")

#set(MKL_INTERFACE "lp64")
#set(ENABLE_BLAS95 true)
#set(ENABLE_LAPACK95 true)
#set(MKL_LINK "static")
#find_package(MKL CONFIG REQUIRED)

set(BLA_STATIC ON)
if(NOT TARGET "BLAS::BLAS")
  find_package("custom-blas" REQUIRED)
endif()

if(NOT TARGET "LAPACK::LAPACK")
  find_package("custom-lapack" REQUIRED)
endif()

set(
  lib-deps
  "OpenMP::OpenMP_Fortran"
  "BLAS::BLAS"
  "LAPACK::LAPACK"
)


# Collect source of the project
set(srcs)
add_subdirectory("src")

# We need the module directory before we finish the configure stage
if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
  make_directory("${PROJECT_BINARY_DIR}/include")
endif()

# Testing library target
add_library(
  "${PROJECT_NAME}-lib"
  "${srcs}"
)
set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
)
target_link_libraries(
  "${PROJECT_NAME}-lib"
  PUBLIC
  "${lib-deps}"
)
target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)
target_compile_definitions(
  "${PROJECT_NAME}-lib"
  PRIVATE
  "-D_DP=1"
  "-D_METIS=0"
  "-D_PARDISO=1"
  "-D_SPAINV=1"
  "-D_VERBOSE=0"
)

# Export targets for other projects
add_library("${PROJECT_NAME}" INTERFACE)
target_link_libraries("${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}-lib")
install(
  TARGETS
  "${PROJECT_NAME}"
  "${PROJECT_NAME}-lib"
  EXPORT
  "${PROJECT_NAME}-targets"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY
  "${PROJECT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)
# Package license files
install(
  FILES
  "LICENSE"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)


# add the testsuite
include(CTest)
if(BUILD_TESTING AND LIBSPARSE_BUILD_TESTING)
  add_subdirectory("test")
endif()
