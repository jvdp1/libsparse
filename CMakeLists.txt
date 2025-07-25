# Based on test-drive
# License test-drive:
# This file is part of test-drive.
# SPDX-Identifier: Apache-2.0 OR MIT
#
# Licensed under either of Apache License, Version 2.0 or MIT license
# at your option; you may not use this file except in compliance with
# the License.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.27.0)

# Include overwrites before setting up the project
set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_CURRENT_SOURCE_DIR}/config/DefaultFlags.cmake)

get_directory_property(is-subproject PARENT_DIRECTORY)

project(
  "libsparse"
  LANGUAGES Fortran C
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

if(NOT TARGET "MKL::MKL")
  list(APPEND CMAKE_PREFIX_PATH "$ENV{MKLROOT}/lib/cmake/mkl")
  set(MKL_INTERFACE "lp64")
  set(MKL_LINK "static")
  set(ENABLE_BLAS95 true)
  set(ENABLE_LAPACK95 true)
  find_package(MKL REQUIRED)
endif()

message(STATUS "MKL found: ${MKL_IMPORTED_TARGETS}")

# Find METIS library
if(NOT DEFINED CMAKE_METIS_LIB)
    find_library(CMAKE_METIS_LIB metis)
endif()

if(CMAKE_METIS_LIB)
  set(WITH_METIS 1 STRING)
  message(STATUS "METIS found: ${CMAKE_METIS_LIB}")
else()
  set(WITH_METIS 0 STRING)
  message(STATUS "METIS not found")
endif()

set(
  lib-deps
  "OpenMP::OpenMP_Fortran"
  "MKL::MKL"
)

if(CMAKE_METIS_LIB)
  list(
    APPEND lib-deps
    "${CMAKE_METIS_LIB}"
  )
endif()

#MKL PARDISO
#if(NOT CMAKE_MKLPARDISO AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
#  find_file(CMAKE_MKLPARDISO "mkl_pardiso.f90" HINTS $ENV{MKLROOT} PATH_SUFFIXES "include")
#endif()

#if(CMAKE_MKLPARDISO OR CMAKE_Fortran_COMPILER_ID MATCHES "^Intel" )
  set(WITH_MKLPARDISO 1 STRING)
#else()
#  set(WITH_MKLPARDISO 0 STRING)
#endif()
message(STATUS "Support of MKL PARDISO: ${WITH_MKLPARDISO}")

#SPAINV
if(NOT DEFINED CMAKE_SPAINV)
  set(WITH_SPAINV 0 CACHE STRING "Support of SPAINV")
else()
  set(WITH_SPAINV ${CMAKE_SPAINV} CACHE STRING "Support of SPAINV")
endif()

#SPAINV
if(NOT DEFINED CMAKE_VERBOSE)
  set(WITH_VERBOSE 0 CACHE STRING "Support of verbosity")
else()
  set(WITH_VERBOSE ${CMAKE_VERBOSE} CACHE STRING "Support of verbosity")
endif()



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
  "$ENV{MKLROOT}/include"  #SHOULD BE OPTIONAL
)
target_compile_definitions(
  "${PROJECT_NAME}-lib"
  PRIVATE
  "-D_DP=1"
  "-D_METIS=${WITH_METIS}"
  "-D_PARDISO=${WITH_MKLPARDISO}"
  "-D_SPAINV=${WITH_SPAINV}"
  "-D_VERBOSE=${WITH_VERBOSE}"
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
