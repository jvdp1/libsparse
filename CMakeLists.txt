cmake_minimum_required(VERSION 3.10)
get_directory_property(is-subproject PARENT_DIRECTORY)

# Include overwrites before setting up the project
set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_CURRENT_SOURCE_DIR}/config/DefaultFlags.cmake)

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

find_package("OpenMP" REQUIRED)

if(NOT TARGET "BLAS::BLAS")
  find_package("custom-blas" REQUIRED)
endif()


#Compiler
#if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
#  include(Compiler-Intel)
#elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
#  include(Compiler-GNU)
#elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
#  include(Compiler-PGI)
#else ()
#  message (FATAL_ERROR "Compiler type (CMAKE_Fortran_COMPILER_ID) not recognized")
#endif ()
#
#set (CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE} "${BUILDCMD_PRE} ${OPENMP_FLAGS} ${STANDARD_CHECK} ${OPTIMIZATION_FLAGS} ${COMPILE_FLAGS} ${PRECISION_FLAGS}")
#set (CMAKE_Fortran_LINK_EXECUTABLE           "${BUILDCMD_PRE} ${CMAKE_Fortran_COMPILER} ${OPENMP_FLAGS}   ${OPTIMIZATION_FLAGS} ${LINKER_FLAGS}")
#
#if (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
#  set (CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE} "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}} ${DEBUG_FLAGS}")
#  set (CMAKE_Fortran_LINK_EXECUTABLE           "${CMAKE_Fortran_LINK_EXECUTABLE} ${DEBUG_FLAGS}")
#endif ()
#
#set (CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}   "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}} ${BUILDCMD_POST}")
#set (CMAKE_Fortran_LINK_EXECUTABLE             "${CMAKE_Fortran_LINK_EXECUTABLE} <OBJECTS> -o <TARGET> <LINK_LIBRARIES> ${BUILDCMD_POST}")
#
#
#message ("Fortran Compiler Flags:\n${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}\n")
#message ("C Compiler Flags:\n${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}\n")
#message ("Fortran Linker Command:\n${CMAKE_Fortran_LINK_EXECUTABLE}\n")


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
target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)
target_compile_definitions(
  "${PROJECT_NAME}-lib"
  PRIVATE
  "WITH_QP=$<BOOL:${WITH_QP}>"
  "-D_DP=1"
  "-D_METIS=0"
  "-D_PARDISO=1"
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



## add the testsuite
#include(CTest)
#if(BUILD_TESTING AND TEST_DRIVE_BUILD_TESTING)
#  add_subdirectory("test")
#endif()
