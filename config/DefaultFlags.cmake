#Based on Fortran stdlib (License: MIT)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(
    CMAKE_Fortran_FLAGS_INIT
    "-cpp -fall-intrinsics"
  )
  set(
    CMAKE_Fortran_FLAGS_RELEASE_INIT
    "-O3"
  )
  set(
    CMAKE_Fortran_FLAGS_DEBUG_INIT
    "-g"
    "-fcheck=all"
    "-fbacktrace"
    "-Wall"
    "-Wextra"
    "-Wimplicit-procedure"
    "-std=f2018"
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^IntelLLVM")
  set(
    CMAKE_Fortran_FLAGS_INIT
    "-fpp -qopt-report=3"
  )
  set(
    CMAKE_Fortran_FLAGS_RELEASE_INIT
    "-O3"
  )
  if(WIN32)
    set(
      CMAKE_Fortran_FLAGS_DEBUG_INIT
      "/stand:f18"
      "/warn:declarations,general,usage,interfaces,unused"
    )
  else()
    set(
      CMAKE_Fortran_FLAGS_DEBUG_INIT
      "-g"
      "-stand f18"
      "-check all,nouninit -traceback -debug extended -debug inline-debug-info -check noarg_temp_created"
      "-warn all"
    )
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set(
    CMAKE_Fortran_FLAGS_INIT
    -fpp -heap-arrays -qopt-report=5
  )
  set(
    CMAKE_Fortran_FLAGS_RELEASE_INIT
    "-O3 -parallel -qopt-matmul"
  )
  if(WIN32)
    set(
      CMAKE_Fortran_FLAGS_DEBUG_INIT
      "/stand:f18"
      "/warn:declarations,general,usage,interfaces,unused"
    )
  else()
    set(
      CMAKE_Fortran_FLAGS_DEBUG_INIT
      "-g"
      "-stand f18"
      "-check all -traceback -debug extended -debug inline-debug-info -check noarg_temp_created"
      "-warn all"
    )
  endif()
else()
  set(
    CMAKE_Fortran_FLAGS_INIT
  )
  set(
    CMAKE_Fortran_FLAGS_RELEASE_INIT
  )
  set(
    CMAKE_Fortran_FLAGS_DEBUG_INIT
  )
endif()
string(REPLACE ";" " " CMAKE_Fortran_FLAGS_INIT "${CMAKE_Fortran_FLAGS_INIT}")
string(REPLACE ";" " " CMAKE_Fortran_FLAGS_RELEASE_INIT "${CMAKE_Fortran_FLAGS_RELEASE_INIT}")
string(REPLACE ";" " " CMAKE_Fortran_FLAGS_DEBUG_INIT "${CMAKE_Fortran_FLAGS_DEBUG_INIT}")
