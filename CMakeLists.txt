cmake_minimum_required(VERSION 3.20)

# Project description
project(uemep LANGUAGES NONE)

# Enable testing
enable_testing()

# Set fortran compiler
if(NOT DEFINED CMAKE_Fortran_COMPILER)
  if(NOT DEFINED ENV{FC})
    set(CMAKE_Fortran_COMPILER "ifort" CACHE FILEPATH "Fortran compiler")
  endif()
endif()
enable_language(Fortran)
message(STATUS "Using Fortran compiler: ${CMAKE_Fortran_COMPILER}")

# Default single-configuration generators to optimized builds unless the user
# explicitly selected another build type.
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Append additional cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Common dialect options; optimization and runtime checks are configuration-specific.
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  add_compile_options(-free -implicitnone)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    add_compile_options(-diag-disable=10448)
  endif()
  set(F90DEBUGFLAGS "-O0 -g -traceback")
  set(F90RELEASEFLAGS "-O3")
  set(F90CHECKFLAGS "-O0 -g -traceback -warn all -warn interfaces -warn declarations -check all -init=snan,arrays -fpe0")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-ffree-form -fimplicit-none -ffree-line-length-none -fdec)
  set(F90DEBUGFLAGS "-O0 -g -fbacktrace")
  set(F90RELEASEFLAGS "-O2")
  set(F90CHECKFLAGS "-O0 -g -fbacktrace -Wall -Wextra -Wimplicit-interface -Wimplicit-procedure -Wcharacter-truncation -Warray-temporaries -Wconversion-extra -fcheck=all -finit-real=snan -ffpe-trap=invalid,zero,overflow")
else()
  message(FATAL_ERROR "Unsupported Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

# Debug supports interactive stepping; Check instruments both the model and tests.
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${F90DEBUGFLAGS}")
set(CMAKE_Fortran_FLAGS_CHECK "${CMAKE_Fortran_FLAGS_CHECK} ${F90CHECKFLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "${F90RELEASEFLAGS}")
add_compile_options("$<$<CONFIG:Release,RelWithDebInfo>:-march=core-avx2>")

# Custom configurations also work with multi-configuration generators.
set(UEMEP_BUILD_TYPES Debug Check Release RelWithDebInfo)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 14)
  list(APPEND UEMEP_BUILD_TYPES Coverage)
  set(CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_COVERAGE} -O0 -g -fbacktrace --coverage -fno-inline")
endif()
if(DEFINED CACHE{CMAKE_BUILD_TYPE})
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${UEMEP_BUILD_TYPES})
endif()
if(CMAKE_CONFIGURATION_TYPES)
  list(APPEND CMAKE_CONFIGURATION_TYPES Check)
  if("Coverage" IN_LIST UEMEP_BUILD_TYPES)
    list(APPEND CMAKE_CONFIGURATION_TYPES Coverage)
  endif()
  list(REMOVE_DUPLICATES CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING "Available build configurations" FORCE)
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" UEMEP_BUILD_TYPE)
if(UEMEP_BUILD_TYPE STREQUAL "COVERAGE" OR "Coverage" IN_LIST CMAKE_CONFIGURATION_TYPES)
  if(NOT "Coverage" IN_LIST UEMEP_BUILD_TYPES)
    message(FATAL_ERROR "The 'Coverage' build type requires GNU Fortran >= 14")
  endif()
endif()

# Import NetCDF
find_package(NetCDF REQUIRED)
include_directories(${NC_INC})

# Add source files
set(SOURCES)
include(src/CMakeLists.txt)

# Build uEMEP as static library
add_library(uemeplib ${SOURCES})

# Build and link test executables
include(test/CMakeLists.txt)

# Build and link main executable
add_executable(uemep src/uEMEP_main.f90)
target_link_libraries(uemep PRIVATE uemeplib ${NC_LIB})
