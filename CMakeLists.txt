cmake_minimum_required(VERSION 3.20)

# Project description
project(uemep)

# Hardcode ifort as the Fortran compiler (gfortran currently does not work)
set(CMAKE_Fortran_COMPILER "ifort")
enable_language(Fortran)

# Append additional cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Define compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(F90DIALECT "") # For future use
  set(F90FLAGS "-O3") # Compiler flags used for release mode
  set(F90DEBUGFLAGS "-O0 -g -debug all -traceback -check all -check bounds") # debug mode
  set(F90DISABLE "-diag-disable=10448") # Disable specific diagnostics 
endif()

# Set compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${F90FLAGS} ${F90DISABLE}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} ${F90DEBUGFLAGS}")

# Import NetCDF
set(NETCDF_F90 YES)
find_package(NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDES})

# Glob source files
file(GLOB_RECURSE SOURCES src/*.f90 src/*.for)

# Create executable and link
add_executable(uemep ${SOURCES})
target_link_libraries(uemep ${NETCDF_LIBRARIES})

# message(STATUS "${NFLIB}")

#file(GLOB_RECURSE SRC src/*.f90 src/*.for)
#add_executable(uEMEP ${SRC})
