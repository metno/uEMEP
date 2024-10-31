cmake_minimum_required(VERSION 3.20)

# Project description
project(uemep)

# Enable testing
enable_testing()

# Hardcode ifort as the Fortran compiler (gfortran currently does not work)
set(CMAKE_Fortran_COMPILER "ifort")
#set(CMAKE_Fortran_COMPILER "gfortran")
enable_language(Fortran)

# Append additional cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Define compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(F90DIALECT "-free -implicitnone")
  set(F90FLAGS "-O3 -march=core-avx2")
  set(F90DEBUGFLAGS "-O0 -g -debug all -traceback -init=snan,arrays -fpe1 -ftrapuv -check all -check bounds -warn all") # debug mode
  set(F90DISABLE "-diag-disable=10448")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(F90DIALECT "-ffree-form -fimplicit-none")
  set(F90FLAGS "-O3 -march=core-avx2 -ffree-line-length-none -fdec")
  set(F90DEBUGFLAGS "-O0 -g -fbacktrace -Wall -Wextra -fcheck=all")
  set(F90DISABLE "")
endif()

# Set compiler flags
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${F90DEBUGFLAGS} ${F90DIALECT} ${F90DISABLE}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${F90FLAGS} ${F90DIALECT} ${F90DISABLE}")

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

