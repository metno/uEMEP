# Append source files
list(APPEND SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/utility_functions.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_definitions.f90
  ${CMAKE_CURRENT_LIST_DIR}/Area_weighted_interpolation_function.f90
  ${CMAKE_CURRENT_LIST_DIR}/lambert_projection.f90
  ${CMAKE_CURRENT_LIST_DIR}/read_esri_ascii_file.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_aggregate_proxy_emission_in_EMEP_grid.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_auto_subgrid.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_calculate_exposure.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_crossreference_grids.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_definitions.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_dispersion_functions.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_grid_roads.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_Kz_functions.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_local_trajectory.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_define_subgrid.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_redistribute_data.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_set_constants.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_set_dispersion_params.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_set_emission_factors.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_set_filenames.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_set_subgrids.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_subgrid_deposition_EMEP.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_subgrid_deposition.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_subgrid_dispersion.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_subgrid_EMEP.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_subgrid_emission_EMEP.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_subgrid_meteo_EMEP.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_tiling_routines.f90
  ${CMAKE_CURRENT_LIST_DIR}/uEMEP_time_functions.f90)

# Include source files in sub-directories
include(${CMAKE_CURRENT_LIST_DIR}/extern/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/global/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/read_data/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/save_data/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/chemistry/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/utilities/CMakeLists.txt)
