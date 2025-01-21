foreach(execid
    utility_functions
    read_namefile_routines
    time_functions
    uemep_logger
    io_functions
    area_interpolation_functions
    )

    add_executable(test_${execid} test/test_${execid}.f90)
    target_link_libraries(test_${execid} PRIVATE uemeplib ${NETCDF_LIBRARIES})
    add_test(NAME test_${execid} COMMAND test_${execid})
endforeach()