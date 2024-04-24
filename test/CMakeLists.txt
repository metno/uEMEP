foreach(execid
    utility_functions
    )

    add_executable(test_${execid} test_${execid}.f90 ${SOURCES})
    target_link_libraries(test_${execid} ${NETCDF_LIBRARIES})
    add_test(NAME test_${execid} COMMAND test_${execid})
endforeach()