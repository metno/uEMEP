program test_utility_functions

    use utility_functions

    implicit none

    logical :: ok = .true. ! Test suite boolean - must be true initially

    logical :: flag_out, flag_test
    character(len=256) :: txtstr_out, txtstr_test
    character(len=2048) :: path, relpath
    integer :: unit_in, iutm_in, isone_in
    real :: lat_in, lon_in, lon0_in, la0_in, utmn_in, utme_in
    real :: utmn_out, utme_out, lat_out, lon_out
    real :: utmn_test, utme_test, lat_test, lon_test
    real :: x0_in, y0_in, x1_in, y1_in, x2_in, y2_in
    real :: xm_out, ym_out, dm_out, wm_out
    real :: xm_test, ym_test, dm_test

    ! Test case 1 (nxtdat): check that file with any remaining uncommented lines returns false
    unit_in = 20
    flag_out = .true.
    flag_test = .false.
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_commented_file.txt", status="old")
    call nxtdat(unit_in, flag_out)
    close(unit_in)
    if (flag_out .neqv. flag_test) then
        ok = .false.
        print "(a)", "Test case 1 failed! Routine: nxtdat"
    end if

    ! Test case 2 (nxtdat): check that read returns the right uncommented line
    unit_in = 20
    txtstr_test = "non-commented line 4"
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_commented_file.txt", status="old")
    call nxtdat(unit_in, flag_out)
    read(unit_in,"(a)") txtstr_out
    close(unit_in)
    if (trim(txtstr_out) /= trim(txtstr_test)) then
        ok = .false.
        print "(a)", "Test case 2 failed! Routine: nxtdat"
    end if

    ! Test case 3 (nxtdat): check that file with only commented lines returns "end of file" = true
    unit_in = 20
    flag_out = .false.
    flag_test = .true.
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"
    open(unit=unit_in, file=trim(relpath)//"test_only_commented_file.txt", status="old")
    call nxtdat(unit_in, flag_out)
    close(unit_in)

    if (flag_out .neqv. flag_test) then
        ok = .false.
        print "(a)", "Test case 3 failed! Routine: nxtdat"
    end if

    ! Test case 4 (ll2utm): Check outputted UTM coordinates against expected values in the Northern hemisphere
    ! 
    ! Note that iutm /= different from 1 is not tested as it is not used!
    iutm_in = 1
    isone_in = 32
    lat_in = 55.0
    lon_in = 12.0
    utmn_out = 0.0
    utme_out = 0.0
    utmn_test = 6098908.0
    utme_test = 691875.6
    call ll2utm(iutm_in, isone_in, lat_in, lon_in, utmn_out, utme_out)
    if (abs(utmn_out - utmn_test) > 1.0e-4 .or. abs(utme_out - utme_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 4 failed! Routine: ll2utm"
    end if

    ! Test case 5 (ll2utm): Check outputted UTM coordinates against expected values in the Southern hemisphere
    ! 
    ! Note that iutm /= different from 1 is not tested as it is not used!
    iutm_in = 1
    isone_in = -32
    lat_in = -30.0
    lon_in = 12.0
    utmn_out = 0.0
    utme_out = 0.0
    utmn_test = 6677424.0
    utme_test = 789409.688
    call ll2utm(iutm_in, isone_in, lat_in, lon_in, utmn_out, utme_out)
    if (abs(utmn_out - utmn_test) > 1.0e-4 .or. abs(utme_out - utme_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 5 failed! Routine: ll2utm"
    end if

    ! Test case 6 (ll2ltm): Check outputted UTM coordinates against expected values in the Northern hemisphere
    ! 
    ! Note that iutm /= different from 1 is not tested as it is not used!
    iutm_in = 1
    lon0_in = 9.0
    lat_in = 55.0
    lon_in = 12.0
    utmn_out = 0.0
    utme_out = 0.0
    utmn_test = 6098908.0
    utme_test = 691875.625
    call ll2ltm(iutm_in, lon0_in, lat_in, lon_in, utmn_out, utme_out)
    if (abs(utmn_out - utmn_test) > 1.0e-4 .or. abs(utme_out - utme_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 6 failed! Routine: ll2ltm"
    end if

    ! Test case 7 (ll2ltm): Check outputted UTM coordinates against expected values in the Southern hemisphere
    ! 
    ! Note that iutm /= different from 1 is not tested as it is not used!
    iutm_in = 1
    lon0_in = 9.0
    lat_in = -30.0
    lon_in = 12.0
    utmn_out = 0.0
    utme_out = 0.0
    utmn_test = 6677424.0
    utme_test = 789409.688
    call ll2ltm(iutm_in, lon0_in, lat_in, lon_in, utmn_out, utme_out)
    if (abs(utmn_out - utmn_test) > 1.0e-4 .or. abs(utme_out - utme_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 7 failed! Routine: ll2ltm"
    end if

    ! Test case 8 (utm2ll): Check outputted lat/lon coordinates against expected values in the Northern hemisphere
    ! 
    ! Note that iutm /= different from 1 is not tested as it is not used!
    iutm_in = 1
    isone_in = 32
    utmn_in = 6098908.0
    utme_in = 691875.625
    lat_out = 0.0
    lon_out = 0.0
    lat_test = 55.0
    lon_test = 12.0
    call utm2ll(iutm_in, isone_in, utmn_in, utme_in, lat_out, lon_out)
    if (abs(lat_out - lat_test) > 1.0e-4 .or. abs(lon_out - lon_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 8 failed! Routine: utm2ll"
    end if

    ! Test case 9 (utm2ll): Check outputted lat/lon coordinates against expected values in the Southern hemisphere
    ! 
    ! Note that iutm /= different from 1 is not tested as it is not used!
    iutm_in = 1
    isone_in = -32
    utmn_in = 6677424.0
    utme_in = 789409.688
    lat_out = 0.0
    lon_out = 0.0
    lat_test = -30.0
    lon_test = 12.0
    call utm2ll(iutm_in, isone_in, utmn_in, utme_in, lat_out, lon_out)
    if (abs(lat_out - lat_test) > 1.0e-4 .or. abs(lon_out - lon_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 9 failed! Routine: utm2ll"
    end if

    ! Test case 10 (ltm2ll): Check outputted lat/lon coordinates against expected values in the Northern hemisphere
    ! 
    ! Note that iutm /= different from 1 is not tested as it is not used!
    iutm_in = 1
    isone_in = 32
    la0_in = 9.0
    utmn_in = 6098908.0
    utme_in = 691875.625
    lat_out = 0.0
    lon_out = 0.0
    lat_test = 55.0
    lon_test = 12.0
    call ltm2ll(iutm_in, isone_in, la0_in, utmn_in, utme_in, lat_out, lon_out)
    if (abs(lat_out - lat_test) > 1.0e-4 .or. abs(lon_out - lon_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 10 failed! Routine: ltm2ll"
    end if

    ! Test case 11 (ltm2ll): Check outputted lat/lon coordinates against expected values in the Southern hemisphere
    ! 
    ! Note that iutm /= different from 1 is not tested as it is not used!
    iutm_in = 1
    isone_in = -32
    la0_in = 9.0
    utmn_in = 6677424.0
    utme_in = 789409.688
    lat_out = 0.0
    lon_out = 0.0
    lat_test = -30.0
    lon_test = 12.0
    call ltm2ll(iutm_in, isone_in, la0_in, utmn_in, utme_in, lat_out, lon_out)
    if (abs(lat_out - lat_test) > 1.0e-4 .or. abs(lon_out - lon_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 11 failed! Routine: ltm2ll"
    end if

    ! Test case 12 (distrl): Check returned distances against expected values when "line is a point"
    ! and when receptor point is on the "line"
    x0_in = 0.0
    y0_in = 0.0
    x1_in = 0.0
    y1_in = 0.0
    x2_in = 0.0
    y2_in = 0.0
    xm_out = -999.0
    ym_out = -999.0
    dm_out = -999.0
    wm_out = 0.0
    xm_test = 0.0
    ym_test = 0.0
    dm_test = 0.0
    call distrl(x0_in, y0_in, x1_in, y1_in, x2_in, y2_in, xm_out, ym_out, dm_out, wm_out)
    if (abs(xm_out - xm_test) > 1.0e-4 .or. abs(ym_out - ym_test) > 1.0e-4 .or. abs(dm_out - dm_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 12 failed! Routine: distrl"
    end if

    ! Test case 13 (distrl): Check returned distances against expected values when "line is a point"
    ! and when receptor point is above the "line"
    x0_in = 2.0
    y0_in = 2.0
    x1_in = 0.0
    y1_in = 0.0
    x2_in = 0.0
    y2_in = 0.0
    xm_out = -999.0
    ym_out = -999.0
    dm_out = -999.0
    wm_out = 0.0
    xm_test = 0.0
    ym_test = 0.0
    dm_test = 2.8284
    call distrl(x0_in, y0_in, x1_in, y1_in, x2_in, y2_in, xm_out, ym_out, dm_out, wm_out)
    if (abs(xm_out - xm_test) > 1.0e-4 .or. abs(ym_out - ym_test) > 1.0e-4 .or. abs(dm_out - dm_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 13 failed! Routine: distrl"
    end if

    ! Test case 14 (distrl): Check returned distances against expected values when line is
    ! "increasing" and receptor point is above the line
    x0_in = 1.0
    y0_in = 4.0
    x1_in = 0.0
    y1_in = 0.0
    x2_in = 3.0
    y2_in = 4.0
    xm_out = -999.0
    ym_out = -999.0
    dm_out = -999.0
    wm_out = 0.0
    xm_test = 2.28
    ym_test = 3.04
    dm_test = 1.6
    call distrl(x0_in, y0_in, x1_in, y1_in, x2_in, y2_in, xm_out, ym_out, dm_out, wm_out)
    if (abs(xm_out - xm_test) > 1.0e-4 .or. abs(ym_out - ym_test) > 1.0e-4 .or. abs(dm_out - dm_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 14 failed! Routine: distrl"
    end if

    ! Test case 15 (distrl): Check returned distances against expected values when line is
    ! "increasing" and receptor point is below the line
    x0_in = 4.0
    y0_in = 1.0
    x1_in = 0.0
    y1_in = 0.0
    x2_in = 3.0
    y2_in = 4.0
    xm_out = -999.0
    ym_out = -999.0
    dm_out = -999.0
    wm_out = 0.0
    xm_test = 1.92
    ym_test = 2.56
    dm_test = 2.6
    call distrl(x0_in, y0_in, x1_in, y1_in, x2_in, y2_in, xm_out, ym_out, dm_out, wm_out)
    if (abs(xm_out - xm_test) > 1.0e-4 .or. abs(ym_out - ym_test) > 1.0e-4 .or. abs(dm_out - dm_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 15 failed! Routine: distrl"
    end if

    ! Test case 16 (distrl): Check returned distances against expected values when line is
    ! "decreasing" and receptor point is above the line
    x0_in = 8.0
    y0_in = 8.0
    x1_in = 1.0
    y1_in = 5.0
    x2_in = 6.0
    y2_in = 2.0
    xm_out = -999.0
    ym_out = -999.0
    dm_out = -999.0
    wm_out = 0.0
    xm_test = 4.8235
    ym_test = 2.7059
    dm_test = 6.1739
    call distrl(x0_in, y0_in, x1_in, y1_in, x2_in, y2_in, xm_out, ym_out, dm_out, wm_out)
    if (abs(xm_out - xm_test) > 1.0e-4 .or. abs(ym_out - ym_test) > 1.0e-4 .or. abs(dm_out - dm_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 16 failed! Routine: distrl"
    end if

    ! Test case 17 (distrl_sqr): Check returned distances against expected values when "line is a point"
    ! and when receptor point is on the "line"
    x0_in = 0.0
    y0_in = 0.0
    x1_in = 0.0
    y1_in = 0.0
    x2_in = 0.0
    y2_in = 0.0
    xm_out = -999.0
    ym_out = -999.0
    dm_out = -999.0
    wm_out = 0.0
    xm_test = 0.0
    ym_test = 0.0
    dm_test = 0.0
    call distrl_sqr(x0_in, y0_in, x1_in, y1_in, x2_in, y2_in, xm_out, ym_out, dm_out, wm_out)
    if (abs(xm_out - xm_test) > 1.0e-4 .or. abs(ym_out - ym_test) > 1.0e-4 .or. abs(dm_out - dm_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 17 failed! Routine: distrl_sqr"
    end if

    ! Test case 18 (distrl_sqr): Check returned distances against expected values when "line is a point"
    ! and when receptor point is above the "line"
    x0_in = 2.0
    y0_in = 2.0
    x1_in = 0.0
    y1_in = 0.0
    x2_in = 0.0
    y2_in = 0.0
    xm_out = -999.0
    ym_out = -999.0
    dm_out = -999.0
    wm_out = 0.0
    xm_test = 0.0
    ym_test = 0.0
    dm_test = 8.0
    call distrl_sqr(x0_in, y0_in, x1_in, y1_in, x2_in, y2_in, xm_out, ym_out, dm_out, wm_out)
    if (abs(xm_out - xm_test) > 1.0e-4 .or. abs(ym_out - ym_test) > 1.0e-4 .or. abs(dm_out - dm_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 18 failed! Routine: distrl_sqr"
    end if

    ! Test case 19 (distrl_sqr): Check returned distances against expected values when line is
    ! "increasing" and receptor point is above the line
    x0_in = 1.0
    y0_in = 4.0
    x1_in = 0.0
    y1_in = 0.0
    x2_in = 3.0
    y2_in = 4.0
    xm_out = -999.0
    ym_out = -999.0
    dm_out = -999.0
    wm_out = 0.0
    xm_test = 2.28
    ym_test = 3.04
    dm_test = 2.56
    call distrl_sqr(x0_in, y0_in, x1_in, y1_in, x2_in, y2_in, xm_out, ym_out, dm_out, wm_out)
    if (abs(xm_out - xm_test) > 1.0e-4 .or. abs(ym_out - ym_test) > 1.0e-4 .or. abs(dm_out - dm_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 19 failed! Routine: distrl_sqr"
    end if

    ! Test case 20 (distrl_sqr): Check returned distances against expected values when line is
    ! "increasing" and receptor point is below the line
    x0_in = 4.0
    y0_in = 1.0
    x1_in = 0.0
    y1_in = 0.0
    x2_in = 3.0
    y2_in = 4.0
    xm_out = -999.0
    ym_out = -999.0
    dm_out = -999.0
    wm_out = 0.0
    xm_test = 1.92
    ym_test = 2.56
    dm_test = 6.76
    call distrl_sqr(x0_in, y0_in, x1_in, y1_in, x2_in, y2_in, xm_out, ym_out, dm_out, wm_out)
    if (abs(xm_out - xm_test) > 1.0e-4 .or. abs(ym_out - ym_test) > 1.0e-4 .or. abs(dm_out - dm_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 20 failed! Routine: distrl_sqr"
    end if

    ! Test case 21 (distrl_sqr): Check returned distances against expected values when line is
    ! "decreasing" and receptor point is above the line
    x0_in = 8.0
    y0_in = 8.0
    x1_in = 1.0
    y1_in = 5.0
    x2_in = 6.0
    y2_in = 2.0
    xm_out = -999.0
    ym_out = -999.0
    dm_out = -999.0
    wm_out = 0.0
    xm_test = 4.8235
    ym_test = 2.7059
    dm_test = 38.1177
    call distrl_sqr(x0_in, y0_in, x1_in, y1_in, x2_in, y2_in, xm_out, ym_out, dm_out, wm_out)
    if (abs(xm_out - xm_test) > 1.0e-4 .or. abs(ym_out - ym_test) > 1.0e-4 .or. abs(dm_out - dm_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 21 failed! Routine: distrl_sqr"
    end if

    ! Return test results summary
    if (ok) then
        print "(a)", "test_utility_functions: All tests passed."
    else
        print "(a)", "test_utility_functions: One ore more tests failed."
        stop 1
    end if

end program test_utility_functions