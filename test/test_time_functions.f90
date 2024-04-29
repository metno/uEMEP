program test_time_functions

    use uemep_constants, only: dp
    use time_functions

    implicit none

    logical :: ok = .true.

    integer :: ref_year
    integer :: date_out(6), date_test(6), date_in(6)
    integer :: julian_out, julian_test
    integer :: wday_out, wday_test
    real(dp) :: date_num_in, date_num_out, date_num_test
    character(len=:), allocatable :: format_in, date_str_in, date_str_test
    character(len=256) :: date_str_out
    logical :: summer_out, summer_test
    real :: lat_in, lon_in
    real :: azimuth_out, zenith_out
    real :: azimuth_test, zenith_test
    
    ! Test case 1 (number_to_date): Check against expected return dates
    ref_year = 1900
    date_num_in = 43582.58
    date_out = [0, 0, 0, 0, 0, 0]
    date_test = [2019, 4, 29, 13, 52, 30]
    call number_to_date(date_num_in, date_out, ref_year)
    if (any(date_out /= date_test)) then
        ok = .false.
        print "(a)", "Test case 1 failed! Routine: number_to_date"
    end if

    ! Test case 2 (number_to_date): Check that January 1st at the reference date is returned at date number 0.0
    date_num_in = 0.0
    ref_year = 1900
    date_out = [0, 0, 0, 0, 0, 0]
    date_test = [1900, 1, 1, 0, 0, 0]
    call number_to_date(date_num_in, date_out, ref_year)
    if (any(date_out /= date_test)) then
        ok = .false.
        print "(a)", "Test case 2 failed! Routine: number_to_date"
    end if

    ! Test case 3 (number_to_date): Check that February 29th is returned in a leap year
    ref_year = 1900
    date_num_in = 45349.5
    date_out = [0, 0, 0, 0, 0, 0]
    date_test = [2024, 2, 29, 12, 0, 0]
    call number_to_date(date_num_in, date_out, ref_year)
    if (any(date_out /= date_test)) then
        ok = .false.
        print "(a)", "Test case 3 failed! Routine: number_to_date"
    end if

    ! Test case 4 (date_to_number): Check against expected return value
    ref_year = 1900
    date_in = [2019, 4, 29, 13, 52, 30]
    date_num_test = 43582.58
    date_num_out = -999.0
    date_num_out = date_to_number(date_in, ref_year)
    if (abs(date_num_out - date_num_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 4 failed! Routine: date_to_number"
    end if

    ! Test case 5 (date_to_number): Check that the reference date returns 0.0
    ref_year = 1900
    date_in = [1900, 1, 1, 0, 0, 0]
    date_num_test = 0.0
    date_num_out = -999.0
    date_num_out = date_to_number(date_in, ref_year)
    if (abs(date_num_out - date_num_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 5 failed! Routine: date_to_number"
    end if

    ! Test case 6 (date_to_number): Check expected value at Feb 29 in a leap year
    ref_year = 1900
    date_in = [2024, 2, 29, 12, 0, 0]
    date_num_test = 45349.5
    date_num_out = -999.0
    date_num_out = date_to_number(date_in, ref_year)
    if (abs(date_num_out - date_num_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 6 failed! Routine: date_to_number"
    end if

    ! Test case 7 (date_to_julian): Check expected value in normal year
    ref_year = 1900
    date_in = [2019, 4, 29, 13, 52, 30]
    julian_test = 119
    julian_out = -999
    julian_out = date_to_julian(date_in, ref_year)
    if (julian_out /= julian_test) then
        ok = .false.
        print "(a)", "Test case 7 failed! Routine: date_to_julian"
    end if

    ! Test case 8 (date_to_julian): Check expected value in leap year
    ref_year = 1900
    date_in = [2024, 12, 31, 12, 0, 0]
    julian_test = 366
    julian_out = -999
    julian_out = date_to_julian(date_in, ref_year)
    if (julian_out /= julian_test) then
        ok = .false.
        print "(a)", "Test case 8 failed! Routine: date_to_julian"
    end if

    ! Test case 9 (datestr_to_date): Check against expected output with yyyy:
    format_in = "yyyy"
    date_str_in = "20190429123005"
    date_test = [2019, 0, 0, 0, 0, 0]
    date_out = [0, 0, 0, 0, 0, 0]
    call datestr_to_date(date_str_in, format_in, date_out)
    if (any(date_out /= date_test)) then
        ok = .false.
        print "(a)", "Test case 9 failed! Routine: datestr_to_date"
    end if
    deallocate(format_in, date_str_in)

    ! Test case 10 (datestr_to_date): Check against expected output with yyyymm:
    format_in = "yyyymm"
    date_str_in = "20190429123005"
    date_test = [2019, 4, 0, 0, 0, 0]
    date_out = [0, 0, 0, 0, 0, 0]
    call datestr_to_date(date_str_in, format_in, date_out)
    if (any(date_out /= date_test)) then
        ok = .false.
        print "(a)", "Test case 10 failed! Routine: datestr_to_date"
    end if
    deallocate(format_in, date_str_in)

    ! Test case 11 (datestr_to_date): Check against expected output with yyyymmdd:
    format_in = "yyyymmdd"
    date_str_in = "20190429123005"
    date_test = [2019, 4, 29, 0, 0, 0]
    date_out = [0, 0, 0, 0, 0, 0]
    call datestr_to_date(date_str_in, format_in, date_out)
    if (any(date_out /= date_test)) then
        ok = .false.
        print "(a)", "Test case 11 failed! Routine: datestr_to_date"
    end if
    deallocate(format_in, date_str_in)

    ! Test case 12 (datestr_to_date): Check against expected output with yyyymmddHH:
    format_in = "yyyymmddHH"
    date_str_in = "20190429123005"
    date_test = [2019, 4, 29, 12, 0, 0]
    date_out = [0, 0, 0, 0, 0, 0]
    call datestr_to_date(date_str_in, format_in, date_out)
    if (any(date_out /= date_test)) then
        ok = .false.
        print "(a)", "Test case 12 failed! Routine: datestr_to_date"
    end if
    deallocate(format_in, date_str_in)

    ! Test case 13 (datestr_to_date): Check against expected output with yyyymmddHHMM:
    format_in = "yyyymmddHHMM"
    date_str_in = "20190429123005"
    date_test = [2019, 4, 29, 12, 30, 0]
    date_out = [0, 0, 0, 0, 0, 0]
    call datestr_to_date(date_str_in, format_in, date_out)
    if (any(date_out /= date_test)) then
        ok = .false.
        print "(a)", "Test case 11 failed3 Routine: datestr_to_date"
    end if
    deallocate(format_in, date_str_in)

    ! Test case 14 (datestr_to_date): Check against expected output with yyyymmddHHMMSS:
    format_in = "yyyymmddHHMMSS"
    date_str_in = "20190429123005"
    date_test = [2019, 4, 29, 12, 30, 5]
    date_out = [0, 0, 0, 0, 0, 0]
    call datestr_to_date(date_str_in, format_in, date_out)
    if (any(date_out /= date_test)) then
        ok = .false.
        print "(a)", "Test case 14 failed! Routine: datestr_to_date"
    end if
    deallocate(format_in, date_str_in)

    ! Test case 15 (date_to_datestr): Check against expected output with yyyymm:
    format_in = "yyyymm"
    date_in = [2019, 4, 29, 12, 30, 5]
    date_str_test = "201904"
    date_str_out = ""
    call date_to_datestr(date_in, format_in, date_str_out)
    if (trim(date_str_out) /= date_str_test) then
        ok = .false.
        print "(a)", "Test case 15 failed! Routine date_to_datestr"
    end if
    deallocate(format_in, date_str_test)

    ! Test case 16 (date_to_datestr): Check against expected output with yyyymmdd:
    format_in = "yyyymmdd"
    date_in = [2019, 4, 29, 12, 30, 5]
    date_str_test = "20190429"
    date_str_out = ""
    call date_to_datestr(date_in, format_in, date_str_out)
    if (trim(date_str_out) /= date_str_test) then
        ok = .false.
        print "(a)", "Test case 16 failed! Routine date_to_datestr"
    end if
    deallocate(format_in, date_str_test)

    ! Test case 17 (date_to_datestr): Check against expected output with yyyymmddHH:
    format_in = "yyyymmddHH"
    date_in = [2019, 4, 29, 12, 30, 5]
    date_str_test = "2019042912"
    date_str_out = ""
    call date_to_datestr(date_in, format_in, date_str_out)
    if (trim(date_str_out) /= date_str_test) then
        ok = .false.
        print "(a)", "Test case 17 failed! Routine date_to_datestr"
    end if
    deallocate(format_in, date_str_test)

    ! Test case 18 (date_to_datestr_bracket): Check against expected output with yyyymm:
    format_in = "<yyyymm>"
    date_in = [2019, 4, 29, 12, 30, 5]
    date_str_test = "201904"
    date_str_out = ""
    call date_to_datestr_bracket(date_in, format_in, date_str_out)
    if (trim(date_str_out) /= date_str_test) then
        ok = .false.
        print "(a)", "Test case 18 failed! Routine date_to_datestr_bracket"
    end if
    deallocate(format_in, date_str_test)

    ! Test case 19 (date_to_datestr_bracket): Check against expected output with yyyymmdd:
    format_in = "<yyyymmdd>"
    date_in = [2019, 4, 29, 12, 30, 5]
    date_str_test = "20190429"
    date_str_out = ""
    call date_to_datestr_bracket(date_in, format_in, date_str_out)
    if (trim(date_str_out) /= date_str_test) then
        ok = .false.
        print "(a)", "Test case 19 failed! Routine date_to_datestr_bracket"
    end if
    deallocate(format_in, date_str_test)

    ! Test case 20 (date_to_datestr_bracket): Check against expected output with yyyymmddHH:
    format_in = "<yyyymmddHH>"
    date_in = [2019, 4, 29, 12, 30, 5]
    date_str_test = "2019042912"
    date_str_out = ""
    call date_to_datestr_bracket(date_in, format_in, date_str_out)
    if (trim(date_str_out) /= date_str_test) then
        ok = .false.
        print "(a)", "Test case 20 failed! Routine date_to_datestr_bracket"
    end if
    deallocate(format_in, date_str_test)

    ! Test case 21 (date_to_datestr_bracket): Check against expected output with yyyymmddHHMM:
    format_in = "<yyyymmddHHMM>"
    date_in = [2019, 4, 29, 12, 30, 5]
    date_str_test = "201904291230"
    date_str_out = ""
    call date_to_datestr_bracket(date_in, format_in, date_str_out)
    if (trim(date_str_out) /= date_str_test) then
        ok = .false.
        print "(a)", "Test case 21 failed! Routine date_to_datestr_bracket"
    end if
    deallocate(format_in, date_str_test)

    ! Test case 22 (date_to_datestr_bracket): Check against expected output with yyyymmddHHMMSS:
    format_in = "<yyyymmddHHMMSS>"
    date_in = [2019, 4, 29, 12, 30, 5]
    date_str_test = "20190429123005"
    date_str_out = ""
    call date_to_datestr_bracket(date_in, format_in, date_str_out)
    if (trim(date_str_out) /= date_str_test) then
        ok = .false.
        print "(a)", "Test case 22 failed! Routine date_to_datestr_bracket"
    end if
    deallocate(format_in, date_str_test)

    ! Test case 23 (date_to_datestr_squarebracket): Check against expected output with yyyymm:
    format_in = "[yyyymm]"
    date_in = [2019, 4, 29, 12, 30, 5]
    date_str_test = "201904"
    date_str_out = ""
    call date_to_datestr_squarebracket(date_in, format_in, date_str_out)
    if (trim(date_str_out) /= date_str_test) then
        ok = .false.
        print "(a)", "Test case 23 failed! Routine date_to_datestr_squarebracket"
    end if
    deallocate(format_in, date_str_test)

    ! Test case 24 (date_to_datestr_squarebracket): Check against expected output with yyyymmdd:
    format_in = "[yyyymmdd]"
    date_in = [2019, 4, 29, 12, 30, 5]
    date_str_test = "20190429"
    date_str_out = ""
    call date_to_datestr_squarebracket(date_in, format_in, date_str_out)
    if (trim(date_str_out) /= date_str_test) then
        ok = .false.
        print "(a)", "Test case 24 failed! Routine date_to_datestr_squarebracket"
    end if
    deallocate(format_in, date_str_test)

    ! Test case 25 (date_to_datestr_squarebracket): Check against expected output with yyyymmddHH:
    format_in = "[yyyymmddHH]"
    date_in = [2019, 4, 29, 12, 30, 5]
    date_str_test = "2019042912"
    date_str_out = ""
    call date_to_datestr_squarebracket(date_in, format_in, date_str_out)
    if (trim(date_str_out) /= date_str_test) then
        ok = .false.
        print "(a)", "Test case 25 failed! Routine date_to_datestr_squarebracket"
    end if
    deallocate(format_in, date_str_test)

    ! Test case 26 (date_to_datestr_bracket): Check against expected output with yyyymmddHHMM:
    format_in = "[yyyymmddHHMM]"
    date_in = [2019, 4, 29, 12, 30, 5]
    date_str_test = "201904291230"
    date_str_out = ""
    call date_to_datestr_squarebracket(date_in, format_in, date_str_out)
    if (trim(date_str_out) /= date_str_test) then
        ok = .false.
        print "(a)", "Test case 26 failed! Routine date_to_datestr_squarebracket"
    end if
    deallocate(format_in, date_str_test)

    ! Test case 27 (date_to_datestr_squarebracket): Check against expected output with yyyymmddHHMMSS:
    format_in = "[yyyymmddHHMMSS]"
    date_in = [2019, 4, 29, 12, 30, 5]
    date_str_test = "20190429123005"
    date_str_out = ""
    call date_to_datestr_squarebracket(date_in, format_in, date_str_out)
    if (trim(date_str_out) /= date_str_test) then
        ok = .false.
        print "(a)", "Test case 27 failed! Routine date_to_datestr_squarebracket"
    end if
    deallocate(format_in, date_str_test)

    ! Test case 28 (day_of_week): Check against expected output on normal day
    date_in = [2019, 5, 3, 12, 30, 5]
    wday_test = 5
    wday_out = -999
    wday_out = day_of_week(date_in)
    if (wday_out /= wday_test) then
        ok = .false.
        print "(a)", "Test case 28 failed! Routine: day_of_week"
    end if

    ! Test case 29 (day_of_week): Check against expected output on Feb 29
    date_in = [2024, 2, 29, 12, 30, 5]
    wday_test = 4
    wday_out = -999
    wday_out = day_of_week(date_in)
    if (wday_out /= wday_test) then
        ok = .false.
        print "(a)", "Test case 28 failed! Routine: day_of_week"
    end if

    ! Test case 30 (summer_time_europe): Check against expected output in summer
    date_in = [2022, 7, 15, 12, 0, 0]
    summer_test = .true.
    summer_out = .false.
    summer_out = summer_time_europe(date_in)
    if (summer_out .neqv. summer_test) then
        ok = .false.
        print "(a)", "Test case 30 failed! Routine: summer_time_europe"
    end if

    ! Test case 31 (summer_time_europe): Check against expected output in winter
    date_in = [2022, 12, 12, 12, 0, 0]
    summer_test = .false.
    summer_out = .true.
    summer_out = summer_time_europe(date_in)
    if (summer_out .neqv. summer_test) then
        ok = .false.
        print "(a)", "Test case 31 failed! Routine: summer_time_europe"
    end if

    ! Test case 32 (global_radiation_sub): Check against expected values Northern hemisphere
    date_in = [1995, 6, 1, 12, 0, 0]
    date_num_in = 51.0
    lat_in = 65.0
    lon_in = 12.0
    azimuth_test = 12.64584
    zenith_test = 92.35858
    azimuth_out = -999.0
    zenith_out = -999.0
    call get_sun_angles(lat_in, lon_in, date_in, date_num_in, 0.0, azimuth_out, zenith_out)
    if (abs(azimuth_out - azimuth_test) > 1.0e-4 .or. abs(zenith_out - zenith_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 32 failed! Routine: global_radiation_sub"
    end if

    ! Test case 33 (global_radiation_sub): Check against expected values Southern hemisphere
    date_in = [1995, 6, 1, 12, 0, 0]
    date_num_in = 51.0
    lat_in = -24.0
    lon_in = 2.0
    azimuth_test = 2.645828
    zenith_test = 176.910
    azimuth_out = -999.0
    zenith_out = -999.0
    call get_sun_angles(lat_in, lon_in, date_in, date_num_in, 0.0, azimuth_out, zenith_out)
    if (abs(azimuth_out - azimuth_test) > 1.0e-4 .or. abs(zenith_out - zenith_test) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 33 failed! Routine: global_radiation_sub"
    end if

    ! Return test results summary
    if (ok) then
        print "(a)", "test_utility_functions: All tests passed."
    else
        print "(a)", "test_utility_functions: One ore more tests failed."
        stop 1
    end if

end program test_time_functions
