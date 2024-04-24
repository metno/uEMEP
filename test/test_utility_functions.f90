program test_utility_functions

    use utility_functions

    implicit none

    logical :: ok = .true.

    ! distrl and distrl_sqr
    real :: xm0, ym0, dm0, wm0, xm1, ym1, dm1, wm1

    ! nxtdat
    character(len=256) :: txtstr0, txtstr1
    character(len=2048) :: path, relpath
    integer :: unit0
    logical :: flag0, flag1

    ! ll2utm and ll2ltm and utm2ll
    integer :: iutm, isone_in
    real :: lat, lon, utmn0, utme0, utmn1, utme1, lon0, utmn_in, utme, lat0, lat1, lon1, la0

    interface
        SUBROUTINE DISTRL(X0,Y0,X1,Y1,X2,Y2,XM,YM,DM,WM)
            REAL :: X0,Y0,X1,Y1,X2,Y2,XM,YM,DM,WM
        END SUBROUTINE DISTRL
        SUBROUTINE DISTRL_SQR(X0,Y0,X1,Y1,X2,Y2,XM,YM,DM_SQR,WM)
            REAL :: X0,Y0,X1,Y1,X2,Y2,XM,YM,DM_SQR,WM
        END SUBROUTINE DISTRL_SQR
    end interface

    interface
        SUBROUTINE NXTDAT(UN,LEOF)
            INTEGER :: UN
            LOGICAL :: LEOF
        END SUBROUTINE NXTDAT
    end interface

    interface
        SUBROUTINE LL2UTM(IUTM,ISONE_IN,LAT,LON,UTMN,UTME)
            INTEGER :: IUTM,ISONE_IN
            REAL :: LAT,LON,UTMN,UTME
        END SUBROUTINE LL2UTM
        SUBROUTINE LL2LTM(IUTM,LON0,LAT,LON,UTMN,UTME)
            INTEGER :: IUTM
            REAL :: LON0,LAT,LON,UTMN,UTME
        END SUBROUTINE LL2LTM
        SUBROUTINE UTM2LL(ISONE_IN,UTMN_IN,UTME,LAT,LON)
            INTEGER :: ISONE_IN
            REAL :: UTMN_IN,UTME,LAT,LON
        END SUBROUTINE UTM2LL
        SUBROUTINE LTM2LL(ISONE_IN,LA0,UTMN_IN,UTME,LAT,LON)
            INTEGER :: ISONE_IN
            REAL :: LA0, UTMN_IN,UTME,LAT,LON
        END SUBROUTINE LTM2LL
    end interface

    ! DISTRL vs distrl_modern

    ! Test case 1: Check that old and new version outputs the same values
    !              when "the line is point"
    call DISTRL(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, xm0, ym0, dm0, wm0)
    call distrl_modern(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, xm1, ym1, dm1, wm1)
    if (abs(xm0 - xm1) > 1.0e-12 .or. abs(ym0 - ym1) > 1.0e-12 .or. &
        abs(dm0 - dm1) > 1.0e-12 .or. abs(wm0 - wm1) > 1.0e-12) then
        ok = .false.
        print "(a)", "Test case 1 failed!"
    end if

    ! Test case 2: Check that old and new version outputs the same values
    !              when the line is increasing and point is below and right of the line
    call DISTRL(1.0, 1.0, 5.0, 5.0, 5.0, 2.0, xm0, ym0, dm0, wm0)
    call distrl_modern(1.0, 1.0, 5.0, 5.0, 5.0, 2.0, xm1, ym1, dm1, wm1)
    if (abs(xm0 - xm1) > 1.0e-12) then
        ok = .false.
        print "(a)", "Test case 2 failed!"
    end if

    ! Test case 3: Check that old and new version outputs the same values
    !              when the line is increasing and point is above and left of the line
    call DISTRL(1.0, 1.0, 5.0, 5.0, 2.0, 5.0, xm0, ym0, dm0, wm0)
    call distrl_modern(1.0, 1.0, 5.0, 5.0, 2.0, 5.0, xm1, ym1, dm1, wm1)
    if (abs(xm0 - xm1) > 1.0e-12) then
        ok = .false.
        print "(a)", "Test case 3 failed!"
    end if

    ! Test case 4: Check that old and new version outputs the same values
    !              when the line is decreasing and point is below the line
    call DISTRL(1.0, 5.0, 5.0, 1.0, 2.0, 2.0, xm0, ym0, dm0, wm0)
    call distrl_modern(1.0, 5.0, 5.0, 1.0, 2.0, 2.0, xm1, ym1, dm1, wm1)
    if (abs(xm0 - xm1) > 1.0e-12) then
        ok = .false.
        print "(a)", "Test case 4 failed!"
    end if

    ! Test case 5: Check that old and new version outputs the same values
    !              when the line is decreasing and point is above the line
    call DISTRL(1.0, 5.0, 5.0, 1.0, 4.0, 4.0, xm0, ym0, dm0, wm0)
    call distrl_modern(1.0, 5.0, 5.0, 1.0, 4.0, 4.0, xm1, ym1, dm1, wm1)
    if (abs(xm0 - xm1) > 1.0e-12) then
        ok = .false.
        print "(a)", "Test case 5 failed!"
    end if

    ! Test case 6: Check that old and new version outputs the same values
    !              when the point is on the line
    call DISTRL(1.0, 1.0, 3.0, 3.0, 2.0, 2.0, xm0, ym0, dm0, wm0)
    call distrl_modern(1.0, 1.0, 3.0, 3.0, 2.0, 2.0, xm1, ym1, dm1, wm1)
    if (abs(xm0 - xm1) > 1.0e-12) then
        ok = .false.
        print "(a)", "Test case 6 failed!"
    end if

    ! DISTRL_SQR vs distrl_sqr_modern

    ! Test case 7: Check that old and new version outputs the same values
    !              when "the line is point"
    call DISTRL_SQR(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, xm0, ym0, dm0, wm0)
    call distrl_sqr_modern(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, xm1, ym1, dm1, wm1)
    if (abs(xm0 - xm1) > 1.0e-12 .or. abs(ym0 - ym1) > 1.0e-12 .or. &
        abs(dm0 - dm1) > 1.0e-12 .or. abs(wm0 - wm1) > 1.0e-12) then
        ok = .false.
        print "(a)", "Test case 7 failed!"
    end if

    ! Test case 8: Check that old and new version outputs the same values
    !              when the line is increasing and point is below and right of the line
    call DISTRL_SQR(1.0, 1.0, 5.0, 5.0, 5.0, 2.0, xm0, ym0, dm0, wm0)
    call distrl_sqr_modern(1.0, 1.0, 5.0, 5.0, 5.0, 2.0, xm1, ym1, dm1, wm1)
    if (abs(xm0 - xm1) > 1.0e-12) then
        ok = .false.
        print "(a)", "Test case 8 failed!"
    end if

    ! Test case 9: Check that old and new version outputs the same values
    !              when the line is increasing and point is above and left of the line
    call DISTRL_SQR(1.0, 1.0, 5.0, 5.0, 2.0, 5.0, xm0, ym0, dm0, wm0)
    call distrl_sqr_modern(1.0, 1.0, 5.0, 5.0, 2.0, 5.0, xm1, ym1, dm1, wm1)
    if (abs(xm0 - xm1) > 1.0e-12) then
        ok = .false.
        print "(a)", "Test case 9 failed!"
    end if

    ! Test case 10: Check that old and new version outputs the same values
    !              when the line is decreasing and point is below the line
    call DISTRL_SQR(1.0, 5.0, 5.0, 1.0, 2.0, 2.0, xm0, ym0, dm0, wm0)
    call distrl_sqr_modern(1.0, 5.0, 5.0, 1.0, 2.0, 2.0, xm1, ym1, dm1, wm1)
    if (abs(xm0 - xm1) > 1.0e-12) then
        ok = .false.
        print "(a)", "Test case 10 failed!"
    end if

    ! Test case 11: Check that old and new version outputs the same values
    !              when the line is decreasing and point is above the line
    call DISTRL_SQR(1.0, 5.0, 5.0, 1.0, 4.0, 4.0, xm0, ym0, dm0, wm0)
    call distrl_sqr_modern(1.0, 5.0, 5.0, 1.0, 4.0, 4.0, xm1, ym1, dm1, wm1)
    if (abs(xm0 - xm1) > 1.0e-12) then
        ok = .false.
        print "(a)", "Test case 11 failed!"
    end if

    ! Test case 12: Check that old and new version outputs the same values
    !              when the point is on the line
    call DISTRL_SQR(1.0, 1.0, 3.0, 3.0, 2.0, 2.0, xm0, ym0, dm0, wm0)
    call distrl_sqr_modern(1.0, 1.0, 3.0, 3.0, 2.0, 2.0, xm1, ym1, dm1, wm1)
    if (abs(xm0 - xm1) > 1.0e-12) then
        ok = .false.
        print "(a)", "Test case 12 failed!"
    end if

    ! NXTDAT vs nextdat_modern

    ! Test case 13: Check that non-positive file unit returns false
    unit0 = 0
    call NXTDAT(unit0, flag0)
    call nxtdat_modern(unit0, flag1)
    if ( .not. flag0 == flag1) then
        ok = .false.
        print "(a)", "Test case 13 failed!"
    end if

    ! Test case 14 + 15: File with commented lines
    unit0 = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"

    open(unit=unit0, file=trim(relpath)//"test_commented_file.txt", status="old")
    call NXTDAT(unit0, flag0)
    read(unit0,"(a)") txtstr0
    close(unit0)

    open(unit=unit0, file=trim(relpath)//"test_commented_file.txt", status="old")
    call nxtdat_modern(unit0, flag1)
    read(unit0,"(a)") txtstr1
    close(unit0)

    if ( .not. flag0 == flag1 .or. flag0 .eqv. .true. .or. flag1 .eqv. .true.) then
        ok = .false.
        print "(a)", "Test case 14 failed!"
    end if

    if ( .not. txtstr0 == txtstr1 .or. trim(txtstr0) /= "non-commented line 4" .or. trim(txtstr1) /= "non-commented line 4") then
        ok = .false.
        print "(a)", "Test case 15 failed!"
    end if

    ! Test case 16: File with only commented lines
    unit0 = 20
    call get_environment_variable("PWD", path)
    relpath = trim(path)//"/../test/test_data/"

    open(unit=unit0, file=trim(relpath)//"test_only_commented_file.txt", status="old")
    call NXTDAT(unit0, flag0)
    close(unit0)

    open(unit=unit0, file=trim(relpath)//"test_only_commented_file.txt", status="old")
    call nxtdat_modern(unit0, flag1)
    close(unit0)

    if ( .not. flag0 == flag1 .or. flag0 .eqv. .false. .or. flag1 .eqv. .false.) then
        ok = .false.
        print "(a)", "Test case 16 failed"
    end if

    ! LL2UTM vs ll2utm_modern

    isone_in = 32
    lat = 55.0
    lon = 12.0

    ! Test case 17-19: UTM WGS84 EUREF 89 (AirQUIS)
    iutm = 1
    call LL2UTM(iutm, isone_in, lat, lon, utmn0, utme0)
    call ll2utm_modern(iutm, isone_in, lat, lon, utmn1, utme1)

    if (abs(utme0 - utme1) > 1.0e-4 .or. abs(utmn0 - utmn1) > 1.0e-4) then
        print *, abs(utme0 - utme1), abs(utmn0 - utmn1)
        ok = .false.
        print "(a)", "Test case 17 failed"
    end if

    if (abs(utmn0 - 6098908.0) > 1.0e-4 .or. abs(utme0 - 691875.6) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 18 failed!"
    end if

    if (abs(utmn1 - 6098908.0) > 1.0e-4 .or. abs(utme1 - 691875.6) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 19 failed!"
    end if

    ! Test case 20-22: UTM WGS84 OLD
    iutm = 2
    call LL2UTM(iutm, isone_in, lat, lon, utmn0, utme0)
    call ll2utm_modern(iutm, isone_in, lat, lon, utmn1, utme1)

    if (abs(utme0 - utme1) > 1.0e-4 .or. abs(utmn0 - utmn1) > 1.0e-4) then
        print *, abs(utme0 - utme1), abs(utmn0 - utmn1)
        ok = .false.
        print "(a)", "Test case 20 failed"
    end if

    if (abs(utmn0 - 6098908.0) > 1.0e-4 .or. abs(utme0 - 691875.6) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 21 failed!"
    end if

    if (abs(utmn1 - 6098908.0) > 1.0e-4 .or. abs(utme1 - 691875.6) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 22 failed!"
    end if

    ! Test case 23-25: UTM ED50
    iutm = 3
    call LL2UTM(iutm, isone_in, lat, lon, utmn0, utme0)
    call ll2utm_modern(iutm, isone_in, lat, lon, utmn1, utme1)

    if (abs(utme0 - utme1) > 1.0e-4 .or. abs(utmn0 - utmn1) > 1.0e-4) then
        print *, abs(utme0 - utme1), abs(utmn0 - utmn1)
        ok = .false.
        print "(a)", "Test case 23 failed"
    end if

    if (abs(utmn0 - 6099040.5) > 1.0e-4 .or. abs(utme0 - 691885.0) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 24 failed!"
    end if

    if (abs(utmn1 - 6099040.5) > 1.0e-4 .or. abs(utme1 - 691885.0) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 25 failed!"
    end if

    ! LL2LTM vs ll2ltm_modern

    ! Test case 26-28: UTM WGS84 EUREF 89 (AirQUIS) with central meridian
    iutm = 1
    lon0 = 3.0
    lat = 55.0
    lon = 12.0
    call LL2LTM(iutm, lon0, lat, lon, utmn0, utme0)
    call ll2ltm_modern(iutm, lon0, lat, lon, utmn1, utme1)
    
    if (abs(utme0 - utme1) > 1.0e-4 .or. abs(utmn0 - utmn1) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 26 failed!"
    end if

    if (abs(utme0 - 1074900.0) > 1.0e-4 .or. abs(utmn0 - 6131905.0) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 27 failed!"
    end if

    if (abs(utme1 - 1074900.0) > 1.0e-4 .or. abs(utmn1 - 6131905.0) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 28 failed!"
    end if

    ! Test case 29-31: UTM WGS84 OLD without central meridian
    iutm = 2
    lon0 = 0.0
    lat = 55.0
    lon = 12.0
    call LL2LTM(iutm, lon0, lat, lon, utmn0, utme0)
    call ll2ltm_modern(iutm, lon0, lat, lon, utmn1, utme1)

    if (abs(utme0 - utme1) > 1.0e-4 .or. abs(utmn0 - utmn1) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 29 failed!"
    end if

    if (abs(utme0 - 1265670.375) > 1.0e-4 .or. abs(utmn0 - 6160873.500) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 30 failed!"
    end if

    if (abs(utme1 - 1265670.375) > 1.0e-4 .or. abs(utmn1 - 6160873.500) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 31 failed!"
    end if

    ! Test case 32-34: UTM ED50 with negative latitude
    iutm = 3
    lon0 = 6.0
    lat = -55.0
    lon = 12.0
    call LL2LTM(iutm, lon0, lat, lon, utmn0, utme0)
    call ll2ltm_modern(iutm, lon0, lat, lon, utmn1, utme1)

    if (abs(utme0 - utme1) > 1.0e-4 .or. abs(utmn0 - utmn1) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 32 failed!"
    end if

    if (abs(utme0 - 883589.375) > 1.0e-4 .or. abs(utmn0 - 3888598.500) > 1.0e-4) then
        print "(2f15.3)", utme0, utmn0
        ok = .false.
        print "(a)", "Test case 33 failed!"
    end if

    if (abs(utme1 - 883589.375) > 1.0e-4 .or. abs(utmn1 - 3888598.500) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 34 failed!"
    end if

    ! UTM2LL vs utm2ll_modern
    
    ! Test case 35-37: UTM WGS84 EUREF 89 (AirQUIS) in Northern Hemisphere
    iutm = 1
    isone_in = 32
    utmn_in = 6098908.0
    utme = 691875.6

    call UTM2LL(isone_in, utmn_in, utme, lat0, lon0)
    call utm2ll_modern(iutm, isone_in, utmn_in, utme, lat1, lon1)

    if (abs(lat0 - lat1) > 1.0e-4 .or. abs(lon0 - lon1) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 35 failed!"
    end if

    if (abs(lat0 - 55.0) > 1.0e-4 .or. abs(lon0 - 12.0) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 36 failed!"
    end if

    if (abs(lat1 - 55.0) > 1.0e-4 .or. abs(lon1 - 12.0) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 37 failed!"
    end if

    ! Test case 38-40: UTM WGS84 EUREF 89 (AirQUIS) in Southern Hemisphere
    iutm = 1
    isone_in = -32
    utmn_in = 3901092.0
    utme = 691875.625

    call UTM2LL(isone_in, utmn_in, utme, lat0, lon0)
    call utm2ll_modern(iutm, isone_in, utmn_in, utme, lat1, lon1)

    if (abs(lat0 - lat1) > 1.0e-4 .or. abs(lon0 - lon1) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 38 failed!"
    end if

    if (abs(lat0 + 55.0) > 1.0e-4 .or. abs(lon0 - 12.0) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 39 failed!"
    end if

    if (abs(lat1 + 55.0) > 1.0e-4 .or. abs(lon1 - 12.0) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 40 failed!"
    end if

    ! LTM2LL vs ltm2ll_modern

    ! Test case 41-43: Local TM with positive ISONE in Northern Hemisphere
    iutm = 1
    isone_in = 32
    la0 = 9.0
    utmn_in = 6098908.0
    utme = 691875.625
    call LTM2LL(isone_in, la0, utmn_in, utme, lat0, lon0)
    call ltm2ll_modern(iutm, isone_in, la0, utmn_in, utme, lat1, lon1)

    if (abs(lat0 - lat1) > 1.0e-4 .or. abs(lon0 - lon1) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 41 failed!"
    end if

    if (abs(lat0 - 55.0) > 1.0e-4 .or. abs(lon0 - 12.0) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 42 failed!"
    end if

    if (abs(lat1 - 55.0) > 1.0e-4 .or. abs(lon1 - 12.0) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 43 failed!"
    end if

    ! Test case 44-46: Local TM with negative ISONE in Southern Hemisphere
    iutm = 1
    isone_in = -32
    la0 = 9.0
    utmn_in = 6677424.0
    utme = 789409.688
    call LTM2LL(isone_in, la0, utmn_in, utme, lat0, lon0)
    call ltm2ll_modern(iutm, isone_in, la0, utmn_in, utme, lat1, lon1)

    if (abs(lat0 - lat1) > 1.0e-4 .or. abs(lon0 - lon1) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 44 failed!"
    end if

    if (abs(lat0 + 30.0) > 1.0e-4 .or. abs(lon0 - 12.0) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 45 failed!"
    end if

    if (abs(lat1 + 30.0) > 1.0e-4 .or. abs(lon1 - 12.0) > 1.0e-4) then
        ok = .false.
        print "(a)", "Test case 46 failed!"
    end if

    if (ok) then
        print "(a)", "test_utility_functions: All tests passed."
    else
        print "(a)", "test_utility_functions: One ore more tests failed."
        stop 1
    end if

    

end program test_utility_functions