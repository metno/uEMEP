module utility_functions

    implicit none
    private

    public :: distrl_modern, distrl_sqr_modern
    public :: nxtdat_modern
    public :: ll2utm_modern, ll2ltm_modern
    public :: utm2ll_modern, ltm2ll_modern

    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: pi = 3.141592653589793

contains

    subroutine distrl_modern(x0, y0, x1, y1, x2, y2, xm, ym, dm, wm)
        !! The subroutine calculates the minimum distance from a given receptor
        !! point to a given line source.
        real, intent(in) :: x0 !! Receptor point x-coordinate
        real, intent(in) :: y0 !! Receptor point y-coordinate
        real, intent(in) :: x1 !! Line source x-coordinate 1
        real, intent(in) :: y1 !! Line source y-coordinate 1
        real, intent(in) :: x2 !! Line source x-coordinate 2
        real, intent(in) :: y2 !! Line source y-coordinate 2
        real, intent(out) :: xm !! Minimum distance x-coordinate 
        real, intent(out) :: ym !! Minimum distance y-coordinate
        real, intent(out) :: dm !! Minimum distance

        ! Local variables
        real :: num, denum, wm

        if (x1 == x2 .and. y1 == y2) then
            wm = 0.5
        else
            num = (x0 - x1)*(x2 - x1) + (y0 - y1)*(y2 - y1)
            denum = (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1)
            wm = num/denum
        end if

        wm = min(wm, 1.0)
        wm = max(wm, 0.0)

        xm = (1.0 - wm)*x1 + wm*x2
        ym = (1.0 - wm)*y1 + wm*y2
        dm = sqrt((x0 - xm)*(x0 - xm) + (y0 - ym)*(y0 - ym))
    end subroutine distrl_modern

    subroutine distrl_sqr_modern(x0, y0, x1, y1, x2, y2, xm, ym, dm_sqr, wm)
        real, intent(in) :: x0 !! Receptor point x-coordinate
        real, intent(in) :: y0 !! Receptor point y-coordinate
        real, intent(in) :: x1 !! Line source x-coordinate 1
        real, intent(in) :: y1 !! Line source y-coordinate 1
        real, intent(in) :: x2 !! Line source x-coordinate 2
        real, intent(in) :: y2 !! Line source y-coordinate 2
        real, intent(out) :: xm !! Minimum distance x-coordinate
        real, intent(out) :: ym !! Minimum distance y-coordinate
        real, intent(out) :: dm_sqr !! Minimum distance

        ! Local variables
        real :: num, denum, wm

        if (x1 == x2 .and. y1 == y2) then
            wm = 0.5
        else
            num = (x0 - x1)*(x2 - x1) + (y0 - y1)*(y2 - y1)
            denum = (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1)
            wm = num/denum
        end if

        wm = min(wm, 1.0)
        wm = max(wm, 0.0)

        xm = (1.0 - wm)*x1 + wm*x2
        ym = (1.0 - wm)*y1 + wm*y2
        dm_sqr = (x0 - xm)*(x0 - xm) + (y0 - ym)*(y0 - ym)
    end subroutine distrl_sqr_modern

    subroutine nxtdat_modern(un, leof)
        !! The subroutine prepares for reading the next uncommented line of data from file
        integer, intent(inout) :: un
        logical, intent(out) :: leof

        ! Local variables
        character(len=256) :: txtstr
        integer :: io

        ! If file unit is non-positive then just return
        if (un .le. 0) return

        leof = .false.

        ! Read lines from file
        do
            read(un,"(a)", iostat=io) txtstr
            if (io /= 0) then
                leof = .true.
                exit
            else
                if (txtstr(1:1) == "*" .or. txtstr(1:1) == "{" .or. txtstr(1:1) == "#") then
                    cycle
                else
                    backspace(un)
                    exit
                end if
            end if
        end do
    end subroutine nxtdat_modern

    subroutine ll2utm_modern(iutm, isone_in, lat, lon, utmn, utme)
        integer, intent(in) :: iutm !! UTM coordinate system indicator
        integer, intent(in) :: isone_in !! UTM zone input?
        real, intent(in) :: lat !! Latitude in decimal degrees
        real, intent(in) :: lon !! Longitude in decimal degress
        real, intent(out) :: utmn !! UTM east-coordinate (y) (meter from west border)
        real, intent(out) :: utme !! UTM north-coordinate (x) (meter from equator)

        ! Local variables and parameters
        real :: isone ! UTM zone
        real(dp), parameter :: deast = 500000.0 ! East movement UTM
        real(dp), parameter :: scale = 0.9996 ! Scale UTM
        real(dp) :: a ! Big semiaxis
        real(dp) :: f ! Flattening
        real(dp) :: latv ! Scaled latitude
        real(dp) :: lonv ! Scaled longitude
        real(dp) :: lon0 ! Central meridian of UTM zone
        real(dp) :: b, e, e2, m, n ! Intermediate values

        select case(iutm)
        case(1) ! UTM WGS84 EUREF 89 (AirQUIS)
            a = 6378137.0
            f = 1.0/298.257222101
        case(2) ! UTM WGS84 OLD
            a = 6378137.0
            f = 1.0/298.257223563
        case(3) ! UTM ED50
            a = 6378388.0
            f = 1.0/297.0
        case default
            print *, "ERROR: Unknown UTM coordinate system indictor, iutm: ", iutm
            stop 1
        end select

        ! Scale coordinates
        isone = abs(isone_in)
        latv = lat*pi/180.0
        lon0 = real(isone - 30.0)*6.0 - 3.0
        lonv = (lon - lon0)*pi/180.0

        ! Calculate some intermediate quantities
        e2 = f*(2.0 - f)
        n = a/dsqrt(1.0 - e2*dsin(latv)*dsin(latv))
        e = dsqrt(e2*dcos(latv)*dcos(latv)/(1.0 - e2))
        m = n/(1.0 + e*e)
        b = ((1.0 - f/2.0 + f*f/16.0 + f*f*f/32.0)*latv &
            - (3.0*f/4.0 - 3.0*f*f*f/128.0)*dsin(2.0*latv) &
            + (15.0*f*f/64.0 + 15.0*f*f*f/128.0)*dsin(4.0*latv) &
            - (35.0*f*f*f/384.0)*dsin(6.0*latv))*a
        
        ! Calculate UTM north coordinate
        utmn = (b + lonv*lonv*n*dsin(latv)*dcos(latv)/2.0 &
            + lonv*lonv*lonv*lonv*n*dsin(latv) &
            * dcos(latv)*dcos(latv)*dcos(latv) &
            * (5.0 - dtan(latv)*dtan(latv) + 9.0*e*e + 4.0*e*e*e*e) &
            / 24.0 + lonv*lonv*lonv*lonv*lonv*lonv*n*dsin(latv) &
            * dcos(latv)*dcos(latv)*dcos(latv)*dcos(latv)*dcos(latv) &
            * (61.0 - 58.0*dtan(latv)*dtan(latv) &
            + dtan(latv)*dtan(latv)*dtan(latv)*dtan(latv)) &
            / 720.0)*scale

        if (lat < 0.0) then
            utmn = utme + 10000000.0
        end if

        ! Calculate UTM east coordinate
        utme = (lonv*n*dcos(latv) &
            + lonv*lonv*lonv*n*dcos(latv)*dcos(latv)*dcos(latv) &
            * (1.0 - dtan(latv)*dtan(latv) + e*e)/6.0 &
            + lonv*lonv*lonv*lonv*lonv*n &
            * dcos(latv)*dcos(latv)*dcos(latv)*dcos(latv)*dcos(latv) &
            * (5.0 - 18.0*dtan(latv)*dtan(latv) &
            + dtan(latv)*dtan(latv)*dtan(latv)*dtan(latv))/120.0) &
            * scale + deast
    end subroutine ll2utm_modern

    subroutine ll2ltm_modern(iutm, lon0, lat, lon, utmn, utme)
        !! Local lon version (of ll2utm) without zone, so just typical Transverse
        !! Mecantor (Local Transverse Mecantor)
        integer, intent(in) :: iutm !! UTM coordinate system indicator
        real, intent(in) :: lon0 !! Central meridian of UTM zone
        real, intent(in) :: lat !! Latitude in decimal degrees
        real, intent(in) :: lon !! Longitude in decimal degrees
        real, intent(out) :: utmn !! UTM north-coordinate (x) (meter from equator)
        real, intent(out) :: utme !! UTM  east-coordinate (y) (meter from west border)

        ! Local variables and paramters
        real(dp), parameter :: deast = 500000.0 ! East movement UTM
        real(dp), parameter :: scale = 0.9996 ! Scale UTM
        real(dp) :: a ! Big semiaxis
        real(dp) :: f ! Flattening
        real(dp) :: latv ! Scaled latitude
        real(dp) :: lonv ! Scaled longitude
        real(dp) :: b, e, e2, m, n ! Intermediate values

        select case(iutm)
        case(1) ! UTM WGS84 EUREF 89 (AirQUIS)
            a = 6378137.0
            f = 1.0/298.257222101
        case(2) ! UTM WGS84 OLD
            a = 6378137.0
            f = 1.0/298.257223563
        case(3) ! UTM ED50
            a = 6378388.0
            f = 1.0/297.0
        case default
            print *, "ERROR: Unknown UTM coordinate system indictor, iutm: ", iutm
            stop 1
        end select

        ! Scale coordinates
        latv = lat*pi/180.0
        lonv = (lon - lon0)*pi/180.0

        ! Calculate some intermediate quantities
        e2 = f*(2.0 - f)
        n = a/dsqrt(1.0 - e2*dsin(latv)*dsin(latv))
        e = dsqrt(e2*dcos(latv)*dcos(latv)/(1.0 - e2))
        m = n/(1.0 + e*e)
        b = ((1.0 - f/2.0 + f*f/16.0 + f*f*f/32.0)*latv &
            - (3.0*f/4.0 - 3.0*f*f*f/128.0)*dsin(2.0*latv) &
            + (15.0*f*f/64.0 + 15.0*f*f*f/128.0)*dsin(4.0*latv) &
            - (35.0*f*f*f/384.0)*dsin(6.0*latv))*a

        ! Calculate UTM north coordinate

        utmn = (b + lonv*lonv*n*dsin(latv)*dcos(latv)/2.0 &
            + lonv*lonv*lonv*lonv*n*dsin(latv) &
            * dcos(latv)*dcos(latv)*dcos(latv) &
            * (5.0 - dtan(latv)*dtan(latv) + 9.0*e*e + 4.0*e*e*e*e) &
            / 24.0 + lonv*lonv*lonv*lonv*lonv*lonv*n*dsin(latv) &
            * dcos(latv)*dcos(latv)*dcos(latv)*dcos(latv)*dcos(latv) &
            * (61.0 - 58.0*dtan(latv)*dtan(latv) &
            + dtan(latv)*dtan(latv)*dtan(latv)*dtan(latv)) &
            / 720.0)*scale
      
        if (lat < 0.0) then
            utmn = utmn + 10000000.0
        end if

        ! Calculate UTM east coordinate
        utme = (lonv*n*dcos(latv) &
            + lonv*lonv*lonv*n*dcos(latv)*dcos(latv)*dcos(latv) &
            * (1.0 - dtan(latv)*dtan(latv) + e*e)/6.0 &
            + lonv*lonv*lonv*lonv*lonv*n &
            * dcos(latv)*dcos(latv)*dcos(latv)*dcos(latv)*dcos(latv) &
            * (5.0 - 18.0*dtan(latv)*dtan(latv) &
            + dtan(latv)*dtan(latv)*dtan(latv)*dtan(latv))/120.0) &
            * scale + deast
    end subroutine ll2ltm_modern

    subroutine utm2ll_modern(iutm, isone_in, utmn_in, utme, lat, lon)
        !! The subroutine converts UTM north- and east-coordinates to latitude
        !! and longitude
        integer, intent(in) :: iutm !! UTM coordinate system indicator
        integer, intent(in) :: isone_in !! UTM zone
        real, intent(in) :: utmn_in !! UTM north-coordinate (X) (meter from equator)
        real, intent(in) :: utme !! UTM  east-coordinate (Y) (meter from west border)
        real, intent(out) :: lat !! Latitude in decimal degrees
        real, intent(out) :: lon !! Longitude in decimal degrees

        ! Local variables and parameters
        integer :: isone
        real :: utmn ! UTM north-coordinate
        real(dp), parameter :: deast = 500000.0
        real(dp), parameter :: scale = 0.9996
        real(dp) :: a ! Big semiaxis
        real(dp) :: f ! Flattening
        real(dp) :: la0 ! Tangeringsmeridian
        real(dp) :: x ! Scaled north-coordinate
        real(dp) :: y ! Scaled east-coordinate
        real(dp) :: bb0, e, e2, fi, m, n ! Intermediate value

        select case(iutm)
        case(1) ! UTM WGS84 EUREF 89 (AirQUIS)
            a = 6378137.0
            f = 1.0/298.257222101
        case(2) ! UTM WGS84 OLD
            a = 6378137.0
            f = 1.0/298.257223563
        case(3) ! UTM ED50
            a = 6378388.0
            f = 1.0/297.0
        case default
            print *, "ERROR: Unknown UTM coordinate system indictor, iutm: ", iutm
            stop 1
        end select

        ! Adjust for southern hemisphere, specified by negative isone_in
        if (isone_in < 0) then
            utmn = utmn_in - 10000000.0
        else
            utmn = utmn_in
        end if
        isone = abs(isone_in)

        ! Scale coordinates
        x = utmn/scale
        y = (utme - deast)/scale
        la0 = real(isone - 30)*6.0 - 3.0

        ! Calculate some intermediate quantities
        e2 = f*(2.0 - f)
        bb0 = (1.0 - f/2.0 + f*f/16.0 + f*f*f/32.0)*a
        fi = x/bb0 &
            + (3.0*f/4.0 + 3.0*f*f/8.0 + 21.0*f*f*f/256.0)*dsin(2.0*x/bb0) &
            + (21.0*f*f/64.0 + 21.0*f*f*f/64.0)*dsin(4.0*x/bb0) &
            + (151.0*f*f*f/768.0)*dsin(6.0*x/bb0)
        n = a/dsqrt(1.0 - e2*dsin(fi)*dsin(fi))
        e = dsqrt(e2*dcos(fi)*dcos(fi)/(1.0 - e2))
        m = n/(1.0 + e*e)

        ! Calculate latitude and longitude in radians
        lat = fi - (y*y*dtan(fi))/(2.0*m*n) &
            + (y*y*y*y*dtan(fi))/(24.0*m*n*n*n) &
            * (5.0 + 3.0*dtan(fi)*dtan(fi) + e*e &
            - 9.0*e*e*dtan(fi)*dtan(fi) - 4.0*e*e*e*e) &
            - (y*y*y*y*y*y*dtan(fi))/(720.0*m*n*n*n*n*n) &
            * (61.0 + 90.0*dtan(fi)*dtan(fi) &
            + 45.0*dtan(fi)*dtan(fi)*dtan(fi)*dtan(fi))

        lon = y/(n*dcos(fi)) &
            - (y*y*y*(1.0 + 2.0*dtan(fi)*dtan(fi) + e*e)) &
            / (6.0*n*n*n*dcos(fi)) &
            + (y*y*y*y*y*(5.0 + 28.0*dtan(fi)*dtan(fi) &
            + 24.0*dtan(fi)*dtan(fi)*dtan(fi)*dtan(fi))) &
            / (120.0*n*n*n*n*n*dcos(fi)) + la0*pi/180.0

        ! Convert from radians to degrees
        lat = lat*180.0/pi
        lon = lon*180.0/pi
    end subroutine utm2ll_modern

    subroutine ltm2ll_modern(iutm, isone_in, la0, utmn_in, utme, lat, lon)
        !! Local Transverse Mecantor version of utm2ll
        integer, intent(in) :: iutm !! UTM coordinate system indicator
        integer, intent(in) :: isone_in !! UTM zone
        real, intent(in) :: la0 !! Tangeringsmeridian
        real, intent(in) :: utmn_in !! UTM north-coordinate (X) (meter from equator)
        real, intent(in) :: utme !! UTM  east-coordinate (Y) (meter from west border)
        real, intent(out) :: lat !! Latitude in decimal degrees
        real, intent(out) :: lon !! Longitude in decimal degrees

        ! Local variables and parameters
        integer :: isone
        real :: utmn ! UTM north-coordinate
        real(dp), parameter :: deast = 500000.0
        real(dp), parameter :: scale = 0.9996
        real(dp) :: a ! Big semiaxis
        real(dp) :: f ! Flattening
        real(dp) :: x ! Scaled north-coordinate
        real(dp) :: y ! Scaled east-coordinate
        real(dp) :: bb0, e, e2, fi, m, n ! Intermediate values

        select case(iutm)
        case(1) ! UTM WGS84 EUREF 89 (AirQUIS)
            a = 6378137.0
            f = 1.0/298.257222101
        case(2) ! UTM WGS84 OLD
            a = 6378137.0
            f = 1.0/298.257223563
        case(3) ! UTM ED50
            a = 6378388.0
            f = 1.0/297.0
        case default
            print *, "ERROR: Unknown UTM coordinate system indictor, iutm: ", iutm
            stop 1
        end select
         
        ! djust for Southern Hemisphere, specified by negative isone_in
        if (isone_in < 0) then
            utmn = utmn_in - 10000000.0
        else
            utmn = utmn_in
        end if
        isone = abs(isone_in)

        ! Scale coordinates
        x = utmn/scale
        y = (utme - deast)/scale

        ! Calculate some intermediate quantities
        e2 = f*(2.0 - f)
        bb0 = (1.0 - f/2.0 + f*f/16.0 + f*f*f/32.0)*a
        fi  = x/bb0 &
            + (3.0*f/4.0 + 3.0*f*f/8.0 + 21.0*f*f*f/256.0)*dsin(2.0*x/bb0) &
            + (21.0*f*f/64.0 + 21.0*f*f*f/64.0)*dsin(4.0*x/bb0) &
            + (151.0*f*f*f/768.0)*dsin(6.0*x/bb0)
        n = a/dsqrt(1.0 - e2*dsin(fi)*dsin(fi))
        e = dsqrt(e2*dcos(fi)*dcos(fi)/(1.0 - e2))
        m = n/(1.0 + e*e)
        
        ! Calculate latitude and longitude in radians
        lat = fi - (y*y*dtan(fi))/(2.0*m*n) &
            + (y*y*y*y*dtan(fi))/(24.0*m*n*n*n) &
            * (5.0 + 3.0*dtan(fi)*dtan(fi) + e*e &
            - 9.0*e*e*dtan(fi)*dtan(fi) - 4.0*e*e*e*e) &
            - (y*y*y*y*y*y*dtan(fi))/(720.0*m*n*n*n*n*n) &
            * (61.0 + 90.0*dtan(fi)*dtan(fi) &
            + 45.0*dtan(fi)*dtan(fi)*dtan(fi)*dtan(fi))

        lon = y/(n*dcos(fi)) &
            - (y*y*y*(1.0 + 2.0*dtan(fi)*dtan(fi) + e*e)) &
            / (6.0*n*n*n*dcos(fi)) &
            + (y*y*y*y*y*(5.0 + 28.0*dtan(fi)*dtan(fi) &
            + 24.0*dtan(fi)*dtan(fi)*dtan(fi)*dtan(fi))) &
            / (120.0*n*n*n*n*n*dcos(fi)) + la0*pi/180.0

        ! Convert from radians to degrees
        lat = lat*180.0/pi
        lon = lon*180.0/PI
    end subroutine ltm2ll_modern

end module utility_functions