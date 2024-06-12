module mod_lambert_projection
    ! Routines for calling the various possible projections for the uEMEP sub-grid to lat lon

    use uemep_constants, only: pi
    use uEMEP_definitions
    use utility_functions, only: ltm2ll, utm2ll, ll2utm, ll2ltm
    use mod_rdm2ll, only: RDM2LL

    implicit none
    private

    public :: lb2lambert2_uEMEP, LL2PS_spherical, PROJ2LL, LL2LAEA, lambert2lb2_uEMEP, &
        lb2lambert_uEMEP

contains

    subroutine testlambert()
        real :: gl, gb, x, y, lon0, lat0, y0, k, F, earth_radius, lat_stand1, GRIDWIDTH_M
        real :: deg2rad
        GRIDWIDTH_M = 2500.0
        lon0 = 15.0
        lat0 = 63.0
        deg2rad = PI/180.0
        earth_radius = 6371000.0
        lat_stand1 = lat0

        k = sin(PI/180.0*lat0)
        F = earth_radius*cos(PI/180.0*lat_stand1) * tan(PI/4 + PI/360.0*lat_stand1)**k/k
        y0 = F*tan(PI/4 - PI/360.0*lat0)**k

        gl = 15.0
        gb = 63.0
        call lb2lambert(x, y, gl, gb, lon0, y0, k, F)
        write(*,*) 'lon = ', gl, 'lat =', gb
        write(*,*) 'give lambert x = ', x, 'y =', y
        write(*,*) 'lambert i = ', (x)/GRIDWIDTH_M,'j =', y/GRIDWIDTH_M
        write(*,*)
        x = -892442.2
        y = 1220678.0
        call lambert2lb(x, y, gl, gb, lon0, y0, k, F)
        write(*,*) 'Lambert x = ', x, 'y =', y
        write(*,*) 'gives lon = ', gl, 'lat =', gb

        call lb2lambert(x, y, gl, gb, lon0, y0, k, F)
        write(*,*) 'and back to Lambert x = ', x, 'y =', y
    end subroutine testlambert

    subroutine lambert2lb(x, y, gl, gb, lon0, y0, k, F)
        real, intent(in) ::x, y, lon0, y0, k, F
        real, intent(out) ::gl, gb
        real :: r, t
        r = sqrt(x*x + (y0 - y)*(y0 - y))
        t = atan(x/(y0 - y))
        gb = 2.0*180.0/PI*atan((F/r)**(1.0/k)) - 90.0
        gl = lon0 + 180.0/PI*t/k
    end subroutine lambert2lb

    subroutine lb2lambert(x, y, gl, gb, lon0, y0, k, F)
        real, intent(in) :: gl, gb, lon0, y0, k, F
        real, intent(out) :: x, y
        real :: r, dr2, dr
        dr = PI/180.0
        dr2 = PI/360.0
        r = F*tan(PI/4 - dr2*gb)**k
        x = r*sin(dr*k*(gl - lon0))
        y = y0 - r*cos(dr*k*(gl - lon0))
    end subroutine lb2lambert

    subroutine lambert2lb_uEMEP(x, y, gl, gb, lon0, lat0)
        real, intent(in) ::x, y, lon0, lat0
        real, intent(out)::gl, gb
        real ::r, t
        real :: earth_radius, k, F, y0

        earth_radius = 6371000.0

        k = sin(PI/180.0*lat0)
        F = earth_radius*cos(PI/180.0*lat0) * (tan(PI/4.0 + PI/360.0*lat0)**k)/k
        y0 = F*tan(PI/4.0 - PI/360.0*lat0)**k
        r = sqrt(x*x + (y0 - y)*(y0 - y))
        t = atan(x/(y0 - y))
        gb = 2.0*180.0/PI*atan((F/r)**(1.0/k)) - 90.0
        gl = lon0 + 180.0/PI*t/k
    end subroutine lambert2lb_uEMEP

    subroutine lambert2lb2_uEMEP(x, y, gl, gb, projection_attr)
        double precision, intent(in) :: projection_attr(10)
        real, intent(in) :: x, y
        real, intent(out):: gl, gb
        real :: r, t
        real :: earth_radius, k, F, y0
        real :: deg2rad, rad2deg, k_lambert, lat0_lambert
        real :: lat0
        real :: lat_stand1_lambert, lat_stand2_lambert, lon0, lat0_in

        lat_stand1_lambert = projection_attr(1)
        lat_stand2_lambert = projection_attr(2)
        lon0 = projection_attr(3)
        lat0_in = projection_attr(4)
        earth_radius = projection_attr(5)

        deg2rad = PI/180.0
        rad2deg = 180.0/PI

        if (lat_stand1_lambert .eq. lat_stand2_lambert) then
            k_lambert = sin(PI/180.0*lat0_in)
        else
            k_lambert = log(cos(deg2rad*lat_stand1_lambert)/cos(deg2rad*lat_stand2_lambert))/ &
                (log(tan(0.25*PI + 0.5*deg2rad*lat_stand2_lambert)/tan(0.25*PI + 0.5*deg2rad*lat_stand1_lambert)))
        end if

        lat0_lambert = rad2deg*asin(k_lambert)
        lat0 = lat0_in
        k = k_lambert
        F = earth_radius*cos(PI/180.0*lat_stand1_lambert)*(tan(PI/4.0+PI/360.0*lat_stand1_lambert)**k)/k
        y0 = F*tan(PI/4.0 - PI/360.0*lat0)**k
        r = sqrt(x*x + (y0 - y)*(y0 - y))
        t = atan(x/(y0 - y))
        gb = 2*180.0/PI*atan((F/r)**(1.0/k)) - 90.0
        gl = lon0 + 180.0/PI*t/k
    end subroutine lambert2lb2_uEMEP

    subroutine lb2lambert_uEMEP(x, y, gl, gb, lon0, lat0)
        real, intent(in) :: gl, gb, lon0, lat0
        real, intent(out):: x, y
        real :: r
        real :: earth_radius, k, F, y0
        real :: rad2deg

        earth_radius = 6371000.0
        rad2deg = PI/180.0

        k = sin(PI/180.0*lat0)
        F = earth_radius*cos(PI/180.0*lat0) * (tan(PI/4.0 + PI/360.0*lat0)**k)/k
        y0 = F*tan(PI/4.0 - PI/360.0*lat0)**k
        r = F*tan(PI/4.0 - PI/360.0*gb)**k
        x = r*sin(PI/180.0*k*(gl - lon0))
        y = y0 - r*cos(PI/180.0*k*(gl - lon0))
    end subroutine lb2lambert_uEMEP

    subroutine lb2lambert2_uEMEP(x, y, gl, gb, projection_attr)
        double precision, intent(in) :: projection_attr(10)
        real, intent(in) :: gl, gb
        real, intent(out):: x, y
        real :: r
        real :: earth_radius, k, F, y0
        real :: deg2rad, rad2deg, k_lambert, lat0_lambert
        real :: lat0
        real :: lat_stand1_lambert, lat_stand2_lambert, lon0, lat0_in

        lat_stand1_lambert = projection_attr(1)
        lat_stand2_lambert = projection_attr(2)
        lon0 = projection_attr(3)
        lat0_in = projection_attr(4)
        earth_radius = projection_attr(5)
        deg2rad = PI/180.0
        rad2deg = 180.0/PI

        if (lat_stand1_lambert .eq. lat_stand2_lambert) then
            k_lambert = sin(PI/180.0*lat0_in)
        else
            k_lambert = log(cos(deg2rad*lat_stand1_lambert)/cos(deg2rad*lat_stand2_lambert))/ &
                (log(tan(0.25*PI + 0.5*deg2rad*lat_stand2_lambert)/tan(0.25*PI + 0.5*deg2rad*lat_stand1_lambert)))
        end if

        lat0_lambert = rad2deg*asin(k_lambert)
        lat0 = lat0_in
        k = k_lambert
        F = earth_radius*cos(PI/180.0*lat_stand1_lambert)*(tan(PI/4.0 + PI/360.0*lat_stand1_lambert)**k)/k
        y0 = F*tan(PI/4.0 - PI/360.0*lat0)**k
        r = F*tan(PI/4.0 - PI/360.0*gb)**k
        x = r*sin(PI/180.0*k*(gl - lon0))
        y = y0 - r*cos(PI/180.0*k*(gl - lon0))
    end subroutine lb2lambert2_uEMEP

    subroutine LL2LAEA_spherical(x, y, lon_in, lat_in, projection_attr)
        ! https://mathworld.wolfram.com/LambertAzimuthalEqual-AreaProjection.html
        ! grid_mapping_name = lambert_azimuthal_equal_area
        ! Map parameters:
        ! longitude_of_projection_origin
        ! latitude_of_projection_origin
        ! false_easting - This parameter is optional (default is 0)
        ! false_northing - This parameter is optional (default is 0)lat_stand1_lambert=projection_attributes(1)
        double precision, intent(in) :: projection_attr(10)
        real, intent(in) :: lon_in, lat_in
        real, intent(out):: x, y

        ! Local variables
        real :: r
        real :: earth_radius
        real :: deg2rad, rad2deg, k_lambert
        real :: lat0, lat0_in, lon0, lon0_in
        real :: false_easting, false_northing
        real :: lon, lat

        lon0_in = projection_attr(1)
        lat0_in = projection_attr(2)
        false_easting = projection_attr(3)
        false_northing = projection_attr(4)
        earth_radius = projection_attr(5)

        deg2rad = PI/180.0
        rad2deg = 180.0/PI
        r = earth_radius

        lat0 = lat0_in*deg2rad
        lon0 = lon0_in*deg2rad
        lon = lon_in*deg2rad
        lat = lat_in*deg2rad

        k_lambert = sqrt(2.0/(1.0 + sin(lat0)*sin(lat) + cos(lat0)*cos(lat)*cos(lon - lon0)))
        x = false_easting + r*k_lambert*cos(lat)*sin(lon - lon0)
        y = false_northing + r*k_lambert*(cos(lat0)*sin(lat) - sin(lat0)*cos(lat)*cos(lon - lon0))
    end subroutine LL2LAEA_spherical

    subroutine LAEA2LL_spherical(x, y, lon, lat, projection_attr)
        ! https://mathworld.wolfram.com/LambertAzimuthalEqual-AreaProjection.html
        ! grid_mapping_name = lambert_azimuthal_equal_area
        ! Map parameters:
        ! longitude_of_projection_origin
        ! latitude_of_projection_origin
        ! false_easting - This parameter is optional (default is 0)
        ! false_northing - This parameter is optional (default is 0)lat_stand1_lambert=projection_attributes(1)
        double precision, intent(in) :: projection_attr(10)
        real, intent(out) :: lon, lat
        real, intent(in):: x, y

        ! Local variables
        real :: r, rho, c
        real :: earth_radius
        real :: deg2rad, rad2deg
        real :: lat0, lat0_in, lon0, lon0_in
        real :: false_easting, false_northing
        
        lon0_in = projection_attr(1)
        lat0_in = projection_attr(2)
        false_easting = projection_attr(3)
        false_northing = projection_attr(4)
        earth_radius = projection_attr(5)

        deg2rad = PI/180.0
        rad2deg = 180.0/PI
        r = earth_radius

        lat0 = lat0_in*deg2rad
        lon0 = lon0_in*deg2rad

        rho = sqrt((x - false_easting)*(x - false_easting) + (y - false_northing)*(y - false_northing))
        c = 2.0*asin(rho*0.5/r)
        lat = asin(cos(c)*sin(lat0)+(y-false_northing)/rho*sin(c)*cos(lat0))
        lon = lon0 + atan((x - false_easting)*sin(c)/(rho*cos(lat0)*cos(c) - (y - false_northing)*sin(lat0)*sin(c)))
        lat = lat*rad2deg
        lon = lon*rad2deg
    end subroutine LAEA2LL_spherical

    subroutine LL2LAEA(x, y, lon_in, lat_in, projection_attr)
        ! https://epsg.io/3035
        ! grid_mapping_name = lambert_azimuthal_equal_area
        ! Map parameters:
        ! longitude_of_projection_origin
        ! latitude_of_projection_origin
        ! false_easting - This parameter is optional (default is 0)
        ! false_northing - This parameter is optional (default is 0)lat_stand1_lambert=projection_attributes(1)
        double precision, intent(in) :: projection_attr(10)
        real, intent(in) :: lon_in, lat_in
        real, intent(out) :: x, y

        ! Local variables
        real :: semi_major_axis
        real :: deg2rad,rad2deg
        real :: lat0, lat0_in, lon0, lon0_in
        real :: false_easting, false_northing
        real :: lon, lat
        real :: inv_f, f, a, e, q_p, q_0, q, beta, beta_0, R_q, D, B
        
        lon0_in = projection_attr(1)
        lat0_in = projection_attr(2)
        false_easting = projection_attr(3)
        false_northing = projection_attr(4)
        semi_major_axis = projection_attr(5)
        inv_f = projection_attr(6) !flattening

        deg2rad = PI/180.0
        rad2deg = 180.0/PI
        a = semi_major_axis
        f = 1.0/inv_f
        e = sqrt(2.0*f - f*f)

        lat0 = lat0_in*deg2rad
        lon0 = lon0_in*deg2rad
        lon = lon_in*deg2rad
        lat = lat_in*deg2rad

        q_p = (1.0 - e*e)*(1.0/(1.0-e*e) - 1.0/(2.0*e)*log((1.0 - e)/(1.0 + e)))
        q_0 = (1.0 - e*e)*(sin(lat0)/(1.0 - e*e*sin(lat0)**2) - 1.0/(2.0*e)*log((1.0 - e*sin(lat0))/(1.0 + e*sin(lat0))))
        q = (1.0 - e*e)*(sin(lat)/(1.0 - e*e*sin(lat)**2) - 1.0/(2.0*e)*log((1.0 - e*sin(lat))/(1.0 + e*sin(lat))))
        beta_0 = asin(q_0/q_p)
        beta = asin(q/q_p)
        R_q = a*sqrt(q_p/2.0)
        D = a*(cos(lat0)/sqrt(1.0 - e*e*sin(lat0)**2)/(R_q*cos(beta_0)))
        B = R_q*sqrt(2.0/(1.0 + sin(beta_0)*sin(beta) + cos(beta_0)*cos(beta)*cos(lon - lon0)))
        x = false_easting + B*D*cos(beta)*sin(lon - lon0)
        y = false_northing + B/D*(cos(beta_0)*sin(beta) - sin(beta_0)*cos(beta)*cos(lon - lon0))
    end subroutine LL2LAEA

    subroutine LAEA2LL(x, y, lon, lat, projection_attr)
        ! www.epsg.org
        ! grid_mapping_name = lambert_azimuthal_equal_area
        ! Map parameters:
        ! longitude_of_projection_origin
        ! latitude_of_projection_origin
        ! false_easting - This parameter is optional (default is 0)
        ! false_northing - This parameter is optional (default is 0)lat_stand1_lambert=projection_attributes(1)
        double precision, intent(in) :: projection_attr(10)
        real, intent(out) :: lon, lat
        real, intent(in) :: x, y
        real :: rho
        real :: semi_major_axis
        real :: deg2rad,rad2deg
        real :: lat0, lat0_in, lon0, lon0_in
        real :: false_easting, false_northing
        real :: inv_f, f, a, e, q_p, q_0, beta_0, R_q, D
        real :: C, beta_d
        
        lon0_in = projection_attr(1)
        lat0_in = projection_attr(2)
        false_easting = projection_attr(3)
        false_northing = projection_attr(4)
        semi_major_axis = projection_attr(5)
        inv_f = projection_attr(6) ! flattening

        deg2rad = PI/180.0
        rad2deg = 180.0/PI
        a = semi_major_axis
        f = 1.0/inv_f
        e = sqrt(2.0*f - f*f)

        lat0 = lat0_in*deg2rad
        lon0 = lon0_in*deg2rad

        q_p = (1.0 - e*e)*(1.0/(1.0 - e*e) - 1.0/(2.0*e)*log((1.0 - e)/(1.0 + e)))
        q_0 = (1.0 - e*e)*(sin(lat0)/(1.0 - e*e*sin(lat0)**2) - 1.0/(2.0*e)*log((1.0 - e*sin(lat0))/(1.0 + e*sin(lat0))))
        beta_0 = asin(q_0/q_p)
        R_q = a*sqrt(q_p/2.0)
        D = a*(cos(lat0)/sqrt(1.0 - e*e*sin(lat0)**2)/(R_q*cos(beta_0)))
        rho = sqrt((x - false_easting)*(x - false_easting)/D/D+D*D*(y - false_northing)*(y - false_northing))
        C = 2.0*asin(rho/2.0/R_q)
        beta_d = asin(cos(C)*sin(beta_0) + D*(y - false_northing)*sin(C)*cos(beta_0)/rho)

        lon = lon0 + atan2((x - false_easting)*sin(C), &
            (D*rho*cos(beta_0)*cos(C) - D*D*(y - false_northing)*sin(beta_0)*sin(C)))
        lat = beta_d + e**2*(1.0/3.0+31.0/180.0*e**2.0 + 517.0/5040.0*e**4)*sin(2.0*beta_d) &
            + e**4.0*(23.0/360.0+251.0/3780.0*e**2)*sin(4.0*beta_d)+(761.0/45360.0*e**6)*sin(6.0*beta_d)

        lat = lat*rad2deg
        lon = lon*rad2deg
    end subroutine LAEA2LL

    subroutine PROJ2LL(x_in, y_in, lon_out, lat_out, projection_attributes_in, projection_type_in)
        ! Definitions only needed for the identification indexes
        double precision, intent(in) :: projection_attributes_in(10)
        real, intent(in) :: x_in, y_in
        integer, intent(in) :: projection_type_in
        real, intent(out) :: lon_out,lat_out
        integer :: utm_zone_in
        real :: ltm_lon0_in

        if (projection_type_in .eq. RDM_projection_index) then
            call RDM2LL(y_in, x_in, lat_out, lon_out)
        else if (projection_type_in .eq. UTM_projection_index) then
            utm_zone_in = floor(projection_attributes_in(1) + 0.5)
            call utm2ll(1, utm_zone_in, y_in, x_in, lat_out, lon_out)
        else if (projection_type_in .eq. LTM_projection_index) then
            utm_zone_in = floor(projection_attributes_in(1) + 0.5)
            ltm_lon0_in = projection_attributes_in(2)
            call ltm2ll(1, utm_zone_in, ltm_lon0_in, y_in, x_in, lat_out, lon_out)
        else if (projection_type_in .eq. LAEA_projection_index) then
            call LAEA2LL(x_in, y_in, lon_out, lat_out, projection_attributes_in)
        else if (projection_type_in .eq. LCC_projection_index) then
            call lambert2lb2_uEMEP(x_in, y_in, lon_out, lat_out, projection_attributes_in)
        else if (projection_type_in .eq. PS_projection_index) then
            call PS2LL_spherical(x_in, y_in, lon_out, lat_out, projection_attributes_in)
        else if (projection_type_in .eq. LL_projection_index) then
            lon_out = x_in
            lat_out = y_in
        else
            write(unit_logfile, '(A,I0)') 'ERROR: This projection type index is not implemented:', projection_type_in
            stop
        end if
    end subroutine PROJ2LL

    subroutine LL2PROJ(lon_in, lat_in, x_out, y_out, projection_attributes_out, projection_type_out)
        ! Definitions only needed for the identification indexes
        double precision, intent(in) :: projection_attributes_out(10)
        real, intent(out) :: x_out, y_out
        integer, intent(in) :: projection_type_out
        real, intent(in) :: lon_in,lat_in
        integer :: utm_zone_out
        real :: ltm_lon0_out

        if (projection_type_out .eq. RDM_projection_index) then
            write(unit_logfile,'(A)') ' ERROR: Conversion from lon-lat to RDM projection is not implemented"'
            stop
        else if (projection_type_out .eq. UTM_projection_index) then
            utm_zone_out = floor(projection_attributes_out(1) + 0.5)
            call ll2utm(1, utm_zone_out, lat_in, lon_in, y_out, x_out)
        else if (projection_type_out .eq. LTM_projection_index) then
            utm_zone_out = floor(projection_attributes_out(1) + 0.5)
            ltm_lon0_out = projection_attributes_out(2)
            call ll2ltm(1, ltm_lon0_out, lat_in, lon_in, y_out, x_out)  !?????????
        else if (projection_type_out .eq. LAEA_projection_index) then
            call LL2LAEA(x_out, y_out, lon_in, lat_in, projection_attributes_out)
        else if (projection_type_out .eq. LCC_projection_index) then
            call lb2lambert2_uEMEP(x_out, y_out, lon_in, lat_in, projection_attributes_out)
        else if (projection_type_out .eq. PS_projection_index) then
            call LL2PS_spherical(x_out, y_out, lon_in, lat_in, projection_attributes_out)
        else if (projection_type_out .eq. LL_projection_index) then
            x_out = lon_in
            y_out = lat_in
        else
            write(unit_logfile, '(A,I0)') 'ERROR: This projection type index is not implemented:', projection_type_in
            stop
        end if
    end subroutine LL2PROJ

    subroutine LL2PS_spherical(x, y, lon_in, lat_in, projection_attr)
        ! https://mathworld.wolfram.com/StereographicProjection.html
        ! grid_mapping_name = Polar_Stereographic
        ! Map parameters:
        ! longitude_of_projection_origin
        ! latitude_of_projection_origin
        ! false_easting - This parameter is optional (default is 0)
        ! false_northing - This parameter is optional (default is 0)lat_stand1_lambert=projection_attributes(1)
        double precision, intent(in) :: projection_attr(10)
        real, intent(in) ::lon_in,lat_in
        real, intent(out)::x,y
        real ::r
        real :: earth_radius
        real deg2rad,rad2deg,k_ps
        real lat0,lat0_in,lon0,lon0_in
        real false_easting,false_northing
        real scaling
        real lon,lat

        lon0_in = projection_attr(1)
        lat0_in = projection_attr(2)
        false_easting = projection_attr(3)
        false_northing = projection_attr(4)
        earth_radius = projection_attr(5)
        scaling = projection_attr(6)

        deg2rad = PI/180.0
        rad2deg = 180.0/PI
        r = earth_radius

        lat0 = lat0_in*deg2rad
        lon0 = lon0_in*deg2rad
        lon = lon_in*deg2rad
        lat = lat_in*deg2rad

        k_ps = 2.0*r*scaling/(1.0 + sin(lat0)*sin(lat) + cos(lat0)*cos(lat)*cos(lon - lon0))
        x = false_easting + k_ps*cos(lat)*sin(lon - lon0)
        y = false_northing + k_ps*(cos(lat0)*sin(lat) - sin(lat0)*cos(lat)*cos(lon - lon0))
    end subroutine LL2PS_spherical

    subroutine PS2LL_spherical(x, y, lon, lat, projection_attr)
        ! https://mathworld.wolfram.com/StereographicProjection.html
        ! grid_mapping_name = Polar_Stereographic
        ! Map parameters:
        ! longitude_of_projection_origin
        ! latitude_of_projection_origin
        ! false_easting - This parameter is optional (default is 0)
        ! false_northing - This parameter is optional (default is 0)
        ! NOTE scale_factor_at_projection_origin =0.5(1+sin(Standard Parallel))
        double precision, intent(in) :: projection_attr(10)
        real, intent(out) :: lon, lat
        real, intent(in) :: x, y
        real :: r, rho, c
        real :: earth_radius
        real :: deg2rad, rad2deg
        real :: lat0, lat0_in, lon0, lon0_in
        real :: false_easting,false_northing
        real :: scaling
        
        lon0_in = projection_attr(1)
        lat0_in = projection_attr(2)
        false_easting = projection_attr(3)
        false_northing = projection_attr(4)
        earth_radius = projection_attr(5)
        scaling = projection_attr(6)

        deg2rad = PI/180.0
        rad2deg = 180.0/PI
        r = earth_radius

        lat0 = lat0_in*deg2rad
        lon0 = lon0_in*deg2rad

        rho = sqrt((x - false_easting)*(x - false_easting) + (y-false_northing)*(y - false_northing))
        c = 2.0*atan(rho*0.5/r/scaling)
        lat = asin(cos(c)*sin(lat0) + (y - false_northing)/rho*sin(c)*cos(lat0))
        lon = lon0 + atan((x - false_easting)*sin(c)/(rho*cos(lat0)*cos(c) - (y - false_northing)*sin(lat0)*sin(c)))
        lat = lat*rad2deg
        lon = lon*rad2deg
    end subroutine PS2LL_spherical

end module mod_lambert_projection

