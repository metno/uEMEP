!compile with
!ifort -r8 Lambert.f90
subroutine testlambert
  implicit none
  real ::gl,gb,x,y,lon0,lat0,y0,k,F,earth_radius,lat_stand1,GRIDWIDTH_M
  real ::PI,x1,y1,k_lambert,F_lambert,deg2rad,lat_stand2
  PI=3.14159265358979323
  GRIDWIDTH_M =2500.0
  lon0=15.0
  lat0=63.0
  deg2rad=PI/180.
  earth_radius = 6371000.
  lat_stand1 = lat0

  k = sin(PI/180.*lat0)
  F = earth_radius*cos(PI/180.*lat_stand1) * tan(PI/4+PI/360.*lat_stand1)**k /k
  y0 = F*tan(PI/4-PI/360.*lat0)**k

  gl =    15.0
  gb = 63.0
  call lb2lambert(x,y,gl,gb,lon0,y0,k,F)
  write(*,*)'lon = ',gl,'lat =',gb
  write(*,*)'give lambert x = ',x,'y =',y
  write(*,*)' lambert i = ',(x)/GRIDWIDTH_M,'j =',y/GRIDWIDTH_M
  write(*,*)
  x=-892442.2
  y=1220678.
  call  lambert2lb(x,y,gl,gb,lon0,y0,k,F)
  write(*,*)'Lambert x = ',x,'y =',y
  write(*,*)'gives lon = ',gl,'lat =',gb

  call lb2lambert(x,y,gl,gb,lon0,y0,k,F)
  write(*,*)'and back to Lambert x = ',x,'y =',y

end subroutine testlambert

    subroutine lambert2lb(x,y,gl,gb,lon0,y0,k,F)
    implicit none
    real, intent(in) ::x,y,lon0,y0,k,F
    real, intent(out)::gl,gb
    real ::r,t
    real ::PI
    PI=3.14159265358979323
    r = sqrt(x*x+(y0-y)*(y0-y))
    t = atan(x/(y0-y))
    gb = 2*180./PI*atan((F/r)**(1.0/k))-90.0
    gl = lon0 + 180./PI*t/k
    end subroutine lambert2lb

    subroutine lb2lambert(x,y,gl,gb,lon0,y0,k,F)
    implicit none
    real, intent(in) ::gl,gb,lon0,y0,k,F
    real, intent(out)::x,y
    real ::r,t,dr2,dr
    real ::PI
    PI=3.14159265358979323
    dr=PI/180.
    dr2=PI/360.
    r = F*tan(PI/4-dr2*gb)**k
    x = r*sin(dr*k*(gl-lon0))
    y = y0 - r*cos(dr*k*(gl-lon0))
    end subroutine lb2lambert

    subroutine lambert2lb_uEMEP(x,y,gl,gb,lon0,lat0)
  
    implicit none
    real, intent(in) ::x,y,lon0,lat0
    real, intent(out)::gl,gb
    real ::r,t
    real ::PI
    real :: earth_radius,k,F,y0
    
    PI=3.14159265358979323
    earth_radius = 6371000.

    k = sin(PI/180.*lat0)
    F = earth_radius*cos(PI/180.*lat0) * (tan(PI/4+PI/360.*lat0)**k) /k
    y0 = F*tan(PI/4-PI/360.*lat0)**k   
    r = sqrt(x*x+(y0-y)*(y0-y))
    t = atan(x/(y0-y))
    gb = 2*180./PI*atan((F/r)**(1.0/k))-90.0
    gl = lon0 + 180./PI*t/k
    
    end subroutine lambert2lb_uEMEP

    subroutine lambert2lb2_uEMEP(x,y,gl,gb,projection_attributes)
  
    implicit none
    double precision, intent(in) :: projection_attributes(10)
    real, intent(in) ::x,y
    real, intent(out)::gl,gb
    real ::r,t
    real ::PI
    real :: earth_radius,k,F,y0
    real deg2rad,rad2deg,k_lambert,lat0_lambert
    real lat0
    real lat_stand1_lambert,lat_stand2_lambert,lon0,lat0_in

    
    lat_stand1_lambert=projection_attributes(1)
    lat_stand2_lambert=projection_attributes(2)
    lon0=projection_attributes(3)
    lat0_in=projection_attributes(4)
    earth_radius = projection_attributes(5)

    PI=3.14159265358979323
    deg2rad=PI/180.
    rad2deg=180./PI

    if (lat_stand1_lambert.eq.lat_stand2_lambert) then
        k_lambert=sin(PI/180.*lat0_in)
    else
        k_lambert = log(cos(deg2rad*lat_stand1_lambert)/cos(deg2rad*lat_stand2_lambert))/&
        (log(tan(0.25*PI+0.5*deg2rad*lat_stand2_lambert)/tan(0.25*PI+0.5*deg2rad*lat_stand1_lambert)))
    endif

    lat0_lambert = rad2deg*asin(k_lambert)  
    lat0=lat0_in
    k=k_lambert

    F = earth_radius*cos(PI/180.*lat_stand1_lambert) * (tan(PI/4.+PI/360.*lat_stand1_lambert)**k) /k
    !k = sin(PI/180.*lat0)
    !F = earth_radius*cos(PI/180.*lat0) * (tan(PI/4+PI/360.*lat0)**k) /k
    y0 = F*tan(PI/4.-PI/360.*lat0)**k   
    r = sqrt(x*x+(y0-y)*(y0-y))
    t = atan(x/(y0-y))
    gb = 2*180./PI*atan((F/r)**(1.0/k))-90.0
    gl = lon0 + 180./PI*t/k
    
    end subroutine lambert2lb2_uEMEP

    subroutine lb2lambert_uEMEP(x,y,gl,gb,lon0,lat0)
    
    implicit none
    real, intent(in) ::gl,gb,lon0,lat0
    real, intent(out)::x,y
    real ::r,t
    real ::PI
    real :: earth_radius,k,F,y0
    real rad2deg
    
    PI=3.14159265358979323
    earth_radius = 6371000.
    rad2deg=PI/180.
    
    k = sin(PI/180.*lat0)
    F = earth_radius*cos(PI/180.*lat0) * (tan(PI/4.+PI/360.*lat0)**k) /k
    y0 = F*tan(PI/4.-PI/360.*lat0)**k
    r = F*tan(PI/4.-PI/360.*gb)**k
    x = r*sin(PI/180.*k*(gl-lon0))
    y = y0 - r*cos(PI/180.*k*(gl-lon0))
    
    end subroutine lb2lambert_uEMEP
    
    subroutine lb2lambert2_uEMEP(x,y,gl,gb,projection_attributes)
    
    implicit none
    double precision, intent(in) :: projection_attributes(10)
    real, intent(in) ::gl,gb
    real, intent(out)::x,y
    real ::r,t
    real ::PI
    real :: earth_radius,k,F,y0
    real deg2rad,rad2deg,k_lambert,lat0_lambert
    real lat0
    real lat_stand1_lambert,lat_stand2_lambert,lon0,lat0_in
    
    lat_stand1_lambert=projection_attributes(1)
    lat_stand2_lambert=projection_attributes(2)
    lon0=projection_attributes(3)
    lat0_in=projection_attributes(4)
    earth_radius = projection_attributes(5)
    PI=3.14159265358979323
    deg2rad=PI/180.
    rad2deg=180./PI
    
    if (lat_stand1_lambert.eq.lat_stand2_lambert) then
        k_lambert=sin(PI/180.*lat0_in)
    else
        k_lambert = log(cos(deg2rad*lat_stand1_lambert)/cos(deg2rad*lat_stand2_lambert))/&
        (log(tan(0.25*PI+0.5*deg2rad*lat_stand2_lambert)/tan(0.25*PI+0.5*deg2rad*lat_stand1_lambert)))
    endif
    
    lat0_lambert = rad2deg*asin(k_lambert)  
    !lat0=lat0_lambert
    lat0=lat0_in
    k=k_lambert
    !write(*,*) lat_stand1_lambert,lat_stand2_lambert,k,lat0_lambert,lat0_in
    !k = sin(PI/180.*lat0)
    F = earth_radius*cos(PI/180.*lat_stand1_lambert) * (tan(PI/4.+PI/360.*lat_stand1_lambert)**k) /k
    y0 = F*tan(PI/4.-PI/360.*lat0)**k
    r = F*tan(PI/4.-PI/360.*gb)**k
    x = r*sin(PI/180.*k*(gl-lon0))
    y = y0 - r*cos(PI/180.*k*(gl-lon0))
    
    end subroutine lb2lambert2_uEMEP
    
    subroutine LL2LAEA(x,y,lon_in,lat_in,projection_attributes)
    !https://mathworld.wolfram.com/LambertAzimuthalEqual-AreaProjection.html
    
    implicit none
    double precision, intent(in) :: projection_attributes(10)
    real, intent(in) ::lon_in,lat_in
    real, intent(out)::x,y
    real ::r
    real ::PI
    real :: earth_radius
    real deg2rad,rad2deg,k_lambert
    real lat0,lat0_in,lon0,lon0_in
    real false_easting,false_northing
    real lon,lat
    
    !grid_mapping_name = lambert_azimuthal_equal_area
    !Map parameters:
    !longitude_of_projection_origin
    !latitude_of_projection_origin
    !false_easting - This parameter is optional (default is 0)
    !false_northing - This parameter is optional (default is 0)lat_stand1_lambert=projection_attributes(1)
    lon0_in=projection_attributes(1)
    lat0_in=projection_attributes(2)
    false_easting=projection_attributes(3)
    false_northing=projection_attributes(4) 
    earth_radius = projection_attributes(5)
    
    PI=3.14159265358979323
    deg2rad=PI/180.
    rad2deg=180./PI
    r=earth_radius
    
    lat0=lat0_in*deg2rad
    lon0=lon0_in*deg2rad
    lon=lon_in*deg2rad
    lat=lat_in*deg2rad
    
    k_lambert=sqrt(2./(1.+sin(lat0)*sin(lat)+cos(lat0)*cos(lat)*cos(lon-lon0)))
    x=false_easting+r*k_lambert*cos(lat)*sin(lon-lon0)
    y=false_northing+r*k_lambert*(cos(lat0)*sin(lat)-sin(lat0)*cos(lat)*cos(lon-lon0))
    
    end subroutine LL2LAEA

    subroutine LAEA2LL(x,y,lon,lat,projection_attributes)
    !https://mathworld.wolfram.com/LambertAzimuthalEqual-AreaProjection.html
     
    implicit none
    double precision, intent(in) :: projection_attributes(10)
    real, intent(out) :: lon,lat
    real, intent(in)::x,y
    real :: r,rho,c
    real :: PI
    real :: earth_radius
    real deg2rad,rad2deg
    real lat0,lat0_in,lon0,lon0_in
    real false_easting,false_northing
    
    !grid_mapping_name = lambert_azimuthal_equal_area
    !Map parameters:
    !longitude_of_projection_origin
    !latitude_of_projection_origin
    !false_easting - This parameter is optional (default is 0)
    !false_northing - This parameter is optional (default is 0)lat_stand1_lambert=projection_attributes(1)
    lon0_in=projection_attributes(1)
    lat0_in=projection_attributes(2)
    false_easting=projection_attributes(3)
    false_northing=projection_attributes(4) 
    earth_radius = projection_attributes(5)
    
    PI=3.14159265358979323
    deg2rad=PI/180.
    rad2deg=180./PI
    r=earth_radius
    
    lat0=lat0_in*deg2rad
    lon0=lon0_in*deg2rad
    
    rho=sqrt((x-false_easting)*(x-false_easting)+(y-false_northing)*(y-false_northing))
    c=2*asin(rho*0.5/r)
    lat=asin(cos(c)*sin(lat0)+(y-false_northing)/rho*sin(c)*cos(lat0))
    lon=lon0+atan((x-false_easting)*sin(c)/(rho*cos(lat0)*cos(c)-(y-false_northing)*sin(lat0)*sin(c)))
    lat=lat*rad2deg
    lon=lon*rad2deg
    
    end subroutine LAEA2LL
    
    !Routine for calling the various possible projections for the uEMEP sub-grid to lat lon
    subroutine PROJ2LL(x_in,y_in,lon_out,lat_out,projection_attributes_in,projection_type_in)

    !Definitions only needed for the identification indexes
    use uEMEP_definitions

    implicit none

    double precision, intent(in) :: projection_attributes_in(10)
    real, intent(in) :: x_in,y_in
    integer, intent(in) :: projection_type_in
    real, intent(out) :: lon_out,lat_out
    integer :: utm_zone_in
    
    if (projection_type.eq.RDM_projection_index) then

        call RDM2LL(y_in,x_in,lat_out,lon_out)

    elseif (projection_type.eq.UTM_projection_index) then

        utm_zone_in=floor(projection_attributes_in(1)+.5)
        call UTM2LL(utm_zone_in,y_in,x_in,lat_out,lon_out)

    elseif (projection_type.eq.LAEA_projection_index) then

        call LAEA2LL(x_in,y_in,lon_out,lat_out,projection_attributes_in)

    elseif (projection_type.eq.LCC_projection_index) then

        call lambert2lb2_uEMEP(x_in,y_in,lon_out,lat_out,projection_attributes_in)

    elseif (projection_type.eq.LL_projection_index) then
    
        lon_out=x_in
        lat_out=y_in
        
    endif

    end subroutine PROJ2LL
