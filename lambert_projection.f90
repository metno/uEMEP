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
    
    !k_lambert = log(cos(deg2rad*lat_stand1_lambert)/cos(deg2rad*lat_stand2_lambert))/&
    !  (log(tan(0.25*PI+0.5*deg2rad*lat_stand2_lambert)/tan(0.25*PI+0.5*deg2rad*lat_stand1_lambert)))

     ! lat0_lambert = rad2deg*asin(k_lambert)  

    k = sin(PI/180.*lat0)
    F = earth_radius*cos(PI/180.*lat0) * (tan(PI/4.+PI/360.*lat0)**k) /k
    y0 = F*tan(PI/4.-PI/360.*lat0)**k
    r = F*tan(PI/4.-PI/360.*gb)**k
    x = r*sin(PI/180.*k*(gl-lon0))
    y = y0 - r*cos(PI/180.*k*(gl-lon0))
    
    end subroutine lb2lambert_uEMEP
    
subroutine lb2lambert2_uEMEP(x,y,gl,gb,lat_stand1_lambert,lat_stand2_lambert,lon0,lat0_in)
    
    implicit none
    real, intent(in) ::gl,gb,lon0,lat0_in,lat_stand1_lambert,lat_stand2_lambert
    real, intent(out)::x,y
    real ::r,t
    real ::PI
    real :: earth_radius,k,F,y0
    real deg2rad,rad2deg,k_lambert,lat0_lambert
    real lat0
    
    PI=3.14159265358979323
    earth_radius = 6370000.
    deg2rad=PI/180.
    rad2deg=180./PI
    
    k_lambert = log(cos(deg2rad*lat_stand1_lambert)/cos(deg2rad*lat_stand2_lambert))/&
      (log(tan(0.25*PI+0.5*deg2rad*lat_stand2_lambert)/tan(0.25*PI+0.5*deg2rad*lat_stand1_lambert)))

    lat0_lambert = rad2deg*asin(k_lambert)  
    lat0=lat0_lambert
    !lat0=lat0_in
    k=k_lambert
    !write(*,*) lat0_lambert,lat0_in
    !k = sin(PI/180.*lat0)
    F = earth_radius*cos(PI/180.*lat_stand1_lambert) * (tan(PI/4.+PI/360.*lat_stand1_lambert)**k) /k
    y0 = F*tan(PI/4.-PI/360.*lat0)**k
    r = F*tan(PI/4.-PI/360.*gb)**k
    x = r*sin(PI/180.*k*(gl-lon0))
    y = y0 - r*cos(PI/180.*k*(gl-lon0))
    
end subroutine lb2lambert2_uEMEP