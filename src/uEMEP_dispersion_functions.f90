module dispersion_functions

    implicit none
    private

    public :: gauss_plume_cartesian_sigma_func, gauss_plume_cartesian_sigma_integral_func, &
        gauss_plume_second_order_rotated_reflected_func, &
        gauss_plume_second_order_rotated_reflected_integral_func, &
        gauss_plume_second_order_rotated_integral_func

contains
    
!==========================================================================
!   uEMEP model gauss_plume_second_order_rotated_reflected_func
!   Rotationally symetric Gaussian plume function to second order
!==========================================================================

    function gauss_plume_second_order_rotated_reflected_func(r,z,ay,by,az,bz,sig_y_0,sig_z_0,z_s,z_pbl)

    implicit none
    real r,z,ay,by,az,bz,sig_y_0,sig_z_0,z_s,z_pbl
    real gauss_plume_second_order_rotated_reflected_func
    real sig_th,sig_z,B,c
    real order_1,order_2
    real z_loop(6)
    integer n_loop,k
    real :: correction=2.
    real pi
    parameter (pi=3.141592)
    
    !Corrected for the B**2 falut in the taylor expansion and for the fact that the integral was only half a circle. 20.08.2019
    r=max(1.,r)
    order_1=1.
    order_2=1.
    
    sig_th=(sig_y_0+ay*(exp(by*log(r))))/r
    sig_z=sig_z_0+az*(exp(bz*log(r)))
    !write(*,*) sig_z,sig_th*r,sig_z_0,sig_y_0
    B=-(sig_th**2)*(bz*(sig_z-sig_z_0)/r/sig_th+by*(r*sig_th-sig_y_0)/sig_z)
    !write(*,*) B
    
        if (z_s.gt.z_pbl.or.z_s+sig_z.lt.z_pbl/3.) then
            n_loop=2
            z_loop(1)=z_s;z_loop(2)=-z_s
        else
            n_loop=5
            z_loop(1)=z_s;z_loop(2)=-z_s;z_loop(3)=2.*z_pbl-z_s;z_loop(4)=2.*z_pbl+z_s;z_loop(5)=-2.*z_pbl+z_s
        endif
        
        if (sig_z.gt.0.9*z_pbl) then
            if (B.gt.-1.) then
                !c=1./(2.*pi*sqrt(2.*pi)*r*sig_z*sqrt(1.+B))*tanh(2/sqrt(pi)*pi/(2.*sqrt(2.))/sig_th*sqrt(1.+B))*(exp((-(z-h)**2)/2./sig_z**2)+exp((-(z+h)**2)/2./sig_z**2))
                c=1./(z_pbl*2.*pi*r*sqrt(1.+B))*erf(correction*pi/(2.*sqrt(2.))/sig_th*sqrt(1.+B))
            else
                c=correction/(2.*sqrt(2*pi)*sig_th*r)*(1-order_1*pi**2*(1.+B)*correction**2/(24*sig_th**2)+order_2*pi**4*((1.+B)**2*correction**4/(640.*sig_th**4)))
            endif
        else
            c=0.
            do k=1,n_loop
            if (B.gt.-1.) then
                !c=1./(2.*pi*sqrt(2.*pi)*r*sig_z*sqrt(1.+B))*tanh(2/sqrt(pi)*pi/(2.*sqrt(2.))/sig_th*sqrt(1.+B))*(exp((-(z-h)**2)/2./sig_z**2)+exp((-(z+h)**2)/2./sig_z**2))
                c=c+1./(2.*pi*sqrt(2.*pi)*r*sig_z*sqrt(1.+B))*erf(correction*pi/(2.*sqrt(2.))/sig_th*sqrt(1.+B))*(exp((-(z-z_loop(k))**2)/2./sig_z**2))
            else
                c=c+correction/(4.*pi*sig_th*r*sig_z)*(1-order_1*pi**2*(1.+B)*correction**2/(24*sig_th**2)+order_2*pi**4*((1.+B)**2*correction**4/(640.*sig_th**4)))*(exp((-(z-z_loop(k))**2)/2./sig_z**2))
                !Perhaps also a correction in the higher orders but must calculate that again
            endif
            !if (r.eq.1.) write(*,*) k,r,c,1./(4.*pi*sig_y_0*sig_z_0)
            enddo
        endif
    !write(*,*) c,z_pbl

        !Original
   ! if (B.gt.-1.) then
    !    !c=1./(2.*pi*sqrt(2.*pi)*r*sig_z*sqrt(1.+B))*tanh(2/sqrt(pi)*pi/(2.*sqrt(2.))/sig_th*sqrt(1.+B))*(exp((-(z-h)**2)/2./sig_z**2)+exp((-(z+h)**2)/2./sig_z**2))
    !    c=1./(2.*pi*sqrt(2.*pi)*r*sig_z*sqrt(1.+B))*erf(pi/(2.*sqrt(2.))/sig_th*sqrt(1.+B))*(exp((-(z-h)**2)/2./sig_z**2)+exp((-(z+h)**2)/2./sig_z**2))
    !else
    !    c=1./(4.*pi*sig_th*r*sig_z)*(1-order_1*pi**2*(1.+B)/(24*sig_th**2)+order_2*pi**4*((1.+B**2)/(640.*sig_th**4)))*(exp((-(z-h)**2)/2./sig_z**2)+exp((-(z+h)**2)/2./sig_z**2))
    !endif
    
    gauss_plume_second_order_rotated_reflected_func=c
    
    end function gauss_plume_second_order_rotated_reflected_func
    
!==========================================================================
!   uEMEP model gauss_plume_second_order_rotated_func
!   Rotationally symetric Gaussian plume function to second order with boundary layer reflection
!==========================================================================

    function gauss_plume_second_order_rotated_func(r,z,ay,by,az,bz,sig_y_0,sig_z_0,h)

    implicit none
    real r,z,ay,by,az,bz,sig_y_0,sig_z_0,h
    real gauss_plume_second_order_rotated_func
    real sig_th,sig_z,B,c
    real order_1,order_2
    real :: correction=2.
    real pi
    parameter (pi=3.141592)
    
    !Corrected for the B**2 falut in the taylor expansion and for the fact that the integral was only half a circle. 20.08.2019
    r=max(0.001,r)
    order_1=1.
    order_2=1.
    
    sig_th=(sig_y_0+ay*(exp(by*log(r))))/r
    sig_z=sig_z_0+az*(exp(bz*log(r)))
    !write(*,*) sig_z,sig_th*r,sig_z_0,sig_y_0
    B=-(sig_th**2)*(bz*(sig_z-sig_z_0)/r/sig_th+by*(r*sig_th-sig_y_0)/sig_z)
    !write(*,*) B
    if (B.gt.-1.) then
        !c=1./(2.*pi*sqrt(2.*pi)*r*sig_z*sqrt(1.+B))*tanh(2/sqrt(pi)*pi/(2.*sqrt(2.))/sig_th*sqrt(1.+B))*(exp((-(z-h)**2)/2./sig_z**2)+exp((-(z+h)**2)/2./sig_z**2))
        c=1./(2.*pi*sqrt(2.*pi)*r*sig_z*sqrt(1.+B))*erf(pi/(2.*sqrt(2.))/sig_th*sqrt(1.+B)*correction)*(exp((-(z-h)**2)/2./sig_z**2)+exp((-(z+h)**2)/2./sig_z**2))
    else
        c=correction/(4.*pi*sig_th*r*sig_z)*(1-order_1*pi**2*(1.+B)*correction**2/(24*sig_th**2)+order_2*pi**4*((1.+B)**2*correction**4/(640.*sig_th**4)))*(exp((-(z-h)**2)/2./sig_z**2)+exp((-(z+h)**2)/2./sig_z**2))
    endif
    
    gauss_plume_second_order_rotated_func=c

    end function gauss_plume_second_order_rotated_func
    
    
!==========================================================================
!   uEMEP model gauss_plume_second_order_rotated_integral_func
!   Rotationally symetric Gaussian plume function to second order vertically integrated from H1 to H2
!==========================================================================

    function gauss_plume_second_order_rotated_integral_func(r,z,ay,by,az,bz,sig_y_0,sig_z_0,h,H1,H2)

    implicit none
    real r,z,ay,by,az,bz,sig_y_0,sig_z_0,h,H1,H2
    real gauss_plume_second_order_rotated_integral_func
    real sig_th,sig_z,B,c_y_int,c_z_int
    real order_1,order_2
    real :: correction=2.
    real pi
    parameter (pi=3.141592)
    
    !Corrected for the B**2 falut in the taylor expansion and for the fact that the integral was only half a circle. 20.08.2019
    !Still need to implement reflections
    r=max(0.001,r)
    order_1=1.
    order_2=1.
    
    sig_th=(sig_y_0+ay*(exp(by*log(r))))/r
    sig_z=sig_z_0+az*(exp(bz*log(r)))

    B=-(sig_th**2)*(bz*(sig_z-sig_z_0)/r/sig_th+by*(r*sig_th-sig_y_0)/sig_z)

    c_z_int=sqrt(pi/2.)*sig_z*(erf((H2-h)/sqrt(2.)/sig_z)-erf((H1-h)/sqrt(2.)/sig_z)+erf((H2+h)/sqrt(2.)/sig_z)-erf((H1+h)/sqrt(2.)/sig_z))/(H2-H1)

    if (B.gt.-1.) then
        !c_int=1./(2.*pi*sqrt(2.*pi)*r*sig_z*sqrt(1.+B))*tanh(2/sqrt(pi)*pi/(2.*sqrt(2.))/sig_th*sqrt(1.+B))
        c_y_int=1./(2.*pi*sqrt(2.*pi)*r*sig_z*sqrt(1.+B))*erf(pi/(2.*sqrt(2.))/sig_th*sqrt(1.+B)*correction)
    else
        c_y_int=correction/(4.*pi*sig_th*r*sig_z)*(1-order_1*pi**2*(1.+B)*correction**2/(24*sig_th**2)+order_2*pi**4*((1.+B)**2*correction**4/(640.*sig_th**4)))
    endif

    gauss_plume_second_order_rotated_integral_func=c_y_int*c_z_int
    
    end function gauss_plume_second_order_rotated_integral_func

!==========================================================================
!   uEMEP model gauss_plume_second_order_rotated_reflected_integral_func
!   Rotationally symetric Gaussian plume function to second order vertically integrated from H1 to H2
!==========================================================================

    function gauss_plume_second_order_rotated_reflected_integral_func(r,z,ay,by,az,bz,sig_y_0,sig_z_0,z_s,z_pbl,H1,H2)

    implicit none
    real r,z,ay,by,az,bz,sig_y_0,sig_z_0,z_s,z_pbl,H1,H2
    real gauss_plume_second_order_rotated_reflected_integral_func
    real sig_th,sig_z,B,c_y_int,c_z_int
    real order_1,order_2
    real :: correction=2.
    integer k,n_loop
    real z_loop(6)
    real pi
    parameter (pi=3.141592)
    
    !Corrected for the B**2 falut in the taylor expansion and for the fact that the integral was only half a circle. 20.08.2019
    !Still need to implement reflections
    r=max(0.001,r)
    order_1=1.
    order_2=1.
    
    sig_th=(sig_y_0+ay*(exp(by*log(r))))/r
    sig_z=sig_z_0+az*(exp(bz*log(r)))
  
    B=-(sig_th**2)*(bz*(sig_z-sig_z_0)/r/sig_th+by*(r*sig_th-sig_y_0)/sig_z)
    
        if (z_s.gt.z_pbl.or.z_s+sig_z.lt.z_pbl/3.) then
            n_loop=2
            z_loop(1)=z_s;z_loop(2)=-z_s
        else
            n_loop=5
            z_loop(1)=z_s;z_loop(2)=-z_s;z_loop(3)=2.*z_pbl-z_s;z_loop(4)=2.*z_pbl+z_s;z_loop(5)=-2.*z_pbl+z_s
        endif

        if (sig_z.gt.0.9*z_pbl) then
            c_z_int=1./z_pbl
        else
            c_z_int=0.
            do k=1,n_loop
                c_z_int=c_z_int+sqrt(pi/2.)*sig_z*(erf((H2-z_loop(k))/sqrt(2.)/sig_z)-erf((H1-z_loop(k))/sqrt(2.)/sig_z))/(H2-H1)
            enddo
        endif       

        if (B.gt.-1.) then
            c_y_int=1./(2.*pi*sqrt(2.*pi)*r*sig_z*sqrt(1.+B))*erf(pi/(2.*sqrt(2.))/sig_th*sqrt(1.+B)*correction)
        else
            c_y_int=correction/(4.*pi*sig_th*r*sig_z)*(1-order_1*pi**2*(1.+B)*correction**2/(24*sig_th**2)+order_2*pi**4*((1.+B)**2*correction**4/(640.*sig_th**4)))
        endif
        !write(*,*) c_y_int,c_z_int

        gauss_plume_second_order_rotated_reflected_integral_func=c_y_int*c_z_int
    
    end function gauss_plume_second_order_rotated_reflected_integral_func

!==========================================================================
!   uEMEP model gauss_plume_cartesian_func
!   Cartesian Gaussian plume function
!==========================================================================

    function gauss_plume_cartesian_func(x_s,y_s,z_s,cos_val,sin_val,x_r,y_r,z_r,ay,by,az,bz,sig_y_0,sig_z_0,delta)

    implicit none
    real x_s,y_s,z_s,x_r,y_r,z_r
    real ay,by,az,bz,sig_y_0,sig_z_0,delta
    real gauss_plume_cartesian_func
    real sig_y,sig_z,x,y
    real cos_val,sin_val
    real pi,sig_limit
    parameter (pi=3.141592,sig_limit=3.)
    
    !r=sqrt((x_s-x_r)**2+(y_s-y_r)**2)
    !if (abs(u_s).lt.001) u_s=0.001
    !th=atan(v_s/u_s)
    !if (u_s.lt.0) th=th+pi
    !cos_val=cos(th)
    !sin_val=sin(th)
    x=(x_r-x_s)*cos_val+(y_r-y_s)*sin_val
    y=-(x_r-x_s)*sin_val+(y_r-y_s)*cos_val
    
    gauss_plume_cartesian_func=0.
    sig_y=sig_y_0+ay*exp(by*log(x))+x*abs(delta)
    if (x.ge.0.and.abs(y).lt.sig_y*sig_limit) then
        sig_z=sig_z_0+az*exp(bz*log(x))

        !write(*,*)sig_y_0,sig_y,sig_z_0,sig_z
        !gauss_plume_cartesian_func=1./(2.*pi*sig_y*sig_z)*exp((-y**2)/2./sig_y**2) &
        !    *(exp((-(z_r-z_s)**2)/2./sig_z**2)+exp((-(z_r+z_s)**2)/2./sig_z**2))
        gauss_plume_cartesian_func=1./(2.*pi*sig_y*sig_z)*exp(-y*y/2./sig_y/sig_y) &
            *(exp(-(z_r-z_s)*(z_r-z_s)/2./sig_z/sig_z)+exp(-(z_r+z_s)*(z_r+z_s)/2./sig_z/sig_z))
        
        
    endif
    
    end function gauss_plume_cartesian_func

!==========================================================================
!   uEMEP model gauss_plume_cartesian_integral_func
!   Cartesian Gaussian plume function
!==========================================================================

    function gauss_plume_cartesian_integral_func(x_s,y_s,z_s,cos_val,sin_val,x_r,y_r,z_r,ay,by,az,bz,sig_y_0,sig_z_0,H1,H2,delta)

    implicit none
    real x_s,y_s,z_s,x_r,y_r,z_r
    real ay,by,az,bz,sig_y_0,sig_z_0,H1,H2,delta
    real gauss_plume_cartesian_integral_func
    real sig_y,sig_z,x,y
    real cos_val,sin_val
    real pi,sig_limit
    parameter (pi=3.141592,sig_limit=3.)
    
    !r=sqrt((x_s-x_r)**2+(y_s-y_r)**2)
    !if (abs(u_s).lt.001) u_s=0.001
    !th=atan(v_s/u_s)
    !if (u_s.lt.0) th=th+pi
    !cos_val=cos(th)
    !sin_val=sin(th)
    x=(x_r-x_s)*cos_val+(y_r-y_s)*sin_val
    y=-(x_r-x_s)*sin_val+(y_r-y_s)*cos_val
    
    gauss_plume_cartesian_integral_func=0.
    sig_y=sig_y_0+ay*exp(by*log(x))+x*abs(delta)
    if (x.ge.0.and.abs(y).lt.sig_y*sig_limit) then
        sig_z=sig_z_0+az*exp(bz*log(x))

        !gauss_plume_cartesian_integral_func=1./(2.*pi*sig_y)*exp((-y**2)/2./sig_y**2) &
        !    *sqrt(pi/2.)*(erf((z_s-H1)/sqrt(2.)/sig_z)-erf((z_s-H2)/sqrt(2.)/sig_z)+erf((z_s+H2)/sqrt(2.)/sig_z)-erf((z_s+H1)/sqrt(2.)/sig_z))/(H2-H1)
        gauss_plume_cartesian_integral_func=1./(2.*pi*sig_y)*exp((-y*y)/2./(sig_y*sig_y)) &
            *sqrt(pi/2.)*(erf((z_s-H1)/sqrt(2.)/sig_z)-erf((z_s-H2)/sqrt(2.)/sig_z)+erf((z_s+H2)/sqrt(2.)/sig_z)-erf((z_s+H1)/sqrt(2.)/sig_z))/(H2-H1)
    endif
    
    end function gauss_plume_cartesian_integral_func

!==========================================================================
!   uEMEP model gauss_plume_cartesian_trajectory_func
!   Cartesian Gaussian plume function that does not calculate direction but uses distance x and perpendicular distance y as input.
!   These are precalculated
!==========================================================================

    !function gauss_plume_cartesian_func(x_s,y_s,z_s,u_s,v_s,x_r,y_r,z_r,ay,by,az,bz,sig_y_0,sig_z_0)
    function gauss_plume_cartesian_trajectory_func(x,y,z_s,z_r,ay,by,az,bz,sig_y_0,sig_z_0,delta)

    implicit none
    real x,y,z_s,z_r
    real ay,by,az,bz,sig_y_0,sig_z_0,delta
    real gauss_plume_cartesian_trajectory_func
    real sig_y,sig_z
    real pi,sig_limit
    parameter (pi=3.141592,sig_limit=3.)
    
    !r=sqrt((x_s-x_r)**2+(y_s-y_r)**2)
    !if (abs(u_s).lt.001) u_s=0.001
    !th=atan(v_s/u_s)
    !if (u_s.lt.0) th=th+pi
    !cos_val=cos(th)
    !sin_val=sin(th)
    !x=(x_r-x_s)*cos_val+(y_r-y_s)*sin_val
    !y=-(x_r-x_s)*sin_val+(y_r-y_s)*cos_val
    
    gauss_plume_cartesian_trajectory_func=0.
    sig_y=sig_y_0+ay*exp(by*log(x))+x*abs(delta)
    if (x.ge.0.and.abs(y).lt.sig_y*sig_limit) then
        sig_z=sig_z_0+az*exp(bz*log(x))

        gauss_plume_cartesian_trajectory_func=1./(2.*pi*sig_y*sig_z)*exp((-y*y)/2./(sig_y*sig_y)) &
            *(exp((-(z_r-z_s)*(z_r-z_s))/2./(sig_z*sig_z))+exp((-(z_r+z_s)*(z_r+z_s))/2./(sig_z*sig_z)))
    endif
    
    end function gauss_plume_cartesian_trajectory_func

!==========================================================================
!   uEMEP model gauss_plume_cartesian_trajectory_integral_func
!   Cartesian Gaussian plume function
!==========================================================================

    function gauss_plume_cartesian_trajectory_integral_func(x,y,z_s,z_r,ay,by,az,bz,sig_y_0,sig_z_0,H1,H2,delta)

    implicit none
    real x,y,z_s,z_r
    real ay,by,az,bz,sig_y_0,sig_z_0,H1,H2,delta
    real gauss_plume_cartesian_trajectory_integral_func
    real sig_y,sig_z
    real pi,sig_limit
    parameter (pi=3.141592,sig_limit=3.)
    
    !r=sqrt((x_s-x_r)**2+(y_s-y_r)**2)
    !if (abs(u_s).lt.001) u_s=0.001
    !th=atan(v_s/u_s)
    !if (u_s.lt.0) th=th+pi
    !cos_val=cos(th)
    !sin_val=sin(th)
    !x=(x_r-x_s)*cos_val+(y_r-y_s)*sin_val
    !y=-(x_r-x_s)*sin_val+(y_r-y_s)*cos_val
    
    gauss_plume_cartesian_trajectory_integral_func=0.
    sig_y=sig_y_0+ay*exp(by*log(x))+x*abs(delta)
    if (x.ge.0.and.abs(y).lt.sig_y*sig_limit) then
        sig_z=sig_z_0+az*exp(bz*log(x))

        !gauss_plume_cartesian_integral_func=1./(2.*pi*sig_y)*exp((-y**2)/2./sig_y**2) &
        !    *sqrt(pi/2.)*(erf((z_s-H1)/sqrt(2.)/sig_z)-erf((z_s-H2)/sqrt(2.)/sig_z)+erf((z_s+H2)/sqrt(2.)/sig_z)-erf((z_s+H1)/sqrt(2.)/sig_z))/(H2-H1)
        gauss_plume_cartesian_trajectory_integral_func=1./(2.*pi*sig_y)*exp((-y*y)/2./(sig_y*sig_y)) &
            *sqrt(pi/2.)*(erf((z_s-H1)/sqrt(2.)/sig_z)-erf((z_s-H2)/sqrt(2.)/sig_z)+erf((z_s+H2)/sqrt(2.)/sig_z)-erf((z_s+H1)/sqrt(2.)/sig_z))/(H2-H1)
    endif
    
    end function gauss_plume_cartesian_trajectory_integral_func

    function gauss_plume_cartesian_sigma_func(x,y,z_s,z_r,sig_z,sig_y,z_pbl,FF)

    implicit none
    real, intent(in) :: x,y,z_s,z_r,sig_y,sig_z,z_pbl,FF
    real gauss_plume_cartesian_sigma_func
    real pi,sig_limit
    parameter (pi=3.141592,sig_limit=3.)
    real z_loop(5)
    real c_z
    integer k,n_loop
    
    gauss_plume_cartesian_sigma_func=0.
    if (x.ge.0.and.abs(y).lt.sig_y*sig_limit) then
        !write(*,*) 'here'
        !If the emission height z_s is greater than the boundary layer height z_pbl then only allow reflection from the surface
        !Also if z_r+z_s<z_pbl/3 then only use surface reflections since the pbl reflection counts for so little
        !Otherwise allow reflection from surface and boundary layer
        if (z_s.gt.z_pbl.or.z_s+sig_z.lt.z_pbl/3.) then
            n_loop=2
            z_loop(1)=z_s;z_loop(2)=-z_s
        else
            n_loop=5
            z_loop(1)=z_s;z_loop(2)=-z_s;z_loop(3)=2.*z_pbl-z_s;z_loop(4)=2.*z_pbl+z_s;z_loop(5)=-2.*z_pbl+z_s
        endif

        !For large sigmaz set the plume to be evenly distributed in the boundary layer
        !A value of 0.9 is used as this is the correct value for surface level releases and within 2% for high level releases
        !Should not be used for above boundary layer releases but is not accounted for
        if (sig_z.gt.0.9*z_pbl) then
            c_z=1./z_pbl
            gauss_plume_cartesian_sigma_func=c_z/(sqrt(2.*pi)*sig_y)*exp(-0.5*(y*y)/(sig_y*sig_y))/FF
        else
            c_z=0.
            do k=1,n_loop
                c_z=c_z+exp(-0.5*((z_r-z_loop(k))/sig_z)*((z_r-z_loop(k))/sig_z))
            enddo
            gauss_plume_cartesian_sigma_func=c_z/(2.*pi*sig_y*sig_z)*exp(-0.5*(y*y)/(sig_y*sig_y))/FF
        endif
        !if (x.eq.0) then
        !    write(*,*) gauss_plume_cartesian_sigma_func*FF,1./(2.*pi*sig_y*sig_z)
        !endif
        

    endif
    
    end function gauss_plume_cartesian_sigma_func
    
    function gauss_plume_cartesian_sigma_integral_func(x,y,z_s,z_r,sig_z,sig_y,z_pbl,FF,H1,H2)

    implicit none
    real, intent(in) :: x,y,z_s,z_r,sig_y,sig_z,z_pbl,FF,H1,H2
    real gauss_plume_cartesian_sigma_integral_func
    real pi,sig_limit
    parameter (pi=3.141592,sig_limit=3.)
    real z_loop(5)
    real c_z
    integer k,n_loop
    
    gauss_plume_cartesian_sigma_integral_func=0.
    if (x.ge.0.and.abs(y).lt.sig_y*sig_limit) then
        
        !If the emission height z_s is greater than the boundary layer height z_pbl then only allow reflection from the surface
        !Also if z_r+z_s<z_pbl/3 then only use surface reflections since the pbl reflection counts for so little
        !Otherwise allow reflection from surface and boundary layer
        if (z_s.gt.z_pbl.or.z_s+sig_z.lt.z_pbl/3.) then
            n_loop=2
            z_loop(1)=z_s;z_loop(2)=-z_s
        else
            n_loop=5
            z_loop(1)=z_s;z_loop(2)=-z_s;z_loop(3)=2.*z_pbl-z_s;z_loop(4)=2.*z_pbl+z_s;z_loop(5)=-2.*z_pbl+z_s
        endif

        !For large sigmaz set the plume to be evenly distributed in the boundary layer
        !A value of 0.9 is used as this is the correct value for surface level releases and within 2% for high level releases
        !Should not be used for above boundary layer releases but is not accounted for
        if (sig_z.gt.0.9*z_pbl) then
            c_z=1./z_pbl
            gauss_plume_cartesian_sigma_integral_func=c_z/(sqrt(2.*pi)*sig_y)*exp(-0.5*(y*y)/(sig_y*sig_y))/FF
        else
            c_z=0.
            do k=1,n_loop
                !c_z=c_z+exp(-0.5*((z_r-z_loop(k))/sig_z)*((z_r-z_loop(k))/sig_z))
                c_z=c_z+sqrt(pi/2.)*sig_z*(erf((z_loop(k)-H1)/sqrt(2.)/sig_z)-erf((z_loop(k)-H2)/sqrt(2.)/sig_z))/(H2-H1)
            enddo
            gauss_plume_cartesian_sigma_integral_func=c_z/(2.*pi*sig_y*sig_z)*exp(-0.5*(y*y)/(sig_y*sig_y))/FF
        endif

        !gauss_plume_cartesian_integral_func=1./(2.*pi*sig_y)*exp((-y*y)/2./(sig_y*sig_y)) &
        !    *sqrt(pi/2.)*(erf((z_s-H1)/sqrt(2.)/sig_z)-erf((z_s-H2)/sqrt(2.)/sig_z)+erf((z_s+H2)/sqrt(2.)/sig_z)-erf((z_s+H1)/sqrt(2.)/sig_z))/(H2-H1)

    endif
    
    end function gauss_plume_cartesian_sigma_integral_func

end module dispersion_functions

