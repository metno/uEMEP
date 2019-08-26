!   Functions for calculating dispersion from Kz and wind profiles
    
!==========================================================================
!   uEMEP_calculate_dispersion_Kz
!   Kz_func
!   z_centremass_gauss_func
!   phi_func
!   u_profile_val_func
!   
!   Calculates dispersion based on Kz and wind profiles
!==========================================================================

    subroutine uEMEP_set_dispersion_sigma_Kz(x_in,sig_z00,sig_y00,sigy_0_subgid_width_scale,sig_z_in,z_emis_loc,h_mix_loc,invL,u_val,z_val,logz0,subgrid_delta,u_star0_in,average_zc_h_in_Kz_flag,n_kz_iterations,sig_z,sig_y,u_zc)

    implicit none
    
    real, intent(in) :: x_in,sig_z00,sig_y00,sigy_0_subgid_width_scale,sig_z_in,z_emis_loc,h_mix_loc,invL,u_val,z_val,logz0,subgrid_delta(2),u_star0_in
    integer, intent(in) :: n_kz_iterations
    logical, intent(in) :: average_zc_h_in_Kz_flag
    real, intent(out) :: sig_z,sig_y,u_zc
    
    integer n_loop,j
    
    real K_z,K_y,L
    real :: z_tau_min=2.,z_tau_max=100.
    real z0,zc,ustar0,u_hmix
    real :: K_min=0.001
    real l_t,f_t
    real u_star0,u_star0_val,tau,zc_start,K_z_start,u_zc_start
    real min_xy,x,h_y
    
    n_loop=n_kz_iterations
    
    !Limit stable L to this positive L value 
    !L=1.e6
    !L=lowest_stable_L
    !if (abs(invL).gt.1./L) L=1./invL
    L=1./invL
    
    z0=exp(logz0)
    min_xy=(subgrid_delta(1)+subgrid_delta(2))/4.
    !min_x=1.
    !x=max(x_in,min_xy)
    !Set x to this value because it simulates that it has already travelled half a grid to get this sig_z0
    !Which is why sg_z is added only to sig_z00   
    x=x_in+min_xy
    
    !Initialise sig_z
    sig_z=sig_z_in
    
    !If the emission is above the boundary layer height then set the initial plume guess to its low turbulence form
    !Same values used in the emulator (1/1000 slender plume)
    if (z_emis_loc/h_mix_loc.ge.1.0) sig_z=sig_z00+0.001*exp(1.0*log(x))
    
    !Set ustar0 for K_z to the value from EMEP
    u_star0=u_star0_in
    
    zc=z_emis_loc
    
    !Set zc and K_z at start of plume based on sig_z0
    !call z_centremass_gauss_func(sig_z0,z_emis_loc,h_mix_loc,zc_start)
    !call Kz_func(h_mix_loc,L,u_star0,zc_start,K_min,K_z_start)
    !call u_profile_val_func(zc_start,L,u_val,z_val,h_mix_loc,z0,u_zc_start,u_star0,u_hmix)
    
    !Put bug back in
    u_zc=u_val
    
    !Calculate the Lagrangian time scale before the iteration loop using a minimum distance for this to make it non zero on the grid
    !l_t=max(x,min_xy)/u_zc
    tau=0.6*max(min(z_tau_max,z_emis_loc),z_tau_min)/u_star0
    !f_t=1.+tau/l_t*(exp(-l_t/tau)-1.)
        !if (x.lt.200) then
        !write(*,'(a,i,5f)') 'S: ',j,tau,l_t,f_t,u_zc,sig_z
        !endif
        
    !All functions commented out 4 secs
    !All functions included 15 secs
    !Without u_profile 11 secs
    !Commenting out all functions below gives 7 secs
    do j=1,n_loop
        
        !Calculate centre of mass for the emission height.
        call z_centremass_gauss_func(sig_z,z_emis_loc,h_mix_loc,zc) !4.5 sec
        
        if (average_zc_h_in_Kz_flag) zc=(z_emis_loc+zc)/2.
        
        !write(*,'(i,4f)') j,sig_z,z_emis_loc,h_mix_loc,zc
        
        !Calculate the wind profile
        !call u_profile_val_func(zc,L,u_val,z_val,h_mix_loc,z0,u_zc,u_star0_val,u_hmix) !7secs. Without phim_func calls 5 secs
        call u_profile_neutral_val_func(zc,u_val,z_val,h_mix_loc,z0,u_zc,u_star0_val) 
        
        !write(*,'(i,9f)') j,zc,L,u_val,z_val,h_mix_loc,z0,u_zc,u_star0,u_hmix
        !u_zc=(u_zc+u_zc_start)/2.
        !write(*,'(i,8f)') j,x,u_star0_in,u_star0,zc,z_val,u_zc,u_zc_start,u_val
        
        !Use the calculated u_star for the dispersion, not the EMEP input
        !u_star0=u_star0_val
        
        !Calculate K_z at the centre of mass
        call Kz_func(h_mix_loc,L,u_star0,zc,K_min,K_z) !5.5 secs
        
        !write(*,'(i,1f)') j,K_z
        !Take average of K_z_start(x=0) and K_z(x)
        !K_z=(K_z+K_z_start)/2.
        
        !calculate l_t based on the centre of mass wind speed
        l_t=max(x,min_xy)/u_zc
        f_t=1.+tau/l_t*(exp(-l_t/tau)-1.)
        
        !Calculate sig_z for the next iteration, using the sig_z0 from the neutral plume approximation
        sig_z=sig_z00+sqrt(2.*K_z*l_t*f_t)
        !write(*,'(i,1f)') j,sig_z
        !write(*,*) K_z,l_t,f_t,sig_z,u_zc
        !write(*,*) j,z_emis_loc,zc,h_mix_loc*0.25,sig_z
    
        !if (x.lt.200) then
        !write(*,'(a,i,5f)') 'S: ',j,tau,l_t,f_t,u_zc,sig_z
        !endif
        
    enddo
    
    !Calculate sigma_y at the maximum K, around the average of the emission height and 0.25 of the boundary layer height
    !THis is new 27.10.2018 and not tested. It will reduce sig_y which is OK
    h_y=(h_mix_loc*0.25+z_emis_loc)/2.
    call Kz_func(h_mix_loc,L,u_star0,h_y,K_min,K_y)
    sig_y=sig_y00+min_xy*sigy_0_subgid_width_scale+sqrt(2.*K_y*l_t*f_t)
    
    !Should change this to what is documented, i.e. 2*sig_z. Need to test
    !Also sigy_0_subgid_width_scale should perhaps be 0.5, not 0.25. Also need to test
    !sig_y=sig_y00+min_xy*sigy_0_subgid_width_scale+(sig_z-sigz00)*2.0
    
    end subroutine uEMEP_set_dispersion_sigma_Kz
    
    subroutine Kz_func(z_pbl,L,u_star0_in,z,K_min,K_z)

    implicit none

    real, intent(in) :: z_pbl,L,u_star0_in,z,K_min
    real, intent(out) :: K_z
    real kappa
    parameter (kappa=0.4)
    real phih,phih_i
    real phih_hs,phih_i_hs
    real phih_hs_p1,phih_i_hs_p1
    real K_zpbl,K_zhs,h_s,delta_K_zhs,u_star0

    h_s=0.04*z_pbl
    
    !ustar0 cannot be 0. Set to a low value
    !u_star0=max(u_star0_in,0.01)
    u_star0=u_star0_in
    
    call phih_func(z,L,phih,phih_i)

    K_zpbl=K_min
    K_z=K_min

    if (L.ge.0) then
        K_z=0.39*u_star0*z*exp(-0.5*(z/0.21/z_pbl)*(z/0.21/z_pbl))  !As in EMEP
        K_z=0.39*u_star0*z*exp(-0.5*(z/0.32/z_pbl)*(z/0.32/z_pbl)) !Adjusted to match under neutral conditions
        K_z=0.39*u_star0*z/phih*exp(-0.5*(z/0.32/z_pbl)*(z/0.32/z_pbl)) !Adjusted but with stability added
    else
        call phih_func(h_s,L,phih_hs,phih_i_hs)
        call phih_func(h_s+1.,L,phih_hs_p1,phih_i_hs_p1)
        K_zhs=u_star0*kappa*h_s/phih_hs
        delta_K_zhs=u_star0*kappa*((h_s+1.)/phih_hs_p1-(h_s)/phih_hs)
        if (z.le.h_s) K_z=u_star0*kappa*z/phih
        if (z.gt.h_s.and.z.le.z_pbl) K_z=K_zpbl+((z_pbl-z)/(z_pbl-h_s))*((z_pbl-z)/(z_pbl-h_s))*(K_zhs-K_zpbl+(z-h_s)*(delta_K_zhs+2.*(K_zhs-K_zpbl)/(z_pbl-h_s)))
    endif

    K_z=max(K_z,K_min)

    end subroutine Kz_func

    subroutine z_centremass_gauss_func(sigma,h,z_pbl,z_c)

    implicit none
    real, intent(in) :: sigma,h,z_pbl
    real, intent(out) :: z_c
    real z_loop(5)
    real H_c
    real c_z,c_av
    integer i_loop,i
    real sqrt_2pi,sqrt_2
    real pi
    parameter (pi=3.141592653589793)

    sqrt_2pi=sqrt(2.*pi)
    sqrt_2=sqrt(2.)

    i_loop=5
    z_loop(1)=h;z_loop(2)=-h;z_loop(3)=2.*z_pbl-h;z_loop(4)=2.*z_pbl+h;z_loop(5)=-2.*z_pbl+h
    z_c=0.
    !c_z=0.;c_av=0.
    H_c=z_pbl

    !If the emission height h is greater than the boundary layer height then only allow reflection from the surface
    !and set the top of the integration H_c to infinity
    if (h.gt.z_pbl) then
        H_c=1.e16
        i_loop=2
    endif
    
    !Reduce the loop size when the reflection from the boundary layer is not important
    if (sigma+h.lt.z_pbl/3.) then
        i_loop=2
    endif

    !Remove this after finished testing
    !i_loop=5
    
    do i=1,i_loop

        z_c=z_c+sigma/sqrt_2pi*(exp(-0.5*(z_loop(i)/sigma)*(z_loop(i)/sigma))-exp(-0.5*((H_c-z_loop(i))/sigma)*((H_c-z_loop(i))/sigma))) &
            +z_loop(i)/2*(erf((H_c-z_loop(i))/sqrt_2/sigma)+erf((z_loop(i))/sqrt_2/sigma))
        
        !c_av=c_av+0.5*(erf((z_pbl-z_loop(i))/sqrt_2/sigma)+erf((z_loop(i))/sqrt_2/sigma))/z_pbl
    
        !c_z=c_z+1./sigma/sqrt_2pi*exp(-0.5*((z-z_loop(i))/sigma)*((z-z_loop(i))/sigma))

    enddo

    !if (sigma.gt.0.9*z_pbl) then
        !c_z=1./z_pbl
        !c_av=1./z_pbl
    !endif

    end subroutine z_centremass_gauss_func
    
    
    subroutine u_profile_val_func(z,L,u_val,z_val_in,z_pbl,z0,u,u_star0,u_pbl)

    implicit none

    real, intent(in) :: z,L,u_val,z_val_in,z_pbl,z0
    real, intent(out) :: u,u_star0,u_pbl
    real a,b,p,q,kappa,pi
    parameter (a=16.,b=5.,p=-0.25,q=-0.5,kappa=0.4,pi=3.141592653589793)
    real z_l,z_val
    real phim,phim_i
    real phim_val,phim_i_val
    !real phim_pbl,phih_pbl,phim_i_pbl,phih_i_pbl

    !If the input height is above the boundary layer then set the height to pbl height and calculate
    if (z_val.ge.z_pbl) then
        z_val=z_pbl
    else
        z_val=z_val_in
    endif

    z_l=0.4*z_pbl

    call phim_func(z,L,phim,phim_i)
    call phim_func(z_val,L,phim_val,phim_i_val)
    
    !call phi_func(z_pbl,L,phim_pbl,phih_pbl,phim_i_pbl,phih_i_pbl)

    if (L.ge.0.) then
        u_star0=u_val*kappa/(log(z_val/z0)-phim_i_val+kappa*z_val/z_l*(1-z_val/2./z_pbl)-z_val/z_pbl*(1+b*z_val/2./L))
        !u_pbl=u_star0/kappa*(log(z_pbl/z0)-phim_i_pbl+kappa*z_pbl/z_l*(1-z_pbl/2./z_pbl)-z_pbl/z_pbl*(1+b*z_pbl/2./L))
        u=u_star0/kappa*(log(z/z0)-phim_i+kappa*z/z_l*(1-z/2./z_pbl)-z/z_pbl*(1+b*z/2./L))
    else
        u_star0=u_val*kappa/(log(z_val/z0)-phim_i_val+kappa*z_val/z_l*(1-z_val/2./z_pbl)-1./z_pbl*((a*z_val-L)*phim_val+L)/a/(p+1))   
        !u_pbl=u_star0/kappa*(log(z_pbl/z0)-phim_i_pbl+kappa*z_pbl/z_l*(1-z_pbl/2./z_pbl)-1./z_pbl*((a*z_pbl-L)*phim_pbl+L)/a/(p+1))  
        u=u_star0/kappa*(log(z/z0)-phim_i+kappa*z/z_l*(1-z/2./z_pbl)-1./z_pbl*((a*z-L)*phim+L)/a/(p+1))   
    endif

    !u_pbl not used and so not calculated here
    u_pbl=u
    if (z.ge.z_pbl) then
        u=u_pbl
    endif

    end subroutine u_profile_val_func
    
    subroutine u_profile_neutral_val_func(z,u_val,z_val_in,z_pbl,z0,u,u_star0)

    implicit none

    real, intent(in) :: z,u_val,z_val_in,z_pbl,z0
    real, intent(out) :: u,u_star0
    real kappa,pi
    parameter (kappa=0.4,pi=3.141592653589793)
    real z_l,z_val

    !If the input height is above the boundary layer then set the height to pbl height and calculate
    z_val=min(z_val_in,z_pbl)

    z_l=0.4*z_pbl

    u_star0=u_val*kappa/(log(z_val/z0)+kappa*z_val/z_l*(1-z_val/2./z_pbl)-z_val/z_pbl)
    u=u_star0/kappa*(log(z/z0)+kappa*z/z_l*(1-z/2./z_pbl)-z/z_pbl)
    
    end subroutine u_profile_neutral_val_func

    subroutine phi_func(z,L,phim,phih,phim_i,phih_i)

    implicit none

    real, intent(in) :: z,L
    real, intent(out) :: phim,phih,phim_i,phih_i
    real a,b,p,q,pi
    parameter (a=16.,b=5.,p=-0.25,q=-0.5,pi=3.141592653589793)
    real eps

    eps=z/L

    if (eps.ge.0) then
        phim=1.+b*eps
        phim_i=-b*eps
        phih=phim
        phih_i=phim_i
    else
        phim=exp(p*log((1.-a*eps)))
        phih=exp(q*log((1.-a*eps)))
        phim_i=2.*log((1.+1./phim)/2.)+log((1.+1./(phim*phim))/2.)-2.*atan(1./phim)+pi/2.
        phih_i=2.*log((1.+1./phih)/2.)
    endif

    end subroutine phi_func
 
    subroutine phim_func(z,L,phim,phim_i)

    implicit none

    real, intent(in) :: z,L
    real, intent(out) :: phim,phim_i
    real a,b,p,q,pi
    parameter (a=16.,b=5.,p=-0.25,q=-0.5,pi=3.141592653589793)
    real eps

    eps=z/L

    if (eps.ge.0) then
        phim=1.+b*eps
        phim_i=-b*eps
    else
        phim=exp(p*log((1.-a*eps)))
        phim_i=2.*log((1.+1./phim)/2.)+log((1.+1./(phim*phim))/2.)-2.*atan(1./phim)+pi/2.
    endif

    end subroutine phim_func
    
    subroutine phih_func(z,L,phih,phih_i)

    implicit none

    real, intent(in) :: z,L
    real, intent(out) :: phih,phih_i
    real a,b,p,q,pi
    parameter (a=16.,b=5.,p=-0.25,q=-0.5,pi=3.141592653589793)
    real eps

    eps=z/L

    if (eps.ge.0) then
        phih=1.+b*eps
        phih_i=-b*eps
    else
        phih=exp(q*log((1.-a*eps)))
        phih_i=2.*log((1.+1./phih)/2.)
    endif

    end subroutine phih_func
    
    
    subroutine z_centremass_gauss_array_func(sig_norm,h_norm,n_array,zc_array)

    implicit none
    real, intent(in) :: sig_norm,h_norm,n_array
    real, intent(out) :: zc_array(n_array)
    real z_loop(5)
    real H_c
    real z_c
    integer i_loop,i,k
    real sqrt_2pi,sqrt_2
    real pi
    parameter (pi=3.141592653589793)

    sqrt_2pi=sqrt(2.*pi)
    sqrt_2=sqrt(2.)

    i_loop=5
    z_loop(1)=h_norm;z_loop(2)=-h_norm;z_loop(3)=2.-h_norm;z_loop(4)=2.+h_norm;z_loop(5)=-2.+h_norm
    H_c=1.

    !If the emission height h is greater than the boundary layer height then only allow reflection from the surface
    !and set the top of the integration H_c to infinity
    if (h_norm.gt.1.) then
        H_c=1.e4
        i_loop=2
    endif
    
    !Reduce the loop size when the reflection from the boundary layer is not important
    !if (sig_norm+h_norm.lt.1./3.) then
    !    i_loop=2
    !endif

    !Remove this after finished testing
    i_loop=5
    
    do k=1,n_array

        z_c=0.
        
        do i=1,i_loop

            z_c=z_c+sig_norm/sqrt_2pi*(exp(-0.5*(z_loop(i)/sig_norm)*(z_loop(i)/sig_norm))-exp(-0.5*((H_c-z_loop(i))/sig_norm)*((H_c-z_loop(i))/sig_norm))) &
                +z_loop(i)/2.*(erf((H_c-z_loop(i))/sqrt_2/sig_norm)+erf((z_loop(i))/sqrt_2/sig_norm))
        
        enddo
        zc_array(k)=z_c
    
    enddo

    end subroutine z_centremass_gauss_array_func
    
