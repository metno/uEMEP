!uEMEP_set_dispersion_params.f90
!Routines for calculating dispersion parameters when sig_z and sig_y are defined as sig=sig0+a*x^b
    
    subroutine uEMEP_set_dispersion_params_simple(source_index,subsource_index)

    use uEMEP_definitions
    
    implicit none

    integer source_index,subsource_index
    
    
    !Set the psedo dispersion parameters here.
    ay(source_index,subsource_index)=0.32
    by(source_index,subsource_index)=0.78
    az(source_index,subsource_index)=0.22
    bz(source_index,subsource_index)=0.78
    !ay(source_index,subsource_index)=0.64!From Liu
    !by(source_index,subsource_index)=0.46!From Liu
    !az(source_index,subsource_index)=0.088!From Liu
    !bz(source_index,subsource_index)=0.72!From Liu 0.72
    !ay(source_index,subsource_index)=0.32
    !by(source_index,subsource_index)=0.78
    az(source_index,subsource_index)=0.2
    bz(source_index,subsource_index)=0.75

    !Alternative to ASME
    !ay(source_index,subsource_index)=0.14
    !by(source_index,subsource_index)=0.9
    !az(source_index,subsource_index)=0.22
    !bz(source_index,subsource_index)=0.85

    !sig_y_0(source_index,subsource_index)=sig_y_00(source_index,subsource_index)
    !sig_z_0(source_index,subsource_index)=sig_z_00(source_index,subsource_index)
    !sig_y_0(source_index,subsource_index)=sig_y_00(source_index,subsource_index)+sqrt(emission_subgrid_delta(x_dim_index,source_index)*emission_subgrid_delta(y_dim_index,source_index))/2.
    !sig_z_0(source_index,subsource_index)=sig_z_00(source_index,subsource_index)+az(source_index,subsource_index)*exp(bz(source_index,subsource_index)*log(sig_y_0(source_index,subsource_index)))
    !h_emis(source_index,subsource_index)=1.
    !z_rec(source_index,subsource_index)=2.
    
    !Exceptions
    !ay=ay*3.
    !az=az*5.
    !bz=0.95
    
    end subroutine uEMEP_set_dispersion_params_simple

    
    subroutine uEMEP_set_dispersion_params_PG(invL,source_index,subsource_index)

    use uEMEP_definitions
    
    implicit none

    integer source_index,subsource_index
    real invL
    integer class_index
    
    real ay_pg(6),by_pg(6),az_pg(6),bz_pg(6)
    real L_class(5),invL_class(5)
    data ay_pg /0.469,0.306,0.230,0.219,0.237,0.237/
    data by_pg /0.903,0.855,0.855,0.764,0.691,0.594/
    data az_pg /0.017,0.072,0.076,0.140,0.217,0.262/
    data bz_pg /1.38,1.021,0.879,0.727,0.610,0.500/
    !data L_class /-10.,-25.,-200.,200.,25./
    data L_class /-10.,-25.,-100.,50.,5./
    
    invL_class=1./L_class
    
    if (invL.le.invL_class(1)) then
        class_index=1
    elseif (invL.gt.invL_class(1).and.invL.le.invL_class(2)) then
        class_index=2
    elseif (invL.gt.invL_class(2).and.invL.le.invL_class(3)) then
        class_index=3
    elseif (invL.gt.invL_class(3).and.invL.le.invL_class(4)) then
        class_index=4
    elseif (invL.gt.invL_class(4).and.invL.le.invL_class(5)) then
        class_index=5
    elseif (invL.gt.invL_class(5)) then
        class_index=6
    else
        class_index=0
        write(*,*) 'No stability class found. Stopping',1./invL
        stop
    endif
        
    !if (class_index.ne.4) write(*,*) invL,class_index
    
    !Set the dispersion parameters here.
    ay(source_index,subsource_index)=ay_pg(class_index)
    by(source_index,subsource_index)=by_pg(class_index)
    az(source_index,subsource_index)=az_pg(class_index)
    bz(source_index,subsource_index)=bz_pg(class_index)

    !if (bz(source_index,subsource_index).ne.0.727) write(*,*) invL,class_index,bz(source_index,subsource_index)
    !ay(source_index,subsource_index)=0.219
    !by(source_index,subsource_index)=0.764
    !az(source_index,subsource_index)=0.114
    !bz(source_index,subsource_index)=0.727


    sig_y_0(source_index,subsource_index)=sig_y_00(source_index,subsource_index)+sqrt(emission_subgrid_delta(x_dim_index,source_index)*emission_subgrid_delta(y_dim_index,source_index))/4.*sigy_0_subgid_width_scale
    sig_z_0(source_index,subsource_index)=sig_z_00(source_index,subsource_index)+az(source_index,subsource_index)*exp(bz(source_index,subsource_index)*log(sig_y_0(source_index,subsource_index)))
    
    
    end subroutine uEMEP_set_dispersion_params_PG
    
    
    subroutine delta_wind_direction (i_cross,j_cross,tt,temp_FF_subgrid,angle_diff)
    
    use uEMEP_definitions
    
    implicit none

    integer, intent(in) :: i_cross,j_cross,tt
    real, intent(in) :: temp_FF_subgrid
    real, intent(out) :: angle_diff
    real :: meandering_degree_max=20.
    real cos_subgrid_loc,sin_subgrid_loc
    
    if (use_last_meteo_in_dispersion) then
        
        cos_subgrid_loc=meteo_subgrid(i_cross,j_cross,tt,cos_subgrid_index)
        sin_subgrid_loc=meteo_subgrid(i_cross,j_cross,tt,sin_subgrid_index)
        angle_diff=abs(asin(meteo_subgrid(i_cross,j_cross,tt,sin_subgrid_index)*last_meteo_subgrid(i_cross,j_cross,cos_subgrid_index) &
                            -meteo_subgrid(i_cross,j_cross,tt,cos_subgrid_index)*last_meteo_subgrid(i_cross,j_cross,sin_subgrid_index) ))/2.
                
        angle_diff=min(angle_diff,3.14159/4.) !Less than 45 degrees
        !write(*,*) i_cross,j_cross,angle_diff(i_cross,j_cross)*180/3.14159
                                        
    else
        
        angle_diff=0.
    
    endif

    if (use_meandering_in_dispersion) then
        angle_diff=angle_diff+(meandering_degree_max*3.14159/180.)*exp(-(temp_FF_subgrid-FF_min_dispersion)/(FF_min_dispersion*2.))
    endif
    
    end subroutine delta_wind_direction
    
    subroutine uEMEP_set_dispersion_sigma_simple(sig_z_00,sig_y_00,sigy_0_subgid_width_scale,subgrid_delta,delta_wind,x,sig_z,sig_y,sig_z_0,sig_y_0)
    
    implicit none

    real, intent(in) :: sig_z_00,sig_y_00,sigy_0_subgid_width_scale,subgrid_delta(2),delta_wind,x
    real, intent(out) :: sig_z,sig_y,sig_z_0,sig_y_0
    
    real ay,by,az,bz
    real min_xy
    
    !Set the psedo dispersion parameters for neutral conditions
    
    !From Klug
    !ay=0.32
    !by=0.78
    !az=0.22
    !bz=0.78
    
    !From Liu
    !ay=0.64
    !by=0.46!From Liu
    !az=0.088!From Liu
    !bz=0.72!From Liu 0.72
    
    !Alternative ASME
    !ay=0.14
    !by=0.9
    !az=0.22
    !bz=0.85

    !Update from K_z fitting
    ay=0.32
    by=0.78
    !az=0.19
    !bz=0.77
    az=0.125
    bz=0.85 !For z0=0.1, corresponding to the same as K_z for wind height at emission source of 1 m
    az=0.245
    bz=0.711 !For z0=0.3, corresponding to the same as K_z for wind height at average of emission source of 1 m and zc
    az=0.21
    bz=0.79 !For z0=0.3, corresponding to the same as K_z for wind height at emission source of 1 m
    !Consistant with uEMEP_set_dispersion_params_PG needs to be fixed. Just one call to the parameters, one calle to the sigma calculation
    az=0.2
    bz=0.75


    min_xy=(subgrid_delta(1)+subgrid_delta(2))/4.
    !Set sig_y_0 to be half of the average x,y grid size
    !Add this here ay*exp(by*log(min_xy)) to be the same as sig_z and the same as the  Kz calculation
    !Does not mean it is correct, just closer to the Kz which is perhaps not so correct
    sig_y_0=sig_y_00+min_xy*sigy_0_subgid_width_scale+ay*exp(by*log(min_xy))
    !Set sig_z_0 to be the size of the plume after travelling half of the grid size
    sig_z_0=sig_z_00+az*exp(bz*log(min_xy))

    !Set sig_y and sig_z = sig_0 + a*x^b +x*delata_wind
    sig_y=sig_y_0+ay*exp(by*log(x))+x*abs(delta_wind)
    sig_z=sig_z_0+az*exp(bz*log(x))

    end subroutine uEMEP_set_dispersion_sigma_simple
    
    subroutine uEMEP_set_dispersion_sigma_PG(invL_in,logz0,sig_z_00,sig_y_00,sigy_0_subgid_width_scale,subgrid_delta,delta_wind,x,sig_z,sig_y,sig_z_0,sig_y_0)
    
    implicit none

    real, intent(in) :: invL_in,logz0,sig_z_00,sig_y_00,sigy_0_subgid_width_scale,subgrid_delta(2),delta_wind,x
    real, intent(out) :: sig_z,sig_y,sig_z_0,sig_y_0
    
    integer i_bot,i_top,i
    real weight
    real ay,by,az,bz
    real min_xy
    real ay_pg(6),by_pg(6),az_pg(6),bz_pg(6)
    real a_invL(6),b_invL(6),invL(6)
    !Use ASME parameters
    data ay_pg /0.40,0.36,0.34,0.32,0.315,0.31/
    data by_pg /0.91,0.86,0.82,0.78,0.75,0.71/
    data az_pg /0.40,0.33,0.27,0.22,0.14,0.06/
    data bz_pg /0.91,0.86,0.82,0.78,0.75,0.71/
    !Conversion of classes to L
    data a_invL /-0.096,-0.037,-0.002,0.0,0.004,0.035/
    data b_invL /0.029,0.029,0.018,0.0,-0.018,-0.036/
    
    invL=a_invL+b_invL*logz0
    min_xy=(subgrid_delta(1)+subgrid_delta(2))/4.

    !Find and interpolate the stability class based on input invL
    if (invL_in.le.invL(1)) then
        i_bot=1
        i_top=1
        weight=0.
    elseif (invL_in.gt.invL(6)) then
        i_bot=6
        i_top=6
        weight=1.
    else
        do i=1,5
            if (invL_in.gt.invL(i).and.invL_in.le.invL(i+1)) then
                i_bot=i
                i_top=i+1
                weight=(invL_in-invL(i))/(invL(i+1)-invL(i))
                exit
            endif
        enddo
    endif
    
    ay=ay_pg(i_bot)*(1.-weight)+ay_pg(i_top)*weight
    by=by_pg(i_bot)*(1.-weight)+by_pg(i_top)*weight
    az=az_pg(i_bot)*(1.-weight)+az_pg(i_top)*weight
    bz=bz_pg(i_bot)*(1.-weight)+bz_pg(i_top)*weight
    
    !Set sig_y_0 to be half of the average x,y grid size
    !Add this here ay*exp(by*log(min_xy)) to be the same as sig_z and the same as the  Kz calculation
    !Does not mean it is correct, just closer to the Kz which is perhaps not so correct
    sig_y_0=sig_y_00+min_xy*sigy_0_subgid_width_scale+ay*exp(by*log(min_xy))
    !Set sig_z_0 to be the size of the plume after travelling half of the grid size
    sig_z_0=sig_z_00+az*exp(bz*log(min_xy))

    !Set sig_y and sig_z = sig_0 + a*x^b +x*delata_wind
    if (x.le.0) then
    sig_y=sig_y_0
    sig_z=sig_z_0
    else
    sig_y=sig_y_0+ay*exp(by*log(x))+x*abs(delta_wind)
    sig_z=sig_z_0+az*exp(bz*log(x))
    endif
    !write(*,'(i,6f)') i,weight,sig_y,sig_z,az,bz,x
    
    end subroutine uEMEP_set_dispersion_sigma_PG
    
    subroutine uEMEP_set_dispersion_sigma_Kz_emulator(z_emis,invL,logz0,z_pbl,sig_z_00,sig_y_00,sigy_0_subgid_width_scale,subgrid_delta,delta_wind,x,sig_z,sig_y,sig_z_0,sig_y_0)
    
    implicit none

    real, intent(in) :: z_emis,invL,logz0,z_pbl,sig_z_00,sig_y_00,sigy_0_subgid_width_scale,subgrid_delta(2),delta_wind,x
    real, intent(out) :: sig_z,sig_y,sig_z_0,sig_y_0
    
    real invL_in,zz_pbl,z0
    real ay,by,az,bz
    real min_xy
    real z0_ref,zz_pbl_ref,zz_pbl_L,logz
    real x_val
    
    !invL_in=1./L
    min_xy=(subgrid_delta(1)+subgrid_delta(2))/4.
    x_val=max(x,min_xy)
    
    !Remove the stable cases as these are not normally done in the full K_z formulation
    !invL_in=min(invL_in,1./100.)
    invL_in=min(invL,0.)
    !invL_in=0.
    
    zz_pbl=min(.95,z_emis/z_pbl) !Cannot be greater than 0.95 as it is outside the emulator region
    zz_pbl_L=z_pbl*invL_in;
    logz=log(z_emis)
    z0=exp(logz0)
    
    !Original
    !az=0.15+0.52*(z0-0.02)-0.15*(1.-EXP(-zz_pbl/0.03))+0.16*SIGN(1.0,invL_in)*(1.-EXP(-ABS(invL_in)/zz_pbl/5.))
    !bz=0.77-0.15*(1.-EXP(-z0/0.3))+0.2*(1.-EXP(-zz_pbl/0.03))-0.4*SIGN(1.0,invL_in)*(1.-EXP(-ABS(invL_in)/zz_pbl/8.))

    !Alternative form
    z0_ref=0.1
    zz_pbl_ref=0.001
    
    !az=0.15+.70*(z0-z0_ref)-0.1*(1.-exp(-(zz_pbl-zz_pbl_ref)*30.))+0.01*sign(1.,invL_in)*(1.-exp(-abs(invL_in)*10.))
    !bz=0.76-0.17*(1.-exp(-(z0-z0_ref)*.3))+0.2*(1.-exp(-(zz_pbl-zz_pbl_ref)*30.))-0.4*sign(1.,invL_in)*(1.-exp(-abs(invL_in)*10.))
    !ay=0.15+.70*(z0-z0_ref)-0.1*(1.-exp(-(zz_pbl)*10.))+.01*sign(1.,invL_in)*(1.-exp(-abs(invL_in)*10.))
    !by=0.76-0.17*(1.-exp(-(z0-z0_ref)*3))+0.70*(1.-exp(-(zz_pbl)*1.))-.4*sign(1.,invL_in)*(1.-exp(-abs(invL_in)*10.))
    
    az=(-tanh((logz-1.2)*3.)*0.08+0.1+(1-exp(-zz_pbl/0.05))*0.06)*(1+(z0-z0_ref)/(z0+0.6)*3.)-.02*sign(1.,invL_in)*(1.-exp(-abs(zz_pbl_L)*10.))
    bz=tanh((logz-1.2)*3)*0.14+0.88-(1-exp(-zz_pbl/0.05))*0.11-log(z0/z0_ref)*0.04-.20*sign(1.,invL_in)*(1.-exp(-abs(invL_in)*10.))
    if (z_emis/z_pbl.gt.1) then
        az=0.001
        bz=1.0
    endif
    
    !Limit to the values explored by the emulator, just in case
    az=min(max(az,0.001),0.7)
    bz=min(max(bz,0.4),1.2)
    
    !Calculate y values, taken as close to the maximum K height of z/z_pbl=0.25
    zz_pbl=0.25
    zz_pbl=(z_pbl*0.25+z_emis)/2./z_pbl

    logz=log(zz_pbl*z_pbl)
    
    !ay=0.15+0.52*(z0-0.02)-0.15*(1.-EXP(-zz_pbl/0.03))+0.16*SIGN(1.0,invL_in)*(1.-EXP(-ABS(invL_in)/zz_pbl/5.))
    !by=0.77-0.15*(1.-EXP(-z0/0.3))+0.2*(1.-EXP(-zz_pbl/0.03))-0.4*SIGN(1.0,invL_in)*(1.-EXP(-ABS(invL_in)/zz_pbl/8.))
    
    !az=0.15+.70*(z0-z0_ref)-0.1*(1.-exp(-(zz_pbl-zz_pbl_ref)*30.))+0.01*sign(1.,invL_in)*(1.-exp(-abs(invL_in)*10.))
    !bz=0.76-0.17*(1.-exp(-(z0-z0_ref)*.3))+0.2*(1.-exp(-(zz_pbl-zz_pbl_ref)*30.))-0.4*sign(1.,invL_in)*(1.-exp(-abs(invL_in)*10.))
    !ay=0.15+.70*(z0-z0_ref)-0.1*(1.-exp(-(zz_pbl)*10.))+.01*sign(1.,invL_in)*(1.-exp(-abs(invL_in)*10.))
    !by=0.76-0.17*(1.-exp(-(z0-z0_ref)*3))+0.70*(1.-exp(-(zz_pbl)*1.))-.4*sign(1.,invL_in)*(1.-exp(-abs(invL_in)*10.))

    ay=(-tanh((logz-1.2)*3.)*0.08+0.1+(1-exp(-zz_pbl/0.05))*0.06)*(1+(z0-z0_ref)/(z0+0.6)*3.)-.02*sign(1.,invL_in)*(1.-exp(-abs(zz_pbl_L)*10.))
    by=tanh((logz-1.2)*3)*0.14+0.88-(1-exp(-zz_pbl/0.05))*0.11-log(z0/z0_ref)*0.04-.20*sign(1.,invL_in)*(1.-exp(-abs(invL_in)*10.))
    if (z_emis/z_pbl.ge.1.) then
        ay=0.001
        by=1.0
    endif

    ay=min(max(ay,0.001),0.7)
    by=min(max(by,0.4),1.2)
    
    !Set sig_y_0 to be half of the average x,y grid size
    !Mulitiply by the scale
    sig_y_0=sig_y_00+min_xy*sigy_0_subgid_width_scale
    
    !Set sig_z_0 to be the size of the plume after travelling half of the grid size
    !sig_z_0=sig_z_00+az*exp(bz*log(min_xy))

    !Set sig_y and sig_z = sig_0 + a*x^b +x*delata_wind
    sig_y=sig_y_0+ay*exp(by*log(x+min_xy))+(x+min_xy)*abs(delta_wind)
    !sig_z=sig_z_0+az*exp(bz*log(x))
    !if (x.lt.10.) write(*,*) x,sig_z_00,sig_z
    sig_z_0=sig_z_00+az*exp(bz*log(min_xy))
    sig_z=sig_z_00+az*exp(bz*log(x+min_xy))
    !if (x.lt.10.) write(*,*) x+min_xy,sig_z_00,sig_z
    !if (x.lt.10.) write(*,*)

    !if (zz_pbl.le.0.or.sig_z.le.0.or.sig_y.le.0) then
    !    write(*,'(7f)') az,bz,ay,by,sig_z,sig_y,x
    !    stop
    !endif
    
    
    end subroutine uEMEP_set_dispersion_sigma_Kz_emulator