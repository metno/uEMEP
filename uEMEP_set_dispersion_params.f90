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
    ay(source_index,subsource_index)=0.32
    by(source_index,subsource_index)=0.78
    az(source_index,subsource_index)=0.19
    bz(source_index,subsource_index)=0.77
    
    sig_y_0(source_index,subsource_index)=sig_y_00(source_index,subsource_index)+sqrt(emission_subgrid_delta(x_dim_index,source_index)*emission_subgrid_delta(y_dim_index,source_index))/2.
    sig_z_0(source_index,subsource_index)=sig_z_00(source_index,subsource_index)+az(source_index,subsource_index)*exp(bz(source_index,subsource_index)*log(sig_y_0(source_index,subsource_index)))
    !h_emis(source_index,subsource_index)=1.
    !z_rec(source_index,subsource_index)=2.
    
    !Exceptions

    
    
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
    data L_class /-10.,-25.,-200.,200.,25./
    
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


    sig_y_0(source_index,subsource_index)=sig_y_00(source_index,subsource_index)+sqrt(emission_subgrid_delta(x_dim_index,source_index)*emission_subgrid_delta(y_dim_index,source_index))/2.
    sig_z_0(source_index,subsource_index)=sig_z_00(source_index,subsource_index)+az(source_index,subsource_index)*exp(bz(source_index,subsource_index)*log(sig_y_0(source_index,subsource_index)))
    
    
    end subroutine uEMEP_set_dispersion_params_PG     
    