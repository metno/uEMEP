
    subroutine uEMEP_chemistry
    !Routine for doing the chemistry calculations in uEMEP
    
    use uEMEP_definitions

    implicit none
    
    integer i,j,k
    real nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature
    real nox_out,no2_out,o3_out,p_bg_out,p_out
    
    integer t,t_start,t_end
    integer i_source,i_subsource,emep_subsource
    integer i_comp,i_file
    character(256) temp_name
    
    integer i_integral,j_integral

    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Calculating chemistry for NO2 (uEMEP_chemistry)'
	write(unit_logfile,'(A)') '================================================================'
    if (no2_chemistry_scheme_flag.eq.0) then
        write(unit_logfile,'(A)') 'No chemistry used'
    elseif (no2_chemistry_scheme_flag.eq.1) then
        write(unit_logfile,'(A)') 'Photostationary state used'
    elseif (no2_chemistry_scheme_flag.eq.2) then
        write(unit_logfile,'(A)') 'Photochemistry with time scale used'
    elseif (no2_chemistry_scheme_flag.eq.3) then
        write(unit_logfile,'(A)') 'Romberg parameterisation used'
    endif
    
    t_start=1
    t_end=subgrid_dim(t_dim_index)
    i_subsource=1
    emep_subsource=1
    comp_subgrid(:,:,:,no2_index)=0
    comp_subgrid(:,:,:,nox_index)=0
    comp_subgrid(:,:,:,o3_index)=0

    nox_bg=0.;no2_bg=0.;o3_bg=0.;nox_loc=0.;f_no2_loc=0.;J_photo=0.;temperature=0.;
        
    !Calculate the weighted travel time from the totals calculated in uEMEP_subgrid_dispersion
    do t=t_start,t_end
        traveltime_subgrid(:,:,t,1,:)=traveltime_subgrid(:,:,t,1,:)/traveltime_subgrid(:,:,t,2,:)
        !Invert it to get the time scale
        !traveltime_subgrid(:,:,t,1)=1./traveltime_subgrid(:,:,t,1)
        !Set none valid to 12 hours (long time)
        where (traveltime_subgrid(:,:,t,2,:).eq.0) traveltime_subgrid(:,:,t,1,:)=3600.*12.
        !write(*,*) t
        !write(*,*) traveltime_subgrid(:,:,t,1)
    enddo
    
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
    !write(*,'(2i4,<subgrid_dim(t_dim_index)>f6.1)') i,j,traveltime_subgrid(i,j,:,1)/60.
    enddo
    enddo
    
    do t=t_start,t_end
    do j=1,subgrid_dim(y_dim_index)
    do i=1,subgrid_dim(x_dim_index)
    if (use_subgrid(i,j,allsource_index)) then
        
        i_integral=crossreference_target_to_integral_subgrid(i,j,x_dim_index)
        j_integral=crossreference_target_to_integral_subgrid(i,j,y_dim_index)

        J_photo=meteo_subgrid(i_integral,j_integral,t,J_subgrid_index)
        temperature=273.
       
        nox_bg=subgrid(i,j,t,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
        !nox_bg=subgrid(i,j,t,emep_subgrid_index,allsource_index,emep_subsource)
        !nox_bg=comp_subgrid(i,j,t,nox_index)*(14.+16.*2.)/14.
        
        o3_bg=comp_EMEP_subgrid(i,j,t,o3_index)
       
        f_no2_loc=0.
        nox_loc=0.
        
        do i_source=1,n_source_index
        if (calculate_source(i_source)) then
            do i_subsource=1,n_subsource(i_source)
            f_no2_loc=f_no2_loc+emission_factor_conversion(no2_index,i_source,i_subsource)/emission_factor_conversion(nox_index,i_source,i_subsource)*subgrid(i,j,t,local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
            nox_loc=nox_loc+subgrid(i,j,t,local_subgrid_index,i_source,pollutant_loop_back_index(nox_nc_index))
            enddo         
        endif
        enddo
        if (nox_loc.ne.0.0) then
            f_no2_loc=f_no2_loc/nox_loc
        else
            f_no2_loc=0.
        endif
        !no2_bg=max(0.,comp_subgrid(i,j,t,no2_index)-f_no2_loc*subgrid(i,j,t,emep_local_subgrid_index,allsource_index,emep_subsource))
        no2_bg=comp_EMEP_subgrid(i,j,t,no2_index)*subgrid(i,j,t,emep_nonlocal_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))/subgrid(i,j,t,emep_subgrid_index,allsource_index,pollutant_loop_back_index(nox_nc_index))
        !no2_bg=comp_subgrid(i,j,t,no2_index)
        
       !if (J_photo.ne.0) then
        !write(*,'(A,3I6,f5.2,6f12.3)') 'IN : ',i,j,t,no2_bg/nox_bg,nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo
       !endif
        
        !if (nox_loc+nox_bg.lt.0.) write(*,*) t,nox_loc,nox_bg,nox_loc+nox_bg
       
        
        if (no2_chemistry_scheme_flag.eq.0) then
            nox_out=nox_bg+nox_loc
            no2_out=no2_bg+nox_loc*f_no2_loc
            o3_out=o3_bg
        elseif (no2_chemistry_scheme_flag.eq.1) then
            call uEMEP_photostationary_NO2(nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,nox_out,no2_out,o3_out,p_bg_out,p_out)
        elseif (no2_chemistry_scheme_flag.eq.2) then
            call uEMEP_phototimescale_NO2(nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,traveltime_subgrid(i,j,t,1,pollutant_loop_back_index(nox_index)),nox_out,no2_out,o3_out,p_bg_out,p_out)
            !write(*,'(7f8.2,f12.2,2f8.2)') nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,traveltime_subgrid(i,j,t,1),no2_out/nox_out,o3_out/o3_bg
        elseif (no2_chemistry_scheme_flag.eq.3) then
            call uEMEP_Romberg_NO2(nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc,nox_out,no2_out,o3_out)        
        endif
        
        !write(*,*) nox_out-subgrid(i,j,t,total_subgrid_index,allsource_index,1)
        
        comp_subgrid(i,j,t,o3_nc_index)=o3_out
        comp_subgrid(i,j,t,no2_nc_index)=no2_out
        comp_subgrid(i,j,t,nox_nc_index)=nox_out
       
        
        !if (J_photo.ne.0) then
        !write(*,'(A,3I6,f5.2,6f12.3)') 'OUT: ',i,j,t,no2_out/nox_out,nox_out,no2_out,o3_out,p_bg_out,p_out
        !write(*,'(A,5I6,f5.2,3f12.3,1es12.2,2f12.3)') 'IN: ',i,j,i_integral,j_integral,t,no2_bg/nox_bg,nox_bg,no2_bg,o3_bg,J_photo,p_bg_out,p_out
       ! write(*,*)
        !endif
        !if (no2_out.gt.nox_out) write(*,*)no2_out,nox_out,nox_bg,nox_loc
    !else
    !    comp_subgrid(i,j,t,:)=0.
    !endif
    
    else
        comp_subgrid(i,j,t,o3_nc_index)=NODATA_value
        comp_subgrid(i,j,t,no2_nc_index)=NODATA_value
        comp_subgrid(i,j,t,nox_nc_index)=NODATA_value
        
    endif
    
    enddo
    enddo
    enddo
    

    
    end subroutine uEMEP_chemistry
   
    
    subroutine uEMEP_photostationary_NO2(nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,nox_out,no2_out,o3_out,p_bg_out,p_out)

    !use uEMEP_definitions
    
    implicit none
    
    real, intent(in) :: nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature
    real, intent(out) :: nox_out,no2_out,o3_out,p_bg_out,p_out
    integer no2_i,no_i,nox_i,o3_i,ox_i,nox_bg_i,no2_bg_i,n_i
    parameter (n_i=7)
    real Na,Na_fac,k1
    real mass(n_i)
    real mmass(n_i)
    real mol(n_i)
    real f_no2,f_ox,Jd,fac_sqrt
    real :: min_nox=1.0e-6

    DATA mmass /46.,30.,46.,48.,47.,46.,46./

    no2_i=1;no_i=2;nox_i=3;o3_i=4;ox_i=5;nox_bg_i=6;no2_bg_i=7
           

    Na=6.022e23        !(molecules/mol)
    Na_fac=Na/1.0e12   !Conversion from ug/m3 to molecules/cm3 included


    k1=1.4e-12*exp(-1310./temperature); !(cm^3/s) and temperature in Kelvin

    !mmass(1:n_i)=(/46.,30.,46.,48.,47.,46.,46.,46./)
    mass(1:n_i)=0.
    mol(1:n_i)=0.

    !Test for 0 NOx. If so leave the routine
    mass(nox_i)=nox_loc+nox_bg   
!    if (mass(nox_i).eq.0.) then
    if (mass(nox_i).le.min_nox) then
        nox_out=0.
        no2_out=0.
        o3_out=o3_bg
        return
    endif
    
    !Check the photostationary assumption for the input data
    mass(nox_i)=nox_bg
    mass(no2_i)=no2_bg
    mass(o3_i)=o3_bg
    mol=mass/mmass*Na_fac !(molecules per cm3)
    mol(ox_i)=mol(o3_i)+mol(no2_i)
    mol(no_i)=max(0.0,mol(nox_i)-mol(no2_i))
    p_bg_out=J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
    !if (J_photo.ne.0.) write(*,*) p_bg_out,J_photo,mol(no2_i),k1,mol(o3_i),mol(no_i)
    
    !Add the local contribution for calculation
    mass(nox_i)=nox_loc+nox_bg
    mass(no2_i)=f_no2_loc*nox_loc+no2_bg
    mass(o3_i)=o3_bg
    !mass(ox_i)=o3_bg+mass(no2_i)
    mol=mass/mmass*Na_fac !(molecules per cm3)

    mol(ox_i)=mol(o3_i)+mol(no2_i)
    mol(no_i)=max(0.0,mol(nox_i)-mol(no2_i))
    
    f_no2=mol(no2_i)/mol(nox_i)
    f_ox=mol(ox_i)/mol(nox_i)
    
    !Test the photostationary state for the input data
    if (J_photo.ne.0.and.mol(no_i).ne.0..and.mol(o3_i).ne.0.) then
        p_out=J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
    else
        p_out=mol(no2_i)/mol(nox_i)
    endif

    !Set the photolysis rate
    Jd=J_photo/k1/mol(nox_i)
    
    !Calculate fraction of NO2 in photostationary state
    fac_sqrt=max(0.0,(1+f_ox+Jd)**2-4*f_ox)
    !if (J_photo.ne.0) then
        f_no2=0.5*((1+f_ox+Jd)-sqrt(fac_sqrt))
    !else
    !    f_no2=min(1.0,f_ox)
    !endif
    
    !write(*,'(A,9ES12.1)') 'MOL:  ',mol(nox_i),mol(no2_i),mol(o3_i),mol(ox_i),f_no2,f_ox,Jd,f_no2,p_bg_out

    !Convert back to mass
    mol(no2_i)=f_no2*mol(nox_i);
    mol(o3_i)=max(0.,mol(ox_i)-mol(no2_i))  !Rounding errors possible
    !mol(no_i)=max(0.,mol(nox_i)-mol(no2_i)) !Rounding errors possible
    mass=mol*mmass/Na_fac !(ug/m3)


    no2_out=mass(no2_i)
    nox_out=mass(nox_i)
    o3_out=mass(o3_i)

    !Check output
    if (J_photo.ne.0.and.mol(no_i).ne.0..and.mol(o3_i).ne.0.) then
        p_out=J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
    else
        p_out=mol(no2_i)/mol(nox_i)
    endif
    

    !write(*,'(A,9ES12.1)') 'MASS: ',mass(nox_i),mass(no2_i),mass(o3_i),mass(ox_i),f_no2,f_ox,Jd,f_no2,p_out

    end subroutine uEMEP_photostationary_NO2
    
    
    subroutine uEMEP_phototimescale_NO2(nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,time_scale,nox_out,no2_out,o3_out,p_bg_out,p_out)

    !use uEMEP_definitions
    
    implicit none
    
    real, intent(in) :: nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,time_scale
    real, intent(out) :: nox_out,no2_out,o3_out,p_bg_out,p_out
    integer no2_i,no_i,nox_i,o3_i,ox_i,nox_bg_i,no2_bg_i,n_i
    parameter (n_i=7)
    real Na,Na_fac,k1
    real mass(n_i)
    real mmass(n_i)
    real mol(n_i)
    real f_no2,f_ox,Jd,Jd_bg,fac_sqrt
    real :: min_nox=1.0e-6
    real c,b,BB,td,f_no2_0,f_no2_ps
    complex(4) AA
    real p_tot_out,f_ox_bg,f_no2_bg_ps,f_no2_bg
    DATA mmass /46.,30.,46.,48.,47.,46.,46./

    no2_i=1;no_i=2;nox_i=3;o3_i=4;ox_i=5;nox_bg_i=6;no2_bg_i=7
           

    Na=6.022e23        !(molecules/mol)
    Na_fac=Na/1.0e12   !Conversion from ug/m3 to molecules/cm3 included


    k1=1.4e-12*exp(-1310./temperature); !(cm^3/s) and temperature in Kelvin

    !mmass(1:n_i)=(/46.,30.,46.,48.,47.,46.,46.,46./)
    mass(1:n_i)=0.
    mol(1:n_i)=0.

    !Test for 0 NOx. If so leave the routine
    mass(nox_i)=nox_loc+nox_bg   
!    if (mass(nox_i).eq.0.) then
    if (mass(nox_i).le.min_nox) then
        nox_out=0.
        no2_out=0.
        o3_out=o3_bg
        return
    endif
    
    !Check the photostationary assumption for the input data
    mass(nox_i)=nox_bg
    mass(no2_i)=no2_bg
    mass(o3_i)=o3_bg
    mol=mass/mmass*Na_fac !(molecules per cm3)

    mol(ox_i)=mol(o3_i)+mol(no2_i)
    mol(no_i)=max(0.0,mol(nox_i)-mol(no2_i))
    f_ox_bg=mol(ox_i)/mol(nox_i)
    Jd_bg=J_photo/k1/mol(nox_i)
    f_no2_bg_ps=0.5*((1+f_ox_bg+Jd_bg)-sqrt((1+f_ox_bg+Jd_bg)**2-4*f_ox_bg))
    f_no2_bg=mol(no2_i)/mol(nox_i)
    p_bg_out=f_no2_bg/f_no2_bg_ps
    !p_bg_out=J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
    !if (J_photo.ne.0.) write(*,*) p_bg_out,J_photo,mol(no2_i),k1,mol(o3_i),mol(no_i)
    
    !Add the local contribution for calculation
    mass(nox_i)=nox_loc+nox_bg
    mass(no2_i)=f_no2_loc*nox_loc+no2_bg
    mass(o3_i)=o3_bg
    !mass(ox_i)=o3_bg+mass(no2_i)
    mol=mass/mmass*Na_fac !(molecules per cm3)
    mol(ox_i)=mol(o3_i)+mol(no2_i)
    mol(no_i)=max(0.0,mol(nox_i)-mol(no2_i))
    
    f_no2=mol(no2_i)/mol(nox_i)
    f_ox=mol(ox_i)/mol(nox_i)

    !Set the photolysis rate
    Jd=J_photo/k1/mol(nox_i)

    !Calculate photostationary for total nox, ox
    f_no2_ps=0.5*((1+f_ox+Jd)-sqrt((1+f_ox+Jd)**2-4*f_ox))
    p_tot_out=f_no2/f_no2_ps

    !Test the photostationary state for the input data
    !if (J_photo.ne.0.and.mol(no_i).ne.0..and.mol(o3_i).ne.0.) then
    !    p_out=J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
    !else
    !    p_out=mol(no2_i)/mol(nox_i)
    !endif

    
    !Calculate fraction of NO2 based on the time scale
    !fac_sqrt=max(0.0,(1+f_ox+Jd)**2-4*f_ox)
    c=f_ox
    b=1+f_ox+Jd
    BB=sqrt(max(0.0,b**2-4.*c))!max avoids roundoff errors
    td=time_scale*k1*mol(nox_i)
    f_no2_0=f_no2

    !if (J_photo.ne.0) then
        !f_no2=0.5*((1+f_ox+Jd)-sqrt(fac_sqrt))
    !if (b.lt.100.)
    AA=clog(cmplx((BB+b-2.*f_no2_0)/(BB-b+2.*f_no2_0)))
    f_no2=real(-BB/2.*((exp(AA+BB*td)-1.)/(exp(AA+BB*td)+1.))+b/2.)
    !if (BB*td.gt.50.) f_no2=-BB/2.+b/2.
    if (isnan(f_no2)) f_no2=-BB/2.+b/2.
    
    f_no2_ps=0.5*((1+f_ox+Jd)-sqrt((1+f_ox+Jd)**2-4*f_ox))
    p_out=f_no2/f_no2_ps
    !write(*,*) p_bg_out,p_tot_out,p_out,AA,BB,b,td,exp(AA+BB*td),f_ox
    !write(*,*) o3_bg,no2_bg,nox_bg
    !write(*,*) c,b,BB,AA,(BB/2.+b/2.-f_no2_0)/(BB/2.-b/2.+f_no2_0),f_no2
    
    !else
    !    f_no2=1.0
    !endif
    
    !write(*,'(A,9ES12.1)') 'MOL:  ',mol(nox_i),mol(no2_i),mol(o3_i),mol(ox_i),f_no2,f_ox,Jd,f_no2,p_bg_out

    !Convert back to mass
    mol(no2_i)=f_no2*mol(nox_i);
    mol(o3_i)=max(0.,mol(ox_i)-mol(no2_i))  !Rounding errors possible
    !mol(no_i)=max(0.,mol(nox_i)-mol(no2_i)) !Rounding errors possible
    mass=mol*mmass/Na_fac !(ug/m3)


    no2_out=mass(no2_i)
    nox_out=mass(nox_i)
    o3_out=mass(o3_i)
    
    if (isnan(no2_out)) then
        write(*,'(8a12)') 'nox_bg','no2_bg','o3_bg','nox_loc','f_no2_loc','J_photo','temperature','time_scale'
        write(*,'(8es12.2)') nox_bg,no2_bg,o3_bg,nox_loc,f_no2_loc,J_photo,temperature,time_scale
        write(*,'(4a12)') 'f_no2','BB','b','b**2-4.*c'
        write(*,'(4es12.2)') f_no2,BB,b,b**2-4.*c
        stop
    endif
    

    !Check output
    !if (J_photo.ne.0.and.mol(no_i).ne.0..and.mol(o3_i).ne.0.) then
    !    p_out=J_photo*mol(no2_i)/k1/mol(o3_i)/mol(no_i)
    !else
    !    p_out=mol(no2_i)/mol(nox_i)
    !endif
    

    !write(*,'(A,9ES12.1)') 'MASS: ',mass(nox_i),mass(no2_i),mass(o3_i),mass(ox_i),f_no2,f_ox,Jd,f_no2,p_out

    end subroutine uEMEP_phototimescale_NO2
    
    subroutine uEMEP_Romberg_NO2(nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc,nox_out,no2_out,o3_out)
    
    implicit none
    
    real, intent(in) :: nox_bg,no2_bg,nox_loc,o3_bg,f_no2_loc
    real, intent(out) :: nox_out,no2_out,o3_out
    
    !From Norwegian obs fit
    real :: a_rom=21
    real :: b_rom=27
    real :: c_rom=0.23
    real ox_init,no2_init
   
    nox_out=nox_bg+nox_loc
    no2_out=a_rom*nox_out/(nox_out+b_rom)+nox_out*c_rom
    no2_init=no2_bg+f_no2_loc*nox_loc
    
    !Small adjustments for molecular weights
    ox_init=no2_init*47./46.+o3_bg*47./48.
    o3_out=ox_init*48./47.-no2_out*48./46.
    
    end subroutine uEMEP_Romberg_NO2
    