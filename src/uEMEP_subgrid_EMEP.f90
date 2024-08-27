module subgrid_emep

    use uEMEP_definitions
    use uemep_constants, only: pi
    use uemep_configuration
    use mod_lambert_projection, only: LL2PROJ, PROJ2LL

    implicit none
    private

    public :: uEMEP_subgrid_EMEP,nlreg_uEMEP_subgrid_EMEP_from_in_region

contains

!==========================================================================
!   uEMEP_subgrid_EMEP
!   Bruce Rolstad Denby
!   MET Norway
!
!   This routine calculates the local and nonlocal contribution of the EMEP
!   grid at each subgrid point using a moving window. The local contribution
!   is then later removed and replaced by the local dispersion or is used
!   to redistribute concentrations. The moving window must take into account
!   edges of the EMEP grid. The nonlocal_correction deals with grid contributions
!   outside the central EMEP grid
!
!   The following options are available:
!   EMEP_grid_interpolation_flag=0 is no interpolation, just uses the EMEP grid it is in
!   EMEP_grid_interpolation_flag=1 is area weighted (same as bilinear interpolation, quickest)
!   EMEP_grid_interpolation_flag=2 is emission subgrid weighted (slowest)
!   EMEP_grid_interpolation_flag=3 is emission aggregated to integral weighted (recommended, faster)
!   EMEP_grid_interpolation_flag=4 is proxy dispersion integral weighted (similar to 3 but needs integral dispersion calculation)
!
!==========================================================================

    subroutine uEMEP_subgrid_EMEP

        integer i,j
        integer ii,jj,tt,iii,jjj
        real, allocatable :: weighting_nc(:,:,:,:),weighting_subgrid(:,:,:,:,:)
        real, allocatable :: total_weighting_nc(:,:,:,:,:),proxy_weighting_nc(:,:,:,:,:)
        real, allocatable :: area_weighting_nc(:,:,:,:,:,:)
        integer i_nc_start,i_nc_end,j_nc_start,j_nc_end
        integer i_start,i_end,j_start,j_end,t_start,t_end
        real xpos_min,xpos_max,ypos_min,ypos_max
        real xpos_area_min,xpos_area_max,ypos_area_min,ypos_area_max
        real xpos_min2,xpos_max2,ypos_min2,ypos_max2
        real xpos_min3,xpos_max3,ypos_min3,ypos_max3
        real xpos_area_min2,xpos_area_max2,ypos_area_min2,ypos_area_max2
        integer i_nc,j_nc
        integer emep_subsource
        real, allocatable :: nonlocal_correction(:,:,:)
        real, allocatable :: nonlocal_correction_average(:,:)
        integer i_source
        integer ii_nc,jj_nc,ii_w,jj_w
        integer :: n_weight=3,ii_w0=2,jj_w0=2
        integer weighting_subgrid_dim(2,n_source_index)
        integer i_cross,j_cross
        integer, allocatable :: crossreference_weighting_to_emep_subgrid(:,:,:,:)
        integer i_w_c,j_w_c
        integer i_nc_c,j_nc_c
        integer count
        integer ii_start,ii_end,jj_start,jj_end
        integer iii_start,iii_end,jjj_start,jjj_end
        real xpos_limit,ypos_limit
        real xpos_limit2,ypos_limit2
        real xpos_subgrid,ypos_subgrid
        real xpos_emission_subgrid,ypos_emission_subgrid
        real xpos_integral_subgrid,ypos_integral_subgrid
        real EMEP_grid_interpolation_size_sqr
        integer :: tt_dim=1
        integer ii_nc_w0,jj_nc_w0,iii_nc_w,jjj_nc_w
        integer ii_nc_w,jj_nc_w

        real weighting_val,weighting_val3
        integer i_pollutant,i_loop
        integer i_sp,ii_sp
        real, allocatable :: EMEP_local_contribution(:,:,:,:)
        real, allocatable :: EMEP_local_contribution_from_in_region(:,:,:,:)

        integer n_weight_nc_x,n_weight_nc_y
        real dgrid_lf_offset_x,dgrid_lf_offset_y
        real amod_temp
        real EMEP_grid_interpolation_size_saved
        real local_fraction_grid_size_scaling_temp
        real weight_check

        real xpos_lf_area_min,xpos_lf_area_max,ypos_lf_area_min,ypos_lf_area_max
        logical :: first_interpolate_lf=.false.
        logical :: set_lf_offset_to_0=.true.

        integer lf_size_index

        real, allocatable :: temp_subgrid(:,:,:,:)
        real, allocatable :: temp_comp_EMEP_subgrid(:,:)
        real, allocatable :: temp_species_EMEP_subgrid(:,:,:)

        real, allocatable :: temp_EMEP(:,:,:,:,:,:)
        real, allocatable :: temp_EMEP_from_in_region(:,:,:,:,:,:)
        real, allocatable :: temp_species_EMEP(:,:,:,:,:)
        real, allocatable :: temp_comp_EMEP(:,:,:,:)


        real distance_grid_x,distance_grid_y
        integer temp_count
        character(len=:), allocatable :: fmt

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Distributing EMEP concentrations to subgrids  (uEMEP_subgrid_EMEP)'
        write(unit_logfile,'(A)') '================================================================'

        !Initialise subgrids to be written to
        subgrid(:,:,:,emep_subgrid_index,:,:)=0
        subgrid(:,:,:,emep_frac_subgrid_index,:,:)=0
        subgrid(:,:,:,emep_local_subgrid_index,:,:)=0
        subgrid(:,:,:,emep_nonlocal_subgrid_index,:,:)=0
        comp_EMEP_subgrid(:,:,:,:)=0
        orig_EMEP_subgrid(:,:,:,:)=0
        if (save_emep_species.or.save_seasalt) species_EMEP_subgrid(:,:,:,:,:)=0

        do i_source=1,n_source_index
            if (calculate_source(i_source).or.calculate_EMEP_source(i_source)) then
                write(unit_logfile,'(2A)') 'Calculating for EMEP source: ',trim(uEMEP_to_EMEP_emis_sector_str(i_source))
            endif
            if (save_EMEP_source(i_source)) then
                write(unit_logfile,'(2A)') 'Saving for EMEP source: ',trim(uEMEP_to_EMEP_emis_sector_str(i_source))
            endif
        enddo


        !Check if the additional EMEP calculation is to be carried out and set parameters
        EMEP_grid_interpolation_size_saved=EMEP_grid_interpolation_size
        lc_local_nc_index=lc_local_nc_loop_index(local_fraction_grid_for_EMEP_grid_interpolation)
        local_fraction_grid_size_scaling_temp=local_fraction_grid_size_scaling
        lf_size_index=1
        if (calculate_EMEP_additional_grid_flag) then
            EMEP_grid_interpolation_size=EMEP_additional_grid_interpolation_size
            local_fraction_grid_size_scaling_temp=local_fraction_additional_grid_size_scaling
            lc_local_nc_index=lc_local_nc_loop_index(local_fraction_grid_for_EMEP_additional_grid_interpolation)
            lf_size_index=2
            write(unit_logfile,'(A,i)') 'Calculating additional EMEP concentrations to subgrids, index:',lc_local_nc_index
        else
            write(unit_logfile,'(A,i)') 'Calculating EMEP concentrations to subgrids, index:',lc_local_nc_index
        endif

        !Set value used later
        EMEP_grid_interpolation_size_sqr=EMEP_grid_interpolation_size*EMEP_grid_interpolation_size

        !Time dimension when external time loop is used
        tt_dim=1

        !Allocate nonlocal_correction
        if (.not.allocated(nonlocal_correction)) allocate (nonlocal_correction(tt_dim,n_source_index,n_pollutant_loop))
        if (.not.allocated(nonlocal_correction_average)) allocate (nonlocal_correction_average(n_source_index,n_pollutant_loop))

        !There are no subsources in EMEP
        emep_subsource=1


        !Initialise a diagnostic check variable for the weighting
        nonlocal_correction_average=0.

        !Save the original EMEP directly to subgrid for comparison and visualisation purposes on the target subgrid
        write(unit_logfile,'(A)')'Calculating original EMEP subgrid using nearest neighbour interpolation'

        do j=1,subgrid_dim(y_dim_index)
            do i=1,subgrid_dim(x_dim_index)
                if (use_subgrid(i,j,allsource_index)) then
                    ii=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
                    jj=crossreference_target_to_emep_subgrid(i,j,y_dim_index)

                    !Nearest neighbour interpolate the EMEP compounds to subgrid
                    if (ii.ge.1.and.ii.le.dim_length_nc(x_dim_nc_index).and.jj.ge.1.and.jj.le.dim_length_nc(y_dim_nc_index)) then

                        do i_pollutant=1,n_emep_pollutant_loop
                            do i_loop=1,n_pollutant_compound_loop(i_pollutant)
                                !write(*,*) trim(pollutant_file_str(pollutant_compound_loop_index(i_pollutant,i_loop)))
                                orig_EMEP_subgrid(i,j,:,pollutant_compound_loop_index(i_pollutant,i_loop))=comp_var3d_nc(ii,jj,:,pollutant_compound_loop_index(i_pollutant,i_loop))
                            enddo
                        enddo

                    endif

                endif
            enddo
        enddo

        !Loop through subgrid and find those subgrids within EMEP grids and allocate concentrations directly from EMEP grids.
        if (EMEP_grid_interpolation_flag.eq.0.or.EMEP_grid_interpolation_flag.eq.5) then

            write(unit_logfile,'(A)')'Calculating EMEP local subgrid contribution using nearest neighbour interpolation'

            jj_start=-1-floor(0.5*(EMEP_grid_interpolation_size-1.))
            ii_start=-1-floor(0.5*(EMEP_grid_interpolation_size-1.))
            jj_end=1+floor(0.5*(EMEP_grid_interpolation_size-1.))
            ii_end=1+floor(0.5*(EMEP_grid_interpolation_size-1.))

            !Window does not extend outside the grid, since it is centred on the EMEP grid
            if (EMEP_grid_interpolation_size.le.1) then
                jj_start=0
                ii_start=0
                jj_end=0
                ii_end=0
            endif

            write(unit_logfile,'(A,4i)') 'LF loop (ii_start,ii_end,jj_start,jj_end): ',ii_start,ii_end,jj_start,jj_end

            ii_nc_w0=xdist_centre_nc
            jj_nc_w0=ydist_centre_nc

            !Set weighting indexes
            n_weight=3+2*floor((EMEP_grid_interpolation_size-1.)*0.5)
            ii_w0=1+floor(n_weight*.5)
            jj_w0=1+floor(n_weight*.5)

            n_weight_nc_x=xdist_centre_nc*2-1
            n_weight_nc_y=ydist_centre_nc*2-1

            !Set the size of the region surounding the target grid that is searched
            xpos_limit=dgrid_nc(lon_nc_index)/2.*EMEP_grid_interpolation_size*local_fraction_grid_size_scaling_temp
            ypos_limit=dgrid_nc(lat_nc_index)/2.*EMEP_grid_interpolation_size*local_fraction_grid_size_scaling_temp

            do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                    !Do the calculation everywhere when EMEP_grid_interpolation_flag=5 since it uses the subgrid values later
                    if (use_subgrid(i,j,allsource_index).or.EMEP_grid_interpolation_flag.eq.5) then

                        i_nc=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
                        j_nc=crossreference_target_to_emep_subgrid(i,j,y_dim_index)

                        if (i_nc.ge.1.and.i_nc.le.dim_length_nc(x_dim_nc_index).and.j_nc.ge.1.and.j_nc.le.dim_length_nc(y_dim_nc_index)) then

                            !Read from the local fraction file
                            subgrid(i,j,:,emep_subgrid_index,:,:)=var3d_nc(i_nc,j_nc,:,conc_nc_index,1:n_source_index,:)
                            !subgrid(i,j,:,emep_local_subgrid_index,:,:)=var3d_nc(i_nc,j_nc,:,local_nc_index,:,:)

                            !Centre of grid
                            xpos_subgrid=var1d_nc(i_nc,lon_nc_index)
                            ypos_subgrid=var1d_nc(j_nc,lat_nc_index)

                            !Set the edges of the search area surounding the target grid
                            xpos_area_min=xpos_subgrid-xpos_limit
                            xpos_area_max=xpos_subgrid+xpos_limit
                            ypos_area_min=ypos_subgrid-ypos_limit
                            ypos_area_max=ypos_subgrid+ypos_limit

                            !Set the offset to the centre of the local fraction grid when using larger local fraction grids
                            !Requires knowledge of the total EMEP grid index so uses dim_start_nc
                            amod_temp=amod(real(dim_start_nc(x_dim_nc_index)-1+i_nc),local_fraction_grid_size_scaling_temp)
                            if (amod_temp.eq.0) amod_temp=local_fraction_grid_size_scaling_temp
                            dgrid_lf_offset_x=0.5-(amod_temp-0.5)/local_fraction_grid_size_scaling_temp
                            amod_temp=amod(real(dim_start_nc(y_dim_nc_index)-1+j_nc),local_fraction_grid_size_scaling_temp)
                            if (amod_temp.eq.0) amod_temp=local_fraction_grid_size_scaling_temp
                            dgrid_lf_offset_y=0.5-(amod_temp-0.5)/local_fraction_grid_size_scaling_temp

                            if (set_lf_offset_to_0) then
                                dgrid_lf_offset_x=0
                                dgrid_lf_offset_y=0
                            endif

                            !dgrid_lf_offset_x=0
                            !dgrid_lf_offset_y=0

                            !Calculate the local fraction contribution from within the moving window, limited on the edges.
                            subgrid(i,j,:,emep_local_subgrid_index,:,:)=0

                            weight_check=0
                            do jj=jj_start,jj_end
                                do ii=ii_start,ii_end

                                    ii_nc=ii+i_nc
                                    jj_nc=jj+j_nc

                                    ii_nc_w=ii+ii_nc_w0
                                    jj_nc_w=jj+jj_nc_w0

                                    ii_w=ii+ii_w0
                                    jj_w=jj+jj_w0

                                    !Put in a limit
                                    if (ii_nc.ge.1.and.ii_nc.le.dim_length_nc(x_dim_nc_index).and.jj_nc.ge.1.and.jj_nc.le.dim_length_nc(y_dim_nc_index) &
                                        .and.ii_nc_w.ge.1.and.ii_nc_w.le.n_weight_nc_x.and.jj_nc_w.ge.1.and.jj_nc_w.le.n_weight_nc_y) then

                                        xpos_lf_area_min=var1d_nc(i_nc,lon_nc_index)+(ii-1/2.+dgrid_lf_offset_x)*dgrid_nc(lon_nc_index)*local_fraction_grid_size_scaling_temp
                                        xpos_lf_area_max=var1d_nc(i_nc,lon_nc_index)+(ii+1/2.+dgrid_lf_offset_x)*dgrid_nc(lon_nc_index)*local_fraction_grid_size_scaling_temp
                                        ypos_lf_area_min=var1d_nc(j_nc,lat_nc_index)+(jj-1/2.+dgrid_lf_offset_y)*dgrid_nc(lat_nc_index)*local_fraction_grid_size_scaling_temp
                                        ypos_lf_area_max=var1d_nc(j_nc,lat_nc_index)+(jj+1/2.+dgrid_lf_offset_y)*dgrid_nc(lat_nc_index)*local_fraction_grid_size_scaling_temp

                                        !Set the edges to an EMEP grid surounding the EMEP grid being assessed
                                        !xpos_min=max(xpos_area_min,var1d_nc(ii_nc,lon_nc_index)-dgrid_nc(lon_nc_index)/2.*local_fraction_grid_size_scaling_temp+dgrid_lf_offset_x*dgrid_nc(lon_nc_index)*local_fraction_grid_size_scaling_temp)
                                        !xpos_max=min(xpos_area_max,var1d_nc(ii_nc,lon_nc_index)+dgrid_nc(lon_nc_index)/2.*local_fraction_grid_size_scaling_temp+dgrid_lf_offset_x*dgrid_nc(lon_nc_index)*local_fraction_grid_size_scaling_temp)
                                        !ypos_min=max(ypos_area_min,var1d_nc(jj_nc,lat_nc_index)-dgrid_nc(lat_nc_index)/2.*local_fraction_grid_size_scaling_temp+dgrid_lf_offset_y*dgrid_nc(lat_nc_index)*local_fraction_grid_size_scaling_temp)
                                        !ypos_max=min(ypos_area_max,var1d_nc(jj_nc,lat_nc_index)+dgrid_nc(lat_nc_index)/2.*local_fraction_grid_size_scaling_temp+dgrid_lf_offset_y*dgrid_nc(lat_nc_index)*local_fraction_grid_size_scaling_temp)
                                        xpos_min=max(xpos_area_min,xpos_lf_area_min)
                                        xpos_max=min(xpos_area_max,xpos_lf_area_max)
                                        ypos_min=max(ypos_area_min,ypos_lf_area_min)
                                        ypos_max=min(ypos_area_max,ypos_lf_area_max)

                                        !Calculate area weighting
                                        if (xpos_max.gt.xpos_min.and.ypos_max.gt.ypos_min) then
                                            weighting_val=(ypos_max-ypos_min)*(xpos_max-xpos_min)/dgrid_nc(lon_nc_index)/dgrid_nc(lat_nc_index)/local_fraction_grid_size_scaling_temp/local_fraction_grid_size_scaling_temp
                                        else
                                            weighting_val=0.
                                        endif

                                        !weighting_val=1.
                                        weight_check=weight_check+weighting_val

                                        !write(*,'(2i,5f12.2,f12.4)') ii,jj,weighting_val,xpos_lf_area_min-xpos_subgrid,xpos_lf_area_max-xpos_subgrid &
                                        !    ,xpos_area_min-xpos_subgrid,xpos_area_max-xpos_subgrid,dgrid_lf_offset_x

                                        subgrid(i,j,:,emep_local_subgrid_index,:,:)=subgrid(i,j,:,emep_local_subgrid_index,:,:) &
                                            +lc_var3d_nc(ii_nc_w,jj_nc_w,i_nc,j_nc,:,lc_local_nc_index,1:n_source_index,:)*weighting_val


                                    endif

                                enddo
                            enddo
                            !write(*,*) 'Check: ',(subgrid(i,j,1,emep_local_subgrid_index,traffic_index,1)),(subgrid(i,j,1,emep_subgrid_index,allsource_index,1))

                            !Interpolate the other EMEP compounds as well to subgrid in the same way. Read from the normal EMEP file
                            !comp_var3d_nc(ii,jj,:,pollutant_compound_loop_index(i_pollutant,i_loop))
                            do i_pollutant=1,n_emep_pollutant_loop
                                do i_loop=1,n_pollutant_compound_loop(i_pollutant)
                                    comp_EMEP_subgrid(i,j,:,pollutant_compound_loop_index(i_pollutant,i_loop))=comp_var3d_nc(i_nc,j_nc,:,pollutant_compound_loop_index(i_pollutant,i_loop))
                                enddo
                            enddo

                            if (save_emep_species.or.save_seasalt) then
                                do i_sp=1,n_species_loop_index
                                    ii_sp=species_loop_index(i_sp)
                                    do i_loop=1,n_pmxx_sp_index
                                        species_EMEP_subgrid(i,j,:,i_loop,i_sp)=species_var3d_nc(i_nc,j_nc,:,i_loop,i_sp)
                                    enddo
                                enddo
                            endif

                        endif

                    endif
                enddo
            enddo


            !Set the non-local for each source individually
            subgrid(:,:,:,emep_nonlocal_subgrid_index,:,:)=subgrid(:,:,:,emep_subgrid_index,:,:)-subgrid(:,:,:,emep_local_subgrid_index,:,:)

            if (calculate_deposition_flag) then
                subgrid(i,j,:,drydepo_nonlocal_subgrid_index,:,:)=subgrid(i,j,:,drydepo_nonlocal_subgrid_index,:,:) &
                    *subgrid(i,j,:,emep_nonlocal_subgrid_index,:,:)/subgrid(i,j,:,emep_subgrid_index,:,:)
                subgrid(i,j,:,wetdepo_nonlocal_subgrid_index,:,:)=subgrid(i,j,:,wetdepo_nonlocal_subgrid_index,:,:) &
                    *subgrid(i,j,:,emep_nonlocal_subgrid_index,:,:)/subgrid(i,j,:,emep_subgrid_index,:,:)
            endif

        endif

        !Loop through subgrid and find those subgrids within EMEP grids and allocate concentrations directly from EMEP grids.
        if (EMEP_grid_interpolation_flag.eq.6) then

            write(unit_logfile,'(A)')'Calculating EMEP local subgrid contribution using method 6. First contribution to grid then area interpolation'

            if (.not.allocated(temp_EMEP)) allocate (temp_EMEP(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index),num_var_nc,n_source_nc_index,n_pollutant_loop))
            if (.not.allocated(temp_comp_EMEP)) allocate (temp_comp_EMEP(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index),n_compound_index))
            temp_EMEP=0
            temp_comp_EMEP=0

            if (save_emep_species.or.save_seasalt) then
                if (.not.allocated(temp_species_EMEP)) allocate (temp_species_EMEP(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),dim_length_nc(time_dim_nc_index),n_pmxx_sp_index,n_species_loop_index))
                temp_species_EMEP=0
            endif

            !Set the extent of the LF grid to be assessed
            jj_start=-1-floor(0.5*(EMEP_grid_interpolation_size-1.))
            ii_start=-1-floor(0.5*(EMEP_grid_interpolation_size-1.))
            jj_end=1+floor(0.5*(EMEP_grid_interpolation_size-1.))
            ii_end=1+floor(0.5*(EMEP_grid_interpolation_size-1.))

            !Window does not extend outside the grid, since it is centred on the EMEP grid
            if (EMEP_grid_interpolation_size.le.1) then
                jj_start=0
                ii_start=0
                jj_end=0
                ii_end=0
            endif

            write(unit_logfile,'(A,4i)') 'LF loop (ii_start,ii_end,jj_start,jj_end): ',ii_start,ii_end,jj_start,jj_end

            ii_nc_w0=xdist_centre_nc
            jj_nc_w0=ydist_centre_nc

            !Set weighting indexes
            n_weight=3+2*floor((EMEP_grid_interpolation_size-1.)*0.5)
            ii_w0=1+floor(n_weight*.5)
            jj_w0=1+floor(n_weight*.5)

            n_weight_nc_x=xdist_centre_nc*2-1
            n_weight_nc_y=ydist_centre_nc*2-1

            !Set the size of the region surounding the target grid that is searched
            xpos_limit=dgrid_nc(lon_nc_index)/2.*EMEP_grid_interpolation_size*local_fraction_grid_size_scaling_temp
            ypos_limit=dgrid_nc(lat_nc_index)/2.*EMEP_grid_interpolation_size*local_fraction_grid_size_scaling_temp

            !Set the size of the region for the area interpoaltion
            xpos_limit2=dgrid_nc(lon_nc_index)/2.
            ypos_limit2=dgrid_nc(lat_nc_index)/2.

            !Limits for the area interpolation
            jjj_start=-1
            iii_start=-1
            jjj_end=1
            iii_end=1

            !Loop through the EMEP grids and create the local contribution in the EMEP grid
            do j_nc=1,dim_length_nc(y_dim_nc_index)
                do i_nc=1,dim_length_nc(x_dim_nc_index)

                    !i_nc=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
                    !j_nc=crossreference_target_to_emep_subgrid(i,j,y_dim_index)

                    !if (i_nc.ge.1.and.i_nc.le.dim_length_nc(x_dim_nc_index).and.j_nc.ge.1.and.j_nc.le.dim_length_nc(y_dim_nc_index)) then

                    !Read from the local fraction file
                    temp_EMEP(i_nc,j_nc,:,emep_subgrid_index,1:n_source_index,:)=var3d_nc(i_nc,j_nc,:,conc_nc_index,1:n_source_index,:)

                    !Centre of EMEP grid in lat lon or local coordinates
                    xpos_subgrid=var1d_nc(i_nc,lon_nc_index)
                    ypos_subgrid=var1d_nc(j_nc,lat_nc_index)

                    !Set the edges of the search area surounding the target grid
                    xpos_area_min=xpos_subgrid-xpos_limit
                    xpos_area_max=xpos_subgrid+xpos_limit
                    ypos_area_min=ypos_subgrid-ypos_limit
                    ypos_area_max=ypos_subgrid+ypos_limit

                    !Set the offset to the centre of the local fraction grid when using larger local fraction grids
                    !Requires knowledge of the total EMEP grid index so uses dim_start_nc
                    amod_temp=amod(real(dim_start_nc(x_dim_nc_index)-1+i_nc),local_fraction_grid_size_scaling_temp)
                    if (amod_temp.eq.0) amod_temp=local_fraction_grid_size_scaling_temp
                    dgrid_lf_offset_x=0.5-(amod_temp-0.5)/local_fraction_grid_size_scaling_temp
                    amod_temp=amod(real(dim_start_nc(y_dim_nc_index)-1+j_nc),local_fraction_grid_size_scaling_temp)
                    if (amod_temp.eq.0) amod_temp=local_fraction_grid_size_scaling_temp
                    dgrid_lf_offset_y=0.5-(amod_temp-0.5)/local_fraction_grid_size_scaling_temp

                    if (set_lf_offset_to_0) then
                        dgrid_lf_offset_x=0
                        dgrid_lf_offset_y=0
                    endif


                    !Calculate the local fraction contribution to the EMEP grid centre, limited on the edges.
                    temp_EMEP(i_nc,j_nc,:,emep_local_subgrid_index,:,:)=0

                    weight_check=0
                    do jj=jj_start,jj_end
                        do ii=ii_start,ii_end

                            ii_nc=ii+i_nc
                            jj_nc=jj+j_nc

                            ii_nc_w=ii+ii_nc_w0
                            jj_nc_w=jj+jj_nc_w0

                            ii_w=ii+ii_w0
                            jj_w=jj+jj_w0

                            !Put in a limit
                            if (ii_nc.ge.1.and.ii_nc.le.dim_length_nc(x_dim_nc_index).and.jj_nc.ge.1.and.jj_nc.le.dim_length_nc(y_dim_nc_index) &
                                .and.ii_nc_w.ge.1.and.ii_nc_w.le.n_weight_nc_x.and.jj_nc_w.ge.1.and.jj_nc_w.le.n_weight_nc_y) then

                                xpos_lf_area_min=var1d_nc(i_nc,lon_nc_index)+(ii-1/2.+dgrid_lf_offset_x)*dgrid_nc(lon_nc_index)*local_fraction_grid_size_scaling_temp
                                xpos_lf_area_max=var1d_nc(i_nc,lon_nc_index)+(ii+1/2.+dgrid_lf_offset_x)*dgrid_nc(lon_nc_index)*local_fraction_grid_size_scaling_temp
                                ypos_lf_area_min=var1d_nc(j_nc,lat_nc_index)+(jj-1/2.+dgrid_lf_offset_y)*dgrid_nc(lat_nc_index)*local_fraction_grid_size_scaling_temp
                                ypos_lf_area_max=var1d_nc(j_nc,lat_nc_index)+(jj+1/2.+dgrid_lf_offset_y)*dgrid_nc(lat_nc_index)*local_fraction_grid_size_scaling_temp

                                !Set the edges to an EMEP grid surounding the EMEP grid being assessed
                                xpos_min=max(xpos_area_min,xpos_lf_area_min)
                                xpos_max=min(xpos_area_max,xpos_lf_area_max)
                                ypos_min=max(ypos_area_min,ypos_lf_area_min)
                                ypos_max=min(ypos_area_max,ypos_lf_area_max)

                                !Calculate area weighting
                                if (xpos_max.gt.xpos_min.and.ypos_max.gt.ypos_min) then
                                    weighting_val=(ypos_max-ypos_min)*(xpos_max-xpos_min)/dgrid_nc(lon_nc_index)/dgrid_nc(lat_nc_index)/local_fraction_grid_size_scaling_temp/local_fraction_grid_size_scaling_temp
                                else
                                    weighting_val=0.
                                endif

                                !weighting_val=1.
                                weight_check=weight_check+weighting_val

                                !write(*,'(2i,5f12.2,f12.4)') ii,jj,weighting_val,xpos_lf_area_min-xpos_subgrid,xpos_lf_area_max-xpos_subgrid &
                                !    ,xpos_area_min-xpos_subgrid,xpos_area_max-xpos_subgrid,dgrid_lf_offset_x

                                temp_EMEP(i_nc,j_nc,:,emep_local_subgrid_index,1:n_source_index,:)=temp_EMEP(i_nc,j_nc,:,emep_local_subgrid_index,1:n_source_index,:) &
                                    +lc_var3d_nc(ii_nc_w,jj_nc_w,i_nc,j_nc,:,lc_local_nc_index,1:n_source_index,:)*weighting_val


                            endif

                        enddo
                    enddo

                    !write(*,*) 'Check: ',(temp_EMEP(i_nc,j_nc,1,emep_local_subgrid_index,traffic_index,1)),(temp_EMEP(i_nc,j_nc,1,emep_subgrid_index,allsource_index,1))

                    !Interpolate the other EMEP compounds as well to subgrid in the same way. Read directly from the normal EMEP file
                    !do i_pollutant=1,n_emep_pollutant_loop
                    !do i_loop=1,n_pollutant_compound_loop(i_pollutant)
                    !    temp_comp_EMEP(i_nc,j_nc,:,pollutant_compound_loop_index(i_pollutant,i_loop))=comp_var3d_nc(i_nc,j_nc,:,pollutant_compound_loop_index(i_pollutant,i_loop))
                    !enddo
                    !enddo
                    temp_comp_EMEP(i_nc,j_nc,:,:)=comp_var3d_nc(i_nc,j_nc,:,:)

                    if (save_emep_species.or.save_seasalt) then
                        temp_species_EMEP(i_nc,j_nc,:,:,:)=species_var3d_nc(i_nc,j_nc,:,:,:)
                        !do i_sp=1,n_species_loop_index
                        !ii_sp=species_loop_index(i_sp)
                        !do i_loop=1,n_pmxx_sp_index
                        !    temp_species_EMEP(i_nc,j_nc,:,i_loop,i_sp)=species_var3d_nc(i_nc,j_nc,:,i_loop,i_sp)
                        !enddo
                        !enddo
                    endif

                    !endif

                    !endif
                enddo
            enddo

            subgrid(:,:,:,emep_local_subgrid_index,:,:)=0
            comp_EMEP_subgrid=0
            species_EMEP_subgrid=0

            !Loop through the subgrids and allocate the temp_EMEP grids using area weighting
            do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                    if (use_subgrid(i,j,allsource_index)) then

                        !Assumes it is never on the edge of the EMEP grid as it is not limitted
                        i_nc=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
                        j_nc=crossreference_target_to_emep_subgrid(i,j,y_dim_index)

                        xpos_subgrid=xproj_subgrid(i,j)
                        ypos_subgrid=yproj_subgrid(i,j)

                        !Are xpos_limit2 in the same coordinates?
                        xpos_area_min2=xpos_subgrid-xpos_limit2
                        xpos_area_max2=xpos_subgrid+xpos_limit2
                        ypos_area_min2=ypos_subgrid-ypos_limit2
                        ypos_area_max2=ypos_subgrid+ypos_limit2

                        !Limit the region. This will still allow an EMEP contribution from half a grid away
                        !Same limit is NOT applied on the emissions in the moving window so inconsistent
                        !write(*,'(2i,4e12.2)') i,j,xpos_area_min,xpos_area_max,ypos_area_min,ypos_area_max
                        if (limit_emep_grid_interpolation_region_to_calculation_region) then
                            xpos_area_min=max(xpos_area_min,subgrid_proj_min(x_dim_index)-xpos_subgrid-xpos_limit2)
                            xpos_area_max=min(xpos_area_max,subgrid_proj_max(x_dim_index)-xpos_subgrid+xpos_limit2)
                            ypos_area_min=max(ypos_area_min,subgrid_proj_min(y_dim_index)-ypos_subgrid-ypos_limit2)
                            ypos_area_max=min(ypos_area_max,subgrid_proj_max(y_dim_index)-ypos_subgrid+ypos_limit2)
                        endif
                        !write(*,'(2i,4e12.2)') i,j,xpos_area_min,xpos_area_max,ypos_area_min,ypos_area_max

                        !Set the offset to the centre of the local fraction grid when using larger local fraction grids
                        amod_temp=amod(real(dim_start_nc(x_dim_nc_index)-1+i_nc),local_fraction_grid_size_scaling_temp)
                        if (amod_temp.eq.0) amod_temp=local_fraction_grid_size_scaling_temp
                        dgrid_lf_offset_x=0.5-(amod_temp-0.5)/local_fraction_grid_size_scaling_temp
                        amod_temp=amod(real(dim_start_nc(y_dim_nc_index)-1+j_nc),local_fraction_grid_size_scaling_temp)
                        if (amod_temp.eq.0) amod_temp=local_fraction_grid_size_scaling_temp
                        dgrid_lf_offset_y=0.5-(amod_temp-0.5)/local_fraction_grid_size_scaling_temp
                        if (set_lf_offset_to_0) then
                            dgrid_lf_offset_x=0
                            dgrid_lf_offset_y=0
                        endif

                        !Loop through the +1, -1 grids
                        do jj=jjj_start,jjj_end
                            do ii=iii_start,iii_end

                                ii_nc=ii+i_nc
                                jj_nc=jj+j_nc

                                ii_w=ii+ii_w0
                                jj_w=jj+jj_w0

                                !Put in a limit
                                if (ii_nc.ge.1.and.ii_nc.le.dim_length_nc(x_dim_nc_index).and.jj_nc.ge.1.and.jj_nc.le.dim_length_nc(y_dim_nc_index)) then

                                    !Set the edges to an EMEP grid surounding the EMEP grid being assessed
                                    xpos_min2=max(xpos_area_min2,var1d_nc(ii_nc,lon_nc_index)-dgrid_nc(lon_nc_index)/2.)
                                    xpos_max2=min(xpos_area_max2,var1d_nc(ii_nc,lon_nc_index)+dgrid_nc(lon_nc_index)/2.)
                                    ypos_min2=max(ypos_area_min2,var1d_nc(jj_nc,lat_nc_index)-dgrid_nc(lat_nc_index)/2.)
                                    ypos_max2=min(ypos_area_max2,var1d_nc(jj_nc,lat_nc_index)+dgrid_nc(lat_nc_index)/2.)


                                    !Calculate area weighting
                                    if (xpos_max2.gt.xpos_min2.and.ypos_max2.gt.ypos_min2) then
                                        weighting_val=(ypos_max2-ypos_min2)*(xpos_max2-xpos_min2)/dgrid_nc(lon_nc_index)/dgrid_nc(lat_nc_index)
                                    else
                                        weighting_val=0.
                                    endif

                                    !write(*,*) i,j,ii_nc,jj_nc,weighting_val,temp_EMEP(ii_nc,jj_nc,1,emep_local_subgrid_index,1,1)

                                    !Area weighting (interpolated) EMEP concentrations, independent of where it comes from
                                    subgrid(i,j,:,emep_subgrid_index,1:n_source_index,:)=subgrid(i,j,:,emep_subgrid_index,1:n_source_index,:)+temp_EMEP(ii_nc,jj_nc,:,emep_subgrid_index,1:n_source_index,:)*weighting_val
                                    subgrid(i,j,:,emep_local_subgrid_index,1:n_source_index,:)=subgrid(i,j,:,emep_local_subgrid_index,1:n_source_index,:)+temp_EMEP(ii_nc,jj_nc,:,emep_local_subgrid_index,1:n_source_index,:)*weighting_val
                                    !var3d_nc(ii_nc,jj_nc,tt,conc_nc_index,1:n_source_index,:)*weighting_val

                                    !do i_pollutant=1,n_emep_pollutant_loop
                                    !do i_loop=1,n_pollutant_compound_loop(i_pollutant)
                                    !    comp_EMEP_subgrid(i,j,:,pollutant_compound_loop_index(i_pollutant,i_loop))=comp_EMEP_subgrid(i,j,:,pollutant_compound_loop_index(i_pollutant,i_loop))+temp_comp_EMEP(ii_nc,jj_nc,:,pollutant_compound_loop_index(i_pollutant,i_loop))*weighting_val
                                    !               ! +comp_var3d_nc(ii_nc,jj_nc,tt,pollutant_compound_loop_index(i_pollutant,i_loop))*weighting_val
                                    !enddo
                                    !enddo
                                    comp_EMEP_subgrid(i,j,:,:)=comp_EMEP_subgrid(i,j,:,:)+temp_comp_EMEP(ii_nc,jj_nc,:,:)*weighting_val

                                    if (save_emep_species.or.save_seasalt) then
                                        species_EMEP_subgrid(i,j,:,:,:)=species_EMEP_subgrid(i,j,:,:,:)+temp_species_EMEP(ii_nc,jj_nc,:,:,:)*weighting_val
                                        !do i_sp=1,n_species_loop_index
                                        !ii_sp=species_loop_index(i_sp)
                                        !do i_loop=1,n_pmxx_sp_index
                                        !    species_EMEP_subgrid(i,j,:,i_loop,i_sp)=species_EMEP_subgrid(i,j,:,i_loop,i_sp)+temp_species_EMEP(ii_nc,jj_nc,:,i_loop,i_sp)*weighting_val
                                        !       ! +species_var3d_nc(ii_nc,jj_nc,tt,i_loop,i_sp)*weighting_val
                                        !enddo
                                        !enddo
                                    endif

                                endif
                            enddo
                        enddo


                    endif
                enddo
            enddo


            !Set the non-local for each source individually
            subgrid(:,:,:,emep_nonlocal_subgrid_index,:,:)=subgrid(:,:,:,emep_subgrid_index,:,:)-subgrid(:,:,:,emep_local_subgrid_index,:,:)

            if (calculate_deposition_flag) then
                subgrid(i,j,:,drydepo_nonlocal_subgrid_index,:,:)=subgrid(i,j,:,drydepo_nonlocal_subgrid_index,:,:) &
                    *subgrid(i,j,:,emep_nonlocal_subgrid_index,:,:)/subgrid(i,j,:,emep_subgrid_index,:,:)
                subgrid(i,j,:,wetdepo_nonlocal_subgrid_index,:,:)=subgrid(i,j,:,wetdepo_nonlocal_subgrid_index,:,:) &
                    *subgrid(i,j,:,emep_nonlocal_subgrid_index,:,:)/subgrid(i,j,:,emep_subgrid_index,:,:)
            endif

            if (allocated(temp_EMEP)) deallocate (temp_EMEP)
            if (allocated(temp_comp_EMEP)) deallocate (temp_comp_EMEP)
            if (allocated(temp_species_EMEP)) deallocate (temp_species_EMEP)
            if (allocated(temp_EMEP_from_in_region)) deallocate (temp_EMEP_from_in_region)

        endif

        !Areal interpolation of the subgrid nearest neighbour calculations EMEP_grid_interpolation_flag=0
        !Do not use, very slow
        if (EMEP_grid_interpolation_flag.eq.5) then

            write(unit_logfile,'(A)')'Interpolating uEMEP local subgrid nearest neighbour contributions using area weighted interpolation'


            if (.not.allocated(temp_subgrid)) allocate (temp_subgrid(subgrid_dim(t_dim_index),n_subgrid_index,n_source_index,n_pollutant_loop))
            if (.not.allocated(temp_comp_EMEP_subgrid)) allocate (temp_comp_EMEP_subgrid(subgrid_dim(t_dim_index),n_compound_index))
            if (.not.allocated(temp_species_EMEP_subgrid)) allocate (temp_species_EMEP_subgrid(subgrid_dim(t_dim_index),n_pmxx_sp_index,n_species_loop_index))

            !Set the loop sizes for the local area interpolation
            jjj_start=-1
            iii_start=-1
            jjj_end=1
            iii_end=1


            do j=1,subgrid_dim(y_dim_index)
                do i=1,subgrid_dim(x_dim_index)
                    !Do this calculation everywhere
                    !if (use_subgrid(i,j,allsource_index)) then

                    if (EMEP_projection_type.eq.LL_projection_index) then
                        distance_grid_x=111000.*dgrid_nc(lon_nc_index)*cos(lat_subgrid(i,j)*pi/180.)
                        distance_grid_y=111000.*dgrid_nc(lat_nc_index)
                    else
                        !Assumed LCC or PS
                        distance_grid_x=dgrid_nc(lon_nc_index)
                        distance_grid_y=dgrid_nc(lat_nc_index)
                    endif
                    xpos_limit=distance_grid_x/2.
                    ypos_limit=distance_grid_y/2.

                    jj_start=-int((ypos_limit)/subgrid_delta(y_dim_index))
                    ii_start=-int((xpos_limit)/subgrid_delta(x_dim_index))
                    jj_end=+int((ypos_limit)/subgrid_delta(y_dim_index))
                    ii_end=+int((xpos_limit)/subgrid_delta(x_dim_index))

                    !Initialise arrays
                    temp_subgrid(:,emep_subgrid_index,:,:)=0
                    temp_subgrid(:,emep_frac_subgrid_index,:,:)=0
                    temp_subgrid(:,emep_local_subgrid_index,:,:)=0
                    temp_subgrid(:,emep_nonlocal_subgrid_index,:,:)=0
                    temp_comp_EMEP_subgrid(:,:)=0
                    temp_species_EMEP_subgrid(:,:,:)=0

                    temp_count=0
                    !write(*,*) i,j,ii_end,jj_end

                    jjj_start=max(1,j+jj_start)
                    iii_start=max(1,i+ii_start)
                    jjj_end=min(subgrid_dim(y_dim_index),j+jj_end)
                    iii_end=min(subgrid_dim(x_dim_index),i+ii_end)

                    temp_subgrid(:,:,:,:)=sum(sum(subgrid(iii_start:iii_end,jjj_start:jjj_end,:,:,:,:),1),1)
                    temp_comp_EMEP_subgrid(:,:)=sum(sum(comp_EMEP_subgrid(iii_start:iii_end,jjj_start:jjj_end,:,:),1),1)
                    temp_species_EMEP_subgrid(:,:,:)=sum(sum(species_EMEP_subgrid(iii_start:iii_end,jjj_start:jjj_end,:,:,:),1),1)
                    temp_count=(iii_end-iii_start+1)*(jjj_end-jjj_start+1)

                    if (temp_count.gt.0) then
                        subgrid(i,j,:,emep_subgrid_index,:,:)=temp_subgrid(:,emep_subgrid_index,:,:)/temp_count
                        subgrid(i,j,:,emep_local_subgrid_index,:,:)=temp_subgrid(:,emep_local_subgrid_index,:,:)/temp_count
                        subgrid(i,j,:,emep_nonlocal_subgrid_index,:,:)=temp_subgrid(:,emep_nonlocal_subgrid_index,:,:)/temp_count
                        comp_EMEP_subgrid(i,j,:,:)=temp_comp_EMEP_subgrid(:,:)/temp_count
                        species_EMEP_subgrid(i,j,:,:,:)=temp_species_EMEP_subgrid(:,:,:)/temp_count
                    else
                        subgrid(i,j,:,emep_subgrid_index,:,:)=0
                        subgrid(i,j,:,emep_local_subgrid_index,:,:)=0
                        subgrid(i,j,:,emep_nonlocal_subgrid_index,:,:)=0
                        species_EMEP_subgrid(i,j,:,:,:)=0
                    endif


                    !endif
                enddo
            enddo


            if (allocated(temp_subgrid)) deallocate (temp_subgrid)
            if (allocated(temp_comp_EMEP_subgrid)) deallocate (temp_comp_EMEP_subgrid)
            if (allocated(temp_species_EMEP_subgrid)) deallocate (temp_species_EMEP_subgrid)

        endif

        !Set the start and end times of the loop
        t_start=1
        t_end=subgrid_dim(t_dim_index)

        !Loop through the time and the subgrids
        do tt=t_start,t_end

            !Quick calculation of area weighting, no edge effects. Does not need to change with time
            !This is done also if there is moving window weighting later as it is used for the nonlocal contribution
            !This is the old version no longer in use
            if (EMEP_grid_interpolation_flag.eq.-1.or.(EMEP_grid_interpolation_flag.gt.1.and.EMEP_grid_interpolation_flag.lt.5)) then

                if (tt.eq.t_start) write(unit_logfile,'(A)')'Calculating EMEP local subgrid contribution using area weighted interpolation (obsolete version)'

                !Set weighting indexes
                n_weight=3+2*floor((EMEP_grid_interpolation_size-1.)*0.5)
                ii_w0=1+floor(n_weight*.5)
                jj_w0=1+floor(n_weight*.5)

                if (tt.eq.t_start) write(unit_logfile,*) 'Weighting grid dimensions and centres: ',n_weight,ii_w0,jj_w0

                if (.not.allocated(weighting_nc)) allocate (weighting_nc(n_weight,n_weight,tt_dim,n_source_index)) !EMEP grid weighting for interpolation. Does not need a source index for area weighting
                if (.not.allocated(area_weighting_nc)) allocate (area_weighting_nc(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),n_weight,n_weight,tt_dim,n_source_index)) !EMEP grid weighting for area interpolation

                !Initialise arrays
                subgrid(:,:,tt,emep_subgrid_index,:,:)=0
                subgrid(:,:,tt,emep_frac_subgrid_index,:,:)=0
                subgrid(:,:,tt,emep_local_subgrid_index,:,:)=0
                subgrid(:,:,tt,emep_nonlocal_subgrid_index,:,:)=0
                comp_EMEP_subgrid(:,:,tt,:)=0
                species_EMEP_subgrid(:,:,tt,:,:)=0

                !Cover the search area necessary for the surounding EMEP grids
                jj_start=-1-floor(0.5*(EMEP_grid_interpolation_size-1.))
                ii_start=-1-floor(0.5*(EMEP_grid_interpolation_size-1.))
                jj_end=1+floor(0.5*(EMEP_grid_interpolation_size-1.))
                ii_end=1+floor(0.5*(EMEP_grid_interpolation_size-1.))
                !jj_start=-1-ceiling(0.5*(EMEP_grid_interpolation_size-1.))
                !ii_start=-1-ceiling(0.5*(EMEP_grid_interpolation_size-1.))
                !jj_end=1+ceiling(0.5*(EMEP_grid_interpolation_size-1.))
                !ii_end=1+ceiling(0.5*(EMEP_grid_interpolation_size-1.))

                !Set the size of the region surounding the target grid that is searched
                xpos_limit=dgrid_nc(lon_nc_index)/2.*EMEP_grid_interpolation_size
                ypos_limit=dgrid_nc(lat_nc_index)/2.*EMEP_grid_interpolation_size

                do j=1,subgrid_dim(y_dim_index)
                    do i=1,subgrid_dim(x_dim_index)
                        if (use_subgrid(i,j,allsource_index)) then

                            !Assumes it is never on the edge of the EMEP grid as it is not limitted
                            i_nc=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
                            j_nc=crossreference_target_to_emep_subgrid(i,j,y_dim_index)

                            !These are the subgrid positions projected to the EMEP projection
                            xpos_subgrid=xproj_subgrid(i,j)
                            ypos_subgrid=yproj_subgrid(i,j)

                            !Set the edges of the search area surounding the target grid
                            xpos_area_min=xpos_subgrid-xpos_limit
                            xpos_area_max=xpos_subgrid+xpos_limit
                            ypos_area_min=ypos_subgrid-ypos_limit
                            ypos_area_max=ypos_subgrid+ypos_limit

                            weighting_nc(:,:,tt_dim,:)=0.

                            do jj=jj_start,jj_end
                                do ii=ii_start,ii_end

                                    ii_nc=ii+i_nc
                                    jj_nc=jj+j_nc

                                    ii_w=ii+ii_w0
                                    jj_w=jj+jj_w0

                                    !Put in a limit
                                    if (ii_nc.ge.1.and.ii_nc.le.dim_length_nc(x_dim_nc_index).and.jj_nc.ge.1.and.jj_nc.le.dim_length_nc(y_dim_nc_index)) then

                                        !Set the edges to an EMEP grid surounding the EMEP grid being assessed
                                        xpos_min=max(xpos_area_min,var1d_nc(ii_nc,lon_nc_index)-dgrid_nc(lon_nc_index)/2.)
                                        xpos_max=min(xpos_area_max,var1d_nc(ii_nc,lon_nc_index)+dgrid_nc(lon_nc_index)/2.)
                                        ypos_min=max(ypos_area_min,var1d_nc(jj_nc,lat_nc_index)-dgrid_nc(lat_nc_index)/2.)
                                        ypos_max=min(ypos_area_max,var1d_nc(jj_nc,lat_nc_index)+dgrid_nc(lat_nc_index)/2.)

                                        !Calculate area weighting
                                        if (xpos_max.gt.xpos_min.and.ypos_max.gt.ypos_min) then
                                            weighting_nc(ii_w,jj_w,tt_dim,:)=(ypos_max-ypos_min)*(xpos_max-xpos_min)/dgrid_nc(lon_nc_index)/dgrid_nc(lat_nc_index)/EMEP_grid_interpolation_size_sqr
                                        else
                                            weighting_nc(ii_w,jj_w,tt_dim,:)=0.
                                        endif

                                        do i_pollutant=1,n_emep_pollutant_loop

                                            subgrid(i,j,tt,emep_local_subgrid_index,:,i_pollutant)=subgrid(i,j,tt,emep_local_subgrid_index,:,i_pollutant) &
                                                +var3d_nc(ii_nc,jj_nc,tt,local_nc_index,1:n_source_index,i_pollutant)*weighting_nc(ii_w,jj_w,tt_dim,:)
                                            subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,i_pollutant)=subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,i_pollutant) &
                                                +(var3d_nc(ii_nc,jj_nc,tt,conc_nc_index,1:n_source_index,i_pollutant)-var3d_nc(ii_nc,jj_nc,tt,local_nc_index,1:n_source_index,i_pollutant))*weighting_nc(ii_w,jj_w,tt_dim,:)

                                        enddo


                                        do i_pollutant=1,n_emep_pollutant_loop
                                            !write(*,*) var3d_nc(ii_nc,jj_nc,tt,local_nc_index,:,i_pollutant)
                                            !Interpolate the other EMEP compounds as well to subgrid
                                            do i_loop=1,n_pollutant_compound_loop(i_pollutant)
                                                comp_EMEP_subgrid(i,j,tt,pollutant_compound_loop_index(i_pollutant,i_loop))=comp_EMEP_subgrid(i,j,tt,pollutant_compound_loop_index(i_pollutant,i_loop)) &
                                                    +comp_var3d_nc(ii_nc,jj_nc,tt,pollutant_compound_loop_index(i_pollutant,i_loop))*weighting_nc(ii_w,jj_w,tt_dim,allsource_index)
                                            enddo

                                        enddo

                                        if (save_emep_species.or.save_seasalt) then
                                            do i_sp=1,n_species_loop_index
                                                ii_sp=species_loop_index(i_sp)
                                                do i_loop=1,n_pmxx_sp_index
                                                    species_EMEP_subgrid(i,j,tt,i_loop,i_sp)=species_EMEP_subgrid(i,j,tt,i_loop,i_sp) &
                                                        +species_var3d_nc(ii_nc,jj_nc,tt,i_loop,i_sp)*weighting_nc(ii_w,jj_w,tt_dim,allsource_index)
                                                enddo
                                            enddo
                                        endif

                                    endif
                                enddo
                            enddo

                            !Calculate the nonlocal correction, weighting for grids beyond the central grid
                            nonlocal_correction(tt_dim,:,:)=0.
                            do jj=jj_start,jj_end
                                do ii=ii_start,ii_end

                                    ii_nc=ii+i_nc
                                    jj_nc=jj+j_nc

                                    ii_w=ii+ii_w0
                                    jj_w=jj+jj_w0

                                    !Put in a limit
                                    if (ii_nc.ge.1.and.ii_nc.le.dim_length_nc(x_dim_nc_index).and.jj_nc.ge.1.and.jj_nc.le.dim_length_nc(y_dim_nc_index)) then

                                        if (jj.ne.0.or.ii.ne.0) then

                                            !First weight is emission, the second is area
                                            do i_pollutant=1,n_emep_pollutant_loop
                                                nonlocal_correction(tt_dim,:,i_pollutant)=nonlocal_correction(tt_dim,:,i_pollutant) &
                                                    -lc_var3d_nc(ii_w0,jj_w0,ii_nc,jj_nc,tt,lc_local_nc_index,1:n_source_index,i_pollutant)*weighting_nc(ii_w,jj_w,tt_dim,:)*weighting_nc(ii_w0,jj_w0,tt_dim,:) &
                                                    -lc_var3d_nc(ii_w,jj_w,i_nc,j_nc,tt,lc_local_nc_index,1:n_source_index,i_pollutant)*weighting_nc(ii_w0,jj_w0,tt_dim,:)*weighting_nc(ii_w,jj_w,tt_dim,:)
                                            enddo
                                        endif

                                    endif

                                enddo
                            enddo


                            !Place the EMEP values in the target subgrid
                            subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,:)=subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,:)+nonlocal_correction(tt_dim,:,:)
                            subgrid(i,j,tt,emep_local_subgrid_index,:,:)=subgrid(i,j,tt,emep_local_subgrid_index,:,:)-nonlocal_correction(tt_dim,:,:)
                            subgrid(i,j,tt,emep_subgrid_index,:,:)=subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,:)+subgrid(i,j,tt,emep_local_subgrid_index,:,:)

                            !Take the already calculated nonlocal depositions to be the fraction of the nonlocal/total EMEP values
                            if (calculate_deposition_flag) then
                                subgrid(i,j,tt,drydepo_nonlocal_subgrid_index,:,:)=subgrid(i,j,tt,drydepo_nonlocal_subgrid_index,:,:) &
                                    *subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,:)/subgrid(i,j,tt,emep_subgrid_index,:,:)
                                subgrid(i,j,tt,wetdepo_nonlocal_subgrid_index,:,:)=subgrid(i,j,tt,wetdepo_nonlocal_subgrid_index,:,:) &
                                    *subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,:)/subgrid(i,j,tt,emep_subgrid_index,:,:)
                            endif

                            !Put the area weighting in the larger array for use later in the emission proxy weighting (if needed)
                            area_weighting_nc(i,j,:,:,tt_dim,:)=weighting_nc(:,:,tt_dim,:)

                            !For diagnostics only
                            nonlocal_correction_average=nonlocal_correction_average+nonlocal_correction(tt_dim,:,:)

                            !write(*,*) subgrid(i,j,tt,emep_nonlocal_subgrid_index,allsource_index,:)

                        endif
                    enddo
                    !write(*,*) 'Subgrid EMEP area interpolation: ',j,' of ',subgrid_dim(2)
                enddo


                if (tt.eq.t_end) then
                    nonlocal_correction_average=nonlocal_correction_average/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
                    do i_pollutant=1,n_emep_pollutant_loop
                        write(fmt,'(A,I0,A)'), '(', n_source_index, 'es12.4)'
                        write(unit_logfile,fmt) 'Nonlocal correction for area weighting ('//trim(pollutant_file_str(pollutant_loop_index(i_pollutant)))//') = ',nonlocal_correction_average(:,i_pollutant)
                    enddo
                endif


            endif



            !This does not work very well for the additional and reigon contributions. Use 6 instead
            if (EMEP_grid_interpolation_flag.eq.1) then

                if (tt.eq.t_start) write(unit_logfile,'(A)')'Calculating EMEP local subgrid contribution using area weighted interpolation v2'

                !Set weighting indexes
                n_weight=3+2*floor((EMEP_grid_interpolation_size-1.)*0.5)
                ii_w0=1+floor(n_weight*.5)
                jj_w0=1+floor(n_weight*.5)

                ii_nc_w0=xdist_centre_nc
                jj_nc_w0=ydist_centre_nc

                n_weight_nc_x=xdist_centre_nc*2-1
                n_weight_nc_y=ydist_centre_nc*2-1
                n_weight_nc_x=dim_length_nc(xdist_dim_nc_index)
                n_weight_nc_y=dim_length_nc(ydist_dim_nc_index)

                if (tt.eq.t_start) write(unit_logfile,'(a,3i)') 'Weighting grid dimensions and centres: ',n_weight,ii_w0,jj_w0
                if (tt.eq.t_start) write(unit_logfile,'(a,4i)') 'EMEP local fraction grid dimensions and centres: ',n_weight_nc_x,n_weight_nc_y,ii_nc_w0,jj_nc_w0

                !if (.not.allocated(EMEP_local_contribution)) allocate (EMEP_local_contribution(n_weight_nc_x,n_weight_nc_y,n_source_index,n_emep_pollutant_loop))
                if (.not.allocated(EMEP_local_contribution)) allocate (EMEP_local_contribution(n_weight_nc_x,n_weight_nc_y,n_source_index,n_pollutant_loop))

                !Initialise arrays
                subgrid(:,:,tt,emep_subgrid_index,:,:)=0
                subgrid(:,:,tt,emep_frac_subgrid_index,:,:)=0
                subgrid(:,:,tt,emep_local_subgrid_index,:,:)=0
                subgrid(:,:,tt,emep_nonlocal_subgrid_index,:,:)=0
                comp_EMEP_subgrid(:,:,tt,:)=0
                species_EMEP_subgrid(:,:,tt,:,:)=0
                EMEP_local_contribution=0

                !Cover the search area necessary for the surounding EMEP grids
                jj_start=-1-floor(0.5*(EMEP_grid_interpolation_size-1.))
                ii_start=-1-floor(0.5*(EMEP_grid_interpolation_size-1.))
                jj_end=1+floor(0.5*(EMEP_grid_interpolation_size-1.))
                ii_end=1+floor(0.5*(EMEP_grid_interpolation_size-1.))

                !Set the loop sizes for the local area interpolation
                jjj_start=-1
                iii_start=-1
                jjj_end=1
                iii_end=1

                !Set the size of the region surounding the target grid that is searched in the LF grid
                xpos_limit=dgrid_nc(lon_nc_index)/2.*EMEP_grid_interpolation_size*local_fraction_grid_size_scaling_temp
                ypos_limit=dgrid_nc(lat_nc_index)/2.*EMEP_grid_interpolation_size*local_fraction_grid_size_scaling_temp
                xpos_limit2=dgrid_nc(lon_nc_index)/2.
                ypos_limit2=dgrid_nc(lat_nc_index)/2.


                !Recheck this!
                if (tt.eq.t_start) write(unit_logfile,'(a,4i)') 'Loop sizes: ',jj_start,jj_end,ii_start,ii_end

                do j=1,subgrid_dim(y_dim_index)
                    do i=1,subgrid_dim(x_dim_index)
                        if (use_subgrid(i,j,allsource_index)) then

                            !Assumes it is never on the edge of the EMEP grid as it is not limitted
                            i_nc=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
                            j_nc=crossreference_target_to_emep_subgrid(i,j,y_dim_index)

                            xpos_subgrid=xproj_subgrid(i,j)
                            ypos_subgrid=yproj_subgrid(i,j)

                            !Set the edges of the search area surounding the target grid
                            xpos_area_min=-xpos_limit
                            xpos_area_max=+xpos_limit
                            ypos_area_min=-ypos_limit
                            ypos_area_max=+ypos_limit

                            xpos_area_min2=xpos_subgrid-xpos_limit2
                            xpos_area_max2=xpos_subgrid+xpos_limit2
                            ypos_area_min2=ypos_subgrid-ypos_limit2
                            ypos_area_max2=ypos_subgrid+ypos_limit2

                            !Limit the region. This will still allow an EMEP contribution from half a grid away
                            !Same limit is NOT applied on the emissions in the moving window so inconsistent
                            !write(*,'(2i,4e12.2)') i,j,xpos_area_min,xpos_area_max,ypos_area_min,ypos_area_max
                            if (limit_emep_grid_interpolation_region_to_calculation_region) then
                                xpos_area_min=max(xpos_area_min,subgrid_proj_min(x_dim_index)-xpos_subgrid-xpos_limit2)
                                xpos_area_max=min(xpos_area_max,subgrid_proj_max(x_dim_index)-xpos_subgrid+xpos_limit2)
                                ypos_area_min=max(ypos_area_min,subgrid_proj_min(y_dim_index)-ypos_subgrid-ypos_limit2)
                                ypos_area_max=min(ypos_area_max,subgrid_proj_max(y_dim_index)-ypos_subgrid+ypos_limit2)
                            endif
                            !write(*,'(2i,4e12.2)') i,j,xpos_area_min,xpos_area_max,ypos_area_min,ypos_area_max

                            !Set the offset to the centre of the local fraction grid when using larger local fraction grids
                            amod_temp=amod(real(dim_start_nc(x_dim_nc_index)-1+i_nc),local_fraction_grid_size_scaling_temp)
                            if (amod_temp.eq.0) amod_temp=local_fraction_grid_size_scaling_temp
                            dgrid_lf_offset_x=0.5-(amod_temp-0.5)/local_fraction_grid_size_scaling_temp
                            amod_temp=amod(real(dim_start_nc(y_dim_nc_index)-1+j_nc),local_fraction_grid_size_scaling_temp)
                            if (amod_temp.eq.0) amod_temp=local_fraction_grid_size_scaling_temp
                            dgrid_lf_offset_y=0.5-(amod_temp-0.5)/local_fraction_grid_size_scaling_temp

                            !First create an interpolated grid around the x,y position for the species and the compounds
                            EMEP_local_contribution=0

                            do jj=jjj_start,jjj_end
                                do ii=iii_start,iii_end

                                    ii_nc=ii+i_nc
                                    jj_nc=jj+j_nc

                                    ii_w=ii+ii_w0
                                    jj_w=jj+jj_w0

                                    !Put in a limit
                                    if (ii_nc.ge.1.and.ii_nc.le.dim_length_nc(x_dim_nc_index).and.jj_nc.ge.1.and.jj_nc.le.dim_length_nc(y_dim_nc_index)) then

                                        !Set the edges to an EMEP grid surounding the EMEP grid being assessed
                                        xpos_min2=max(xpos_area_min2,var1d_nc(ii_nc,lon_nc_index)-dgrid_nc(lon_nc_index)/2.)
                                        xpos_max2=min(xpos_area_max2,var1d_nc(ii_nc,lon_nc_index)+dgrid_nc(lon_nc_index)/2.)
                                        ypos_min2=max(ypos_area_min2,var1d_nc(jj_nc,lat_nc_index)-dgrid_nc(lat_nc_index)/2.)
                                        ypos_max2=min(ypos_area_max2,var1d_nc(jj_nc,lat_nc_index)+dgrid_nc(lat_nc_index)/2.)


                                        !Calculate area weighting
                                        if (xpos_max2.gt.xpos_min2.and.ypos_max2.gt.ypos_min2) then
                                            weighting_val=(ypos_max2-ypos_min2)*(xpos_max2-xpos_min2)/dgrid_nc(lon_nc_index)/dgrid_nc(lat_nc_index)
                                        else
                                            weighting_val=0.
                                        endif

                                        !write(*,*) ii,jj,weighting_val,(ypos_max2-ypos_min2),(xpos_max2-xpos_min2)

                                        !Area weighting (interpolated) EMEP concentrations, independent of where it comes from
                                        subgrid(i,j,tt,emep_subgrid_index,:,:)=subgrid(i,j,tt,emep_subgrid_index,1:n_source_index,:)+var3d_nc(ii_nc,jj_nc,tt,conc_nc_index,1:n_source_index,:)*weighting_val

                                        do i_pollutant=1,n_emep_pollutant_loop
                                            do i_loop=1,n_pollutant_compound_loop(i_pollutant)
                                                comp_EMEP_subgrid(i,j,tt,pollutant_compound_loop_index(i_pollutant,i_loop))=comp_EMEP_subgrid(i,j,tt,pollutant_compound_loop_index(i_pollutant,i_loop)) &
                                                    +comp_var3d_nc(ii_nc,jj_nc,tt,pollutant_compound_loop_index(i_pollutant,i_loop))*weighting_val
                                            enddo
                                        enddo

                                        if (save_emep_species.or.save_seasalt) then
                                            do i_sp=1,n_species_loop_index
                                                ii_sp=species_loop_index(i_sp)
                                                do i_loop=1,n_pmxx_sp_index
                                                    species_EMEP_subgrid(i,j,tt,i_loop,i_sp)=species_EMEP_subgrid(i,j,tt,i_loop,i_sp) &
                                                        +species_var3d_nc(ii_nc,jj_nc,tt,i_loop,i_sp)*weighting_val
                                                enddo
                                            enddo
                                        endif

                                        if (first_interpolate_lf) then
                                            !Set the area surrounding the multi LF grid, centred on 0
                                            xpos_lf_area_min=(ii-1/2.+dgrid_lf_offset_x)*dgrid_nc(lon_nc_index)*local_fraction_grid_size_scaling_temp
                                            xpos_lf_area_max=(ii+1/2.+dgrid_lf_offset_x)*dgrid_nc(lon_nc_index)*local_fraction_grid_size_scaling_temp
                                            ypos_lf_area_min=(jj-1/2.+dgrid_lf_offset_y)*dgrid_nc(lat_nc_index)*local_fraction_grid_size_scaling_temp
                                            ypos_lf_area_max=(jj+1/2.+dgrid_lf_offset_y)*dgrid_nc(lat_nc_index)*local_fraction_grid_size_scaling_temp

                                            xpos_min3=max(xpos_area_min,xpos_lf_area_min)
                                            xpos_max3=min(xpos_area_max,xpos_lf_area_max)
                                            ypos_min3=max(ypos_area_min,ypos_lf_area_min)
                                            ypos_max3=min(ypos_area_max,ypos_lf_area_max)

                                            !Calculate area weighting LF grid
                                            if (xpos_max3.gt.xpos_min3.and.ypos_max3.gt.ypos_min3) then
                                                weighting_val3=(ypos_max3-ypos_min3)*(xpos_max3-xpos_min3)/dgrid_nc(lon_nc_index)/dgrid_nc(lat_nc_index)/local_fraction_grid_size_scaling_temp/local_fraction_grid_size_scaling_temp
                                            else
                                                weighting_val3=0.
                                            endif

                                        else
                                            weighting_val3=weighting_val
                                        endif

                                        !write(*,*) ii,jj,ii_nc,jj_nc,weighting_val3

                                        EMEP_local_contribution(:,:,:,:)=EMEP_local_contribution(:,:,:,:)+lc_var3d_nc(:,:,ii_nc,jj_nc,tt,lc_local_nc_index,1:n_source_index,:)*weighting_val3

                                    endif
                                enddo
                            enddo


                            !Set the offset to the centre of the local fraction grid when using larger local fraction grids
                            !Requires knowledge of the total EMEP grid index so uses dim_start_nc
                            amod_temp=amod(real(dim_start_nc(x_dim_nc_index)-1+i_nc),local_fraction_grid_size_scaling_temp)
                            if (amod_temp.eq.0) amod_temp=local_fraction_grid_size_scaling_temp
                            dgrid_lf_offset_x=0.5-(amod_temp-0.5)/local_fraction_grid_size_scaling_temp
                            amod_temp=amod(real(dim_start_nc(y_dim_nc_index)-1+j_nc),local_fraction_grid_size_scaling_temp)
                            if (amod_temp.eq.0) amod_temp=local_fraction_grid_size_scaling_temp
                            dgrid_lf_offset_y=0.5-(amod_temp-0.5)/local_fraction_grid_size_scaling_temp
                            !write(*,*) dim_start_nc(x_dim_nc_index)-1+i_nc,dim_start_nc(y_dim_nc_index)-1+j_nc,amod_temp,dgrid_lf_offset_x,dgrid_lf_offset_y

                            if (first_interpolate_lf) then
                                dgrid_lf_offset_x=0
                                dgrid_lf_offset_y=0
                            endif
                            if (set_lf_offset_to_0) then
                                dgrid_lf_offset_x=0
                                dgrid_lf_offset_y=0
                            endif

                            !Still need to change the dispersion routines for distance calculated and the size of the domain read in by EMEP and the other emissions

                            !Calculate the local contribution within the moving window area
                            do jjj=jj_start,jj_end
                                do iii=ii_start,ii_end

                                    ii_nc=ii+i_nc
                                    jj_nc=jj+j_nc

                                    ii_w=ii+ii_w0
                                    jj_w=jj+jj_w0

                                    iii_nc_w=iii+ii_nc_w0
                                    jjj_nc_w=jjj+jj_nc_w0

                                    !Put in a limit
                                    if (iii_nc_w.ge.1.and.iii_nc_w.le.n_weight_nc_x.and.jjj_nc_w.ge.1.and.jjj_nc_w.le.n_weight_nc_y) then

                                        !Set the edges to an EMEP grid surounding the EMEP grid being assessed
                                        xpos_min=max(xpos_area_min,(iii+dgrid_lf_offset_x-1/2.)*dgrid_nc(lon_nc_index)*local_fraction_grid_size_scaling_temp)
                                        xpos_max=min(xpos_area_max,(iii+dgrid_lf_offset_x+1/2.)*dgrid_nc(lon_nc_index)*local_fraction_grid_size_scaling_temp)
                                        ypos_min=max(ypos_area_min,(jjj+dgrid_lf_offset_y-1/2.)*dgrid_nc(lat_nc_index)*local_fraction_grid_size_scaling_temp)
                                        ypos_max=min(ypos_area_max,(jjj+dgrid_lf_offset_y+1/2.)*dgrid_nc(lat_nc_index)*local_fraction_grid_size_scaling_temp)

                                        !Calculate area weighting
                                        if (xpos_max.gt.xpos_min.and.ypos_max.gt.ypos_min) then
                                            weighting_val=(ypos_max-ypos_min)*(xpos_max-xpos_min)/dgrid_nc(lon_nc_index)/dgrid_nc(lat_nc_index)/local_fraction_grid_size_scaling_temp/local_fraction_grid_size_scaling_temp
                                        else
                                            weighting_val=0.
                                        endif


                                        !write(*,*) iii_nc_w,jjj_nc_w,weighting_val
                                        subgrid(i,j,tt,emep_local_subgrid_index,:,:)=subgrid(i,j,tt,emep_local_subgrid_index,:,:)+EMEP_local_contribution(iii_nc_w,jjj_nc_w,1:n_source_index,:)*weighting_val

                                        !write(*,*) iii,jjj,weighting_val,EMEP_local_contribution(iii_nc_w,jjj_nc_w,traffic_nc_index,allsource_index)

                                    endif
                                enddo
                            enddo

                            !Place the EMEP values in the target subgrid
                            !subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,:)=subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,:)+nonlocal_correction(tt_dim,:,:)
                            subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,:)=subgrid(i,j,tt,emep_subgrid_index,:,:)-subgrid(i,j,tt,emep_local_subgrid_index,:,:)
                            !subgrid(i,j,tt,emep_subgrid_index,:,:)=subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,:)+subgrid(i,j,tt,emep_local_subgrid_index,:,:)
                            !write(*,*) i,j,sum(subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,:))
                            !Take the already calculated nonlocal depositions to be the fraction of the nonlocal/total EMEP values
                            if (calculate_deposition_flag) then
                                subgrid(i,j,tt,drydepo_nonlocal_subgrid_index,:,:)=subgrid(i,j,tt,drydepo_nonlocal_subgrid_index,:,:) &
                                    *subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,:)/subgrid(i,j,tt,emep_subgrid_index,:,:)
                                subgrid(i,j,tt,wetdepo_nonlocal_subgrid_index,:,:)=subgrid(i,j,tt,wetdepo_nonlocal_subgrid_index,:,:) &
                                    *subgrid(i,j,tt,emep_nonlocal_subgrid_index,:,:)/subgrid(i,j,tt,emep_subgrid_index,:,:)
                            endif


                        endif
                    enddo
                enddo


            endif

            !Loop through subgrid and carry out a subgrid weighted moving window interpolation using emissions. Not used anymore
            if (EMEP_grid_interpolation_flag.gt.1.and.EMEP_grid_interpolation_flag.lt.5) then

                nonlocal_correction_average=0.

                !n_weight is already set in the previous call

                if (tt.eq.t_start) write(unit_logfile,'(A,2i)')'Calculating EMEP local subgrid contribution using moving window interpolation method ',EMEP_grid_interpolation_flag

                if (.not.allocated(total_weighting_nc)) allocate (total_weighting_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),tt_dim,n_source_index,n_pollutant_loop)) !EMEP grid weighting for interpolation
                if (.not.allocated(proxy_weighting_nc)) allocate (proxy_weighting_nc(n_weight,n_weight,tt_dim,n_source_index,n_pollutant_loop)) !EMEP grid weighting for interpolation

                !Set the index offset for the local contribution
                i_w_c=1+floor(n_weight*.5)
                j_w_c=1+floor(n_weight*.5)

                subgrid(:,:,tt,emep_local_subgrid_index,:,:)=0

                !Emission moving window only works properly when the emission and subgrids are the same
                !Emission moving window adds all the subsource emissions since EMEP does not understand subsources
                !This means that combine_emission_subsources_during_dispersion must be set to true whenever using the redistribution of EMEP

                !Emission weighting
                if (EMEP_grid_interpolation_flag.eq.2) then

                    if (.not.allocated(weighting_subgrid)) allocate (weighting_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),tt_dim,n_source_index,n_pollutant_loop))
                    if (.not.allocated(crossreference_weighting_to_emep_subgrid)) allocate (crossreference_weighting_to_emep_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),2,n_source_index))
                    weighting_subgrid(:,:,tt_dim,:,:)=emission_subgrid(:,:,tt,:,:)
                    weighting_subgrid_dim(:,:)=emission_subgrid_dim(1:2,:)
                    crossreference_weighting_to_emep_subgrid=crossreference_emission_to_emep_subgrid

                endif

                if (EMEP_grid_interpolation_flag.eq.3) then

                    !Aggregated emissions first on integral grid to increase speed
                    if (.not.allocated(weighting_subgrid)) allocate (weighting_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),tt_dim,n_source_index,n_pollutant_loop))
                    if (.not.allocated(crossreference_weighting_to_emep_subgrid)) allocate (crossreference_weighting_to_emep_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),2,n_source_index))
                    do i_source=1,n_source_index
                        if (calculate_source(i_source)) then
                            crossreference_weighting_to_emep_subgrid(:,:,:,i_source)=crossreference_integral_to_emep_subgrid
                            weighting_subgrid_dim(:,i_source)=integral_subgrid_dim(1:2)
                        endif
                    enddo
                    weighting_subgrid(:,:,tt_dim,:,:)=0.
                    do i_source=1,n_source_index
                        if (calculate_source(i_source)) then
                            do j=1,emission_subgrid_dim(y_dim_index,i_source)
                                do i=1,emission_subgrid_dim(x_dim_index,i_source)
                                    i_cross=crossreference_emission_to_integral_subgrid(i,j,x_dim_index,i_source)
                                    j_cross=crossreference_emission_to_integral_subgrid(i,j,y_dim_index,i_source)
                                    weighting_subgrid(i_cross,j_cross,tt_dim,i_source,:)=weighting_subgrid(i_cross,j_cross,tt_dim,i_source,:)+emission_subgrid(i,j,tt,i_source,:)
                                enddo
                            enddo
                        endif
                    enddo

                endif

                !Integral proxy weighting
                if (EMEP_grid_interpolation_flag.eq.4) then
                    !Aggregated proxy on integral grid to increase speed
                    if (.not.allocated(weighting_subgrid)) allocate (weighting_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),integral_subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop))
                    if (.not.allocated(crossreference_weighting_to_emep_subgrid)) allocate (crossreference_weighting_to_emep_subgrid(integral_subgrid_dim(x_dim_index),integral_subgrid_dim(y_dim_index),2,n_source_index))
                    do i_source=1,n_source_index
                        if (calculate_source(i_source)) then
                            crossreference_weighting_to_emep_subgrid(:,:,:,i_source)=crossreference_integral_to_emep_subgrid
                            weighting_subgrid_dim(:,i_source)=integral_subgrid_dim(1:2)
                        endif
                    enddo
                    !Set the weighting subgrid to the sum of all subsource integral emissions
                    weighting_subgrid(:,:,tt_dim,:,:)=integral_subgrid(:,:,tt,hsurf_integral_subgrid_index,:,:)

                endif

                !Calculate weighting sum for each EMEP grid.
                total_weighting_nc=0.
                do i_source=1,n_source_index
                    if (calculate_source(i_source)) then
                        do j=1,weighting_subgrid_dim(y_dim_index,i_source)
                            do i=1,weighting_subgrid_dim(x_dim_index,i_source)
                                i_nc=crossreference_weighting_to_emep_subgrid(i,j,x_dim_index,i_source)
                                j_nc=crossreference_weighting_to_emep_subgrid(i,j,y_dim_index,i_source)
                                total_weighting_nc(i_nc,j_nc,tt_dim,i_source,:)=total_weighting_nc(i_nc,j_nc,tt_dim,i_source,:)+weighting_subgrid(i,j,tt_dim,i_source,:)
                                !write(*,*) i_source,i,j,i_nc,j_nc,weighting_subgrid(i,j,:,i_source)
                            enddo
                        enddo
                    endif
                enddo

                nonlocal_correction_average=0.

                xpos_limit=dgrid_nc(lon_nc_index)/2.*EMEP_grid_interpolation_size
                ypos_limit=dgrid_nc(lat_nc_index)/2.*EMEP_grid_interpolation_size


                do i_source=1,n_source_index
                    if (calculate_source(i_source)) then

                        !Calculate the proxy weighting in the nearest emep grids for each subgrid

                        do j=1,subgrid_dim(y_dim_index)
                            do i=1,subgrid_dim(x_dim_index)
                                if (use_subgrid(i,j,allsource_index)) then

                                    proxy_weighting_nc=0.

                                    xpos_subgrid=xproj_subgrid(i,j)
                                    ypos_subgrid=yproj_subgrid(i,j)

                                    !Set the edges of the search area surounding the target grid
                                    xpos_area_min=xpos_subgrid-xpos_limit
                                    xpos_area_max=xpos_subgrid+xpos_limit
                                    ypos_area_min=ypos_subgrid-ypos_limit
                                    ypos_area_max=ypos_subgrid+ypos_limit

                                    if (EMEP_grid_interpolation_flag.eq.2) then
                                        i_cross=crossreference_target_to_emission_subgrid(i,j,x_dim_index,i_source)
                                        j_cross=crossreference_target_to_emission_subgrid(i,j,y_dim_index,i_source)
                                        i_nc_c=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
                                        j_nc_c=crossreference_target_to_emep_subgrid(i,j,y_dim_index)
                                        !Limit the loop so that it doesn't go over more than necessary subgrids and does not go outside the domain
                                        i_start=max(1,i_cross-emission_subgrid_loop_index(x_dim_index,i_source))
                                        i_end=min(emission_subgrid_dim(x_dim_index,i_source),i_cross+emission_subgrid_loop_index(x_dim_index,i_source))
                                        j_start=max(1,j_cross-emission_subgrid_loop_index(y_dim_index,i_source))
                                        j_end=min(emission_subgrid_dim(y_dim_index,i_source),j_cross+emission_subgrid_loop_index(y_dim_index,i_source))

                                        do jj=j_start,j_end
                                            do ii=i_start,i_end

                                                xpos_emission_subgrid=xproj_emission_subgrid(ii,jj,i_source)
                                                ypos_emission_subgrid=yproj_emission_subgrid(ii,jj,i_source)

                                                if (abs(xpos_subgrid-xpos_emission_subgrid).le.xpos_limit &
                                                    .and.abs(ypos_subgrid-ypos_emission_subgrid).le.ypos_limit) then
                                                    i_nc=crossreference_emission_to_emep_subgrid(ii,jj,x_dim_index,i_source)
                                                    j_nc=crossreference_emission_to_emep_subgrid(ii,jj,y_dim_index,i_source)
                                                    proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt_dim,i_source,:)=proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt_dim,i_source,:)+weighting_subgrid(ii,jj,tt_dim,i_source,:)
                                                    !write(*,*) tt, proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt_dim,i_source), weighting_subgrid(ii,jj,tt_dim,i_source)

                                                endif
                                            enddo
                                        enddo

                                        !endif

                                    elseif (EMEP_grid_interpolation_flag.eq.3.or.EMEP_grid_interpolation_flag.eq.4)  then
                                        !Find the cross reference to the integral grid from the target grid

                                        i_cross=crossreference_target_to_integral_subgrid(i,j,x_dim_index)
                                        j_cross=crossreference_target_to_integral_subgrid(i,j,y_dim_index)
                                        i_nc_c=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
                                        j_nc_c=crossreference_target_to_emep_subgrid(i,j,y_dim_index)
                                        i_start=max(1,i_cross-integral_subgrid_loop_index(x_dim_index))
                                        i_end=min(integral_subgrid_dim(x_dim_index),i_cross+integral_subgrid_loop_index(x_dim_index))
                                        j_start=max(1,j_cross-integral_subgrid_loop_index(y_dim_index))
                                        j_end=min(integral_subgrid_dim(y_dim_index),j_cross+integral_subgrid_loop_index(y_dim_index))

                                        do jj=j_start,j_end
                                            do ii=i_start,i_end

                                                xpos_integral_subgrid=xproj_integral_subgrid(ii,jj)
                                                ypos_integral_subgrid=yproj_integral_subgrid(ii,jj)

                                                if (abs(xpos_subgrid-xpos_integral_subgrid).le.xpos_limit &
                                                    .and.abs(ypos_subgrid-ypos_integral_subgrid).le.ypos_limit) then
                                                    i_nc=crossreference_integral_to_emep_subgrid(ii,jj,x_dim_index)
                                                    j_nc=crossreference_integral_to_emep_subgrid(ii,jj,y_dim_index)
                                                    proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt_dim,i_source,:)=proxy_weighting_nc(i_nc-i_nc_c+i_w_c,j_nc-j_nc_c+j_w_c,tt_dim,i_source,:)+weighting_subgrid(ii,jj,tt_dim,i_source,:)
                                                endif
                                            enddo
                                        enddo

                                    endif

                                    i_nc=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
                                    j_nc=crossreference_target_to_emep_subgrid(i,j,y_dim_index)

                                    i_nc_start=max(1+i_nc-i_w_c,i_nc-1-floor((EMEP_grid_interpolation_size-1.)*0.5))
                                    i_nc_end=min(dim_length_nc(x_dim_nc_index)+i_nc-i_w_c,i_nc+1+floor((EMEP_grid_interpolation_size-1.)*0.5))
                                    j_nc_start=max(1+j_nc-j_w_c,j_nc-1-floor((EMEP_grid_interpolation_size-1.)*0.5))
                                    j_nc_end=min(dim_length_nc(y_dim_nc_index)+j_nc-j_w_c,j_nc+1+floor((EMEP_grid_interpolation_size-1.)*0.5))
                                    i_nc_start=max(1+i_nc-i_w_c,i_nc-1-ceiling((EMEP_grid_interpolation_size-1.)*0.5))
                                    i_nc_end=min(dim_length_nc(x_dim_nc_index)+i_nc-i_w_c,i_nc+1+ceiling((EMEP_grid_interpolation_size-1.)*0.5))
                                    j_nc_start=max(1+j_nc-j_w_c,j_nc-1-ceiling((EMEP_grid_interpolation_size-1.)*0.5))
                                    j_nc_end=min(dim_length_nc(y_dim_nc_index)+j_nc-j_w_c,j_nc+1+ceiling((EMEP_grid_interpolation_size-1.)*0.5))


                                    do jj=j_nc_start,j_nc_end
                                        do ii=i_nc_start,i_nc_end

                                            do i_pollutant=1,n_emep_pollutant_loop
                                                if (total_weighting_nc(ii,jj,tt_dim,i_source,i_pollutant).ne.0.) then
                                                    proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,tt_dim,i_source,i_pollutant)=proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,tt_dim,i_source,i_pollutant)/total_weighting_nc(ii,jj,tt_dim,i_source,i_pollutant)/EMEP_grid_interpolation_size_sqr
                                                    !write(*,*)  tt,proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,tt_dim,i_source)
                                                else
                                                    proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,tt_dim,i_source,i_pollutant)=0.
                                                endif
                                                if (proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,tt_dim,i_source,i_pollutant).gt.1) write(*,'(A,8i6,f12.2)')'WEIGHTING>1: ',i_pollutant,tt,i,j,ii,jj,ii-i_nc,jj-j_nc,proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,tt_dim,i_source,i_pollutant)
                                            enddo

                                        enddo
                                    enddo

                                    !Add up the contributing weights
                                    do jj=j_nc_start,j_nc_end
                                        do ii=i_nc_start,i_nc_end

                                            subgrid(i,j,tt,emep_local_subgrid_index,i_source,:)=subgrid(i,j,tt,emep_local_subgrid_index,i_source,:) &
                                                +var3d_nc(ii,jj,tt,local_nc_index,i_source,:)*proxy_weighting_nc(ii-i_nc+i_w_c,jj-j_nc+j_w_c,tt_dim,i_source,:)

                                        enddo
                                    enddo

                                    !Subtract the additional local emissions from the nonlocal using the new scheme
                                    nonlocal_correction(tt_dim,i_source,:)=0.

                                    !           do jj=-1-floor((EMEP_grid_interpolation_size-1.)*0.5),+1+floor((EMEP_grid_interpolation_size-1.)*0.5)
                                    !           do ii=-1-floor((EMEP_grid_interpolation_size-1.)*0.5),+1+floor((EMEP_grid_interpolation_size-1.)*0.5)
                                    do jj=-1-ceiling((EMEP_grid_interpolation_size-1.)*0.5),+1+ceiling((EMEP_grid_interpolation_size-1.)*0.5)
                                        do ii=-1-ceiling((EMEP_grid_interpolation_size-1.)*0.5),+1+ceiling((EMEP_grid_interpolation_size-1.)*0.5)

                                            ii_nc=ii+i_nc
                                            jj_nc=jj+j_nc
                                            ii_w=ii+ii_w0
                                            jj_w=jj+jj_w0

                                            if (jj.ne.0.or.ii.ne.0) then
                                                !First weight is emission, the second is area
                                                nonlocal_correction(tt_dim,i_source,:)=nonlocal_correction(tt_dim,i_source,:) &
                                                    -lc_var3d_nc(ii_w0,jj_w0,ii_nc,jj_nc,tt,lc_local_nc_index,i_source,:)*proxy_weighting_nc(ii+i_w_c,jj+j_w_c,tt_dim,i_source,:)*area_weighting_nc(i,j,ii_w0,jj_w0,tt_dim,i_source) &
                                                    -lc_var3d_nc(ii_w,jj_w,i_nc,j_nc,tt,lc_local_nc_index,i_source,:)*proxy_weighting_nc(ii+i_w_c,jj+j_w_c,tt_dim,i_source,:)*area_weighting_nc(i,j,ii_w,jj_w,tt_dim,i_source)
                                            endif

                                        enddo
                                    enddo

                                    subgrid(i,j,tt,emep_nonlocal_subgrid_index,i_source,:)=subgrid(i,j,tt,emep_nonlocal_subgrid_index,i_source,:)+nonlocal_correction(tt_dim,i_source,:)
                                    subgrid(i,j,tt,emep_local_subgrid_index,i_source,:)=subgrid(i,j,tt,emep_local_subgrid_index,i_source,:)-nonlocal_correction(tt_dim,i_source,:)
                                    subgrid(i,j,tt,emep_subgrid_index,i_source,:)=subgrid(i,j,tt,emep_nonlocal_subgrid_index,i_source,:)+subgrid(i,j,tt,emep_local_subgrid_index,i_source,:)

                                    !Averaged over time for diagnostic purposes only
                                    nonlocal_correction_average(i_source,:)=nonlocal_correction_average(i_source,:)+nonlocal_correction(tt_dim,i_source,:)

                                endif !use subgrid

                            enddo
                        enddo


                    endif !End if calculate_source
                enddo !End source loop

                if (tt.eq.t_end) then
                    nonlocal_correction_average=nonlocal_correction_average/subgrid_dim(x_dim_index)/subgrid_dim(y_dim_index)/subgrid_dim(t_dim_index)
                    do i_pollutant=1,n_emep_pollutant_loop
                        write(fmt,'(A,I0,A)') '(', n_source_index, 'es12.4)'
                        write(unit_logfile,fmt) 'Nonlocal correction for proxy weighting ('//trim(pollutant_file_str(pollutant_loop_index(i_pollutant)))//') = ',nonlocal_correction_average(:,i_pollutant)
                    enddo
                endif

            endif

            if (mod(j,1).eq.0) write(*,'(a,i5,a,i5)') 'Gridding EMEP for hour ',tt,' of ',subgrid_dim(t_dim_index)

        enddo   !End time loop

        !Create the all source version of the local and nonlocal contribution after calculating all the source contributions
        !The nonlocal contribution uses the difference between the local and total, here the total is based on the area interpolation. Is this correct?
        subgrid(:,:,:,emep_local_subgrid_index,allsource_index,:)=0.
        subgrid(:,:,:,emep_subgrid_index,allsource_index,:)=0.
        subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,:)=0. !-subgrid(:,:,:,emep_subgrid_index,allsource_index,:)

        count=0
        do i_source=1,n_source_index
            !do i_source=1,n_source_calculate_index
            if (calculate_source(i_source).or.calculate_EMEP_source(i_source).and.i_source.ne.allsource_index) then

                !Check values for local and totals for each source
                !write(*,*) trim(source_file_str(i_source))
                do i_pollutant=1,n_emep_pollutant_loop
                    if (minval(subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,i_pollutant)).lt.0.0) then
                        write(unit_logfile,'(A,A,f12.4,A)') 'WARNING: Min nonlocal source less than 0 for ',trim(source_file_str(i_source))//' '//trim(pollutant_file_str(pollutant_loop_index(i_pollutant))),minval(subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,emep_subsource)),' Setting to 0 and adding to local'
                    endif
                enddo

                !Set any negative nonlocal to 0 and add the value back into the local. Indicates a problem with the moving window method
                do j=1,subgrid_dim(y_dim_index)
                    do i=1,subgrid_dim(x_dim_index)
                        where (subgrid(i,j,:,emep_nonlocal_subgrid_index,i_source,:).lt.0.)
                            subgrid(i,j,:,emep_local_subgrid_index,i_source,:)=subgrid(i,j,:,emep_local_subgrid_index,i_source,:)-subgrid(i,j,:,emep_nonlocal_subgrid_index,i_source,:)
                            subgrid(i,j,:,emep_nonlocal_subgrid_index,i_source,:)=0.
                        endwhere
                    enddo
                enddo

                !Add the local subgrid sources together to get an allsource local contribution
                subgrid(:,:,:,emep_local_subgrid_index,allsource_index,:)=subgrid(:,:,:,emep_local_subgrid_index,allsource_index,:)+subgrid(:,:,:,emep_local_subgrid_index,i_source,:)
                subgrid(:,:,:,emep_subgrid_index,i_source,:)=subgrid(:,:,:,emep_nonlocal_subgrid_index,i_source,:)+subgrid(:,:,:,emep_local_subgrid_index,i_source,:)
                do i_pollutant=1,n_pollutant_loop
                    !write(*,'(a,2i,3f12.1)')'Sum  (subgrid_index,i_source,comp_EMEP,original): ',i_source,i_pollutant,sum(subgrid(:,:,:,emep_subgrid_index,i_source,i_pollutant)),sum(comp_EMEP_subgrid(:,:,:,pollutant_loop_index(i_pollutant))),sum(orig_EMEP_subgrid(:,:,:,pollutant_loop_index(i_pollutant)))
                enddo

                !Add up the total EMEP for all source (will be averaged with count)
                subgrid(:,:,:,emep_subgrid_index,allsource_index,:)=subgrid(:,:,:,emep_subgrid_index,allsource_index,:)+subgrid(:,:,:,emep_subgrid_index,i_source,:)
                count=count+1
                do i_pollutant=1,n_pollutant_loop
                    if (minval(subgrid(:,:,:,emep_subgrid_index,i_source,i_pollutant)).lt.0.0) then
                        write(unit_logfile,'(A,A,f12.4)') 'ERROR: Min total source less than 0 for ',trim(source_file_str(i_source))//' '//trim(pollutant_file_str(pollutant_loop_index(pollutant_loop_index(i_pollutant)))),minval(subgrid(:,:,:,emep_subgrid_index,i_source,i_pollutant))
                        stop
                    endif
                enddo

            endif
        enddo

        !Set the allsource nonlocal value to the average of the remainder. This can be negative
        if (count.gt.0) then
            subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,:)=(subgrid(:,:,:,emep_subgrid_index,allsource_index,:)/real(count)-subgrid(:,:,:,emep_local_subgrid_index,allsource_index,:))
            !write(*,*) calculate_EMEP_additional_grid_flag,sum(subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,:))
            !stop
        endif

        do i_pollutant=1,n_emep_pollutant_loop
            if (minval(subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,i_pollutant)).lt.0.0) then
                write(unit_logfile,'(A,f12.4,A)') 'WARNING: Min nonlocal allsource less than 0 with '//trim(source_file_str(allsource_index))//' '//trim(pollutant_file_str(pollutant_loop_index(i_pollutant))),minval(subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,i_pollutant)),' Setting to 0'
            endif
        enddo

        !Remove any negative values.
        do j=1,subgrid_dim(y_dim_index)
            do i=1,subgrid_dim(x_dim_index)
                where (subgrid(i,j,:,emep_nonlocal_subgrid_index,allsource_index,:).lt.0.)
                    subgrid(i,j,:,emep_nonlocal_subgrid_index,allsource_index,:)=0.
                endwhere
            enddo
        enddo

        !write(*,'(a,2f12.1)')'Sum (local,nonlocal): ',sum(subgrid(:,:,:,emep_local_subgrid_index,allsource_index,:)), sum(subgrid(:,:,:,emep_nonlocal_subgrid_index,allsource_index,:))

        !Add up the sources and calculate fractions
        do i_pollutant=1,n_emep_pollutant_loop
            !write(*,'(a,i,f12.1)')'Sum before  (emep_subgrid allsource): ',i_pollutant,sum(subgrid(:,:,:,emep_subgrid_index,allsource_index,i_pollutant))
        enddo
        subgrid(:,:,:,emep_subgrid_index,:,:)=subgrid(:,:,:,emep_nonlocal_subgrid_index,:,:)+subgrid(:,:,:,emep_local_subgrid_index,:,:)
        do i_pollutant=1,n_emep_pollutant_loop
            !write(*,'(a,i,f12.1)')'Sum after (emep_subgrid allsource): ',i_pollutant,sum(subgrid(:,:,:,emep_subgrid_index,allsource_index,i_pollutant))
        enddo
        subgrid(:,:,:,emep_frac_subgrid_index,:,:)=subgrid(:,:,:,emep_local_subgrid_index,:,:)/subgrid(:,:,:,emep_subgrid_index,:,:)
        do i_pollutant=1,n_emep_pollutant_loop
            if (minval(subgrid(:,:,:,emep_subgrid_index,allsource_index,:)).lt.0.0) then
                write(unit_logfile,'(A,f12.4)')'ERROR: Minimum total allsource less than 0 with '//trim(source_file_str(allsource_index))//' '//trim(pollutant_file_str(pollutant_loop_index(i_pollutant))),minval(subgrid(:,:,:,emep_subgrid_index,allsource_index,i_pollutant))
                stop
            endif
        enddo

        !Check if the additional EMEP calculation is to be carried out and set parameters
        EMEP_grid_interpolation_size=EMEP_grid_interpolation_size_saved
        if (calculate_EMEP_additional_grid_flag) then
            subgrid(:,:,:,emep_additional_local_subgrid_index,:,:)=subgrid(:,:,:,emep_local_subgrid_index,:,:)
            subgrid(:,:,:,emep_additional_nonlocal_subgrid_index,:,:)=subgrid(:,:,:,emep_nonlocal_subgrid_index,:,:)
        endif

        if (allocated(weighting_nc)) deallocate(weighting_nc)
        if (allocated(area_weighting_nc)) deallocate(area_weighting_nc)
        if (allocated(total_weighting_nc)) deallocate(total_weighting_nc)
        if (allocated(proxy_weighting_nc)) deallocate(proxy_weighting_nc)
        if (allocated(weighting_subgrid)) deallocate(weighting_subgrid)
        if (allocated(crossreference_weighting_to_emep_subgrid)) deallocate(crossreference_weighting_to_emep_subgrid)
        if (allocated(nonlocal_correction)) deallocate(nonlocal_correction)
        if (allocated(nonlocal_correction_average)) deallocate(nonlocal_correction_average)
        if (allocated(EMEP_local_contribution)) deallocate (EMEP_local_contribution)
        if (allocated(EMEP_local_contribution_from_in_region)) deallocate (EMEP_local_contribution_from_in_region)



    end subroutine uEMEP_subgrid_EMEP

    subroutine nlreg_uEMEP_subgrid_EMEP_from_in_region

        ! contribution to centre of each EMEP cell from in-region, within MW and outside MW, for each region
        real, allocatable :: EMEP_local_from_in_region(:, :, :, :, :, :) ! (x,y,region,t,source,pollutant)
        real, allocatable :: EMEP_semilocal_from_in_region(:, :, :, :, :, :) ! (x,y,region,t,source,pollutant)
        ! arrays for accumulating data for a single EMEP cell of the above arrays (is this actually faster than accumulating directly in the big arrays????)
        real, allocatable :: temp_EMEP_local_from_in_region(:, :, :, :)
        real, allocatable :: temp_EMEP_semilocal_from_in_region(:, :, :, :)
        ! additional increment for the centre of each EMEP cell, for each region
        real, allocatable :: EMEP_additional_increment_from_in_region(:, :, :, :, :, :)  ! (x,y,region,t,source,pollutant)
        ! additional increment for a single EMEP cell, for each region, to be accumulated
        real, allocatable :: temp_EMEP_additional_increment_from_in_region(:, :, :, :)  ! (region,t,source,pollutant)
        ! additional increment of a single big LF grid to a single receptor EMEP grid (before weighing by region)
        real, allocatable :: EMEP_additional_increment_current_lfgrid(:, :, :)  ! (t,source,pollutant)
        ! weighting of the additional increment, for each region
        real, allocatable :: weights_EMEP_additional_increment_current_lfgrid(:)  ! (region)
        ! weighting of an LF cell for local and semilocal contribution
        real weighting_value_local,weighting_value_semilocal
        ! indexers for looping
        integer i,j             ! target subgrids
        integer ii,jj           ! EMEP grids
        integer i_dist,j_dist   ! LF dimensions
        integer iiii,jjjj       ! small LF grids within a big LF grid
        integer i_sub,j_sub     ! subsamples of an EMEP grid
        integer i_region        ! region dimension
        integer i_source        ! source dimension
        ! indexers for determining positioning of additional LF grids
        integer ii_start,jj_start
        ! additinoal indexers
        integer iii,jjj
        integer iii_nc,jjj_nc
        ! displacement distances in LF grids (whole grids)
        integer x_dist,y_dist
        integer xdist_big,ydist_big
        integer xdist_small,ydist_small
        integer xdist_small_first,ydist_small_first
        integer idist_small,jdist_small
        integer max_x_dist,max_y_dist
        ! indexers for the location of a 1x1 LF grid cell in the EMEP grid
        integer iiii_nc,jjjj_nc
        ! indexers for the location of a 1x1 LF grid cell in the extended EMEP grid
        integer iiii_extended,jjjj_extended
        ! indexers for finding local contributions in LF array
        integer lc_index,lc_additional_index
        ! counters
        integer counter,counter_local,counter_semilocal
        ! variables to hold a region index
        integer current_region_index
        ! fractional position of uEMEP subgrid within an EMEP grid
        real ii_frac_target,jj_frac_target
        ! location of subgrid in EMEP's coordinate system
        real x_temp,y_temp
        ! distance between target subgrid and an EMEP grid centre or EMEP subsample location
        real x_dist_target,y_dist_target
        ! half-size of the moving window
        real n_EMEP_grids_to_edge_of_moving_window
        ! fraction of an EMEP grid that is in the correct region
        real current_EMEP_region_fraction
        ! weighting for area interpolation
        real weighting_val,weight_check

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '====================================================================================================='
        write(unit_logfile,'(A)') 'Calculation local and semilocal from-in-region EMEP contributions to subgrids (nlreg_uEMEP_subgrid_EMEP_from_in_region)'
        write(unit_logfile,'(A)') '====================================================================================================='

        write(unit_logfile,'(A,8I8)') 'dims: (x_emep,y_emep,reg,time,source,pollutant,xdist,ydist) = ',dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),nlreg_n_regions,subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop,dim_length_nc(xdist_dim_nc_index),dim_length_nc(ydist_dim_nc_index)

        ! ************************************************
        ! PART 1: Calculate contributions to the EMEP grid
        ! ************************************************

        ! Indices in lc_var3d_nc to find the normal and additional local contributions
        lc_index = lc_local_nc_loop_index(local_fraction_grid_for_EMEP_grid_interpolation)
        lc_additional_index = lc_local_nc_loop_index(local_fraction_grid_for_EMEP_additional_grid_interpolation)

        ! distance in x or y (EMEP grid) from receptor subgrid to edge of moving window, in units of EMEP grids
        ! (normally whole number, but not if EMEP_grid_interpolation_size is odd number)
        n_EMEP_grids_to_edge_of_moving_window = EMEP_grid_interpolation_size*0.5
        write(unit_logfile,'(A,f8.1)') 'n_EMEP_grids_to_edge_of_moving_window=', n_EMEP_grids_to_edge_of_moving_window

        ! Calculate local and semilocal contributions to centre of all EMEP grids
        write(unit_logfile,'(A)') 'Allocating arrays for calculating in-region local and semilocal contribution to the EMEP grid'

        ! Allocate arrays
        allocate(EMEP_local_from_in_region(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),nlreg_n_regions,subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop))
        EMEP_local_from_in_region = 0.0
        allocate(EMEP_semilocal_from_in_region(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),nlreg_n_regions,subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop))
        EMEP_semilocal_from_in_region = 0.0
        allocate(temp_EMEP_local_from_in_region(nlreg_n_regions,subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop))
        allocate(temp_EMEP_semilocal_from_in_region(nlreg_n_regions,subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop))

        ! Loop over all EMEP cells and calculate local and semilocal contributions when moving window is centered at the centre of the cell
        do ii = 1, dim_length_nc(x_dim_nc_index)
            do jj = 1, dim_length_nc(y_dim_nc_index)
                temp_EMEP_local_from_in_region = 0.0
                temp_EMEP_semilocal_from_in_region = 0.0
                ! Loop over all (small) LF contribution cells to this EMEP cell
                do i_dist = 1, dim_length_nc(xdist_dim_nc_index)
                    do j_dist = 1, dim_length_nc(ydist_dim_nc_index)
                        ! number of grids displaced relative to the receptor grid cell
                        x_dist = i_dist - xdist_centre_nc
                        y_dist = j_dist - ydist_centre_nc
                        ! index in the netcdf file of this cell
                        iiii_nc = ii + x_dist
                        jjjj_nc = jj + y_dist
                        ! index in the extended EMEP grid of this cell
                        iiii_extended = iiii_nc + nlreg_ngrid_extended_margin
                        jjjj_extended = jjjj_nc + nlreg_ngrid_extended_margin

                        ! Check if this cell is covered by the extended grid. If not, assume it is completely outside all regions
                        if (.not. (iiii_extended >= 1 .and. iiii_extended <= nlreg_nx_EMEP_extended .and. jjjj_extended >= 1 .and. jjjj_extended <= nlreg_ny_EMEP_extended)) then
                            cycle
                        end if

                        ! Loop over all regions
                        do i_region = 1, nlreg_n_regions

                            ! get fraction of the LF grid that is in the region
                            current_region_index = nlreg_regionindex_loop_index(i_region)
                            current_EMEP_region_fraction = nlreg_regionfraction_per_EMEP_extended_grid(iiii_extended, jjjj_extended, i_region)

                            ! Determine the weights to use for this LF cell
                            ! for local, the weight is the area fraction that is in the region and inside the moving window
                            ! for semilocal, the weight is the area fraction that is in the region but outside the moving window
                            if (current_EMEP_region_fraction <= 0) then
                                ! No part of this EMEP grid is within the region
                                cycle
                            else if (abs(x_dist) <= n_EMEP_grids_to_edge_of_moving_window-1 .and. abs(y_dist) <= n_EMEP_grids_to_edge_of_moving_window-1) then
                                ! this LF grid is sure to be completely within the moving window
                                weighting_value_local = current_EMEP_region_fraction
                                weighting_value_semilocal = 0.0
                            else if (abs(x_dist) >= n_EMEP_grids_to_edge_of_moving_window+1 .or. abs(y_dist) >= n_EMEP_grids_to_edge_of_moving_window+1) then
                                ! this LF grid is sure to be completely outside the moving window
                                weighting_value_local = 0.0
                                weighting_value_semilocal = current_EMEP_region_fraction
                            else
                                ! this LF grid might be partly covered by the moving window
                                ! -> we must go through all subsamples of that grid to find the overlap between region and moving window
                                counter_local = 0       ! count the subsample grids inside-moving-window & in-region
                                counter_semilocal = 0   ! count the subsample grids outside-moving-window & in-region
                                do i_sub = 1, nlreg_n_subsamples_per_EMEP_grid
                                    do j_sub = 1, nlreg_n_subsamples_per_EMEP_grid
                                        ! first check if this subsample is in the region
                                        if (nlreg_EMEP_extended_subsample_region_id(i_sub,j_sub,iiii_extended,jjjj_extended) == current_region_index) then
                                            ! deduce x- and y-distance (in number of EMEP grids) from this subsample location to the midpoint of the receptor grid-cell
                                            x_dist_target = x_dist + (i_sub-0.5)/nlreg_n_subsamples_per_EMEP_grid - 0.5
                                            y_dist_target = y_dist + (j_sub-0.5)/nlreg_n_subsamples_per_EMEP_grid - 0.5
                                            ! use these distances to determine whether the subsample is inside the moving window
                                            if (abs(x_dist_target) <= n_EMEP_grids_to_edge_of_moving_window .and. abs(y_dist_target) <= n_EMEP_grids_to_edge_of_moving_window) then
                                                ! this subsample is inside the moving window
                                                counter_local = counter_local + 1
                                            else
                                                ! this subsample is outside the moving window
                                                counter_semilocal = counter_semilocal + 1
                                            end if
                                        end if
                                    end do
                                end do
                                ! Divide counters by total number of subsamples per EMEP grid to get the area fraction
                                weighting_value_local = counter_local*1.0/nlreg_n_subsamples_per_EMEP_grid**2
                                weighting_value_semilocal = counter_semilocal*1.0/nlreg_n_subsamples_per_EMEP_grid**2
                            end if
                            temp_EMEP_local_from_in_region(i_region,:,:,:) = temp_EMEP_local_from_in_region(i_region,:,:,:) + lc_var3d_nc(i_dist,j_dist,ii,jj,:,lc_index,1:n_source_index,:)*weighting_value_local
                            temp_EMEP_semilocal_from_in_region(i_region,:,:,:) = temp_EMEP_semilocal_from_in_region(i_region,:,:,:) + lc_var3d_nc(i_dist,j_dist,ii,jj,:,lc_index,1:n_source_index,:)*weighting_value_semilocal
                        end do ! i_region = 1, nlreg_n_regions

                    end do
                end do
                EMEP_local_from_in_region(ii,jj,:,:,:,:) = temp_EMEP_local_from_in_region
                EMEP_semilocal_from_in_region(ii,jj,:,:,:,:) = temp_EMEP_semilocal_from_in_region
            end do
        end do

        ! Deallocate temporary arrays
        deallocate(temp_EMEP_local_from_in_region)
        deallocate(temp_EMEP_semilocal_from_in_region)

        ! Calculate the additional incremental contribution to each EMEP grid for each region
        ! NB: we let the interpolation domain extend to the very edge of the LF data, not limited by the EMEP_additional_grid_interpolation_size (this is only used to check >0 to determine whether we should calculate it or not)
        if (EMEP_additional_grid_interpolation_size > 0.0) then

            write(unit_logfile,'(A)') 'Allocating arrays for additional increment calculation'

            allocate(EMEP_additional_increment_from_in_region(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index),nlreg_n_regions,subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop))
            allocate(temp_EMEP_additional_increment_from_in_region(nlreg_n_regions,subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop))
            allocate(EMEP_additional_increment_current_lfgrid(subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop))
            allocate(weights_EMEP_additional_increment_current_lfgrid(nlreg_n_regions))

            write(unit_logfile,'(A)') 'Calculating additional increment to all EMEP grids'

            ! Use the starting position of the read in EMEP file to initialise the starting point
            ii_start = mod(dim_start_EMEP_nc(x_dim_nc_index)-1,local_fraction_grid_size(2))
            jj_start = mod(dim_start_EMEP_nc(y_dim_nc_index)-1,local_fraction_grid_size(2))

            ! Deduce the max distance (+/-) we have LF data for
            max_x_dist = (dim_length_nc(xdist_dim_nc_index)-1)/2
            max_y_dist = (dim_length_nc(ydist_dim_nc_index)-1)/2

            ! Loop over all EMEP grids
            EMEP_additional_increment_from_in_region = 0.0
            do ii = 1, dim_length_nc(x_dim_nc_index)
                do jj = 1, dim_length_nc(y_dim_nc_index)
                    ! Initialize the additional increment of this EMEP cell to zero
                    temp_EMEP_additional_increment_from_in_region = 0.0
                    ! EMEP grid index of bottom-left-corner-cell of the additional grid associated with that EMEP grid
                    iii = int((ii-1+ii_start)/local_fraction_grid_size(2))*local_fraction_grid_size(2) + 1 - ii_start
                    jjj = int((jj-1+jj_start)/local_fraction_grid_size(2))*local_fraction_grid_size(2) + 1 - jj_start

                    ! Loop over all big (additional) LF grids that give contributions to the EMEP grid
                    do i_dist = 1, dim_length_nc(xdist_dim_nc_index)
                        do j_dist = 1, dim_length_nc(ydist_dim_nc_index)

                            ! Initialize additional increment to the total additional contribution from EMEP from this cell (all times, sources and pollutants)
                            EMEP_additional_increment_current_lfgrid = lc_var3d_nc(i_dist,j_dist,ii,jj,:,lc_additional_index,1:n_source_index,:)
                            
                            ! the x_dist and y_dist of this big LF grid (..., -1, 0, 1, ...)
                            xdist_big = i_dist - xdist_centre_nc
                            ydist_big = j_dist - ydist_centre_nc
                            ! deduce what is the xdist and ydist in the 1x1 LF grid of the lower-left EMEP cell falling within this big LF grid
                            xdist_small_first = iii - ii + xdist_big*local_fraction_grid_size(2)
                            ydist_small_first = jjj - jj + ydist_big*local_fraction_grid_size(2)
                            ! loop over all EMEP grids contained within this big LF grid
                            ! Reset weights and counter
                            weights_EMEP_additional_increment_current_lfgrid = 0.0
                            counter = 0  ! count the grids not covered by the small LF grid
                            do iiii = 1, local_fraction_grid_size(2)
                                do jjjj = 1, local_fraction_grid_size(2)
                                    ! Deduce the xdist and ydist of this EMEP grid in the small LF domain
                                    xdist_small = xdist_small_first - 1 + iiii
                                    ydist_small = ydist_small_first - 1 + jjjj
                                    ! and corresponding index in the LF array
                                    idist_small = xdist_small + xdist_centre_nc
                                    jdist_small = ydist_small + ydist_centre_nc
                                    if (abs(xdist_small) > max_x_dist .or. abs(ydist_small) > max_y_dist) then
                                        ! This grid is NOT covered by 1x1 LF data, so add its region coverage to the weight
                                        ! find index in the normal EMEP grid
                                        iiii_nc = ii + xdist_small
                                        jjjj_nc = jj + ydist_small
                                        ! find corresponding index in the extended EMEP grid of region fractions
                                        iiii_extended = iiii_nc + nlreg_ngrid_extended_margin
                                        jjjj_extended = jjjj_nc + nlreg_ngrid_extended_margin
                                        ! check if this is within the extended region mask grid
                                        if (iiii_extended >= 1 .and. iiii_extended <= nlreg_nx_EMEP_extended .and. jjjj_extended >= 1 .and. jjjj_extended <= nlreg_ny_EMEP_extended) then
                                            weights_EMEP_additional_increment_current_lfgrid = weights_EMEP_additional_increment_current_lfgrid + nlreg_regionfraction_per_EMEP_extended_grid(iiii_extended, jjjj_extended, :)
                                        end if
                                        ! NB: If the cell is is outside the extended EMEP grid, it is considered to be outside all the regions. If the extended grid size is set properly, this should never be the case for the receptor EMEP grids (ii,jj) that overlap with the target grid.
                                        counter = counter + 1
                                    else
                                        ! The grid is covered by 1x1 LF data: subtract the 1x1 LF from the additional increment
                                        EMEP_additional_increment_current_lfgrid=EMEP_additional_increment_current_lfgrid-lc_var3d_nc(idist_small,jdist_small,ii,jj,:,lc_index,1:n_source_index,:)
                                    end if
                                end do
                            end do
                            ! ensure the additional increment is not smaller than zero
                            where (EMEP_additional_increment_current_lfgrid < 0.0) EMEP_additional_increment_current_lfgrid = 0.0
                            ! normalize the weights by the number of grids
                            ! NB: if counter = 0, then the weights are zero so we can go to next LF source grid
                            if (counter > 0) then
                                ! normalize weights by the number of cells summed over
                                weights_EMEP_additional_increment_current_lfgrid = weights_EMEP_additional_increment_current_lfgrid/counter
                                ! For each region, multiply the additional increment by the weight calculated for that region
                                ! and accumulate this in the array for the total additional increment (to be accumulated over all the big LF cells)
                                do i_region = 1, nlreg_n_regions
                                    temp_EMEP_additional_increment_from_in_region(i_region,:,:,:)=temp_EMEP_additional_increment_from_in_region(i_region,:,:,:)+EMEP_additional_increment_current_lfgrid*weights_EMEP_additional_increment_current_lfgrid(i_region)
                                end do
                            end if
                        end do
                    end do
                    ! The additional increment has now been accumulated over all xdist and ydist source grids and weighted for each region
                    ! So now it can be inserted into the main array
                    EMEP_additional_increment_from_in_region(ii,jj,:,:,:,:) = temp_EMEP_additional_increment_from_in_region
                end do
            end do

            EMEP_semilocal_from_in_region = EMEP_semilocal_from_in_region + EMEP_additional_increment_from_in_region

            deallocate(temp_EMEP_additional_increment_from_in_region)
            deallocate(EMEP_additional_increment_current_lfgrid)
            deallocate(weights_EMEP_additional_increment_current_lfgrid)
            deallocate(EMEP_additional_increment_from_in_region)

        end if

        ! *****************************************
        ! PART 2: Interpolate to the target subgrid
        ! *****************************************

        write(unit_logfile,'(A)') 'Allocating arrays for in-region local and semilocal contribution to the target grid'

        ! Allocate the arrays for holding the results and initialize them to zero
        if (allocated(nlreg_subgrid_EMEP_local_from_in_region)) deallocate(nlreg_subgrid_EMEP_local_from_in_region)
        allocate(nlreg_subgrid_EMEP_local_from_in_region(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop))
        nlreg_subgrid_EMEP_local_from_in_region = 0.0
        if (allocated(nlreg_subgrid_EMEP_semilocal_from_in_region)) deallocate(nlreg_subgrid_EMEP_semilocal_from_in_region)
        allocate(nlreg_subgrid_EMEP_semilocal_from_in_region(subgrid_dim(x_dim_index),subgrid_dim(y_dim_index),subgrid_dim(t_dim_index),n_source_index,n_pollutant_loop))
        nlreg_subgrid_EMEP_semilocal_from_in_region = 0.0

        write(unit_logfile,'(A)') 'Interpolating local and semilocal contribution from-in-region to the target subgrid'

        ! go through all target subgrids and interpolate the EMEP contributions calculated above
        do i = 1, subgrid_dim(x_dim_index)
            do j = 1, subgrid_dim(y_dim_index)

                ! Find position along the region dimension of this region index
                if (nlreg_subgrid_region_index(i,j) > 0) then
                    i_region = nlreg_regionindex_loop_back_index(nlreg_subgrid_region_index(i,j))
                else
                    ! this subgrid is not in any region, so keep it as 0
                    cycle
                end if

                ! Find which EMEP grid the current subgrid is in
                ii=crossreference_target_to_emep_subgrid(i,j,x_dim_index)
                jj=crossreference_target_to_emep_subgrid(i,j,y_dim_index)

                ! Find out where in this EMEP grid we are in EMEP grid coordinates
                ! E.g. (0,0) for lower-left corner and (1,1) for upper right corner
                call LL2PROJ(lon_subgrid(i,j),lat_subgrid(i,j),x_temp,y_temp,EMEP_projection_attributes,EMEP_projection_type)
                ! fractional position inside the EMEP grid (center is (0,0), lower-left corner is (-0.5,-0.5), upper-right corner is (+0.5,+0.5))
                ii_frac_target = (x_temp-var1d_nc(ii,x_dim_nc_index))/dgrid_nc(x_dim_nc_index)
                jj_frac_target = (y_temp-var1d_nc(jj,y_dim_nc_index))/dgrid_nc(y_dim_nc_index)
                ! verify that this is within 0-1 (if not, crossreference has gone wrong...)
                if (ii_frac_target < -0.5 .or. ii_frac_target > 0.5 .or. jj_frac_target < -0.5 .or. jj_frac_target > 0.5) then
                    write(unit_logfile,'(A,2I12)') 'Something went wrong with locating target subgrid within EMEP grid!',ii_frac_target,jj_frac_target
                    stop
                end if

                ! Interpolate to the target grid
                if (EMEP_grid_interpolation_flag == 0) then
                    ! No interpolation: just pick the value of the EMEP cell we are in
                    nlreg_subgrid_EMEP_local_from_in_region(i,j,:,:,:) = EMEP_local_from_in_region(ii,jj,i_region,:,:,:)
                    nlreg_subgrid_EMEP_semilocal_from_in_region(i,j,:,:,:) = EMEP_semilocal_from_in_region(ii,jj,i_region,:,:,:)
                else if (EMEP_grid_interpolation_flag == 6) then
                    ! Use area-weighted interpolation
                    weight_check = 0.0
                    ! loop over a 3x3 EMEP cell domain centered at the closest EMEP cell to the target subgrid
                    do iii = -1,1
                        do jjj = -1,1
                            iii_nc = ii + iii
                            jjj_nc = jj + jjj
                            ! verify the EMEP grid covers this index
                            if (.not. (iii_nc >= 1 .and. iii_nc <= dim_length_nc(x_dim_nc_index) .and. jjj_nc >= 1 .and. jjj_nc <= dim_length_nc(y_dim_nc_index))) then
                                write(unit_logfile,'(A)') 'ERROR: EMEP grid did not go far enough out to allow area interpolation to the target grid!'
                                stop
                            end if
                            ! calculate the weighting value as the area fraction of an EMEP cell centered at the target subgrid that falls within this EMEP cell
                            x_dist_target = iii - ii_frac_target
                            y_dist_target = jjj - jj_frac_target
                            weighting_val = max(0.0, (1.0 - abs(x_dist_target))) * max(0.0, (1.0 - abs(y_dist_target)))
                            weight_check = weight_check + weighting_val
                            ! use this weighting for the data at this EMEP cell
                            if (weighting_val > 0) then
                                nlreg_subgrid_EMEP_local_from_in_region(i,j,:,:,:) = nlreg_subgrid_EMEP_local_from_in_region(i,j,:,:,:) + EMEP_local_from_in_region(iii_nc,jjj_nc,i_region,:,:,:)*weighting_val
                                nlreg_subgrid_EMEP_semilocal_from_in_region(i,j,:,:,:) = nlreg_subgrid_EMEP_semilocal_from_in_region(i,j,:,:,:) + EMEP_semilocal_from_in_region(iii_nc,jjj_nc,i_region,:,:,:)*weighting_val
                            end if
                        end do
                    end do
                else
                    write(unit_logfile,'(A,I0)') 'ERROR: nlreg_uEMEP_subgrid_EMEP_from_in_region is not implemented for EMEP_grid_interpolation_flag =', EMEP_grid_interpolation_flag
                    stop
                end if
            end do
        end do

        deallocate(EMEP_local_from_in_region)
        deallocate(EMEP_semilocal_from_in_region)

        ! Set allsources to be the sum of only the sources we calculate for
        nlreg_subgrid_EMEP_local_from_in_region(:,:,:,allsource_index,:)=0
        nlreg_subgrid_EMEP_semilocal_from_in_region(:,:,:,allsource_index,:)=0
        do i_source = 1, n_source_index
            if (calculate_source(i_source).or.calculate_EMEP_source(i_source).and.i_source.ne.allsource_index) then
                nlreg_subgrid_EMEP_local_from_in_region(:,:,:,allsource_index,:)=nlreg_subgrid_EMEP_local_from_in_region(:,:,:,allsource_index,:)+nlreg_subgrid_EMEP_local_from_in_region(:,:,:,i_source,:)
                nlreg_subgrid_EMEP_semilocal_from_in_region(:,:,:,allsource_index,:)=nlreg_subgrid_EMEP_semilocal_from_in_region(:,:,:,allsource_index,:)+nlreg_subgrid_EMEP_semilocal_from_in_region(:,:,:,i_source,:)
            end if
        end do

    end subroutine nlreg_uEMEP_subgrid_EMEP_from_in_region

end module subgrid_emep

