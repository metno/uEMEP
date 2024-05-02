module subgrid_meteo_emep

    use uemep_constants, only: pi
    use uemep_configuration
    
    implicit none
    private

    public :: uEMEP_subgrid_meteo_EMEP

contains

!==========================================================================
!   uEMEP_subgrid_meteo_EMEP
!   Bruce Rolstad Denby
!   MET Norway
!
!   This routine interpolates EMEP meteo data to the integral subgrid
!   using either nearest neighbour EMEP_meteo_grid_interpolation_flag=0
!   or area weighted interpolation EMEP_meteo_grid_interpolation_flag=1
!   It also generates cos and sin wind field data for dispersion modelling
!==========================================================================

    subroutine uEMEP_subgrid_meteo_EMEP

        use uEMEP_definitions

        implicit none

        integer i,j
        integer ii,jj,iii,jjj
        real, allocatable :: weighting_nc(:,:)
        integer t_start,t_end
        real xpos_min,xpos_max,ypos_min,ypos_max
        real xpos_area_min,xpos_area_max,ypos_area_min,ypos_area_max
        integer i_nc,j_nc
        real angle_utm,angle_lcc,angle_utm2
        real, allocatable :: u_utm(:),v_utm(:),th(:),ff(:)
        real dlatx,dlaty
        ! real xpos_subgrid,ypos_subgrid
        real xpos_limit,ypos_limit
        real xpos_integral_subgrid,ypos_integral_subgrid

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Interpolating meteo data to subgrids (uEMEP_subgrid_meteo_EMEP)'
        write(unit_logfile,'(A)') '================================================================'

        if (.not.allocated(u_utm)) allocate (u_utm(dim_length_nc(time_dim_nc_index)))
        if (.not.allocated(v_utm)) allocate (v_utm(dim_length_nc(time_dim_nc_index)))
        if (.not.allocated(th)) allocate (th(dim_length_nc(time_dim_nc_index)))
        if (.not.allocated(ff)) allocate (ff(dim_length_nc(time_dim_nc_index)))

        !Set the last meteo data to be the previous for the external time loop before updating
        if (use_single_time_loop_flag) then
            if (t_loop.gt.start_time_loop_index) then
                last_meteo_subgrid(:,:,:)=meteo_subgrid(:,:,1,:)
            endif
        endif

        !Initialise all meteo subgrid fields
        meteo_subgrid=0.

        !Set the time dimensions for transfering the alternative meteorology which has a time index that starts at 0
        t_start=1
        t_end=subgrid_dim(t_dim_index)


        write(unit_logfile,'(A,I4)')'Setting EMEP subgrid meteo data using method ',EMEP_meteo_grid_interpolation_flag

        !Loop through the integral subgrid and find those subgrids within EMEP grids and allocate values directly from EMEP grids. Nearest neighbour
        if (EMEP_meteo_grid_interpolation_flag.eq.0) then
            do j=1,integral_subgrid_dim(y_dim_index)
                do i=1,integral_subgrid_dim(x_dim_index)

                    !Assumes it is never on the edge of the EMEP grid, not limitted
                    i_nc=crossreference_integral_to_emep_subgrid(i,j,x_dim_index)
                    j_nc=crossreference_integral_to_emep_subgrid(i,j,y_dim_index)

                    meteo_subgrid(i,j,:,ugrid_subgrid_index)=var4d_nc(i_nc,j_nc,surface_level_nc,:,ugrid_nc_index,allsource_index,meteo_p_loop_index)
                    meteo_subgrid(i,j,:,vgrid_subgrid_index)=var4d_nc(i_nc,j_nc,surface_level_nc,:,vgrid_nc_index,allsource_index,meteo_p_loop_index)
                    meteo_subgrid(i,j,:,FFgrid_subgrid_index)=var4d_nc(i_nc,j_nc,surface_level_nc,:,FFgrid_nc_index,allsource_index,meteo_p_loop_index)
                    meteo_subgrid(i,j,:,u10_subgrid_index)=var3d_nc(i_nc,j_nc,:,u10_nc_index,allsource_index,meteo_p_loop_index)
                    meteo_subgrid(i,j,:,v10_subgrid_index)=var3d_nc(i_nc,j_nc,:,v10_nc_index,allsource_index,meteo_p_loop_index)
                    meteo_subgrid(i,j,:,FF10_subgrid_index)=var3d_nc(i_nc,j_nc,:,FF10_nc_index,allsource_index,meteo_p_loop_index)
                    meteo_subgrid(i,j,:,kz_subgrid_index)=var4d_nc(i_nc,j_nc,surface_level_nc,:,kz_nc_index,allsource_index,meteo_p_loop_index)
                    meteo_subgrid(i,j,:,hmix_subgrid_index)=var3d_nc(i_nc,j_nc,:,hmix_nc_index,allsource_index,meteo_p_loop_index)
                    meteo_subgrid(i,j,:,logz0_subgrid_index)=var3d_nc(i_nc,j_nc,:,logz0_nc_index,allsource_index,meteo_p_loop_index)
                    meteo_subgrid(i,j,:,invL_subgrid_index)=var3d_nc(i_nc,j_nc,:,invL_nc_index,allsource_index,meteo_p_loop_index)
                    meteo_subgrid(i,j,:,inv_FFgrid_subgrid_index)=var4d_nc(i_nc,j_nc,surface_level_nc,:,inv_FFgrid_nc_index,allsource_index,meteo_p_loop_index)
                    meteo_subgrid(i,j,:,inv_FF10_subgrid_index)=var3d_nc(i_nc,j_nc,:,inv_FF10_nc_index,allsource_index,meteo_p_loop_index)
                    meteo_subgrid(i,j,:,ustar_subgrid_index)=var3d_nc(i_nc,j_nc,:,ustar_nc_index,allsource_index,meteo_p_loop_index)
                    meteo_subgrid(i,j,:,J_subgrid_index)=var4d_nc(i_nc,j_nc,surface_level_nc,:,J_nc_index,allsource_index,meteo_p_loop_index)
                    meteo_subgrid(i,j,:,t2m_subgrid_index)=var3d_nc(i_nc,j_nc,:,t2m_nc_index,allsource_index,meteo_p_loop_index)
                    !meteo_subgrid(i,j,:,precip_subgrid_index)=var4d_nc(i_nc,j_nc,surface_level_nc,:,precip_nc_index,allsource_index,meteo_p_loop_index)
                    meteo_subgrid(i,j,:,precip_subgrid_index)=var3d_nc(i_nc,j_nc,:,precip_nc_index,allsource_index,meteo_p_loop_index)

                    if (use_alternative_meteorology_flag) then
                        i_nc=crossreference_integral_to_meteo_nc_subgrid(i,j,x_dim_index)
                        j_nc=crossreference_integral_to_meteo_nc_subgrid(i,j,y_dim_index)
                        meteo_subgrid(i,j,t_start:t_end,ugrid_subgrid_index)=meteo_var3d_nc(i_nc,j_nc,t_start:t_end,ugrid_nc_index)
                        meteo_subgrid(i,j,t_start:t_end,vgrid_subgrid_index)=meteo_var3d_nc(i_nc,j_nc,t_start:t_end,vgrid_nc_index)
                        meteo_subgrid(i,j,t_start:t_end,u10_subgrid_index)=meteo_var3d_nc(i_nc,j_nc,t_start:t_end,u10_nc_index)
                        meteo_subgrid(i,j,t_start:t_end,v10_subgrid_index)=meteo_var3d_nc(i_nc,j_nc,t_start:t_end,v10_nc_index)
                        meteo_subgrid(i,j,t_start:t_end,FFgrid_subgrid_index)=meteo_var3d_nc(i_nc,j_nc,t_start:t_end,FFgrid_nc_index)
                        meteo_subgrid(i,j,t_start:t_end,FF10_subgrid_index)=meteo_var3d_nc(i_nc,j_nc,t_start:t_end,FF10_nc_index)
                        meteo_subgrid(i,j,t_start:t_end,hmix_subgrid_index)=meteo_var3d_nc(i_nc,j_nc,t_start:t_end,hmix_nc_index)
                        meteo_subgrid(i,j,t_start:t_end,logz0_subgrid_index)=meteo_var3d_nc(i_nc,j_nc,t_start:t_end,logz0_nc_index)
                        meteo_subgrid(i,j,t_start:t_end,invL_subgrid_index)=meteo_var3d_nc(i_nc,j_nc,t_start:t_end,invL_nc_index)
                        !meteo_subgrid(i,j,t_start:t_end,inv_FFgrid_subgrid_index)=meteo_var3d_nc(i_nc,j_nc,,t_start:t_end,inv_FFgrid_nc_index)
                        !meteo_subgrid(i,j,t_start:t_end,inv_FF10_subgrid_index)=meteo_var3d_nc(i_nc,j_nc,t_start:t_end,inv_FF10_nc_index)
                        meteo_subgrid(i,j,t_start:t_end,ustar_subgrid_index)=meteo_var3d_nc(i_nc,j_nc,t_start:t_end,ustar_nc_index)
                        meteo_subgrid(i,j,t_start:t_end,t2m_subgrid_index)=meteo_var3d_nc(i_nc,j_nc,t_start:t_end,t2m_nc_index)
                        meteo_subgrid(i,j,t_start:t_end,precip_subgrid_index)=meteo_var3d_nc(i_nc,j_nc,t_start:t_end,precip_nc_index)
                    endif

                    !Not properly implemented. Always false
                    if (use_alternative_z0_flag) then
                        meteo_subgrid(i,j,t_start:t_end,logz0_subgrid_index)=meteo_var3d_nc(i_nc,j_nc,t_start:t_end,logz0_nc_index)
                    endif

                    !        write(*,*) meteo_subgrid(i,j,t_start:t_end,precip_subgrid_index),var4d_nc(i_nc,j_nc,surface_level_nc,:,precip_nc_index,allsource_index,meteo_p_loop_index)
                enddo
            enddo

            !write(*,*) maxval(meteo_subgrid(:,:,t_start:t_end,precip_subgrid_index))
            !write(*,*) maxval(var4d_nc(:,:,surface_level_nc,:,precip_nc_index,allsource_index,meteo_p_loop_index))
            !stop
        endif

        !Area weighted interpolation of meteorology to integral grid
        if (EMEP_meteo_grid_interpolation_flag.ge.1) then

            xpos_limit=dgrid_nc(lon_nc_index)/2.
            ypos_limit=dgrid_nc(lat_nc_index)/2.

            allocate (weighting_nc(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))

            do j=1,integral_subgrid_dim(y_dim_index)
                do i=1,integral_subgrid_dim(x_dim_index)

                    !Assumes it is never on the edge of the EMEP grid. Something wrong if it fails
                    i_nc=crossreference_integral_to_emep_subgrid(i,j,x_dim_index)
                    j_nc=crossreference_integral_to_emep_subgrid(i,j,y_dim_index)

                    xpos_integral_subgrid=xproj_integral_subgrid(i,j)
                    ypos_integral_subgrid=yproj_integral_subgrid(i,j)

                    xpos_area_max=xpos_integral_subgrid+xpos_limit
                    xpos_area_min=xpos_integral_subgrid-xpos_limit
                    ypos_area_max=ypos_integral_subgrid+ypos_limit
                    ypos_area_min=ypos_integral_subgrid-ypos_limit

                    do jjj=j_nc-1,j_nc+1
                        do iii=i_nc-1,i_nc+1

                            ii=max(1,iii);ii=min(dim_length_nc(x_dim_nc_index),iii)
                            jj=max(1,jjj);jj=min(dim_length_nc(y_dim_nc_index),jjj)
                            !write(*,*) ii,jj,iii,jjj
                            xpos_min=max(xpos_area_min,var1d_nc(ii,lon_nc_index)-dgrid_nc(lon_nc_index)/2.)
                            xpos_max=min(xpos_area_max,var1d_nc(ii,lon_nc_index)+dgrid_nc(lon_nc_index)/2.)
                            ypos_min=max(ypos_area_min,var1d_nc(jj,lat_nc_index)-dgrid_nc(lat_nc_index)/2.)
                            ypos_max=min(ypos_area_max,var1d_nc(jj,lat_nc_index)+dgrid_nc(lat_nc_index)/2.)

                            !Determine the area intersection of the EMEP grid and an EMEP grid size centred on the integral subgrid
                            if (xpos_max.gt.xpos_min.and.ypos_max.gt.ypos_min) then
                                weighting_nc(ii,jj)=(ypos_max-ypos_min)*(xpos_max-xpos_min)/dgrid_nc(lon_nc_index)/dgrid_nc(lat_nc_index)
                            else
                                weighting_nc(ii,jj)=0.
                            endif

                            meteo_subgrid(i,j,:,ugrid_subgrid_index)=meteo_subgrid(i,j,:,ugrid_subgrid_index)+var4d_nc(ii,jj,surface_level_nc,:,ugrid_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            meteo_subgrid(i,j,:,vgrid_subgrid_index)=meteo_subgrid(i,j,:,vgrid_subgrid_index)+var4d_nc(ii,jj,surface_level_nc,:,vgrid_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            meteo_subgrid(i,j,:,FFgrid_subgrid_index)=meteo_subgrid(i,j,:,FFgrid_subgrid_index)+var4d_nc(ii,jj,surface_level_nc,:,FFgrid_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            meteo_subgrid(i,j,:,kz_subgrid_index)=meteo_subgrid(i,j,:,kz_subgrid_index)+var4d_nc(ii,jj,surface_level_nc,:,kz_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            meteo_subgrid(i,j,:,hmix_subgrid_index)=meteo_subgrid(i,j,:,hmix_subgrid_index)+var3d_nc(ii,jj,:,hmix_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            meteo_subgrid(i,j,:,u10_subgrid_index)=meteo_subgrid(i,j,:,u10_subgrid_index)+var3d_nc(ii,jj,:,u10_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            meteo_subgrid(i,j,:,v10_subgrid_index)=meteo_subgrid(i,j,:,v10_subgrid_index)+var3d_nc(ii,jj,:,v10_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            meteo_subgrid(i,j,:,FF10_subgrid_index)=meteo_subgrid(i,j,:,FF10_subgrid_index)+var3d_nc(ii,jj,:,FF10_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            meteo_subgrid(i,j,:,logz0_subgrid_index)=meteo_subgrid(i,j,:,logz0_subgrid_index)+var3d_nc(ii,jj,:,logz0_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            meteo_subgrid(i,j,:,invL_subgrid_index)=meteo_subgrid(i,j,:,invL_subgrid_index)+var3d_nc(ii,jj,:,invL_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            meteo_subgrid(i,j,:,inv_FFgrid_subgrid_index)=meteo_subgrid(i,j,:,inv_FFgrid_subgrid_index)+var4d_nc(ii,jj,surface_level_nc,:,inv_FFgrid_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            meteo_subgrid(i,j,:,inv_FF10_subgrid_index)=meteo_subgrid(i,j,:,inv_FF10_subgrid_index)+var3d_nc(ii,jj,:,inv_FF10_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            meteo_subgrid(i,j,:,ustar_subgrid_index)=meteo_subgrid(i,j,:,ustar_subgrid_index)+var3d_nc(ii,jj,:,ustar_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            meteo_subgrid(i,j,:,J_subgrid_index)=meteo_subgrid(i,j,:,J_subgrid_index)+var4d_nc(ii,jj,surface_level_nc,:,J_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            meteo_subgrid(i,j,:,t2m_subgrid_index)=meteo_subgrid(i,j,:,t2m_subgrid_index)+var3d_nc(ii,jj,:,t2m_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            !meteo_subgrid(i,j,:,precip_subgrid_index)=meteo_subgrid(i,j,:,precip_subgrid_index)+var4d_nc(ii,jj,surface_level_nc,:,precip_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            meteo_subgrid(i,j,:,precip_subgrid_index)=meteo_subgrid(i,j,:,precip_subgrid_index)+var3d_nc(ii,jj,:,precip_nc_index,allsource_index,meteo_p_loop_index)*weighting_nc(ii,jj)
                            !if (var4d_nc(ii,jj,surface_level_nc,1,precip_nc_index,allsource_index,meteo_p_loop_index).gt.0) write(*,*) var4d_nc(ii,jj,surface_level_nc,1,precip_nc_index,allsource_index,meteo_p_loop_index)
                        enddo
                    enddo
                    !if (meteo_subgrid(i,j,1,precip_subgrid_index).gt.0.) write(*,*) i,j,meteo_subgrid(i,j,:,precip_subgrid_index)

                enddo
            enddo
            !write(*,*) maxval(meteo_subgrid(:,:,t_start:t_end,precip_subgrid_index))
            !write(*,*) maxval(var4d_nc(:,:,surface_level_nc,:,precip_nc_index,allsource_index,meteo_p_loop_index))
            !stop

            if (use_alternative_meteorology_flag.or.use_alternative_z0_flag) then

                if (allocated(weighting_nc)) deallocate(weighting_nc)
                allocate (weighting_nc(dim_length_meteo_nc(x_dim_nc_index),dim_length_meteo_nc(y_dim_nc_index)))

                xpos_limit=meteo_dgrid_nc(lon_nc_index)/2.
                ypos_limit=meteo_dgrid_nc(lat_nc_index)/2.

                if (use_alternative_meteorology_flag) then
                    meteo_subgrid(:,:,:,ugrid_subgrid_index)=0.
                    meteo_subgrid(:,:,:,vgrid_subgrid_index)=0.
                    meteo_subgrid(:,:,:,FFgrid_subgrid_index)=0.
                    meteo_subgrid(:,:,:,hmix_subgrid_index)=0.
                    meteo_subgrid(:,:,:,FF10_subgrid_index)=0.
                    meteo_subgrid(:,:,:,logz0_subgrid_index)=0.
                    meteo_subgrid(:,:,:,invL_subgrid_index)=0.
                    !meteo_subgrid(:,:,:,inv_FFgrid_subgrid_index)=0.
                    !meteo_subgrid(:,:,:,inv_FF10_subgrid_index)=0.
                    meteo_subgrid(:,:,:,ustar_subgrid_index)=0.
                    meteo_subgrid(:,:,:,t2m_subgrid_index)=0.
                    meteo_subgrid(:,:,:,precip_subgrid_index)=0.
                endif

                if (use_alternative_z0_flag) then
                    meteo_subgrid(:,:,:,logz0_subgrid_index)=0.
                endif


                do j=1,integral_subgrid_dim(y_dim_index)
                    do i=1,integral_subgrid_dim(x_dim_index)

                        !Assumes it is never on the edge of the EMEP grid, not limitted
                        i_nc=crossreference_integral_to_meteo_nc_subgrid(i,j,x_dim_index)
                        j_nc=crossreference_integral_to_meteo_nc_subgrid(i,j,y_dim_index)

                        !write(*,*) i,i_nc,j,j_nc
                        xpos_integral_subgrid=meteo_nc_xproj_integral_subgrid(i,j)
                        ypos_integral_subgrid=meteo_nc_yproj_integral_subgrid(i,j)

                        xpos_area_max=xpos_integral_subgrid+xpos_limit
                        xpos_area_min=xpos_integral_subgrid-xpos_limit
                        ypos_area_max=ypos_integral_subgrid+ypos_limit
                        ypos_area_min=ypos_integral_subgrid-ypos_limit

                        do jj=j_nc-1,j_nc+1
                            do ii=i_nc-1,i_nc+1

                                xpos_min=max(xpos_area_min,meteo_var1d_nc(ii,lon_nc_index)-dgrid_nc(lon_nc_index)/2.)
                                xpos_max=min(xpos_area_max,meteo_var1d_nc(ii,lon_nc_index)+dgrid_nc(lon_nc_index)/2.)
                                ypos_min=max(ypos_area_min,meteo_var1d_nc(jj,lat_nc_index)-dgrid_nc(lat_nc_index)/2.)
                                ypos_max=min(ypos_area_max,meteo_var1d_nc(jj,lat_nc_index)+dgrid_nc(lat_nc_index)/2.)

                                !Determine the area intersection of the EMEP grid and an EMEP grid size centred on the integral subgrid
                                if (xpos_max.gt.xpos_min.and.ypos_max.gt.ypos_min) then
                                    weighting_nc(ii,jj)=(ypos_max-ypos_min)*(xpos_max-xpos_min)/dgrid_nc(lon_nc_index)/dgrid_nc(lat_nc_index)
                                else
                                    weighting_nc(ii,jj)=0.
                                endif

                                if (use_alternative_meteorology_flag) then
                                    meteo_subgrid(i,j,t_start:t_end,ugrid_subgrid_index)=meteo_subgrid(i,j,t_start:t_end,ugrid_subgrid_index)+meteo_var3d_nc(ii,jj,t_start:t_end,ugrid_nc_index)*weighting_nc(ii,jj)
                                    meteo_subgrid(i,j,t_start:t_end,vgrid_subgrid_index)=meteo_subgrid(i,j,t_start:t_end,vgrid_subgrid_index)+meteo_var3d_nc(ii,jj,t_start:t_end,vgrid_nc_index)*weighting_nc(ii,jj)
                                    meteo_subgrid(i,j,t_start:t_end,FFgrid_subgrid_index)=meteo_subgrid(i,j,t_start:t_end,FFgrid_subgrid_index)+meteo_var3d_nc(ii,jj,t_start:t_end,FFgrid_nc_index)*weighting_nc(ii,jj)
                                    meteo_subgrid(i,j,t_start:t_end,hmix_subgrid_index)=meteo_subgrid(i,j,t_start:t_end,hmix_subgrid_index)+meteo_var3d_nc(ii,jj,t_start:t_end,hmix_nc_index)*weighting_nc(ii,jj)
                                    meteo_subgrid(i,j,t_start:t_end,u10_subgrid_index)=meteo_subgrid(i,j,t_start:t_end,u10_subgrid_index)+meteo_var3d_nc(ii,jj,t_start:t_end,ugrid_nc_index)*weighting_nc(ii,jj)
                                    meteo_subgrid(i,j,t_start:t_end,v10_subgrid_index)=meteo_subgrid(i,j,t_start:t_end,v10_subgrid_index)+meteo_var3d_nc(ii,jj,t_start:t_end,vgrid_nc_index)*weighting_nc(ii,jj)
                                    meteo_subgrid(i,j,t_start:t_end,FF10_subgrid_index)=meteo_subgrid(i,j,t_start:t_end,FF10_subgrid_index)+meteo_var3d_nc(ii,jj,t_start:t_end,FF10_nc_index)*weighting_nc(ii,jj)
                                    meteo_subgrid(i,j,t_start:t_end,logz0_subgrid_index)=meteo_subgrid(i,j,t_start:t_end,logz0_subgrid_index)+meteo_var3d_nc(ii,jj,t_start:t_end,logz0_nc_index)*weighting_nc(ii,jj)
                                    meteo_subgrid(i,j,t_start:t_end,invL_subgrid_index)=meteo_subgrid(i,j,t_start:t_end,invL_subgrid_index)+meteo_var3d_nc(ii,jj,t_start:t_end,invL_nc_index)*weighting_nc(ii,jj)
                                    !meteo_subgrid(i,j,t_start:t_end,inv_FFgrid_subgrid_index)=meteo_subgrid(i,j,t_start:t_end,inv_FFgrid_subgrid_index)+meteo_var3d_nc(ii,jj,t_start:t_end,inv_FFgrid_nc_index)*weighting_nc(ii,jj)
                                    !meteo_subgrid(i,j,t_start:t_end,inv_FF10_subgrid_index)=meteo_subgrid(i,j,t_start:t_end,inv_FF10_subgrid_index)+meteo_var3d_nc(ii,jj,t_start:t_end,inv_FF10_nc_index)*weighting_nc(ii,jj)
                                    meteo_subgrid(i,j,t_start:t_end,ustar_subgrid_index)=meteo_subgrid(i,j,t_start:t_end,ustar_subgrid_index)+meteo_var3d_nc(ii,jj,t_start:t_end,ustar_nc_index)*weighting_nc(ii,jj)
                                    meteo_subgrid(i,j,t_start:t_end,t2m_subgrid_index)=meteo_subgrid(i,j,t_start:t_end,t2m_subgrid_index)+meteo_var3d_nc(ii,jj,t_start:t_end,t2m_nc_index)*weighting_nc(ii,jj)
                                    meteo_subgrid(i,j,t_start:t_end,precip_subgrid_index)=meteo_subgrid(i,j,t_start:t_end,precip_subgrid_index)+meteo_var3d_nc(ii,jj,t_start:t_end,precip_nc_index)*weighting_nc(ii,jj)
                                endif

                                if (use_alternative_z0_flag) then
                                    meteo_subgrid(i,j,t_start:t_end,logz0_subgrid_index)=meteo_subgrid(i,j,t_start:t_end,logz0_subgrid_index)+meteo_var3d_nc(ii,jj,t_start:t_end,logz0_nc_index)*weighting_nc(ii,jj)
                                endif

                            enddo
                        enddo

                    enddo
                enddo
            endif

        endif

        !Rotate wind fields if necessary and define cos and sin of the wind direction
        do j=1,integral_subgrid_dim(y_dim_index)
            do i=1,integral_subgrid_dim(x_dim_index)

                !Assumes it is never on the edge of the EMEP grid, not limitted
                i_nc=crossreference_integral_to_emep_subgrid(i,j,x_dim_index)
                j_nc=crossreference_integral_to_emep_subgrid(i,j,y_dim_index)

                !Adjust wind direction to utm projection.
                !First determine rotation for grids that are not lat lon
                !Will fail for 90 degree rotations though this should never be the case
                if (EMEP_projection_type.ne.LL_projection_index) then
                    dlatx=var2d_nc(i_nc+1,j_nc,lat_nc_index)-var2d_nc(i_nc-1,j_nc,lat_nc_index)
                    dlaty=var2d_nc(i_nc,j_nc+1,lat_nc_index)-var2d_nc(i_nc,j_nc-1,lat_nc_index)
                    angle_lcc=atan(dlatx/dlaty)
                else
                    angle_lcc=0.
                endif

                if (use_alternative_meteorology_flag) then
                    !Assumes it is never on the edge of the EMEP grid, not limitted
                    i_nc=crossreference_integral_to_meteo_nc_subgrid(i,j,x_dim_index)
                    j_nc=crossreference_integral_to_meteo_nc_subgrid(i,j,y_dim_index)

                    !Adjust wind direction to utm projection.
                    !First determine rotation for grids that are not lat lon
                    !Will fail for 90 degree rotations though this should never be the case
                    if (meteo_nc_projection_type.ne.LL_projection_index) then
                        dlatx=meteo_var2d_nc(i_nc+1,j_nc,lat_nc_index)-meteo_var2d_nc(i_nc-1,j_nc,lat_nc_index)
                        dlaty=meteo_var2d_nc(i_nc,j_nc+1,lat_nc_index)-meteo_var2d_nc(i_nc,j_nc-1,lat_nc_index)
                        angle_lcc=atan(dlatx/dlaty)
                    else
                        angle_lcc=0.
                    endif
                endif

                !Rotation from lat lon grid to UTM grid. No alternatives
                if (projection_type.eq.UTM_projection_index) then
                    angle_utm2 = atan(tan((lon_integral_subgrid(i,j)-utm_lon0)/180.*pi)*sin(lat_integral_subgrid(i,j)/180.*pi))
                elseif (projection_type.eq.LTM_projection_index) then
                    angle_utm2 = atan(tan((lon_integral_subgrid(i,j)-ltm_lon0)/180.*pi)*sin(lat_integral_subgrid(i,j)/180.*pi))
                else
                    angle_utm2 = 0.
                endif

                !Sum the angles (Triple check this is correct in regard to signs)
                angle_utm=angle_utm2+angle_lcc

                !write(*,*) i,j,angle_lcc*180./pi,angle_utm2*180./pi,angle_utm*180./pi

                !Rotate
                u_utm = meteo_subgrid(i,j,:,ugrid_subgrid_index)*cos(angle_utm)+meteo_subgrid(i,j,:,vgrid_subgrid_index)*sin(angle_utm)
                v_utm =-meteo_subgrid(i,j,:,ugrid_subgrid_index)*sin(angle_utm)+meteo_subgrid(i,j,:,vgrid_subgrid_index)*cos(angle_utm)
                meteo_subgrid(i,j,:,ugrid_subgrid_index)=u_utm
                meteo_subgrid(i,j,:,vgrid_subgrid_index)=v_utm

                !Create cos and sin's of the lowest level wind direction for efficient use in the dispersion equations
                if (wind_vectors_10m_available) then
                    ff=sqrt(meteo_subgrid(i,j,:,v10_subgrid_index)*meteo_subgrid(i,j,:,v10_subgrid_index)+meteo_subgrid(i,j,:,u10_subgrid_index)*meteo_subgrid(i,j,:,u10_subgrid_index))
                    meteo_subgrid(i,j,:,sin_subgrid_index)=meteo_subgrid(i,j,:,v10_subgrid_index)/ff
                    meteo_subgrid(i,j,:,cos_subgrid_index)=meteo_subgrid(i,j,:,u10_subgrid_index)/ff

                else

                    ff=sqrt(meteo_subgrid(i,j,:,vgrid_subgrid_index)*meteo_subgrid(i,j,:,vgrid_subgrid_index)+meteo_subgrid(i,j,:,ugrid_subgrid_index)*meteo_subgrid(i,j,:,ugrid_subgrid_index))
                    meteo_subgrid(i,j,:,sin_subgrid_index)=meteo_subgrid(i,j,:,vgrid_subgrid_index)/ff
                    meteo_subgrid(i,j,:,cos_subgrid_index)=meteo_subgrid(i,j,:,ugrid_subgrid_index)/ff
                endif

                !In case where no wind. Hopefully this never happens
                where (ff.eq.0.)
                    meteo_subgrid(i,j,:,sin_subgrid_index)=0.
                    meteo_subgrid(i,j,:,cos_subgrid_index)=1.
                endwhere

            enddo
        enddo

        !Set the last meteo data to be the current for the first value of the external time loop
        if (use_single_time_loop_flag) then
            if (t_loop.eq.start_time_loop_index) then
                last_meteo_subgrid(:,:,:)=meteo_subgrid(:,:,1,:)
            endif
        endif

        if (allocated(u_utm)) deallocate(u_utm)
        if (allocated(v_utm)) deallocate(v_utm)
        if (allocated(th)) deallocate(th)
        if (allocated(ff)) deallocate(ff)
        if (allocated(weighting_nc)) deallocate(weighting_nc)

    end subroutine uEMEP_subgrid_meteo_EMEP

end module subgrid_meteo_emep

