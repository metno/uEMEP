module read_landuse_data

    use uemep_configuration
    use uEMEP_definitions
    use crossreference_grids, only: uEMEP_crossreference_grids
    use mod_read_esri_ascii_file, only: read_esri_ascii_file, read_esri_ascii_header
    use mod_lambert_projection, only: PROJ2LL, lb2lambert2_uEMEP, LL2PS_spherical

    implicit none

    private

    public :: uEMEP_set_landuse_classes, uEMEP_read_netcdf_landuse_latlon

    integer Continuous_urban_fabric_value,Discontinuous_urban_fabric_value,Industrial_or_commercial_units_value,Road_and_rail_networks_and_associated_land_value,Port_areas_value
    integer Airports_value,Mineral_extraction_sites_value,Dump_sites_value,Construction_sites_value,Green_urban_areas_value
    integer Sport_and_leisure_facilities_value,Non_irrigated_arable_land_value,Permanently_irrigated_land_value,Rice_fields_value,Vineyards_value
    integer Fruit_trees_and_berry_plantations_value,Olive_groves_value,Pastures_value,Annual_crops_associated_with_permanent_crops_value,Complex_cultivation_patterns_value
    integer Land_principally_occupied_by_agriculture_value,Agro_forestry_areas_value,Broad_leaved_forest_value,Coniferous_forest_value,Mixed_forest_value
    integer Natural_grasslands_value,Moors_and_heathland_value,Sclerophyllous_vegetation_value,Transitional_woodland_shrub_value,Beaches_dunes_sands_value
    integer Bare_rocks_value,Sparsely_vegetated_areas_value,Burnt_areas_value,Glaciers_and_perpetual_snow_value,Inland_marshes_value
    integer Peat_bogs_value,Salt_marshes_value,Salines_value,Intertidal_flats_value,Water_courses_value,Water_bodies_value,Coastal_lagoons_value,Estuaries_value,Sea_and_ocean_value
    integer NODATA_clc_value
    integer n_corine_landuse_index
    parameter (n_corine_landuse_index=48)
    integer Corine_to_EMEP_landuse(n_corine_landuse_index)

    double precision :: scale_factor_nc, add_offset_nc

contains

    subroutine uEMEP_read_netcdf_landuse_latlon

        use uEMEP_definitions
        use netcdf

        implicit none
        integer status_nc
        integer i,j
        integer i_dim,id_nc
        character(256) var_name_nc_temp,dimname_temp
        integer var_id_nc
        integer i_landuse_index,j_landuse_index
        real delta_landuse_nc(num_dims_landuse_nc)
        integer dim_id_nc(num_dims_landuse_nc)
        logical reduce_landuse_region_flag
        real temp_lon(4),temp_lat(4),temp_x(4),temp_y(4)
        real temp_x_min,temp_x_max,temp_y_min,temp_y_max
        integer i_temp_min,i_temp_max,j_temp_min,j_temp_max
        real temp_delta(num_dims_landuse_nc)
        real correct_lon(2)
        real temp_scale
        integer i_source,i_landuse
        real buffer_delta
        logical :: exists

        !Temporary reading variables
        real, allocatable :: landuse_nc_dp(:,:)
        double precision, allocatable :: var2d_nc_dp(:,:)
        double precision, allocatable :: temp_var2d_nc_dp(:,:)

        !Functions
        !real area_weighted_extended_vectorgrid_interpolation_function

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Reading landuse data (uEMEP_read_netcdf_landuse_latlon)'
        write(unit_logfile,'(A)') '================================================================'

        !Set the filename
        pathfilename_landuse=trim(pathname_landuse)//trim(filename_landuse)

        !Test existence. If does not exist then stop
        inquire(file=trim(pathfilename_landuse),exist=exists)
        if (.not.exists) then
            write(unit_logfile,'(A,A)') ' ERROR: Netcdf file does not exist: ', trim(pathfilename_landuse)
            write(unit_logfile,'(A)') '  STOPPING'
            stop
        endif

        !Open the netcdf file for reading
        write(unit_logfile,'(2A)') ' Opening netcdf file: ',trim(pathfilename_landuse)
        status_nc = NF90_OPEN (pathfilename_landuse, nf90_nowrite, id_nc)
        if (status_nc .NE. NF90_NOERR) then
            write(unit_logfile,'(A,I)') 'ERROR opening netcdf file. Stopping: ',status_nc
            stop
        endif

        !Find the (lon,lat) dimensions of the file. Use the meteo id's as these are x and y
        do i_dim=1,num_dims_landuse_nc
            status_nc = NF90_INQ_DIMID (id_nc,dim_name_landuse_nc(i_dim),dim_id_nc(i_dim))
            status_nc = NF90_INQUIRE_DIMENSION (id_nc,dim_id_nc(i_dim),dimname_temp,dim_length_landuse_nc(i_dim))
            if (status_nc .NE. NF90_NOERR) then
                write(unit_logfile,'(A,A,A,I)') 'No dimension information available for ',trim(dim_name_landuse_nc(i_dim)),' Setting to 1 with status: ',status_nc
                dim_length_landuse_nc(i_dim)=1
            endif
        enddo

        write(unit_logfile,'(A,6I)') ' Size of landuse dimensions (lon,lat): ',dim_length_landuse_nc


        !Reduce the size of the grid to the heating emission grid size
        reduce_landuse_region_flag=.true.
        if (reduce_landuse_region_flag) then
            write(unit_logfile,'(A)') 'Reducing landuse domain for reading'
            !Determine the LL cordinates of the target grid
            !if (EMEP_projection_type.eq.LCC_projection_index) then
            !Retrieve the four corners of the target grid in lat and lon
            buffer_delta=10
            call PROJ2LL(landuse_subgrid_min(x_dim_index)-buffer_delta*landuse_subgrid_delta(x_dim_index),landuse_subgrid_min(y_dim_index)-buffer_delta*landuse_subgrid_delta(y_dim_index),temp_lon(1),temp_lat(1),projection_attributes,projection_type)
            call PROJ2LL(landuse_subgrid_max(x_dim_index)+buffer_delta*landuse_subgrid_delta(x_dim_index),landuse_subgrid_max(y_dim_index)+buffer_delta*landuse_subgrid_delta(y_dim_index),temp_lon(2),temp_lat(2),projection_attributes,projection_type)
            call PROJ2LL(landuse_subgrid_min(x_dim_index)-buffer_delta*landuse_subgrid_delta(x_dim_index),landuse_subgrid_max(y_dim_index)+buffer_delta*landuse_subgrid_delta(y_dim_index),temp_lon(3),temp_lat(3),projection_attributes,projection_type)
            call PROJ2LL(landuse_subgrid_max(x_dim_index)+buffer_delta*landuse_subgrid_delta(x_dim_index),landuse_subgrid_min(y_dim_index)-buffer_delta*landuse_subgrid_delta(y_dim_index),temp_lon(4),temp_lat(4),projection_attributes,projection_type)


            temp_x_min=1.e32;temp_y_min=1.e32
            temp_x_max=-1.e32;temp_y_max=-1.e32

            temp_x=temp_lon;temp_y=temp_lat
            do i=1,4
                write(*,*) i,temp_x(i),temp_y(i)
                if (temp_x(i).lt.temp_x_min) temp_x_min=temp_x(i)
                if (temp_y(i).lt.temp_y_min) temp_y_min=temp_y(i)
                if (temp_x(i).gt.temp_x_max) temp_x_max=temp_x(i)
                if (temp_y(i).gt.temp_y_max) temp_y_max=temp_y(i)
            enddo
            write(unit_logfile,'(A,2f12.2)') 'Min: ',temp_x_min,temp_y_min
            write(unit_logfile,'(A,2f12.2)') 'Max: ',temp_x_max,temp_y_max


            !Read the lon and lat values to get the delta and size. Put in temporary array

            !Allocate the temporary arrays for lat,lon and population
            if (.not.allocated(temp_var2d_nc_dp)) allocate (temp_var2d_nc_dp(max(dim_length_landuse_nc(x_dim_nc_index),dim_length_landuse_nc(y_dim_nc_index)),num_dims_landuse_nc)) !x and y

            dim_start_landuse_nc=1
            do i=1,num_dims_landuse_nc
                !Identify the variable name and ID in the nc file and read it
                var_name_nc_temp=dim_name_landuse_nc(i)
                status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                if (status_nc .EQ. NF90_NOERR) then
                    status_nc = nf90_get_att(id_nc, var_id_nc, 'scale_factor', scale_factor_nc)
                    if (status_nc /= nf90_noerr) scale_factor_nc = 1.0d0
                    status_nc = nf90_get_att(id_nc, var_id_nc, 'add_offset', add_offset_nc)
                    if (status_nc /= nf90_noerr) add_offset_nc = 0.0d0
                    status_nc = NF90_GET_VAR(id_nc, var_id_nc, temp_var2d_nc_dp(1:dim_length_landuse_nc(i),i), start=[dim_start_landuse_nc(i)], count=[dim_length_landuse_nc(i)])
                    temp_var2d_nc_dp(1:dim_length_landuse_nc(i),i) = temp_var2d_nc_dp(1:dim_length_landuse_nc(i),i) * scale_factor_nc + add_offset_nc
                else
                    write(unit_logfile,'(A,A,A,I)') 'No information available for ',trim(var_name_nc_temp),' Status: ',status_nc
                endif
            enddo

            delta_landuse_nc=temp_var2d_nc_dp(2,:)-temp_var2d_nc_dp(1,:)
            write(unit_logfile,'(A,2f12.6)') 'Landuse grid delta (degrees): ',delta_landuse_nc

            !write(*,*) temp_var1d_nc_dp
            temp_delta(1)=delta_landuse_nc(1)
            temp_delta(2)=delta_landuse_nc(2)

            !write(*,*) temp_delta
            !Find grid position of the max and min coordinates and add2 grids*EMEP_grid_interpolation_size
            i_temp_min=1+floor((temp_x_min-temp_var2d_nc_dp(1,1))/temp_delta(1)+0.5)
            i_temp_max=1+floor((temp_x_max-temp_var2d_nc_dp(1,1))/temp_delta(1)+0.5)
            j_temp_min=1+floor((temp_y_min-temp_var2d_nc_dp(1,2))/temp_delta(2)+0.5)
            j_temp_max=1+floor((temp_y_max-temp_var2d_nc_dp(1,2))/temp_delta(2)+0.5)
            !write(unit_logfile,'(A,2I)') ' Reading EMEP i grids: ',i_temp_min,i_temp_max
            !write(unit_logfile,'(A,2I)') ' Reading EMEP j grids: ',j_temp_min,j_temp_max
            !Increase the region by 5 grids to be certain
            i_temp_min=max(1,i_temp_min-10)
            i_temp_max=min(dim_length_landuse_nc(x_dim_nc_index),i_temp_max+10)
            j_temp_min=max(1,j_temp_min-10)
            j_temp_max=min(dim_length_landuse_nc(y_dim_nc_index),j_temp_max+10)
            dim_length_landuse_nc(x_dim_nc_index)=i_temp_max-i_temp_min+1
            dim_length_landuse_nc(y_dim_nc_index)=j_temp_max-j_temp_min+1
            dim_start_landuse_nc(x_dim_nc_index)=i_temp_min
            dim_start_landuse_nc(y_dim_nc_index)=j_temp_min
            write(unit_logfile,'(A,3I)') ' Reading landuse i grids: ',i_temp_min,i_temp_max,dim_length_landuse_nc(x_dim_nc_index)
            write(unit_logfile,'(A,3I)') ' Reading landuse j grids: ',j_temp_min,j_temp_max,dim_length_landuse_nc(y_dim_nc_index)
            write(unit_logfile,'(A,2f12.2)') ' Reading landuse lon grids (min,max): ',temp_var2d_nc_dp(i_temp_min,x_dim_nc_index),temp_var2d_nc_dp(i_temp_max,x_dim_nc_index)
            write(unit_logfile,'(A,2f12.2)') ' Reading landuse lat grids (min,max): ',temp_var2d_nc_dp(j_temp_min,y_dim_nc_index),temp_var2d_nc_dp(j_temp_max,y_dim_nc_index)
            !endif

        endif

        if (i_temp_min.ge.i_temp_max.or.j_temp_min.ge.j_temp_max) then
            !No population data available
            write(unit_logfile,'(A)') ' WARNING: No landuse data available in this region. Setting to 0'
            landuse_subgrid(:,:,clc_index)=0

        else

            if (.not.allocated(landuse_nc_dp)) allocate (landuse_nc_dp(dim_length_landuse_nc(x_dim_nc_index),dim_length_landuse_nc(y_dim_nc_index))) !Lat and lon
            if (.not.allocated(var2d_nc_dp)) allocate (var2d_nc_dp(max(dim_length_landuse_nc(x_dim_nc_index),dim_length_landuse_nc(y_dim_nc_index)),num_dims_landuse_nc)) !x and y

            !Read the lon and lat values to get the delta
            do i=1,num_dims_landuse_nc
                !Identify the variable name and ID in the nc file and read it
                var_name_nc_temp=dim_name_landuse_nc(i)
                status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                if (status_nc .EQ. NF90_NOERR) then
                    status_nc = nf90_get_att(id_nc, var_id_nc, 'scale_factor', scale_factor_nc)
                    if (status_nc /= nf90_noerr) scale_factor_nc = 1.0d0
                    status_nc = nf90_get_att(id_nc, var_id_nc, 'add_offset', add_offset_nc)
                    if (status_nc /= nf90_noerr) add_offset_nc = 0.0d0
                    status_nc = NF90_GET_VAR(id_nc, var_id_nc, var2d_nc_dp(1:dim_length_landuse_nc(i),i), start=[dim_start_landuse_nc(i)], count=[dim_length_landuse_nc(i)])
                    var2d_nc_dp(1:dim_length_landuse_nc(i),i) = var2d_nc_dp(1:dim_length_landuse_nc(i),i) * scale_factor_nc + add_offset_nc
                else
                    write(unit_logfile,'(A,A,A,I)') 'No information available for ',trim(var_name_nc_temp),' Status: ',status_nc
                endif
            enddo
            delta_landuse_nc=var2d_nc_dp(2,:)-var2d_nc_dp(1,:)
            write(unit_logfile,'(A,2f12.6)') 'Landuse grid delta (degrees): ',delta_landuse_nc
            !write(*,*) var2d_nc_dp(1,1),var2d_nc_dp(dim_length_population_nc(x_dim_nc_index),1)
            !write(*,*) var2d_nc_dp(1,2),var2d_nc_dp(dim_length_population_nc(y_dim_nc_index),2)

            !Read the landuse data
            i=1
            !Uses the population_nc_index as index, =1, but not logical
            !Identify the variable name and ID in the nc file and read it
            var_name_nc_temp=var_name_landuse_nc(i)
            status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
            if (status_nc .EQ. NF90_NOERR) then
                status_nc = nf90_get_att(id_nc, var_id_nc, 'scale_factor', scale_factor_nc)
                if (status_nc /= nf90_noerr) scale_factor_nc = 1.0d0
                status_nc = nf90_get_att(id_nc, var_id_nc, 'add_offset', add_offset_nc)
                if (status_nc /= nf90_noerr) add_offset_nc = 0.0d0
                status_nc = NF90_GET_VAR(id_nc, var_id_nc, landuse_nc_dp(:,:), &
                    start=[dim_start_landuse_nc(x_dim_nc_index),dim_start_landuse_nc(y_dim_nc_index)], &
                    count=[dim_length_landuse_nc(x_dim_nc_index),dim_length_landuse_nc(y_dim_nc_index)])
                landuse_nc_dp(:,:) = landuse_nc_dp(:,:) * scale_factor_nc + add_offset_nc
                write(unit_logfile,'(2a,2f12.2)') 'Landuse variable min and max: ',trim(var_name_nc_temp),minval(landuse_nc_dp(:,:)),maxval(landuse_nc_dp(:,:))
            else
                write(unit_logfile,'(A,A,A,I)') 'No information available for ',trim(var_name_nc_temp),' Status: ',status_nc
            endif
            !enddo
            !write(*,*) 'Finished reading landuse data'

            !Loop through the landuse data and put it in the landuse grid
            !Converting from lat lon to the subgrid coordinates and then finding the nearest neighbour
            landuse_subgrid(:,:,clc_index)=0
            where (landuse_nc_dp.lt.0) landuse_nc_dp=0.
            write(unit_logfile,'(2a,2f12.2)') 'Landuse min and max: ',trim(var_name_nc_temp),minval(landuse_nc_dp(:,:)),maxval(landuse_nc_dp(:,:))
            !stop


            do j=1,landuse_subgrid_dim(y_dim_nc_index)
                do i=1,landuse_subgrid_dim(x_dim_nc_index)

                    !Project the centre position to lat lon
                    call PROJ2LL(x_landuse_subgrid(i,j),y_landuse_subgrid(i,j),temp_lon(1),temp_lat(1),projection_attributes,projection_type)
                    !Project both sides to get the delta
                    call PROJ2LL(x_landuse_subgrid(i,j)-landuse_subgrid_delta(x_dim_index)/2.,y_landuse_subgrid(i,j),temp_lon(2),temp_lat(2),projection_attributes,projection_type)
                    call PROJ2LL(x_landuse_subgrid(i,j)+landuse_subgrid_delta(x_dim_index)/2.,y_landuse_subgrid(i,j),temp_lon(3),temp_lat(3),projection_attributes,projection_type)
                    temp_delta(x_dim_index)=temp_lon(3)-temp_lon(2)
                    call PROJ2LL(x_landuse_subgrid(i,j),y_landuse_subgrid(i,j)-landuse_subgrid_delta(y_dim_index)/2.,temp_lon(2),temp_lat(2),projection_attributes,projection_type)
                    call PROJ2LL(x_landuse_subgrid(i,j),y_landuse_subgrid(i,j)+landuse_subgrid_delta(y_dim_index)/2.,temp_lon(3),temp_lat(3),projection_attributes,projection_type)
                    temp_delta(y_dim_index)=temp_lat(3)-temp_lat(2)

                    !Make a local correction to lon so it is essentially in the same units as lat so area averaging is correct
                    correct_lon(1)=cos(3.14159/180.*temp_lat(1))
                    correct_lon(2)=1.

                    !Take the nearest instead
                    i_landuse_index=1+floor((temp_lon(1)-var2d_nc_dp(1,x_dim_nc_index))/delta_landuse_nc(1)+0.5)
                    j_landuse_index=1+floor((temp_lat(1)-var2d_nc_dp(1,y_dim_nc_index))/delta_landuse_nc(2)+0.5)
                    !write(*,*) i,j,temp_lon(1),temp_lat(1)


                    landuse_subgrid(i,j,clc_index)=landuse_nc_dp(i_landuse_index,j_landuse_index)

                    !Place the clc landuse in the EMEP landuse
                    if (landuse_subgrid(i,j,clc_index).gt.0) then
                        landuse_subgrid(i,j,Corine_to_EMEP_landuse(landuse_subgrid(i,j,clc_index)))=1
                    else
                        landuse_subgrid(i,j,Corine_to_EMEP_landuse(NODATA_clc_value))=1
                    endif

                    !Do the interpolation on the same grid then scale afterwards. Equivalent to interpolating density then rescaling with grid size
                    !landuse_subgrid(i,j,clc_index)=area_weighted_extended_vectorgrid_interpolation_function( &
                    !    real(var2d_nc_dp(1:dim_length_landuse_nc(x_dim_nc_index),x_dim_nc_index))*correct_lon(1),real(var2d_nc_dp(1:dim_length_landuse_nc(y_dim_nc_index),y_dim_nc_index)) &
                    !    ,landuse_nc_dp(:,:,landuse_nc_index),dim_length_landuse_nc(x_dim_nc_index),dim_length_landuse_nc(y_dim_nc_index) &
                    !    ,delta_landuse_nc*correct_lon,temp_lon(1)*correct_lon(1),temp_lat(1),delta_landuse_nc*correct_lon)

                    !temp_scale=(temp_delta(1)*correct_lon(1)*temp_delta(2)*correct_lon(2))/(delta_landuse_nc(1)*correct_lon(1)*delta_landuse_nc(2)*correct_lon(2))
                    !write(*,*) temp_scale
                    !landuse_subgrid(i,j,clc_index)=landuse_subgrid(i,j,clc_index)*temp_scale

                    if (isnan(landuse_subgrid(i,j,clc_index))) then
                        write(*,*) 'Stopping, nan in landuse_subgrid'
                        write(*,*) temp_scale,correct_lon,delta_landuse_nc,temp_delta,temp_lon
                        stop
                    endif
                    if (landuse_subgrid(i,j,clc_index).lt.0.) then
                        write(*,*) 'Stopping, negative value in landuse_subgrid'
                        write(*,*) temp_scale,correct_lon,delta_landuse_nc,temp_delta,temp_lon
                        stop
                    endif

                enddo
            enddo

            write(unit_logfile,'(A,2f12.2)') 'Max and min landuse in read domain: ',maxval(landuse_nc_dp(:,:)),minval(landuse_nc_dp(:,:))
            write(unit_logfile,'(A,2f12.2)') 'Max and min landuse in subgrid domain: ',maxval(landuse_subgrid(:,:,clc_index)),minval(landuse_subgrid(:,:,clc_index))

            if (use_landuse_as_proxy) then
                !Place the landuse as a proxy emission with the appropriate weights
                !loop through all sources and landuses
                do i_source=1,n_source_index
                    !If source is to be downscaled and at least one landuse is selected then calculate the proxy emission weighting
                    !write(*,*) i_source,calculate_source(i_source),sum(landuse_proxy_weighting(i_source,:))
                    if (calculate_source(i_source).and.sum(landuse_proxy_weighting(i_source,:)).gt.0) then
                        proxy_emission_subgrid(:,:,i_source,:)=0.
                        do i_landuse=1,n_clc_landuse_index
                            !write(*,*) i_source,i_landuse,landuse_proxy_weighting(i_source,i_landuse)
                            if (landuse_proxy_weighting(i_source,i_landuse).gt.0) then
                                write(unit_logfile,'(A,i4,A,A,a,i4)') 'Distributing landuse index ',i_landuse,' to uEMEP sector "',trim(source_file_str(i_source)),'" and GNFR sector ',uEMEP_to_EMEP_sector(i_source)

                                do j=1,emission_subgrid_dim(y_dim_nc_index,i_source)
                                    do i=1,emission_subgrid_dim(x_dim_nc_index,i_source)

                                        i_landuse_index=crossreference_emission_to_landuse_subgrid(i,j,x_dim_index,i_source)
                                        j_landuse_index=crossreference_emission_to_landuse_subgrid(i,j,y_dim_index,i_source)

                                        if (int(landuse_subgrid(i_landuse_index,j_landuse_index,clc_index)).eq.i_landuse) then
                                            proxy_emission_subgrid(i,j,i_source,:)=proxy_emission_subgrid(i,j,i_source,:)+landuse_proxy_weighting(i_source,i_landuse)
                                            !write(*,'(6i,2f12.2)') i,j,i_landuse_index,j_landuse_index,i_source,i_landuse,landuse_proxy_weighting(i_source,i_landuse),proxy_emission_subgrid(i,j,i_source,1)
                                        endif

                                        !If there is no data (0 or greater than the maximum number of landuse categories) then distribute emissions evenly on the EMEP grid
                                        if (int(landuse_subgrid(i_landuse_index,j_landuse_index,clc_index)).gt.n_clc_landuse_index.or.int(landuse_subgrid(i_landuse_index,j_landuse_index,clc_index)).lt.1) then
                                            proxy_emission_subgrid(i,j,i_source,:)=1.
                                        endif

                                    enddo
                                enddo
                            endif
                        enddo
                    endif
                enddo
            endif

        endif !No landuse available

        if (allocated(landuse_nc_dp)) deallocate (landuse_nc_dp)
        if (allocated(var2d_nc_dp)) deallocate (var2d_nc_dp)
        if (allocated(temp_var2d_nc_dp)) deallocate (temp_var2d_nc_dp)

    end subroutine uEMEP_read_netcdf_landuse_latlon

    subroutine uEMEP_set_landuse_classes

        use uEMEP_definitions

        implicit none


        !test
        !landuse_proxy_weighting(heating_index,Continuous_urban_fabric_value)=1.
        !landuse_proxy_weighting(heating_index,Discontinuous_urban_fabric_value)=0.5

        Continuous_urban_fabric_value=1
        Discontinuous_urban_fabric_value=2
        Industrial_or_commercial_units_value=3
        Road_and_rail_networks_and_associated_land_value=4
        Port_areas_value=5
        Airports_value=6
        Mineral_extraction_sites_value=7
        Dump_sites_value=8
        Construction_sites_value=9
        Green_urban_areas_value=10
        Sport_and_leisure_facilities_value=11
        Non_irrigated_arable_land_value=12
        Permanently_irrigated_land_value=13
        Rice_fields_value=14
        Vineyards_value=15
        Fruit_trees_and_berry_plantations_value=16
        Olive_groves_value=17
        Pastures_value=18
        Annual_crops_associated_with_permanent_crops_value=19
        Complex_cultivation_patterns_value=20
        Land_principally_occupied_by_agriculture_value=21
        Agro_forestry_areas_value=22
        Broad_leaved_forest_value=23
        Coniferous_forest_value=24
        Mixed_forest_value=25
        Natural_grasslands_value=26
        Moors_and_heathland_value=27
        Sclerophyllous_vegetation_value=28
        Transitional_woodland_shrub_value=29
        Beaches_dunes_sands_value=30
        Bare_rocks_value=31
        Sparsely_vegetated_areas_value=32
        Burnt_areas_value=33
        Glaciers_and_perpetual_snow_value=34
        Inland_marshes_value=35
        Peat_bogs_value=36
        Salt_marshes_value=37
        Salines_value=38
        Intertidal_flats_value=39
        Water_courses_value=40
        Water_bodies_value=41
        Coastal_lagoons_value=42
        Estuaries_value=43
        Sea_and_ocean_value=44
        NODATA_clc_value=48

        Corine_to_EMEP_landuse(Continuous_urban_fabric_value)=urban_index
        Corine_to_EMEP_landuse(Discontinuous_urban_fabric_value)=urban_index
        Corine_to_EMEP_landuse(Industrial_or_commercial_units_value)=urban_index
        Corine_to_EMEP_landuse(Road_and_rail_networks_and_associated_land_value)=urban_index
        Corine_to_EMEP_landuse(Port_areas_value)=urban_index
        Corine_to_EMEP_landuse(Airports_value)=urban_index
        Corine_to_EMEP_landuse(Mineral_extraction_sites_value)=urban_index
        Corine_to_EMEP_landuse(Dump_sites_value)=urban_index
        Corine_to_EMEP_landuse(Construction_sites_value)=urban_index
        Corine_to_EMEP_landuse(Green_urban_areas_value)=grass_index
        Corine_to_EMEP_landuse(Sport_and_leisure_facilities_value)=urban_index
        Corine_to_EMEP_landuse(Non_irrigated_arable_land_value)=grass_index
        Corine_to_EMEP_landuse(Permanently_irrigated_land_value)=grass_index
        Corine_to_EMEP_landuse(Rice_fields_value)=wetlands_index
        Corine_to_EMEP_landuse(Vineyards_value)=med_crop_index
        Corine_to_EMEP_landuse(Fruit_trees_and_berry_plantations_value)=med_crop_index
        Corine_to_EMEP_landuse(Olive_groves_value)=med_crop_index
        Corine_to_EMEP_landuse(Pastures_value)=grass_index
        Corine_to_EMEP_landuse(Annual_crops_associated_with_permanent_crops_value)=temp_crop_index
        Corine_to_EMEP_landuse(Complex_cultivation_patterns_value)=temp_crop_index
        Corine_to_EMEP_landuse(Land_principally_occupied_by_agriculture_value)=temp_crop_index
        Corine_to_EMEP_landuse(Agro_forestry_areas_value)=temp_decid_index
        Corine_to_EMEP_landuse(Broad_leaved_forest_value)=temp_decid_index
        Corine_to_EMEP_landuse(Coniferous_forest_value)=med_needle_index
        Corine_to_EMEP_landuse(Mixed_forest_value)=med_broadleaf_index
        Corine_to_EMEP_landuse(Natural_grasslands_value)=grass_index
        Corine_to_EMEP_landuse(Moors_and_heathland_value)=moorland_index
        Corine_to_EMEP_landuse(Sclerophyllous_vegetation_value)=medscrub_index
        Corine_to_EMEP_landuse(Transitional_woodland_shrub_value)=moorland_index
        Corine_to_EMEP_landuse(Beaches_dunes_sands_value)=desert_index
        Corine_to_EMEP_landuse(Bare_rocks_value)=urban_index
        Corine_to_EMEP_landuse(Sparsely_vegetated_areas_value)=medscrub_index
        Corine_to_EMEP_landuse(Burnt_areas_value)=medscrub_index
        Corine_to_EMEP_landuse(Glaciers_and_perpetual_snow_value)=ice_index
        Corine_to_EMEP_landuse(Inland_marshes_value)=wetlands_index
        Corine_to_EMEP_landuse(Peat_bogs_value)=wetlands_index
        Corine_to_EMEP_landuse(Salt_marshes_value)=wetlands_index
        Corine_to_EMEP_landuse(Salines_value)=wetlands_index
        Corine_to_EMEP_landuse(Intertidal_flats_value)=wetlands_index
        Corine_to_EMEP_landuse(Water_courses_value)=wetlands_index
        Corine_to_EMEP_landuse(Water_bodies_value)=water_index
        Corine_to_EMEP_landuse(Coastal_lagoons_value)=water_index
        Corine_to_EMEP_landuse(Estuaries_value)=water_index
        Corine_to_EMEP_landuse(Sea_and_ocean_value)=water_index
        Corine_to_EMEP_landuse(NODATA_clc_value)=grid_index


    end subroutine uEMEP_set_landuse_classes

end module read_landuse_data

