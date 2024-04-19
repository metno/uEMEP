module read_ssb_data

    use mod_lambert_projection, only: LL2LAEA, PROJ2LL
    use mod_area_interpolation, only: area_weighted_extended_vectorgrid_interpolation_function

    implicit none
    private

    public :: uEMEP_read_netcdf_population, uEMEP_read_SSB_data, uEMEP_read_netcdf_population_latlon

contains

!uEMEP_read_SSB_data.f90
    
    
    subroutine uEMEP_read_SSB_data
 
    use uEMEP_definitions
    
    implicit none
    
    integer i,j,k
    character(256) temp_name
    character(256) temp_str,temp_str1,temp_str2
    integer unit_in
    integer exists
    integer count,index_val
    integer temp_int
    integer*8 ssb_id
    real dwe_todw,dwe_mult
    real pop_tot,emp_tot
    integer i_ssb_index,j_ssb_index
    integer source_index,subsource_index
    integer t
    integer, allocatable :: count_subgrid(:,:)
    real, allocatable :: temp1_subgrid(:,:),temp2_subgrid(:,:),temp3_subgrid(:,:)
   
    real x_ssb,y_ssb
    real :: f_easting=2.e6
    integer SSB_file_index
    real :: ssb_dx=250.,ssb_dy=250.
    real heating_proxy
    integer :: use_region=0

    character(256) region_number_str
    integer n_search
    parameter (n_search=5)
    character(16) search_str(n_search)
    real search_delta(n_search)
    integer temp_search
    
    data search_str /'1000m','500m','250m','100m','50m'/
    data search_delta /1000.,500.,250.,100.,50./

    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading SSB data  (uEMEP_read_SSB_data)'
	write(unit_logfile,'(A)') '================================================================'
    
    source_index=heating_index
    n_subsource(source_index)=1
    t=1

    !Initialise the use_grid array to false if population is to be used for the auto subgridding
    if (use_population_positions_for_auto_subgrid_flag) then
        use_subgrid=.false.
    endif
    
    !If dwellings are read then allocate the emission heating arrays. Otherwise allocate the population arrays
    if (SSB_data_type.eq.dwelling_index) then
        proxy_emission_subgrid(:,:,source_index,:)=0.
        allocate (count_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index)))
        allocate (temp1_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index)))
        allocate (temp2_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index)))
        allocate (temp3_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index)))
    else
        allocate (count_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index)))
        allocate (temp1_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index)))
        allocate (temp2_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index)))
        allocate (temp3_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index)))
    endif
    
    count_subgrid=0
        
    SSB_file_index=SSB_data_type
    
    if (SSB_data_type.eq.dwelling_index) then
        pathfilename_heating(SSB_file_index)=trim(pathname_heating(SSB_file_index))//trim(filename_heating(SSB_file_index))
 
        !Test existence of the heating filename. If does not exist then stop
        inquire(file=trim(pathfilename_heating(SSB_file_index)),exist=exists)
        if (.not.exists) then
            write(unit_logfile,'(A,A)') ' ERROR: SSB file does not exist: ', trim(pathfilename_heating(SSB_file_index))
            stop
        endif
        
        temp_name=pathfilename_heating(SSB_file_index)
        
    elseif (.not.use_region_select_and_mask_flag) then
        pathfilename_population(SSB_file_index)=trim(pathname_population(SSB_file_index))//trim(filename_population(SSB_file_index))
        
        !Test existence of the heating filename. If does not exist then stop
        inquire(file=trim(pathfilename_population(SSB_file_index)),exist=exists)
        if (.not.exists) then
            write(unit_logfile,'(A,A)') ' ERROR: SSB file does not exist: ', trim(pathfilename_population(SSB_file_index))
            stop
        endif
        
        temp_name=pathfilename_population(SSB_file_index)
       
    elseif (use_region_select_and_mask_flag.and.SSB_data_type.eq.population_index) then
        
        region_number_str=''
        write(region_number_str,*) region_index
        region_number_str=trim(region_number_str)//'_'
        pathfilename_population(SSB_file_index)=trim(pathname_population(SSB_file_index))//trim(adjustl(region_number_str))//trim(filename_population(SSB_file_index))

        !Test existence of the heating filename. If does not exist then use default
        inquire(file=trim(pathfilename_population(SSB_file_index)),exist=exists)
        if (.not.exists) then
            write(unit_logfile,'(A,A)') ' ERROR: SSB file does not exist: ', trim(pathfilename_population(SSB_file_index))
            stop
        endif
        
        !Search file name to define the grid size
        ssb_dx=0.;ssb_dy=0.
        do k=1,n_search
            temp_search=index(filename_population(SSB_file_index),trim(adjustl(search_str(k))))
            if (temp_search.ne.0) then
                ssb_dx=search_delta(k)
                ssb_dy=search_delta(k)
                write(unit_logfile,'(i,A)') temp_search,' Reading municipality population data with resolution '//trim(adjustl(search_str(k)))
                limit_population_delta=search_delta(k)
            endif
        enddo
    
        if (ssb_dx.eq.0) then
            write(unit_logfile,'(A)') 'Cannot find a valid SSB grid size. Stopping. '//trim(filename_population(SSB_file_index))
            stop
        endif
        
        temp_name=pathfilename_population(SSB_file_index)
    
    endif
    

    
    !Open the file for reading
    unit_in=20
    open(unit_in,file=temp_name,access='sequential',status='old',readonly)  
    write(unit_logfile,'(a)') ' Opening SSB file '//trim(temp_name)
    
    rewind(unit_in)

    subsource_index=1
    
    !Read header SSBID0250M;dwe_todw;dwe_det;dwe_2dw;dwe_row;dwe_mult;dwe_com;dwe_oth;dwe_area
    read(unit_in,'(A)') temp_str
    write(unit_logfile,'(A)') 'Header: '//trim(temp_str)
    !read(unit_in,'(A)') temp_str
    !write(*,*) trim(temp_str)
    count=0
    do while(.not.eof(unit_in))
        ssb_id=0;dwe_todw=0;dwe_mult=0;pop_tot=0;emp_tot=0
        if (SSB_data_type.eq.dwelling_index) then
            
            !Read in file string    
            read(unit_in,'(A)') temp_str
            !Extract the ssb id for the coordinates
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) ssb_id
            !Extract the total number of dwellings
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) dwe_todw

            !Skip over some values not to be used
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) temp_int
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) temp_int
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) temp_int
        
            !Extract the multiple dwellings number
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) dwe_mult
            
        endif

        if (SSB_data_type.eq.population_index) then

            !Read in file string    
            read(unit_in,'(A)') temp_str
            !write(*,*) trim(temp_str)
            !Extract the ssb id for the coordinates
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) ssb_id
            !write(*,*) trim(temp_str)
            !index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) pop_tot
            read(temp_str,*) pop_tot
            !write (*,*) ssb_id,pop_tot,index_val
        endif
    
        if (SSB_data_type.eq.establishment_index) then

            !Read in file string    
            read(unit_in,'(A)') temp_str
            !write(*,*) trim(temp_str)
            !Extract the ssb id for the coordinates
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) ssb_id
            !write(*,*) trim(temp_str)
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) temp_int
            read(temp_str,*) pop_tot
            !write (*,*) ssb_id,pop_tot,index_val
        endif

        if (SSB_data_type.eq.kindergaten_index.or.SSB_data_type.eq.school_index) then

            !Read in file string    
            read(unit_in,'(A)') temp_str
            !write(*,'(a)') trim(temp_str)
            !Extract the ssb id for the coordinates
            index_val=index(temp_str,',',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) x_ssb
            index_val=index(temp_str,',',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) y_ssb
            index_val=index(temp_str,',',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) temp_str2
            index_val=index(temp_str,',',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) temp_int
            !write(*,'(a)') trim(temp_str)
            index_val=index(temp_str,',',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) pop_tot
            !read(temp_str,*) pop_tot
            !write (*,*) x_ssb,y_ssb,pop_tot
        endif

        if (SSB_data_type.eq.home_index) then

            !Read in file string    
            read(unit_in,'(A)') temp_str
            !write(*,'(a)') trim(temp_str)
            !Extract the ssb id for the coordinates
            index_val=index(temp_str,',',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) temp_str2
            index_val=index(temp_str,',',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) y_ssb
            index_val=index(temp_str,',',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) x_ssb
            !index_val=index(temp_str,',',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) pop_tot
            read(temp_str,*) pop_tot
            !write (*,*) trim(temp_str2),y_32,x_32,pop_tot
            
            !Convert from UTM32 to 33
            !call UTM2LL(utm_zone-1,y_32,x_32,lat_32,lon_32)
            !write(*,*) lat_32,lon_32
            !call LL2UTM(1,utm_zone,lat_32,lon_32,y_ssb,x_ssb)
            !write(*,*) y_ssb,x_ssb
        endif

        count=count+1
        !if (mod(count,100000).eq.0) write(*,*) count,ssb_id,dwe_todw,dwe_mult,pop_tot
        
        if  (dwe_todw.gt.0.or.pop_tot.gt.0) then
            
            !Convert id to grid centre coordinates that are already in UTM33 for SSB data
            if (SSB_data_type.eq.dwelling_index.or.SSB_data_type.eq.establishment_index.or.SSB_data_type.eq.population_index) then
                x_ssb=ssb_id/10000000-f_easting+ssb_dx/2.
                y_ssb=mod(ssb_id,10000000)+ssb_dy/2.
            endif

            !Convert lat lon to utm coords
            !call LL2UTM(1,utm_zone,ddlatitude,ddlongitude,y_ship,x_ship)
        
            !Add to heating emission proxy subgrid       
            if (SSB_data_type.eq.dwelling_index) then
                
                !Find the grid index it belongs to
                i_ssb_index=1+floor((x_ssb-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
                j_ssb_index=1+floor((y_ssb-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))

                if (i_ssb_index.ge.1.and.i_ssb_index.le.emission_subgrid_dim(x_dim_index,source_index) &
                    .and.j_ssb_index.ge.1.and.j_ssb_index.le.emission_subgrid_dim(y_dim_index,source_index)) then

                    !write(*,*) x_ssb,y_ssb,emission_subgrid_delta(x_dim_index,source_index),i_ssb_index,j_ssb_index

                    !Reduce the number of dwellings when they are in a multiple dwelling by factor of 5. i.e. the proxy is reduced in blocks with the assumption that only 1 in 5 use their wood heater
                    heating_proxy=dwe_todw
                    heating_proxy=max(0.,dwe_todw-dwe_mult)+dwe_mult/5.
                    proxy_emission_subgrid(i_ssb_index,j_ssb_index,source_index,:)=proxy_emission_subgrid(i_ssb_index,j_ssb_index,source_index,:)+heating_proxy
                    count_subgrid(i_ssb_index,j_ssb_index)=count_subgrid(i_ssb_index,j_ssb_index)+1
                    !write(*,*) count,proxy_emission_subgrid(i_ssb_index,j_ssb_index,source_index,subsource_index)
                endif
                    
            else
                
                !Find the grid index it belongs to in the population grid
                i_ssb_index=1+floor((x_ssb-population_subgrid_min(x_dim_index))/population_subgrid_delta(x_dim_index))
                j_ssb_index=1+floor((y_ssb-population_subgrid_min(y_dim_index))/population_subgrid_delta(y_dim_index))

                if (i_ssb_index.ge.1.and.i_ssb_index.le.population_subgrid_dim(x_dim_index) &
                    .and.j_ssb_index.ge.1.and.j_ssb_index.le.population_subgrid_dim(y_dim_index).and.pop_tot.gt.0) then

                    population_subgrid(i_ssb_index,j_ssb_index,SSB_data_type)=population_subgrid(i_ssb_index,j_ssb_index,SSB_data_type)+pop_tot
                    count_subgrid(i_ssb_index,j_ssb_index)=count_subgrid(i_ssb_index,j_ssb_index)+1
                    !write(*,*) count,proxy_emission_subgrid(i_ssb_index,j_ssb_index,source_index,subsource_index)
                    
                endif

                if (use_population_positions_for_auto_subgrid_flag) then
                    !Cover the grids when target grids are smaller than population grids
                    if (SSB_data_type.eq.population_index) then
                        use_region=floor(population_subgrid_delta(x_dim_index)/subgrid_delta(x_dim_index)/2.)
                    endif
                    !Find the grid index it belongs to in the target grid
                    i_ssb_index=1+floor((x_ssb-subgrid_min(x_dim_index))/subgrid_delta(x_dim_index))
                    j_ssb_index=1+floor((y_ssb-subgrid_min(y_dim_index))/subgrid_delta(y_dim_index))
                    if (i_ssb_index-use_region.ge.1.and.i_ssb_index+use_region.le.subgrid_dim(x_dim_index) &
                        .and.j_ssb_index-use_region.ge.1.and.j_ssb_index+use_region.le.subgrid_dim(y_dim_index).and.pop_tot.gt.0) then
                         use_subgrid(i_ssb_index-use_region:i_ssb_index+use_region,j_ssb_index-use_region:j_ssb_index+use_region,:)=.true.
                    endif
                        
                endif

            endif
            
        endif
             
    enddo
    
    if (SSB_data_type.eq.dwelling_index) then
        write(unit_logfile,'(A,I)') 'Dwelling counts = ',count
        write(unit_logfile,'(A,es12.3)') 'Total dwellings = ',sum(proxy_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index,1))
        write(unit_logfile,'(A,I,a,i,a)') 'Number of grid placements = ',sum(count_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index))),' of ',emission_subgrid_dim(x_dim_index,source_index)*emission_subgrid_dim(y_dim_index,source_index),' grids'
    else
        write(unit_logfile,'(A,I)') 'Population type index = ',SSB_data_type
        write(unit_logfile,'(A,I)') 'Population counts = ',count
        write(unit_logfile,'(A,es12.3)') 'Total population = ',sum(population_subgrid(:,:,SSB_data_type))
        write(unit_logfile,'(A,I,a,i,a)') 'Number of grid placements = ',sum(count_subgrid),' of ',subgrid_dim(x_dim_index)*subgrid_dim(y_dim_index),' grids'
    endif
    
    close(unit_in)
    
    
    deallocate (count_subgrid)

    !Find the number of subgrids to be used
    if (use_population_positions_for_auto_subgrid_flag.and.SSB_data_type.ne.dwelling_index) then
        count=0
        do j=1,subgrid_dim(y_dim_index)
        do i=1,subgrid_dim(x_dim_index)
            if (use_subgrid(i,j,allsource_index)) count=count+1
        enddo
        enddo
        write(unit_logfile,'(a,i,a,i)') ' Using population for subgrids. Number of subgrids to be calculated based on population = ', count,' of ',subgrid_dim(y_dim_index)*subgrid_dim(x_dim_index)
    endif
    
   
    deallocate (temp1_subgrid,temp2_subgrid,temp3_subgrid)
    
    end subroutine uEMEP_read_SSB_data
    

    !Read in the population data in netcdf format
    !This is particularly used for reading in the global population dataset
    !This is not used any longer, replaced by uEMEP_read_netcdf_population_latlon
    !Attempt to read in local but still not working
    subroutine uEMEP_read_netcdf_population
 
    use uEMEP_definitions
    use netcdf
    
    implicit none
    integer status_nc,exists
    integer i_split,j_split,n_delta_split
    integer i,j
    integer i_dim,id_nc
    character(256) var_name_nc_temp,dimname_temp
    integer var_id_nc,var_id_nc_temp
    real x_ssb,y_ssb
    integer i_ssb_index,j_ssb_index
    real delta_pop_nc(num_dims_population_nc)
    integer dim_id_nc(num_dims_population_nc)
    integer dim_length_population_nc(num_dims_population_nc)
    real y_pop,x_pop
    integer source_index

    
    !Temporary reading rvariables
    real, allocatable :: population_nc_dp(:,:,:)
    double precision, allocatable :: var2d_nc_dp(:,:)
        
    !If the data type is dwelling then this means it is for heating
    if (SSB_data_type.eq.dwelling_index) source_index=heating_index
    
        !Set the filename
        pathfilename_population(SSB_data_type)=trim(pathname_population(SSB_data_type))//trim(filename_population(SSB_data_type))
     
        !Test existence. If does not exist then stop
        inquire(file=trim(pathfilename_population(SSB_data_type)),exist=exists)
        if (.not.exists) then
            write(unit_logfile,'(A,A)') ' ERROR: Netcdf file does not exist: ', trim(pathfilename_population(SSB_data_type))
            write(unit_logfile,'(A)') '  STOPPING'
            stop
        endif

        !Open the netcdf file for reading
        write(unit_logfile,'(2A)') ' Opening netcdf file: ',trim(pathfilename_population(SSB_data_type))
        status_nc = NF90_OPEN (pathfilename_population(SSB_data_type), nf90_nowrite, id_nc)
        if (status_nc .NE. NF90_NOERR) then
            write(unit_logfile,'(A,I)') 'ERROR opening netcdf file. Stopping: ',status_nc
            stop
        endif
        
        !Find the projection. If no projection then in lat lon coordinates
        !status_nc = NF90_INQ_VARID (id_nc,'projection_lambert',var_id_nc)
        !status_nc = NF90_INQ_VARID (id_nc,'mollweide',var_id_nc_temp)
        !if (status_nc.eq.NF90_NOERR) then
        !    population_nc_projection_type=mollweide_projection_index
        !    var_id_nc=var_id_nc_temp
        !endif
        status_nc = NF90_INQ_VARID (id_nc,'transverse_mercator',var_id_nc_temp)
        if (status_nc.eq.NF90_NOERR) then
            population_nc_projection_type=UTM_projection_index
            var_id_nc=var_id_nc_temp
        endif
        status_nc = NF90_INQ_VARID (id_nc,'projection_utm',var_id_nc_temp)
        if (status_nc.eq.NF90_NOERR) then
            population_nc_projection_type=UTM_projection_index
             var_id_nc=var_id_nc_temp
        endif
       status_nc = NF90_INQ_VARID (id_nc,'projection_lambert',var_id_nc_temp)
        if (status_nc.eq.NF90_NOERR) then
            population_nc_projection_type=LCC_projection_index
             var_id_nc=var_id_nc_temp
        endif
       status_nc = NF90_INQ_VARID (id_nc,'projection_ETRS89_LAEA',var_id_nc_temp)
        if (status_nc.eq.NF90_NOERR) then
            population_nc_projection_type=LAEA_projection_index
             var_id_nc=var_id_nc_temp
        endif
       
        if (population_nc_projection_type.ne.LL_projection_index) then
            !If there is a projection then read in the attributes. All these are doubles
            !status_nc = nf90_inquire_variable(id_nc, var_id_nc, natts = numAtts_projection)
                status_nc = nf90_get_att(id_nc, var_id_nc, 'Central_Meridian', population_nc_projection_attributes(1))
                status_nc = nf90_get_att(id_nc, var_id_nc, 'False_Easting', population_nc_projection_attributes(2))
                status_nc = nf90_get_att(id_nc, var_id_nc, 'False_Northing', population_nc_projection_attributes(3))
                status_nc = nf90_get_att(id_nc, var_id_nc, 'longitude_of_prime_meridian', population_nc_projection_attributes(4))
                status_nc = nf90_get_att(id_nc, var_id_nc, 'semi_major_axis', population_nc_projection_attributes(5))
                status_nc = nf90_get_att(id_nc, var_id_nc, 'inverse_flattening', population_nc_projection_attributes(6))
                        
                write(unit_logfile,'(A,4i,6f12.2)') 'Reading projection (index, params) ',population_nc_projection_type,population_nc_projection_attributes(1:6)

        endif
 
        !Even if projection is Mollweide still use the lat lon positions
        !population_nc_projection_type=LL_projection_index             

        !Find the (x,y) dimensions of the file. Use the meteo id's as these are x and y
        do i_dim=1,num_dims_population_nc
            status_nc = NF90_INQ_DIMID (id_nc,dim_name_population_nc(i_dim),dim_id_nc(i_dim))
            status_nc = NF90_INQUIRE_DIMENSION (id_nc,dim_id_nc(i_dim),dimname_temp,dim_length_population_nc(i_dim))
            if (status_nc .NE. NF90_NOERR) then
                write(unit_logfile,'(A,A,A,I)') 'No dimension information available for ',trim(dim_name_population_nc(i_dim)),' Setting to 1 with status: ',status_nc
                dim_length_population_nc(i_dim)=1
            endif
        enddo
                       
        write(unit_logfile,'(A,6I)') ' Size of population dimensions (x,y): ',dim_length_population_nc
        
        
        !Allocate the temporary arrays for lat,lon and population
        if (.not.allocated(population_nc_dp)) allocate (population_nc_dp(dim_length_population_nc(x_dim_nc_index),dim_length_population_nc(y_dim_nc_index),num_var_population_nc)) !Lat and lon
        if (.not.allocated(var2d_nc_dp)) allocate (var2d_nc_dp(max(dim_length_population_nc(x_dim_nc_index),dim_length_population_nc(y_dim_nc_index)),num_dims_population_nc)) !x and y

       !Read the x and y values to get the delta
        do i=1,num_dims_population_nc
            !Identify the variable name and ID in the nc file and read it
            var_name_nc_temp=dim_name_population_nc(i)
            status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
            if (status_nc .EQ. NF90_NOERR) then
                !status_nc = nf90_get_att(id_nc, var_id_nc, "units", unit_dim_meteo_nc(i))
                status_nc = NF90_GET_VAR (id_nc, var_id_nc, var2d_nc_dp(1:dim_length_population_nc(i),i),start=(/1/),count=(/dim_length_population_nc(i)/))
            else
                write(unit_logfile,'(A,A,A,I)') 'No information available for ',trim(var_name_nc_temp),' Status: ',status_nc
            endif            
        enddo
        
        delta_pop_nc=var2d_nc_dp(2,:)-var2d_nc_dp(1,:)
        write(unit_logfile,'(A,2f12.2)') 'Population grid delta (m): ',delta_pop_nc    
        
        do i=1,num_var_population_nc
            !Identify the variable name and ID in the nc file and read it
            var_name_nc_temp=var_name_population_nc(i)
            status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
            if (status_nc .EQ. NF90_NOERR) then
                !status_nc = nf90_get_att(id_nc, var_id_nc, "units", unit_dim_meteo_nc(i))
                status_nc = NF90_GET_VAR (id_nc, var_id_nc, population_nc_dp(:,:,i),start=(/1,1/),count=(/dim_length_population_nc(x_dim_nc_index),dim_length_population_nc(y_dim_nc_index)/))
                write(unit_logfile,'(2a,2f12.2)') 'Population variable min and max: ',trim(var_name_nc_temp),minval(population_nc_dp(:,:,i)),maxval(population_nc_dp(:,:,i))
            else
                write(unit_logfile,'(A,A,A,I)') 'No information available for ',trim(var_name_nc_temp),' Status: ',status_nc
            endif            
        enddo
        
        !Loop through the population data and put it in the population grid
        !Converting from lat lon to the subgrid coordinates and then finding the nearest neighbour

        !If population_subgrid_delta<delta_pop_nc then split the input grid into subgrids that are smaller than subgrid
        !Split any way into four even if they are the same grids for better interpolation
        n_delta_split=floor(sqrt(delta_pop_nc(x_dim_nc_index)**2+delta_pop_nc(y_dim_nc_index)**2)/sqrt(population_subgrid_delta(x_dim_index)**2+population_subgrid_delta(y_dim_index)**2)+.001)
        n_delta_split=10*max(n_delta_split,1)
        write(unit_logfile,'(A,i)') 'Population grids split into this many segments for numerical nearest neighbour interpolation: ',n_delta_split    
        
        !Already defined as 0. Do it again
        population_subgrid(:,:,SSB_data_type)=0
        
        do j=1,dim_length_population_nc(y_dim_nc_index)
        do i=1,dim_length_population_nc(x_dim_nc_index)
            if (population_nc_dp(i,j,population_nc_index).gt.0) then
                
                if (projection_type.eq.UTM_projection_index) then

                    call LL2UTM(1,utm_zone,population_nc_dp(i,j,lat_nc_index),population_nc_dp(i,j,lon_nc_index),y_pop,x_pop)
                    !write(*,*) population_nc_dp(i,j,lat_nc_index),population_nc_dp(i,j,lon_nc_index),y_pop,x_pop

                elseif (projection_type.eq.LTM_projection_index) then

                    call LL2LTM(1,ltm_lon0,population_nc_dp(i,j,lat_nc_index),population_nc_dp(i,j,lon_nc_index),y_pop,x_pop)
                
                elseif (projection_type.eq.LAEA_projection_index) then

                    call LL2LAEA(x_pop,y_pop,population_nc_dp(i,j,lon_nc_index),population_nc_dp(i,j,lat_nc_index),projection_attributes)
                    !call LAEA2LL(x_pop,y_pop,x_ssb,y_ssb,projection_attributes)
                    !write(*,*) population_nc_dp(i,j,lon_nc_index)-x_ssb,population_nc_dp(i,j,lat_nc_index)-y_ssb

                else
                    write(unit_logfile,'(A)') 'No supported projection available for population. Stopping '
                    stop
                endif
                
                !Loop through the n_delta_split
                do i_split=1,n_delta_split
                do j_split=1,n_delta_split
                    
                    x_ssb=x_pop-delta_pop_nc(x_dim_nc_index)*0.5+(i_split-0.5)*delta_pop_nc(x_dim_nc_index)/n_delta_split
                    y_ssb=y_pop-delta_pop_nc(y_dim_nc_index)*0.5+(j_split-0.5)*delta_pop_nc(y_dim_nc_index)/n_delta_split
                
                    if (SSB_data_type.eq.population_index) then
                        
                        i_ssb_index=1+floor((x_ssb-population_subgrid_min(x_dim_index))/population_subgrid_delta(x_dim_index))
                        j_ssb_index=1+floor((y_ssb-population_subgrid_min(y_dim_index))/population_subgrid_delta(y_dim_index))

                        if (i_ssb_index.ge.1.and.i_ssb_index.le.population_subgrid_dim(x_dim_index) &
                        .and.j_ssb_index.ge.1.and.j_ssb_index.le.population_subgrid_dim(y_dim_index)) then
                        
                            population_subgrid(i_ssb_index,j_ssb_index,SSB_data_type)=population_subgrid(i_ssb_index,j_ssb_index,SSB_data_type) &
                            +population_nc_dp(i,j,population_nc_index)/(n_delta_split**2)
                        
                            !write(*,'(2i,6f12.1)') i_split,j_split,x_pop,y_pop,x_ssb,y_ssb,(i_split-0.5)*delta_pop_nc(x_dim_nc_index)/n_delta_split,(j_split-0.5)*delta_pop_nc(y_dim_nc_index)/n_delta_split
                            !write(*,*) i_split,j_split,i_ssb_index,j_ssb_index
                                
                            !write(*,*) i_ssb_index,j_ssb_index,SSB_data_type,population_subgrid(i_ssb_index,j_ssb_index,SSB_data_type)
                        endif
                        
                    endif
                    
                    if (SSB_data_type.eq.dwelling_index) then
                        
                        i_ssb_index=1+floor((x_ssb-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
                        j_ssb_index=1+floor((y_ssb-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))
 
                        if (i_ssb_index.ge.1.and.i_ssb_index.le.emission_subgrid_dim(x_dim_index,source_index) &
                        .and.j_ssb_index.ge.1.and.j_ssb_index.le.emission_subgrid_dim(y_dim_index,source_index)) then

                            proxy_emission_subgrid(i_ssb_index,j_ssb_index,source_index,:)=proxy_emission_subgrid(i_ssb_index,j_ssb_index,source_index,:) &
                            +population_nc_dp(i,j,population_nc_index)/(n_delta_split**2)
                        
                        endif 

                    endif

                enddo
                enddo
               
            endif
            
        enddo
        enddo
        
        if (allocated(population_nc_dp)) deallocate (population_nc_dp)
        if (allocated(var2d_nc_dp)) deallocate (var2d_nc_dp)
    
    end subroutine uEMEP_read_netcdf_population

    !Read in the population data in netcdf format in latlon grid
    !This is particularly used for reading in the global population dataset
    subroutine uEMEP_read_netcdf_population_latlon
 
    use uEMEP_definitions
    use netcdf
    
    implicit none
    integer status_nc,exists
    integer i,j
    integer i_dim,id_nc
    character(256) var_name_nc_temp,dimname_temp
    integer var_id_nc
    real delta_pop_nc(num_dims_population_nc)
    integer dim_id_nc(num_dims_population_nc)
    integer dim_length_population_nc(num_dims_population_nc)
    integer dim_start_population_nc(num_dims_population_nc)
    integer source_index
    logical reduce_population_region_flag
    real temp_lon(4),temp_lat(4),temp_x(4),temp_y(4)
    real temp_x_min,temp_x_max,temp_y_min,temp_y_max
    integer i_temp_min,i_temp_max,j_temp_min,j_temp_max
    real temp_delta(num_dims_population_nc)
    real correct_lon(2)
    real temp_scale
    integer :: name_index=0
    
    !Temporary reading rvariables
    real, allocatable :: population_nc_dp(:,:,:)
    !double precision, allocatable :: population_nc_dp(:,:,:)
    double precision, allocatable :: var2d_nc_dp(:,:)
    double precision, allocatable :: temp_var2d_nc_dp(:,:)
    
    !If the data type is dwelling then this means it is for heating
    !Always set to this since emissions will be the largest domain, when reducing the reading domain
    !source_index=heating_index
    if (SSB_data_type.eq.dwelling_index) source_index=heating_index
    if (SSB_data_type.eq.dwelling_index) name_index=dwelling_nc_index
    if (SSB_data_type.eq.population_index) name_index=population_nc_index
    
        !Set the filename
        pathfilename_population(SSB_data_type)=trim(pathname_population(SSB_data_type))//trim(filename_population(SSB_data_type))
     
        !Test existence. If does not exist then stop
        inquire(file=trim(pathfilename_population(SSB_data_type)),exist=exists)
        if (.not.exists) then
            write(unit_logfile,'(A,A)') ' ERROR: Netcdf file does not exist: ', trim(pathfilename_population(SSB_data_type))
            write(unit_logfile,'(A)') '  STOPPING'
            stop
        endif

        !Open the netcdf file for reading
        write(unit_logfile,'(2A)') ' Opening netcdf file: ',trim(pathfilename_population(SSB_data_type))
        status_nc = NF90_OPEN (pathfilename_population(SSB_data_type), nf90_nowrite, id_nc)
        if (status_nc .NE. NF90_NOERR) then
            write(unit_logfile,'(A,I)') 'ERROR opening netcdf file. Stopping: ',status_nc
            stop
        endif
        
        !Find the (lon,lat) dimensions of the file. Use the meteo id's as these are x and y
        do i_dim=1,num_dims_population_nc
            status_nc = NF90_INQ_DIMID (id_nc,dim_name_population_nc(i_dim),dim_id_nc(i_dim))
            status_nc = NF90_INQUIRE_DIMENSION (id_nc,dim_id_nc(i_dim),dimname_temp,dim_length_population_nc(i_dim))
            if (status_nc .NE. NF90_NOERR) then
                write(unit_logfile,'(A,A,A,I)') 'No dimension information available for ',trim(dim_name_population_nc(i_dim)),' Setting to 1 with status: ',status_nc
                dim_length_population_nc(i_dim)=1
            endif
        enddo
                       
        write(unit_logfile,'(A,6I)') ' Size of population dimensions (lon,lat): ',dim_length_population_nc
        
        
        !Reduce the size of the grid to the heating emission grid size
        reduce_population_region_flag=.true.
        if (reduce_population_region_flag) then
            write(unit_logfile,'(A)') 'Reducing population domain for reading'
            !Determine the LL cordinates of the target grid
            !if (EMEP_projection_type.eq.LCC_projection_index) then
                !Retrieve the four corners of the target grid in lat and lon
            if (SSB_data_type.eq.dwelling_index) then
            call PROJ2LL(emission_subgrid_min(x_dim_index,source_index),emission_subgrid_min(y_dim_index,source_index),temp_lon(1),temp_lat(1),projection_attributes,projection_type)
            call PROJ2LL(emission_subgrid_max(x_dim_index,source_index),emission_subgrid_max(y_dim_index,source_index),temp_lon(2),temp_lat(2),projection_attributes,projection_type)
            call PROJ2LL(emission_subgrid_min(x_dim_index,source_index),emission_subgrid_max(y_dim_index,source_index),temp_lon(3),temp_lat(3),projection_attributes,projection_type)
            call PROJ2LL(emission_subgrid_max(x_dim_index,source_index),emission_subgrid_min(y_dim_index,source_index),temp_lon(4),temp_lat(4),projection_attributes,projection_type)
            endif
            if (SSB_data_type.eq.population_index) then
            call PROJ2LL(population_subgrid_min(x_dim_index),population_subgrid_min(y_dim_index),temp_lon(1),temp_lat(1),projection_attributes,projection_type)
            call PROJ2LL(population_subgrid_max(x_dim_index),population_subgrid_max(y_dim_index),temp_lon(2),temp_lat(2),projection_attributes,projection_type)
            call PROJ2LL(population_subgrid_min(x_dim_index),population_subgrid_max(y_dim_index),temp_lon(3),temp_lat(3),projection_attributes,projection_type)
            call PROJ2LL(population_subgrid_max(x_dim_index),population_subgrid_min(y_dim_index),temp_lon(4),temp_lat(4),projection_attributes,projection_type)
            endif
            
                            
                temp_x_min=1.e32;temp_y_min=1.e32
                temp_x_max=-1.e32;temp_y_max=-1.e32
                
                temp_x=temp_lon;temp_y=temp_lat    
                do i=1,4
                    !write(*,*) i,temp_x(i),temp_y(i)
                    if (temp_x(i).lt.temp_x_min) temp_x_min=temp_x(i)
                    if (temp_y(i).lt.temp_y_min) temp_y_min=temp_y(i)
                    if (temp_x(i).gt.temp_x_max) temp_x_max=temp_x(i)
                    if (temp_y(i).gt.temp_y_max) temp_y_max=temp_y(i)
                enddo
                write(unit_logfile,'(A,2f12.2)') 'Min: ',temp_x_min,temp_y_min
                write(unit_logfile,'(A,2f12.2)') 'Max: ',temp_x_max,temp_y_max
                
            
                !Read the lon and lat values to get the delta and size. Put in temporary array
                
                !Allocate the temporary arrays for lat,lon and population
                if (.not.allocated(temp_var2d_nc_dp)) allocate (temp_var2d_nc_dp(max(dim_length_population_nc(x_dim_nc_index),dim_length_population_nc(y_dim_nc_index)),num_dims_population_nc)) !x and y

                dim_start_population_nc=1
                do i=1,num_dims_population_nc
                    !Identify the variable name and ID in the nc file and read it
                    var_name_nc_temp=dim_name_population_nc(i)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc .EQ. NF90_NOERR) then
                        !status_nc = nf90_get_att(id_nc, var_id_nc, "units", unit_dim_meteo_nc(i))
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, temp_var2d_nc_dp(1:dim_length_population_nc(i),i),start=(/dim_start_population_nc(i)/),count=(/dim_length_population_nc(i)/))
                    else
                        write(unit_logfile,'(A,A,A,I)') 'No information available for ',trim(var_name_nc_temp),' Status: ',status_nc
                    endif            
                enddo
        
                delta_pop_nc=temp_var2d_nc_dp(2,:)-temp_var2d_nc_dp(1,:)
                write(unit_logfile,'(A,2f12.6)') 'Population grid delta (degrees): ',delta_pop_nc    
               
                !write(*,*) temp_var1d_nc_dp
                temp_delta(1)=delta_pop_nc(1)
                temp_delta(2)=delta_pop_nc(2)

                !write(*,*) temp_delta
                !Find grid position of the max and min coordinates and add2 grids*EMEP_grid_interpolation_size
                i_temp_min=1+floor((temp_x_min-temp_var2d_nc_dp(1,1))/temp_delta(1)+0.5)
                i_temp_max=1+floor((temp_x_max-temp_var2d_nc_dp(1,1))/temp_delta(1)+0.5)
                j_temp_min=1+floor((temp_y_min-temp_var2d_nc_dp(1,2))/temp_delta(2)+0.5)
                j_temp_max=1+floor((temp_y_max-temp_var2d_nc_dp(1,2))/temp_delta(2)+0.5)
                !write(unit_logfile,'(A,2I)') ' Reading EMEP i grids: ',i_temp_min,i_temp_max
                !write(unit_logfile,'(A,2I)') ' Reading EMEP j grids: ',j_temp_min,j_temp_max
                !Increase the region by 5 grids to be certain
                i_temp_min=max(1,i_temp_min-5)
                i_temp_max=min(dim_length_population_nc(x_dim_nc_index),i_temp_max+5)
                j_temp_min=max(1,j_temp_min-5)
                j_temp_max=min(dim_length_population_nc(y_dim_nc_index),j_temp_max+5)
                dim_length_population_nc(x_dim_nc_index)=i_temp_max-i_temp_min+1
                dim_length_population_nc(y_dim_nc_index)=j_temp_max-j_temp_min+1
                dim_start_population_nc(x_dim_nc_index)=i_temp_min
                dim_start_population_nc(y_dim_nc_index)=j_temp_min
                write(unit_logfile,'(A,3I)') ' Reading population i grids: ',i_temp_min,i_temp_max,dim_length_population_nc(x_dim_nc_index)
                write(unit_logfile,'(A,3I)') ' Reading population j grids: ',j_temp_min,j_temp_max,dim_length_population_nc(y_dim_nc_index)
                write(unit_logfile,'(A,2f12.2)') ' Reading population lon grids (min,max): ',temp_var2d_nc_dp(i_temp_min,x_dim_nc_index),temp_var2d_nc_dp(i_temp_max,x_dim_nc_index)
                write(unit_logfile,'(A,2f12.2)') ' Reading population lat grids (min,max): ',temp_var2d_nc_dp(j_temp_min,y_dim_nc_index),temp_var2d_nc_dp(j_temp_max,y_dim_nc_index)
            !endif

        endif
        
        if (i_temp_min.ge.i_temp_max.or.j_temp_min.ge.j_temp_max) then
            !No population data available
            write(unit_logfile,'(A)') ' WARNING: No population data available in this region. Setting to 0'
            proxy_emission_subgrid(:,:,source_index,:)=0.
            population_subgrid(:,:,SSB_data_type)=0
            
        else
            
        if (.not.allocated(population_nc_dp)) allocate (population_nc_dp(dim_length_population_nc(x_dim_nc_index),dim_length_population_nc(y_dim_nc_index),num_var_population_nc)) !Lat and lon
        if (.not.allocated(var2d_nc_dp)) allocate (var2d_nc_dp(max(dim_length_population_nc(x_dim_nc_index),dim_length_population_nc(y_dim_nc_index)),num_dims_population_nc)) !x and y

        !Read the lon and lat values to get the delta
        do i=1,num_dims_population_nc
            !Identify the variable name and ID in the nc file and read it
            var_name_nc_temp=dim_name_population_nc(i)
            status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
            if (status_nc .EQ. NF90_NOERR) then
                !status_nc = nf90_get_att(id_nc, var_id_nc, "units", unit_dim_meteo_nc(i))
                status_nc = NF90_GET_VAR (id_nc, var_id_nc, var2d_nc_dp(1:dim_length_population_nc(i),i),start=(/dim_start_population_nc(i)/),count=(/dim_length_population_nc(i)/))
            else
                write(unit_logfile,'(A,A,A,I)') 'No information available for ',trim(var_name_nc_temp),' Status: ',status_nc
            endif            
        enddo
        delta_pop_nc=var2d_nc_dp(2,:)-var2d_nc_dp(1,:)
        write(unit_logfile,'(A,2f12.6)') 'Population grid delta (degrees): ',delta_pop_nc    
       !write(*,*) var2d_nc_dp(1,1),var2d_nc_dp(dim_length_population_nc(x_dim_nc_index),1)
       !write(*,*) var2d_nc_dp(1,2),var2d_nc_dp(dim_length_population_nc(y_dim_nc_index),2)

        !Read the population data 
        !write(*,*) 'Reading population data'
        !do i=1,num_var_population_nc
        i=name_index
        !Uses the population_nc_index as index, =1, but not logical
            !Identify the variable name and ID in the nc file and read it
            var_name_nc_temp=var_name_population_nc(i)
            status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
            if (status_nc .EQ. NF90_NOERR) then
                !status_nc = nf90_get_att(id_nc, var_id_nc, "units", unit_dim_meteo_nc(i))
                status_nc = NF90_GET_VAR (id_nc, var_id_nc, population_nc_dp(:,:,population_nc_index),start=(/dim_start_population_nc(x_dim_nc_index),dim_start_population_nc(y_dim_nc_index)/),count=(/dim_length_population_nc(x_dim_nc_index),dim_length_population_nc(y_dim_nc_index)/))
                write(unit_logfile,'(2a,2f12.2)') 'Population variable min and max: ',trim(var_name_nc_temp),minval(population_nc_dp(:,:,population_nc_index)),maxval(population_nc_dp(:,:,population_nc_index))
            else
                write(unit_logfile,'(A,A,A,I)') 'No information available for ',trim(var_name_nc_temp),' Status: ',status_nc
            endif            
        !enddo
        !write(*,*) 'Finished reading population data'
        
        !Loop through the population data and put it in the population grid
        !Converting from lat lon to the subgrid coordinates and then finding the nearest neighbour
        !Interpolate to the population grid in lat lon coordinates
        population_subgrid(:,:,SSB_data_type)=0.
        !where (population_nc_dp.lt.0.0D00) population_nc_dp=0.0D00
        where (population_nc_dp.lt.0.0) population_nc_dp=0.0
        write(unit_logfile,'(2a,2f12.2)') 'Population min and max: ',trim(var_name_nc_temp),minval(population_nc_dp(:,:,population_nc_index)),maxval(population_nc_dp(:,:,population_nc_index))
        !stop
        
        if (SSB_data_type.eq.population_index) then

        do j=1,population_subgrid_dim(y_dim_nc_index)
        do i=1,population_subgrid_dim(x_dim_nc_index)
            
            !Project the centre position to lat lon
            call PROJ2LL(x_population_subgrid(i,j),y_population_subgrid(i,j),temp_lon(1),temp_lat(1),projection_attributes,projection_type)
            !Project both sides to get the delta
            call PROJ2LL(x_population_subgrid(i,j)-population_subgrid_delta(x_dim_index)/2.,y_population_subgrid(i,j),temp_lon(2),temp_lat(2),projection_attributes,projection_type)
            call PROJ2LL(x_population_subgrid(i,j)+population_subgrid_delta(x_dim_index)/2.,y_population_subgrid(i,j),temp_lon(3),temp_lat(3),projection_attributes,projection_type)            
            temp_delta(x_dim_index)=temp_lon(3)-temp_lon(2)
            call PROJ2LL(x_population_subgrid(i,j),y_population_subgrid(i,j)-population_subgrid_delta(y_dim_index)/2.,temp_lon(2),temp_lat(2),projection_attributes,projection_type)
            call PROJ2LL(x_population_subgrid(i,j),y_population_subgrid(i,j)+population_subgrid_delta(y_dim_index)/2.,temp_lon(3),temp_lat(3),projection_attributes,projection_type)            
            temp_delta(y_dim_index)=temp_lat(3)-temp_lat(2)
            
            !Make a local correction to lon so it is essentially in the same units as lat so area averaging is correct
            correct_lon(1)=cos(3.14159/180.*temp_lat(1))
            correct_lon(2)=1.
            !write(*,*) correct_lon
           ! population_subgrid(i,j,SSB_data_type)=area_weighted_extended_vectorgrid_interpolation_function( &
           !     real(var2d_nc_dp(1:dim_length_population_nc(x_dim_nc_index),x_dim_nc_index))*correct_lon(1),real(var2d_nc_dp(1:dim_length_population_nc(y_dim_nc_index),y_dim_nc_index)) &
           !     ,population_nc_dp(:,:,population_nc_index),dim_length_population_nc(x_dim_nc_index),dim_length_population_nc(y_dim_nc_index) &
           !     ,delta_pop_nc*correct_lon,temp_lon(1)*correct_lon(1),temp_lat(1),temp_delta*correct_lon) 
           ! write(*,*) temp_lon(1),temp_lat(1),var2d_nc_dp(int(dim_length_population_nc(x_dim_nc_index)/2),x_dim_nc_index),var2d_nc_dp(int(dim_length_population_nc(y_dim_nc_index)/2),y_dim_nc_index),population_subgrid(i,j,SSB_data_type)
        
            !Take the nearest instead for a check
            !i_ssb_index=1+floor((temp_lon(1)-var2d_nc_dp(1,x_dim_nc_index))/delta_pop_nc(1)+0.5)
            !j_ssb_index=1+floor((temp_lat(1)-var2d_nc_dp(1,y_dim_nc_index))/delta_pop_nc(2)+0.5)
            !population_subgrid(i,j,SSB_data_type)=population_nc_dp(i_ssb_index,j_ssb_index,population_nc_index)
            
            !Do the interpolation on the same grid then scale afterwards. Equivalent to interpolating density then rescaling with grid size
            population_subgrid(i,j,SSB_data_type)=area_weighted_extended_vectorgrid_interpolation_function( &
                real(var2d_nc_dp(1:dim_length_population_nc(x_dim_nc_index),x_dim_nc_index))*correct_lon(1),real(var2d_nc_dp(1:dim_length_population_nc(y_dim_nc_index),y_dim_nc_index)) &
                ,population_nc_dp(:,:,population_nc_index),dim_length_population_nc(x_dim_nc_index),dim_length_population_nc(y_dim_nc_index) &
                ,delta_pop_nc*correct_lon,temp_lon(1)*correct_lon(1),temp_lat(1),delta_pop_nc*correct_lon) 
            
            temp_scale=(temp_delta(1)*correct_lon(1)*temp_delta(2)*correct_lon(2))/(delta_pop_nc(1)*correct_lon(1)*delta_pop_nc(2)*correct_lon(2))
            !write(*,*) temp_scale
            population_subgrid(i,j,SSB_data_type)=population_subgrid(i,j,SSB_data_type)*temp_scale
        
            if (isnan(population_subgrid(i,j,SSB_data_type))) then
            write(*,*) 'Stopping, nan in population_subgrid'
            write(*,*) temp_scale,correct_lon,delta_pop_nc,temp_delta,temp_lon
            stop
            endif
            if (population_subgrid(i,j,SSB_data_type).lt.0.) then
            write(*,*) 'Stopping, negative value in population_subgrid'
            write(*,*) temp_scale,correct_lon,delta_pop_nc,temp_delta,temp_lon
            stop
            endif

        enddo
        enddo
        
                write(unit_logfile,'(A,f12.2)') 'Total population in read domain: ',sum(population_nc_dp(:,:,population_nc_index))
                write(unit_logfile,'(A,f12.2)') 'Total population in subgrid domain: ',sum(population_subgrid(:,:,SSB_data_type))

        endif

        if (SSB_data_type.eq.dwelling_index) then
        proxy_emission_subgrid(:,:,source_index,:)=0.
        do j=1,emission_subgrid_dim(y_dim_nc_index,source_index)
        do i=1,emission_subgrid_dim(x_dim_nc_index,source_index)
            !Project the centre position to lat lon
            call PROJ2LL(x_emission_subgrid(i,j,source_index),y_emission_subgrid(i,j,source_index),temp_lon(1),temp_lat(1),projection_attributes,projection_type)
            !Project both sides to get the delta
            call PROJ2LL(x_emission_subgrid(i,j,source_index)-emission_subgrid_delta(x_dim_index,source_index)/2.,y_emission_subgrid(i,j,source_index),temp_lon(2),temp_lat(2),projection_attributes,projection_type)
            call PROJ2LL(x_emission_subgrid(i,j,source_index)+emission_subgrid_delta(x_dim_index,source_index)/2.,y_emission_subgrid(i,j,source_index),temp_lon(3),temp_lat(3),projection_attributes,projection_type)            
            temp_delta(x_dim_index)=temp_lon(3)-temp_lon(2)
            call PROJ2LL(x_emission_subgrid(i,j,source_index),y_emission_subgrid(i,j,source_index)-emission_subgrid_delta(y_dim_index,source_index)/2.,temp_lon(2),temp_lat(2),projection_attributes,projection_type)
            call PROJ2LL(x_emission_subgrid(i,j,source_index),y_emission_subgrid(i,j,source_index)+emission_subgrid_delta(y_dim_index,source_index)/2.,temp_lon(3),temp_lat(3),projection_attributes,projection_type)            
            temp_delta(y_dim_index)=temp_lat(3)-temp_lat(2)
            
            !Make a local correction to lon so it is essentially in the same units as lat so area averaging is correct
            correct_lon(1)=cos(3.14159/180.*temp_lat(1))
            correct_lon(2)=1.

            !Interpolate on same grid then scale, equivalent to interpolating density and then recalculating
            proxy_emission_subgrid(i,j,source_index,:)=area_weighted_extended_vectorgrid_interpolation_function( &
                real(var2d_nc_dp(1:dim_length_population_nc(x_dim_nc_index),x_dim_nc_index))*correct_lon(1),real(var2d_nc_dp(1:dim_length_population_nc(y_dim_nc_index),y_dim_nc_index)) &
                ,population_nc_dp(:,:,population_nc_index),dim_length_population_nc(x_dim_nc_index),dim_length_population_nc(y_dim_nc_index) &
                ,delta_pop_nc*correct_lon,temp_lon(1)*correct_lon(1),temp_lat(1),delta_pop_nc*correct_lon)
            
            temp_scale=(temp_delta(1)*correct_lon(1)*temp_delta(2)*correct_lon(2))/(delta_pop_nc(1)*correct_lon(1)*delta_pop_nc(2)*correct_lon(2))
            proxy_emission_subgrid(i,j,source_index,:)=proxy_emission_subgrid(i,j,source_index,:)*temp_scale
            
            if (isnan(proxy_emission_subgrid(i,j,source_index,1))) then
            write(*,*) 'Stopping, nan in proxy_emission_subgrid'
            write(*,*) temp_scale,correct_lon,delta_pop_nc,temp_delta,temp_lon
            stop
            endif
            if (proxy_emission_subgrid(i,j,source_index,1).lt.0.) then
            write(*,*) 'Stopping, negative value in proxy_emission_subgrid'
            write(*,*) temp_scale,correct_lon,delta_pop_nc,temp_delta,temp_lon
            stop
            endif

        enddo
        enddo
        
        !Apply the power law scaling for population to reduce the distribution in higher population areas
        write(unit_logfile,'(A,f12.2)') 'Power scaling of population using: ',population_power_scale
        proxy_emission_subgrid(:,:,source_index,:)=proxy_emission_subgrid(:,:,source_index,:)**population_power_scale
    
        endif

        endif !No population available
       
        if (allocated(population_nc_dp)) deallocate (population_nc_dp)
        if (allocated(var2d_nc_dp)) deallocate (var2d_nc_dp)
        if (allocated(temp_var2d_nc_dp)) deallocate (temp_var2d_nc_dp)
    
    end subroutine uEMEP_read_netcdf_population_latlon

end module read_ssb_data

