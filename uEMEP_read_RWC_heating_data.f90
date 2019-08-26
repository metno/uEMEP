!uEMEP_read_RWC_heating_data.f90
!Reads in MetVed data in SSB format at 250 m
!Reads in HDD cdf file in same file
!Allocates emissions to 250 m grid
    
    subroutine uEMEP_read_RWC_heating_data
    
    use uEMEP_definitions
    use netcdf

    implicit none
    
    logical exists
    logical nxtdat_flag
    integer source_index
    character(256) header_str(5)
    integer count,count_grid
    integer unit_in
    integer RWC_pm25_index,RWC_pm10_index,RWC_nox_index
    parameter (RWC_pm25_index=1,RWC_pm10_index=2,RWC_nox_index=3)
    integer RWC_HDD11_index,RWC_HDD15_index
    parameter (RWC_HDD11_index=1,RWC_HDD15_index=2)
    integer*8 ssb_id
    real x_ssb,y_ssb,lon_ssb,lat_ssb
    integer i_ssb_index,j_ssb_index
    integer :: threshold_index=0
    real :: f_easting=2.e6
    real :: ssb_dx=250.,ssb_dy=250.
    integer RWC_compound_index
    real sum_RWC_grid_emission(3)
    integer :: subsource_index=1
    real emission_scaling(3)
    logical :: read_file_with_nox_and_kommune_number=.true.
    integer :: n_region_max,n_region_back_max
    parameter (n_region_max=1000,n_region_back_max=10000)
    integer :: region_scaling_id(n_region_max)=0
    real total_emissions(n_compound_nc_index)
    integer region_scaling_id_back_index(n_region_back_max)
    real, allocatable :: region_heating_scaling(:,:)
    integer n_region,k,k_index
    
    
    integer i_pollutant
    
    
    if (filename_heating(RWC_heating_index).eq.'') then
        write(unit_logfile,'(A)') 'WARNING: No RWC heating data file available. Will not use RWC data'
        return
    endif
    
    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading RWC heating data  (uEMEP_read_RWC_heating_data)'
	write(unit_logfile,'(A)') '================================================================'
    
    !threshold_index=RWC_HDD15_index
    if (HDD_threshold_value.eq.15) then
        threshold_index=RWC_HDD15_index
    elseif (HDD_threshold_value.eq.11) then
        threshold_index=RWC_HDD11_index
    else
        write(unit_logfile,'(A,f12.1)') 'HDD_threshold_value is not valid. Stopping. ',HDD_threshold_value
        stop
    endif
     
    subsource_index=1
    
    !Emission scaling for nox compared to pm25.
    emission_scaling=1.
    emission_scaling(RWC_nox_index)=emission_factor(nox_index,heating_index,subsource_index)/emission_factor(pm25_index,heating_index,subsource_index)
    write(unit_logfile,'(A,f12.4,2es12.2)') 'NOx/PM2.5 emission ratio =',emission_scaling(RWC_nox_index),emission_factor(nox_index,heating_index,subsource_index),emission_factor(pm25_index,heating_index,subsource_index)
    
    source_index=heating_index
    n_subsource(source_index)=1

    proxy_emission_subgrid(:,:,source_index,:)=0.

    !Read in the data in the first g_loop and t_loop
    if (g_loop.eq.1) then
        pathfilename_heating(RWC_heating_index)=trim(pathname_heating(RWC_heating_index))//trim(filename_heating(RWC_heating_index))
 
        !Test existence of the heating filename. If does not exist then stop
        inquire(file=trim(pathfilename_heating(RWC_heating_index)),exist=exists)
        if (.not.exists) then
            write(unit_logfile,'(A,A)') ' ERROR: RWC file does not exist: ', trim(pathfilename_heating(RWC_heating_index))
            stop
        endif

   
       !Open the file for reading
        unit_in=20
        open(unit_in,file=pathfilename_heating(RWC_heating_index),access='sequential',status='old',readonly)  
        write(unit_logfile,'(a)') ' Opening RWC heating file '//trim(pathfilename_heating(RWC_heating_index))
    
        rewind(unit_in)
        !read (unit_in,'(a)') header_str(1)
       ! read (unit_in,'(a)') header_str(2)
        !read (unit_in,'(a)') header_str(3)
        !write(*,*) header_str(1:3)
    
        call NXTDAT(unit_in,nxtdat_flag)
        !Read number of grids
        read(unit_in,'(i)') n_RWC_grids
        write(unit_logfile,'(A,i)') 'Number of RWC grids =',n_RWC_grids
    
        !Allocate the arrays in the first g_loop and t_loop
        allocate (RWC_grid_emission(n_RWC_grids,3))
        allocate (RWC_grid_HDD(n_RWC_grids,2))
        allocate (RWC_grid_id(n_RWC_grids))
        allocate (RWC_region_id(n_RWC_grids))
    
        !Read header      SSBID    PM25_2016    PM10_2016       HDD_11       HDD_15
        read(unit_in,*) header_str
       ! write(unit_logfile,'(6A24)') 'Headers: ',trim(header_str(1)),trim(header_str(2)),trim(header_str(3)),trim(header_str(4)),trim(header_str(5))

        count=0
        RWC_grid_emission=0.
        if (read_file_with_nox_and_kommune_number) then
        !Read new version
        do while(.not.eof(unit_in))
            count=count+1
            read(unit_in,*) RWC_grid_id(count),RWC_grid_emission(count,RWC_pm25_index),RWC_grid_emission(count,RWC_pm10_index),RWC_grid_emission(count,RWC_nox_index),RWC_grid_HDD(count,RWC_HDD11_index),RWC_grid_HDD(count,RWC_HDD15_index),RWC_region_id(count)
            !write(*,'(i,5es,i)') RWC_grid_id(count),RWC_grid_emission(count,RWC_pm25_index),RWC_grid_emission(count,RWC_pm10_index),RWC_grid_emission(count,RWC_nox_index),RWC_grid_HDD(count,RWC_HDD11_index),RWC_grid_HDD(count,RWC_HDD15_index),RWC_region_id(count)
        enddo
        else
        !Read old version
        do while(.not.eof(unit_in))
            count=count+1
            !read(unit_in,'(i,4es)') RWC_grid_id(count),RWC_grid_val(count,1:4)
            read(unit_in,*) RWC_grid_id(count),RWC_grid_emission(count,RWC_pm25_index),RWC_grid_emission(count,RWC_pm10_index),RWC_grid_HDD(count,RWC_HDD11_index),RWC_grid_HDD(count,RWC_HDD15_index)
            !write(*,'(2i,4es)') count,RWC_grid_id(count),RWC_grid_val(count,1:4)
            RWC_grid_emission(count,RWC_nox_index)=RWC_grid_emission(count,RWC_pm25_index)*emission_scaling(RWC_nox_index)
            RWC_region_id(count)=1
        enddo        
        endif
    
        if (count.ne.n_RWC_grids) then
            write(unit_logfile,'(A,2i)') 'ERROR: Total number of RWC grids in file is not the same as given. Stopping: ',n_RWC_grids,count
            stop
        endif
        
        write(unit_logfile,'(A,i)') 'Total number of RWC grids read =',count
        write(unit_logfile,'(A,3f12.2)') 'Total emissions PM2.5, PM10 and NOX (tonne/year) =',sum(RWC_grid_emission(:,RWC_pm25_index))/1.e6,sum(RWC_grid_emission(:,RWC_pm10_index))/1.e6,sum(RWC_grid_emission(:,RWC_nox_index))/1.e6
   
    endif
    
    !Read in regional scaling data if available
    pathfilename_region_heating_scaling=trim(inpath_region_heating_scaling)//trim(infile_region_heating_scaling)

    !Test existence of the filename. If does not exist then use default
    inquire(file=trim(pathfilename_region_heating_scaling),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' WARNING: Regional heating scaling data file does not exist and no scaling will be applied: ', trim(pathfilename_region_heating_scaling)
        n_region=1
        region_scaling_id_back_index=n_region
        allocate (region_heating_scaling(n_region,n_compound_nc_index))
        region_heating_scaling=1.
    else
    
        unit_in=20
        open(unit_in,file=pathfilename_region_heating_scaling,access='sequential',status='old',readonly)  
        write(unit_logfile,'(a)') ' Opening regional scaling heating file: '//trim(pathfilename_region_heating_scaling)
    
        !Skip over header lines starting with *
        rewind(unit_in)
        call NXTDAT(unit_in,nxtdat_flag)
       
        !Read the data
        read(unit_in,*,ERR=10) n_region
        write(unit_logfile,'(a,i)') ' Number of regions read: ',n_region
        
        if (allocated(region_heating_scaling)) deallocate (region_heating_scaling)
        allocate (region_heating_scaling(n_region,n_compound_nc_index))
        
        call NXTDAT(unit_in,nxtdat_flag)
        do k=1,n_region          
            read(unit_in,*,ERR=10) &
                k_index,region_scaling_id(k), &
                total_emissions(pm25_index),total_emissions(pm10_index),total_emissions(nox_index), &
                region_heating_scaling(k,pm25_nc_index),region_heating_scaling(k,pm10_nc_index),region_heating_scaling(k,nox_nc_index)
            !write(*,'(2i,3es10.2,3f10.3)')  &
            !    k_index,region_scaling_id(k), &
            !    total_emissions(pm25_index),total_emissions(pm10_index),total_emissions(nox_index), &
            !    region_heating_scaling(k,pm25_nc_index),region_heating_scaling(k,pm10_nc_index),region_heating_scaling(k,nox_nc_index)
            
            !Set a back referencing value for speed. Note the region_scaling_id cannot be larger than the hardcoded dimensions
            region_scaling_id_back_index(region_scaling_id(k))=k
            
        enddo
         
        close(unit_in,status='keep')

    endif

    
    !Put the data into the subgrid for all g_loop calls
    count_grid=0
    sum_RWC_grid_emission=0
    
    do count=1,n_RWC_grids

        ssb_id=RWC_grid_id(count)
        x_ssb=ssb_id/10000000-f_easting+ssb_dx/2.
        y_ssb=mod(ssb_id,10000000)+ssb_dy/2.
        
        !Convert to EMEP coordinates
        if (save_emissions_for_EMEP(heating_index)) then
            call UTM2LL(utm_zone,y_ssb,x_ssb,lat_ssb,lon_ssb)
            if (projection_type.eq.LL_projection_index) then
                x_ssb=lon_ssb
                y_ssb=lat_ssb
            elseif (projection_type.eq.LCC_projection_index) then
                call lb2lambert2_uEMEP(x_ssb,y_ssb,lon_ssb,lat_ssb,EMEP_projection_attributes)
            endif
        endif
                
        k=region_scaling_id_back_index(RWC_region_id(count))
        if (k.lt.1.or.k.gt.n_region) then
            write(unit_logfile,'(A,3i)') ' ERROR: Region index out of bounds, stopping: ',count,RWC_region_id(count),k
            stop
        endif
            
        
        !Find the grid index it belongs to
        i_ssb_index=1+floor((x_ssb-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
        j_ssb_index=1+floor((y_ssb-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))

        if (i_ssb_index.ge.1.and.i_ssb_index.le.emission_subgrid_dim(x_dim_index,source_index) &
            .and.j_ssb_index.ge.1.and.j_ssb_index.le.emission_subgrid_dim(y_dim_index,source_index)) then

            !write(*,*) x_ssb,y_ssb,emission_subgrid_delta(x_dim_index,source_index),i_ssb_index,j_ssb_index
            !Set the proxy emssion subgrid. This will be multiplied by the hdd later in the read_time_profiles routine
            do i_pollutant=1,n_pollutant_loop
                
                if (pollutant_loop_index(i_pollutant).eq.pm10_nc_index) then
                    RWC_compound_index=RWC_pm10_index
                    proxy_emission_subgrid(i_ssb_index,j_ssb_index,source_index,i_pollutant)=proxy_emission_subgrid(i_ssb_index,j_ssb_index,source_index,i_pollutant) &
                    +RWC_grid_emission(count,RWC_compound_index)/RWC_grid_HDD(count,threshold_index)*region_heating_scaling(k,pm10_nc_index)
                    sum_RWC_grid_emission(RWC_compound_index)=sum_RWC_grid_emission(RWC_compound_index)+RWC_grid_emission(count,RWC_compound_index)
                elseif (pollutant_loop_index(i_pollutant).eq.pm25_nc_index) then
                    RWC_compound_index=RWC_pm25_index
                    proxy_emission_subgrid(i_ssb_index,j_ssb_index,source_index,i_pollutant)=proxy_emission_subgrid(i_ssb_index,j_ssb_index,source_index,i_pollutant) &
                    +RWC_grid_emission(count,RWC_compound_index)/RWC_grid_HDD(count,threshold_index)*region_heating_scaling(k,pm25_nc_index)
                    sum_RWC_grid_emission(RWC_compound_index)=sum_RWC_grid_emission(RWC_compound_index)+RWC_grid_emission(count,RWC_compound_index)
                elseif (pollutant_loop_index(i_pollutant).eq.nox_nc_index) then
                    RWC_compound_index=RWC_nox_index
                    proxy_emission_subgrid(i_ssb_index,j_ssb_index,source_index,i_pollutant)=proxy_emission_subgrid(i_ssb_index,j_ssb_index,source_index,i_pollutant) &
                    +RWC_grid_emission(count,RWC_compound_index)/RWC_grid_HDD(count,threshold_index)*region_heating_scaling(k,nox_nc_index)
                    sum_RWC_grid_emission(RWC_compound_index)=sum_RWC_grid_emission(RWC_compound_index)+RWC_grid_emission(count,RWC_compound_index)
                endif
            
            enddo
           
            !count_subgrid(i_ssb_index,j_ssb_index)=count_subgrid(i_ssb_index,j_ssb_index)+1
            count_grid=count_grid+1
            !write(*,*) count,proxy_emission_subgrid(i_ssb_index,j_ssb_index,source_index,subsource_index)
        endif
                    
             
    enddo
    
    write(unit_logfile,'(A,i)') 'Number of RWC grid placements in emission subgrid =',count_grid
    
    do i_pollutant=1,n_pollutant_loop
        if (pollutant_loop_index(i_pollutant).eq.pm10_nc_index) then
            RWC_compound_index=RWC_pm10_index
            write(unit_logfile,'(A,f12.2)') 'Total subgrid RWC emissions for '//trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//' (tonne/year) =',sum_RWC_grid_emission(RWC_compound_index)/1.e6
        elseif (pollutant_loop_index(i_pollutant).eq.pm25_nc_index) then
            RWC_compound_index=RWC_pm25_index
            write(unit_logfile,'(A,f12.2)') 'Total subgrid RWC emissions for '//trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//' (tonne/year) =',sum_RWC_grid_emission(RWC_compound_index)/1.e6
        elseif (pollutant_loop_index(i_pollutant).eq.nox_nc_index) then
            RWC_compound_index=RWC_nox_index
            write(unit_logfile,'(A,f12.2)') 'Total subgrid RWC emissions for '//trim(var_name_nc(conc_nc_index,pollutant_loop_index(i_pollutant),allsource_index))//' (tonne/year) =',sum_RWC_grid_emission(RWC_compound_index)/1.e6
        endif
    enddo
             
    if (allocated(region_heating_scaling)) deallocate(region_heating_scaling)
    
    return
10  write(unit_logfile,'(A)') 'ERROR reading scaling heating file'
    stop 

    end subroutine uEMEP_read_RWC_heating_data