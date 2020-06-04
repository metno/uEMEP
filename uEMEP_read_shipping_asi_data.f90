!uEMEP_read_shipping_asi_data.f90
    
    subroutine uEMEP_read_weekly_shipping_asi_data
 
    use uEMEP_definitions
    
    implicit none
    
    integer i,j,k
    character(256) temp_name
    character(2048) temp_str
    character(256) temp_str1,temp_str2
    real temp_val
    integer unit_in
    integer exists
    integer count,index_val
    integer temp_int
    real totalnoxemission,totalparticulatematteremission
    real y_ship,x_ship
    integer i_ship_index,j_ship_index
    integer source_index,subsource_index
    integer t
    integer, allocatable :: count_subgrid(:,:,:)
   
    integer i_pollutant
    
    integer a(6)
    character(256) format_temp,week_of_year_str
    double precision date_num,date_num_start
    integer week_of_year,ship_week,ship_counts
    logical nxtdat_flag
    real ship_delta_x,ship_delta_y
    real lat_ship,lon_ship
    
    double precision date_to_number
    
    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading weekly shipping asi data  (uEMEP_read_weekly_shipping_asi_data)'
	write(unit_logfile,'(A)') '================================================================'
    
    source_index=shipping_index
    n_subsource(source_index)=1
    proxy_emission_subgrid(:,:,source_index,:)=0.
    t=1

    allocate (count_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index),n_pollutant_loop))
    count_subgrid=0
    
    !Determine week of year approximately using date string. Needs proper function    
    format_temp='yyyymmdd'
    call datestr_to_date(config_date_str,format_temp,a)
    date_num=date_to_number(a,ref_year_meteo)
    a(2)=1;a(3)=1;a(4)=1;
    date_num_start=date_to_number(a,ref_year_meteo)
    week_of_year=1+int((date_num-date_num_start)/7.)
    week_of_year=max(min(week_of_year,52),1)
    write(*,*) week_of_year
    write(week_of_year_str,'(i2)') week_of_year
    write(unit_logfile,'(a)') 'Week of year: '//trim(ADJUSTL(week_of_year_str))
    
    pathfilename_ship(1)=trim(pathname_ship(1))//trim(filename_ship(1))//'_week_'//trim(ADJUSTL(week_of_year_str))//'.txt'

    !Test existence of the shipping filename. If does not exist then use default
    inquire(file=trim(pathfilename_ship(1)),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Shipping file does not exist: ', trim(pathfilename_ship(1))
        stop
    endif

    !Open the file for reading
    unit_in=20
    open(unit_in,file=pathfilename_ship(1),access='sequential',status='old',readonly)  
    write(unit_logfile,'(a)') ' Opening shipping file '//trim(pathfilename_ship(1))
    
    rewind(unit_in)

    subsource_index=1
    
    !Skip over lines starting with #
    call NXTDAT(unit_in,nxtdat_flag)
    
    !Read data
    read(unit_in,*) temp_str1,ship_delta_x
    read(unit_in,*) temp_str1,ship_delta_y
    read(unit_in,*) temp_str1,ship_week
    read(unit_in,*) temp_str1,ship_counts
    
    !Skip header
    read(unit_in,*) temp_str1
    do i=1,ship_counts
        read(unit_in,*) x_ship,y_ship,totalnoxemission,totalparticulatematteremission
 
        !Special case when saving emissions, convert to either latlon or lambert
        if (save_emissions_for_EMEP(shipping_index)) then
            call PROJ2LL(x_ship,y_ship,lon_ship,lat_ship,projection_attributes,projection_type)
            !call UTM2LL(utm_zone,y_ship,x_ship,lat_ship,lon_ship)
            if (EMEP_projection_type.eq.LL_projection_index) then
                x_ship=lon_ship
                y_ship=lat_ship       
            elseif (EMEP_projection_type.eq.LCC_projection_index) then
                call lb2lambert2_uEMEP(x_ship,y_ship,lon_ship,lat_ship,EMEP_projection_attributes)
            endif
        endif

        i_ship_index=1+floor((x_ship-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
        j_ship_index=1+floor((y_ship-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))

        !Add to subgrid 
        if (i_ship_index.ge.1.and.i_ship_index.le.emission_subgrid_dim(x_dim_index,source_index) &
            .and.j_ship_index.ge.1.and.j_ship_index.le.emission_subgrid_dim(y_dim_index,source_index)) then
            do i_pollutant=1,n_pollutant_loop
                if  (totalnoxemission.gt.0.and.pollutant_loop_index(i_pollutant).eq.nox_nc_index) then
                    proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,i_pollutant)=proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,i_pollutant)+totalnoxemission
                    count_subgrid(i_ship_index,j_ship_index,i_pollutant)=count_subgrid(i_ship_index,j_ship_index,i_pollutant)+1
                elseif (totalparticulatematteremission.gt.0.and.(pollutant_loop_index(i_pollutant).eq.pm25_nc_index.or.pollutant_loop_index(i_pollutant).eq.pm10_nc_index)) then
                    proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,i_pollutant)=proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,i_pollutant)+totalparticulatematteremission
                    count_subgrid(i_ship_index,j_ship_index,i_pollutant)=count_subgrid(i_ship_index,j_ship_index,i_pollutant)+1
                endif
            enddo
            endif
 
    enddo

    write(unit_logfile,'(A,I)') 'Shipping counts = ',ship_counts
    do i_pollutant=1,n_pollutant_loop   
    write(unit_logfile,'(A,es12.3)') 'Total emission (g/hr) '//trim(pollutant_file_str(pollutant_loop_index(i_pollutant)))//' = ',sum(proxy_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index,i_pollutant))
    enddo
        
    close(unit_in)
    
    deallocate (count_subgrid)
    
    end subroutine uEMEP_read_weekly_shipping_asi_data

    subroutine uEMEP_read_monthly_and_daily_shipping_asi_data
 
    use uEMEP_definitions
    
    implicit none
    
    integer i,j,k
    character(256) temp_name
    character(2048) temp_str
    character(256) temp_str1,temp_str2
    real temp_val
    integer unit_in
    integer exists
    integer count,index_val
    integer temp_int
    real totalnoxemission,totalparticulatematteremission
    real y_ship,x_ship
    integer i_ship_index,j_ship_index
    integer source_index,subsource_index
    integer t,tt
    integer, allocatable :: count_subgrid(:,:,:)
   
    integer i_pollutant
    
    integer a(6)
    character(256) format_temp,month_of_year_str
    double precision date_num,date_num_start
    integer month_of_year,ship_month,ship_counts
    logical nxtdat_flag
    real ship_delta_x,ship_delta_y
    real lat_ship,lon_ship
    real daily_cycle(24)
    integer i_ship_range,j_ship_range
    integer date_array(6)
    double precision date_num_temp
    
    double precision date_to_number
    
    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading monthly shipping asi data  (uEMEP_read_monthly_and_daily_shipping_asi_data)'
	write(unit_logfile,'(A)') '================================================================'
    
    source_index=shipping_index
    n_subsource(source_index)=1
    proxy_emission_subgrid(:,:,source_index,:)=0.
    daily_cycle=1.
    t=1

    allocate (count_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index),n_pollutant_loop))
    count_subgrid=0
    
    !Determine week of year approximately using date string. Needs proper function    
    format_temp='yyyymmdd'
    call datestr_to_date(config_date_str,format_temp,a)
    !date_num=date_to_number(a,ref_year_meteo)
    !a(2)=1;a(3)=1;a(4)=1;
    !date_num_start=date_to_number(a,ref_year_meteo)
    !week_of_year=1+int((date_num-date_num_start)/7.)
    month_of_year=a(2)
   ! write(*,*) month_of_year
    write(month_of_year_str,'(i0.2)') month_of_year
    write(unit_logfile,'(a)') 'Month of year: '//trim(ADJUSTL(month_of_year_str))
    
    pathfilename_ship(1)=trim(pathname_ship(1))//trim(filename_ship(1))//'_'//trim(ADJUSTL(month_of_year_str))//'.txt'

    !Test existence of the shipping filename. If does not exist then use default
    inquire(file=trim(pathfilename_ship(1)),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Shipping file does not exist: ', trim(pathfilename_ship(1))
        stop
    endif

    !Open the file for reading
    unit_in=20
    open(unit_in,file=pathfilename_ship(1),access='sequential',status='old',readonly)  
    write(unit_logfile,'(a)') ' Opening shipping file '//trim(pathfilename_ship(1))
    
    rewind(unit_in)

    subsource_index=1
    
    !Skip over lines starting with #
    call NXTDAT(unit_in,nxtdat_flag)
    
    !Read data
    read(unit_in,*) temp_str1,ship_delta_x
    read(unit_in,*) temp_str1,ship_delta_y
    read(unit_in,*) temp_str1,ship_month
    read(unit_in,*) temp_str1,ship_counts
    
    !Skip header
    read(unit_in,*) temp_str1
    do i=1,ship_counts
        read(unit_in,*) x_ship,y_ship,totalnoxemission,totalparticulatematteremission
 
        !Special case when saving emissions, convert to either latlon or lambert
        if (save_emissions_for_EMEP(shipping_index)) then
            call PROJ2LL(x_ship,y_ship,lon_ship,lat_ship,projection_attributes,projection_type)
            !call UTM2LL(utm_zone,y_ship,x_ship,lat_ship,lon_ship)
            if (projection_type.eq.LL_projection_index) then
                x_ship=lon_ship
                y_ship=lat_ship       
            elseif (projection_type.eq.LCC_projection_index) then
                call lb2lambert2_uEMEP(x_ship,y_ship,lon_ship,lat_ship,EMEP_projection_attributes)
            endif
        endif

        i_ship_index=1+floor((x_ship-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
        j_ship_index=1+floor((y_ship-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))

        !Add to subgrid 
        if (i_ship_index.ge.1.and.i_ship_index.le.emission_subgrid_dim(x_dim_index,source_index) &
            .and.j_ship_index.ge.1.and.j_ship_index.le.emission_subgrid_dim(y_dim_index,source_index)) then
            do i_pollutant=1,n_pollutant_loop
                if  (totalnoxemission.gt.0.and.pollutant_loop_index(i_pollutant).eq.nox_nc_index) then
                    proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,i_pollutant)=proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,i_pollutant)+totalnoxemission
                    count_subgrid(i_ship_index,j_ship_index,i_pollutant)=count_subgrid(i_ship_index,j_ship_index,i_pollutant)+1
                elseif (totalparticulatematteremission.gt.0.and.(pollutant_loop_index(i_pollutant).eq.pm25_nc_index.or.pollutant_loop_index(i_pollutant).eq.pm10_nc_index)) then
                    proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,i_pollutant)=proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,i_pollutant)+totalparticulatematteremission
                    count_subgrid(i_ship_index,j_ship_index,i_pollutant)=count_subgrid(i_ship_index,j_ship_index,i_pollutant)+1
                endif
            enddo
            endif
 
    enddo

    write(unit_logfile,'(A,I)') 'Shipping counts for monthly mean = ',ship_counts
    do i_pollutant=1,n_pollutant_loop   
    write(unit_logfile,'(A,es12.3)') 'Total emission (g/hr) '//trim(pollutant_file_str(pollutant_loop_index(i_pollutant)))//' = ',sum(proxy_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index,i_pollutant))
    enddo
        
    close(unit_in)
    
    !Now read in the daily cycle data for the same month and applyt it to the emission grid
    !Read in as utm33
    pathfilename_ship(2)=trim(pathname_ship(2))//trim(filename_ship(2))//'_'//trim(ADJUSTL(month_of_year_str))//'.txt'

    !Test existence of the shipping filename. If does not exist then use default
    inquire(file=trim(pathfilename_ship(2)),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Shipping file does not exist: ', trim(pathfilename_ship(2))
        stop
    endif

    count_subgrid=0
    
    !Set the emission time profile to the default of 1
    emission_time_profile_subgrid(:,:,:,source_index,:)=1.

    !Open the file for reading
    unit_in=20
    open(unit_in,file=pathfilename_ship(2),access='sequential',status='old',readonly)  
    write(unit_logfile,'(a)') ' Opening shipping file '//trim(pathfilename_ship(2))
    
    rewind(unit_in)

    subsource_index=1
    
    !Skip over lines starting with #
    call NXTDAT(unit_in,nxtdat_flag)
    
    !Read data
    read(unit_in,*) temp_str1,ship_delta_x
    read(unit_in,*) temp_str1,ship_delta_y
    read(unit_in,*) temp_str1,ship_month
    read(unit_in,*) temp_str1,ship_counts
    
    !Skip header
    read(unit_in,*) temp_str1
    do i=1,ship_counts
        !read(unit_in,'(2f16.1,24f8.2)') x_ship,y_ship,(daily_cycle(t),t=1,24)
        read(unit_in,*) x_ship,y_ship,(daily_cycle(t),t=1,24)
        !write(*,'(i,2f16.1,24f8.2)') i,x_ship,y_ship,(daily_cycle(t),t=1,24)
        !Convert to EMEP coordinates if it is to be saved. emission grids are already in the EMEP coordinate system
        if (save_emissions_for_EMEP(shipping_index)) then
            call PROJ2LL(x_ship,y_ship,lon_ship,lat_ship,projection_attributes,projection_type)
            !call UTM2LL(utm_zone,y_ship,x_ship,lat_ship,lon_ship)
            call lb2lambert2_uEMEP(x_ship,y_ship,lon_ship,lat_ship,EMEP_projection_attributes)
        endif

        i_ship_index=1+floor((x_ship-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
        j_ship_index=1+floor((y_ship-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))
        i_ship_range=floor(ship_delta_x/emission_subgrid_delta(x_dim_index,source_index)/2.)
        j_ship_range=floor(ship_delta_y/emission_subgrid_delta(x_dim_index,source_index)/2.)
        
        !Add to subgrid 
        if (i_ship_index.ge.1+i_ship_range.and.i_ship_index.le.emission_subgrid_dim(x_dim_index,source_index)-i_ship_range &
            .and.j_ship_index.ge.1+j_ship_range.and.j_ship_index.le.emission_subgrid_dim(y_dim_index,source_index)-j_ship_range) then
            !do i_pollutant=1,n_pollutant_loop
                do t=1,dim_length_nc(time_dim_nc_index)
                    date_num_temp=val_dim_nc(t,time_dim_nc_index)        
                    call number_to_date(date_num_temp,date_array,ref_year_EMEP)
                    tt=date_array(4)
                    if (tt.eq.0) tt=24
                    !write(*,'(4i,f6.2)') i,t,date_array(4),tt,daily_cycle(tt)
                    emission_time_profile_subgrid(i_ship_index-i_ship_range:i_ship_index+i_ship_range,j_ship_index-j_ship_range:j_ship_index+j_ship_range,t,source_index,:)=daily_cycle(tt)
                enddo
                count_subgrid(i_ship_index,j_ship_index,:)=count_subgrid(i_ship_index,j_ship_index,:)+1
            !enddo
        endif
 
    enddo

    write(unit_logfile,'(A,I)') 'Shipping counts for daily cycle = ',ship_counts
    write(unit_logfile,'(A,I)') 'Shipping grids found for daily cycle = ',sum(count_subgrid(:,:,1))
        
    close(unit_in)

    deallocate (count_subgrid)
    
    
    
    end subroutine uEMEP_read_monthly_and_daily_shipping_asi_data

    
    subroutine uEMEP_read_shipping_asi_data
    !Reads in the original ais raw data
    use uEMEP_definitions
    
    implicit none
    
    integer i,j,k
    character(256) temp_name
    character(2048) temp_str
    character(256) temp_str1,temp_str2
    real temp_val
    integer unit_in
    integer exists
    integer count,index_val
    integer temp_int
    real ddlatitude,ddlongitude,totalnoxemission,totalparticulatematteremission
    real y_ship,x_ship
    integer i_ship_index,j_ship_index
    integer source_index,subsource_index
    integer t
    integer, allocatable :: count_subgrid(:,:,:)
    real, allocatable :: temp1_subgrid(:,:),temp2_subgrid(:,:),temp3_subgrid(:,:)
   
    integer i_pollutant
    
    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading shipping asi data  (uEMEP_read_shipping_asi_data)'
	write(unit_logfile,'(A)') '================================================================'
    
    source_index=shipping_index
    n_subsource(source_index)=1
    proxy_emission_subgrid(:,:,source_index,:)=0.
    t=1

    allocate (count_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index),n_pollutant_loop))
    count_subgrid=0
    allocate (temp1_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index)))
    allocate (temp2_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index)))
    allocate (temp3_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index)))
    
    
    pathfilename_ship(1)=trim(pathname_ship(1))//trim(filename_ship(1))
    if (use_aggregated_shipping_emissions_flag) pathfilename_ship(1)=trim(pathname_ship(1))//'Aggregated_'//trim(filename_ship(1))

    !Test existence of the shipping filename. If does not exist then use default
    inquire(file=trim(pathfilename_ship(1)),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Shipping file does not exist: ', trim(pathfilename_ship(1))
        stop
    endif

    !Open the file for reading
    unit_in=20
    open(unit_in,file=pathfilename_ship(1),access='sequential',status='old',readonly)  
    write(unit_logfile,'(a)') ' Opening shipping file '//trim(pathfilename_ship(1))
    
    rewind(unit_in)

    subsource_index=1
    
    !Read header ddlatitude;ddlongitude;totalnoxemission;totalparticulatematteremission;fk_vessellloydstype;fk_ais_norwegianmainvesselcategory;date;time
    read(unit_in,'(A)') temp_str
    !write(*,*) trim(temp_str)
    count=0
    do while(.not.eof(unit_in))
        read(unit_in,'(A)') temp_str
        
        ddlatitude=0.;ddlongitude=0.;totalnoxemission=0.;totalparticulatematteremission=0.
        !Extract the values in the temp_str
        index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
        read(temp_str1,*) ddlatitude
        !write (*,*) ddlatitude
        index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
        read(temp_str1,*) ddlongitude
        !write (*,*) ddlongitude
        index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
        !write(*,*) index_val,trim(temp_str1),trim(temp_str)
        if (index_val.gt.1) read(temp_str1,*) totalnoxemission
        !write (*,*) totalnoxemission
        index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
        if (index_val.gt.1) read(temp_str1,*) totalparticulatematteremission
        
        !index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
        !if (index_val.gt.1) read(temp_str1,*) temp_int
        !index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
        !if (index_val.gt.1) read(temp_str1,*) temp_int
        !index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
        !if (index_val.gt.1) read(temp_str1,*) temp_str2
        !write (*,*) trim(temp_str1)
        !temp_str1=temp_str
        !if (len(temp_str1).gt.0) read(temp_str1,*) temp_str2
        !write (*,*) trim(temp_str1)
        
        !write(*,*) count,ddlatitude,ddlongitude,totalnoxemission,totalparticulatematteremission
        count=count+1
        !if (mod(count,100000).eq.0) write(*,*) count,ddlatitude,ddlongitude,totalnoxemission,totalparticulatematteremission
        
        if  (totalnoxemission.gt.0.or.totalparticulatematteremission.gt.0) then
            
        !Convert to EMEP coordinates if it is to be saved. emission grids are already in the EMEP coordinate system
        !This will not work for lat lon as it is now written but will never be called either
        if (save_emissions_for_EMEP(shipping_index)) then
            if (EMEP_projection_type.eq.LCC_projection_index) then
                call lb2lambert2_uEMEP(x_ship,y_ship,ddlongitude,ddlatitude,EMEP_projection_attributes)
            !elseif (EMEP_projection_type.eq.LL_projection_index) then
                !lon_ship=ddlongitude
                !lat_ship=ddlatitude
            endif           
        else
            !Convert lat lon to utm coords
            if (projection_type.eq.UTM_projection_index) then
                call LL2UTM(1,utm_zone,ddlatitude,ddlongitude,y_ship,x_ship)
            elseif (projection_type.eq.LAEA_projection_index) then
                call LL2LAEA(x_ship,y_ship,ddlongitude,ddlatitude,projection_attributes)
            endif
       endif        

    !Find the grid index it belongs to
        i_ship_index=1+floor((x_ship-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
        j_ship_index=1+floor((y_ship-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))
        !(x_subgrid(i,j)-subgrid_min(1))/+subgrid_delta(1)+1=i
        
        !Add to subgrid       
        if (i_ship_index.ge.1.and.i_ship_index.le.emission_subgrid_dim(x_dim_index,source_index) &
            .and.j_ship_index.ge.1.and.j_ship_index.le.emission_subgrid_dim(y_dim_index,source_index)) then
            do i_pollutant=1,n_pollutant_loop
                if  (totalnoxemission.gt.0.and.pollutant_loop_index(i_pollutant).eq.nox_nc_index) then
                    proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,i_pollutant)=proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,i_pollutant)+totalnoxemission
                    count_subgrid(i_ship_index,j_ship_index,i_pollutant)=count_subgrid(i_ship_index,j_ship_index,i_pollutant)+1
                elseif (totalparticulatematteremission.gt.0.and.(pollutant_loop_index(i_pollutant).eq.pm25_nc_index.or.pollutant_loop_index(i_pollutant).eq.pm10_nc_index)) then
                    proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,i_pollutant)=proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,i_pollutant)+totalparticulatematteremission
                    count_subgrid(i_ship_index,j_ship_index,i_pollutant)=count_subgrid(i_ship_index,j_ship_index,i_pollutant)+1
                endif
            enddo
        endif
            
        endif
        
    enddo
    write(unit_logfile,'(A,I)') 'Shipping counts = ',count
    do i_pollutant=1,n_pollutant_loop   
    write(unit_logfile,'(A,es12.3)') 'Total emission '//trim(pollutant_file_str(pollutant_loop_index(i_pollutant)))//' = ',sum(proxy_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index,i_pollutant))
    enddo
        
    close(unit_in)
    
    
    deallocate (count_subgrid)
    deallocate (temp1_subgrid,temp2_subgrid,temp3_subgrid)
    
    end subroutine uEMEP_read_shipping_asi_data
    
    !Read in the shipping data in netcdf format in latlon grid
    !This is particularly used for reading in the global shipping dataset
    subroutine uEMEP_read_netcdf_shipping_latlon
 
    use uEMEP_definitions
    use netcdf
    
    implicit none
    integer status_nc,exists
    integer i_split,j_split,n_delta_split
    integer i,j
    integer i_dim,id_nc
    character(256) var_name_nc_temp,dimname_temp
    integer var_id_nc
    real x_ssb,y_ssb
    integer i_ssb_index,j_ssb_index
    real delta_shipping_nc(num_dims_shipping_nc)
    integer dim_id_nc(num_dims_shipping_nc)
    integer dim_length_shipping_nc(num_dims_shipping_nc)
    integer dim_start_shipping_nc(num_dims_shipping_nc)
    real y_pop,x_pop
    integer source_index
    logical reduce_shipping_region_flag
    real temp_lon(4),temp_lat(4),temp_x(4),temp_y(4)
    real temp_x_min,temp_x_max,temp_y_min,temp_y_max
    integer i_temp_min,i_temp_max,j_temp_min,j_temp_max
    real temp_delta(num_dims_shipping_nc)
    real correct_lon(2)
    real temp_scale
    integer i_ship
    
    !Temporary reading rvariables
    real, allocatable :: shipping_nc_dp(:,:,:)
    double precision, allocatable :: var2d_nc_dp(:,:)
    double precision, allocatable :: temp_var2d_nc_dp(:,:)
    
    !Functions
    real area_weighted_extended_vectorgrid_interpolation_function
        
    source_index=shipping_nc_index
    
    !Set the filename
    pathfilename_ship(1)=trim(pathname_ship(1))//trim(filename_ship(1))
     
    !Test existence. If does not exist then stop
    inquire(file=trim(pathfilename_ship(1)),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Netcdf file does not exist: ', trim(pathfilename_ship(1))
        write(unit_logfile,'(A)') '  STOPPING'
        stop
    endif

    !Open the netcdf file for reading
    write(unit_logfile,'(2A)') ' Opening netcdf file: ',trim(pathfilename_ship(1))
    status_nc = NF90_OPEN (pathfilename_ship(1), nf90_nowrite, id_nc)
    if (status_nc .NE. NF90_NOERR) then
        write(unit_logfile,'(A,I)') 'ERROR opening netcdf file. Stopping: ',status_nc
        stop
    endif
        
        !Find the (lon,lat) dimensions of the file. Use the meteo id's as these are x and y
        do i_dim=1,num_dims_shipping_nc
            status_nc = NF90_INQ_DIMID (id_nc,dim_name_shipping_nc(i_dim),dim_id_nc(i_dim))
            status_nc = NF90_INQUIRE_DIMENSION (id_nc,dim_id_nc(i_dim),dimname_temp,dim_length_shipping_nc(i_dim))
            if (status_nc .NE. NF90_NOERR) then
                write(unit_logfile,'(A,A,A,I)') 'No dimension information available for ',trim(dim_name_shipping_nc(i_dim)),' Setting to 1 with status: ',status_nc
                dim_length_shipping_nc(i_dim)=1
            endif
        enddo
                       
        write(unit_logfile,'(A,6I)') ' Size of shipping dimensions (lon,lat): ',dim_length_shipping_nc
        
        
        !Reduce the size of the grid to the heating emission grid size
        reduce_shipping_region_flag=.true.
        if (reduce_shipping_region_flag) then
            write(unit_logfile,'(A)') 'Reducing shipping domain for reading'
            !Determine the LL cordinates of the target grid
            !if (EMEP_projection_type.eq.LCC_projection_index) then
                !Retrieve the four corners of the target grid in lat and lon
            call PROJ2LL(emission_subgrid_min(x_dim_index,source_index),emission_subgrid_min(y_dim_index,source_index),temp_lon(1),temp_lat(1),projection_attributes,projection_type)
            call PROJ2LL(emission_subgrid_max(x_dim_index,source_index),emission_subgrid_max(y_dim_index,source_index),temp_lon(2),temp_lat(2),projection_attributes,projection_type)
            call PROJ2LL(emission_subgrid_min(x_dim_index,source_index),emission_subgrid_max(y_dim_index,source_index),temp_lon(3),temp_lat(3),projection_attributes,projection_type)
            call PROJ2LL(emission_subgrid_max(x_dim_index,source_index),emission_subgrid_min(y_dim_index,source_index),temp_lon(4),temp_lat(4),projection_attributes,projection_type)
            
                            
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
                
                !Allocate the temporary arrays for lat,lon and shipping
                if (.not.allocated(temp_var2d_nc_dp)) allocate (temp_var2d_nc_dp(max(dim_length_shipping_nc(x_dim_nc_index),dim_length_shipping_nc(y_dim_nc_index)),num_dims_shipping_nc)) !x and y

                dim_start_shipping_nc=1
                do i=1,num_dims_shipping_nc
                    !Identify the variable name and ID in the nc file and read it
                    var_name_nc_temp=dim_name_shipping_nc(i)
                    status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
                    if (status_nc .EQ. NF90_NOERR) then
                        !status_nc = nf90_get_att(id_nc, var_id_nc, "units", unit_dim_meteo_nc(i))
                        status_nc = NF90_GET_VAR (id_nc, var_id_nc, temp_var2d_nc_dp(1:dim_length_shipping_nc(i),i),start=(/dim_start_shipping_nc(i)/),count=(/dim_length_shipping_nc(i)/))
                    else
                        write(unit_logfile,'(A,A,A,I)') 'No information available for ',trim(var_name_nc_temp),' Status: ',status_nc
                    endif            
                enddo
        
                delta_shipping_nc=temp_var2d_nc_dp(2,:)-temp_var2d_nc_dp(1,:)
                write(unit_logfile,'(A,2f12.6)') 'Shipping grid delta (degrees): ',delta_shipping_nc    
               
                !write(*,*) temp_var1d_nc_dp
                temp_delta(1)=delta_shipping_nc(1)
                temp_delta(2)=delta_shipping_nc(2)

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
                i_temp_max=min(dim_length_shipping_nc(x_dim_nc_index),i_temp_max+5)
                j_temp_min=max(1,j_temp_min-5)
                j_temp_max=min(dim_length_shipping_nc(y_dim_nc_index),j_temp_max+5)
                dim_length_shipping_nc(x_dim_nc_index)=i_temp_max-i_temp_min+1
                dim_length_shipping_nc(y_dim_nc_index)=j_temp_max-j_temp_min+1
                dim_start_shipping_nc(x_dim_nc_index)=i_temp_min
                dim_start_shipping_nc(y_dim_nc_index)=j_temp_min
                write(unit_logfile,'(A,3I)') ' Reading shipping i grids: ',i_temp_min,i_temp_max,dim_length_shipping_nc(x_dim_nc_index)
                write(unit_logfile,'(A,3I)') ' Reading shipping j grids: ',j_temp_min,j_temp_max,dim_length_shipping_nc(y_dim_nc_index)
                write(unit_logfile,'(A,2f12.2)') ' Reading shipping lon grids (min,max): ',temp_var2d_nc_dp(i_temp_min,x_dim_nc_index),temp_var2d_nc_dp(i_temp_max,x_dim_nc_index)
                write(unit_logfile,'(A,2f12.2)') ' Reading shipping lat grids (min,max): ',temp_var2d_nc_dp(j_temp_min,y_dim_nc_index),temp_var2d_nc_dp(j_temp_max,y_dim_nc_index)
            !endif

        endif
        
        if (.not.allocated(shipping_nc_dp)) allocate (shipping_nc_dp(dim_length_shipping_nc(x_dim_nc_index),dim_length_shipping_nc(y_dim_nc_index),num_var_shipping_nc)) !Lat and lon
        if (.not.allocated(var2d_nc_dp)) allocate (var2d_nc_dp(max(dim_length_shipping_nc(x_dim_nc_index),dim_length_shipping_nc(y_dim_nc_index)),num_dims_shipping_nc)) !x and y

        !Read the lon and lat values to get the delta
        do i=1,num_dims_shipping_nc
            !Identify the variable name and ID in the nc file and read it
            var_name_nc_temp=dim_name_shipping_nc(i)
            status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
            if (status_nc .EQ. NF90_NOERR) then
                !status_nc = nf90_get_att(id_nc, var_id_nc, "units", unit_dim_meteo_nc(i))
                status_nc = NF90_GET_VAR (id_nc, var_id_nc, var2d_nc_dp(1:dim_length_shipping_nc(i),i),start=(/dim_start_shipping_nc(i)/),count=(/dim_length_shipping_nc(i)/))
            else
                write(unit_logfile,'(A,A,A,I)') 'No information available for ',trim(var_name_nc_temp),' Status: ',status_nc
            endif            
        enddo
        delta_shipping_nc=var2d_nc_dp(2,:)-var2d_nc_dp(1,:)
        write(unit_logfile,'(A,2f12.6)') 'Shipping grid delta (degrees): ',delta_shipping_nc    
       !write(*,*) var2d_nc_dp(1,1),var2d_nc_dp(dim_length_shipping_nc(x_dim_nc_index),1)
       !write(*,*) var2d_nc_dp(1,2),var2d_nc_dp(dim_length_shipping_nc(y_dim_nc_index),2)

        !Read the shipping data 
        do i_ship=1,num_var_shipping_nc
            !Identify the variable name and ID in the nc file and read it
            var_name_nc_temp=var_name_shipping_nc(i_ship)
            status_nc = NF90_INQ_VARID (id_nc, trim(var_name_nc_temp), var_id_nc)
            if (status_nc .EQ. NF90_NOERR) then
                !status_nc = nf90_get_att(id_nc, var_id_nc, "units", unit_dim_meteo_nc(i))
                status_nc = NF90_GET_VAR (id_nc, var_id_nc, shipping_nc_dp(:,:,i_ship),start=(/dim_start_shipping_nc(x_dim_nc_index),dim_start_shipping_nc(y_dim_nc_index)/),count=(/dim_length_shipping_nc(x_dim_nc_index),dim_length_shipping_nc(y_dim_nc_index)/))
                write(unit_logfile,'(2a,2f12.2)') 'Shipping variable min and max: ',trim(var_name_nc_temp),minval(shipping_nc_dp(:,:,i_ship)),maxval(shipping_nc_dp(:,:,i_ship))
            else
                write(unit_logfile,'(A,A,A,I)') 'No information available for ',trim(var_name_nc_temp),' Status: ',status_nc
            endif            
        enddo
        
        !Loop through the shipping data and put it in the shipping emission grid
        !Interpolate to the shipping grid in lat lon coordinates
        !Temporary emissions cutoff
        i_ship=1
        where (shipping_nc_dp.lt.min_proxy_emission_shipping_value) shipping_nc_dp=0.
        write(unit_logfile,'(2a,2f12.2)') 'Shipping min and max: ',trim(var_name_nc_temp),minval(shipping_nc_dp(:,:,i_ship)),maxval(shipping_nc_dp(:,:,i_ship))


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
            correct_lon(1)=1./cos(3.14159/180.*temp_lat(1))
            correct_lon(2)=1.

            !Interpolate on same grid then scale, equivalent to interpolating density and then recalculating
            proxy_emission_subgrid(i,j,source_index,:)=area_weighted_extended_vectorgrid_interpolation_function( &
                real(var2d_nc_dp(1:dim_length_shipping_nc(x_dim_nc_index),x_dim_nc_index))*correct_lon(1),real(var2d_nc_dp(1:dim_length_shipping_nc(y_dim_nc_index),y_dim_nc_index)) &
                ,shipping_nc_dp(:,:,i_ship),dim_length_shipping_nc(x_dim_nc_index),dim_length_shipping_nc(y_dim_nc_index) &
                ,delta_shipping_nc*correct_lon,temp_lon(1)*correct_lon(1),temp_lat(1),delta_shipping_nc*correct_lon)
            
            temp_scale=(temp_delta(1)*correct_lon(1)*temp_delta(2)*correct_lon(2))/(delta_shipping_nc(1)*correct_lon(1)*delta_shipping_nc(2)*correct_lon(2))
            proxy_emission_subgrid(i,j,source_index,:)=proxy_emission_subgrid(i,j,source_index,:)*temp_scale
            
  
            if (isnan(proxy_emission_subgrid(i,j,source_index,1))) then
            write(*,*) 'Stopping, nan in proxy_emission_subgrid'
            write(*,*) temp_scale,correct_lon,delta_shipping_nc,temp_delta,temp_lon
            stop
            endif
            if (proxy_emission_subgrid(i,j,source_index,1).lt.0.) then
            write(*,*) 'Stopping, negative value in proxy_emission_subgrid'
            write(*,*) temp_scale,correct_lon,delta_shipping_nc,temp_delta,temp_lon
            stop
            endif

        enddo
        enddo
        

       
        if (allocated(shipping_nc_dp)) deallocate (shipping_nc_dp)
        if (allocated(var2d_nc_dp)) deallocate (var2d_nc_dp)
        if (allocated(temp_var2d_nc_dp)) deallocate (temp_var2d_nc_dp)
    
    end subroutine uEMEP_read_netcdf_shipping_latlon
    
    !----------------------------------------------------------------------
    subroutine match_string_multi_val(match_str,unit_in,unit_output,val,n_val)
    !Finds a leading string and returns all the integer variables that follows it
    !Tab delimitted before and free format after
    implicit none
    
    integer n_val
    real val(n_val)
    character (*) match_str
    character(256) temp_str1,temp_str2,temp_str
    integer unit_in,unit_output
    integer index_val
    
    val=-999.
    temp_str1=''
    temp_str2='Not available'
    rewind(unit_in)
    do while (index(temp_str1,match_str).eq.0)
        read(unit_in,'(a)',end=10) temp_str
        index_val=index(temp_str,achar(9))
        temp_str1=temp_str(1:index_val-1)
        temp_str=temp_str(index_val+1:)
        index_val=index(temp_str,achar(9))
        !if (index_val.gt.0) then
        !    temp_str2=temp_str(1:index_val-1)
        !else
        !    temp_str2=temp_str
        !endif
    end do
    if (LEN(trim(temp_str)).gt.0) then
        read(temp_str,*) val(1:n_val)
    else
        goto 15
    endif
    
    if (unit_output.ge.0) then
        write(unit_output,'(A40,A3,<n_val>es10.2)') trim(match_str),' = ',val
    endif
    return
    
10  write(unit_output,*) 'WARNING: No match found to "'//trim(match_str)//'" in input files. Set to -999'
    return
15	write(unit_output,*) 'WARNING: No values for "'//trim(match_str)//'" in input files'

    end subroutine match_string_multi_val
!----------------------------------------------------------------------

    subroutine uEMEP_preaggregate_shipping_asi_data
 
    !This routine aggregates shipping data in space and time
    !Reads in ASI data from standard files and aggregates in UTM33 100 m grids
    !Writes the data out again in standard ASI format
    !This routine is called if the flag 'preaggregate_shipping_asi_data_flag' is set to true
    use uEMEP_definitions
    
    implicit none
    
    character(256) temp_name
    character(1024) temp_str
    character(256) temp_str1,temp_str2
    real temp_val
    integer unit_in
    integer exists
    integer count,index_val
    integer temp_int
    real ddlatitude,ddlongitude,totalnoxemission,totalparticulatematteremission
    real y_ship,x_ship
    integer i_ship_index,j_ship_index
    integer source_index,subsource_index
    integer t
    integer ship_i_dim_index,ship_j_dim_index,ship_x_dim_index,ship_y_dim_index,ship_lat_dim_index,ship_lon_dim_index,ship_count_dim_index,ship_pm_dim_index,ship_nox_dim_index
    parameter (ship_i_dim_index=1,ship_j_dim_index=2,ship_count_dim_index=3)
    parameter (ship_x_dim_index=1,ship_y_dim_index=2,ship_lat_dim_index=3,ship_lon_dim_index=4,ship_pm_dim_index=5,ship_nox_dim_index=6)
    integer max_ship_dim_index
    parameter (max_ship_dim_index=500000)
    integer ship_index(max_ship_dim_index,3)
    real ship_value(max_ship_dim_index,6)
    integer ship_index_count,i_count
    logical found_index
    
    real :: ship_delta=250.
    integer :: i_ship_min=1000000,i_ship_max=-1000000,j_ship_min=1000000,j_ship_max=-1000000
    integer i_ship_dim_min,i_ship_dim_max,j_ship_dim_min,j_ship_dim_max
    parameter (i_ship_dim_min=-400,i_ship_dim_max=4500,j_ship_dim_min=25000,j_ship_dim_max=32000)
    integer ship_array_index(i_ship_dim_min:i_ship_dim_max,j_ship_dim_min:j_ship_dim_max)
    logical :: havbase_data_type=.false.
    
    if (.not.calculate_aggregated_shipping_emissions_flag) return

    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Aggregating shipping asi data  (uEMEP_preaggregate_shipping_asi_data)'
	write(unit_logfile,'(A)') '================================================================'
    
    source_index=shipping_index
    n_subsource(source_index)=1
    proxy_emission_subgrid(:,:,source_index,:)=0.
    t=1

    
    pathfilename_ship(1)=trim(pathname_ship(1))//trim(filename_ship(1))
    
    !Test existence of the shipping filename. If does not exist then use default
    inquire(file=trim(pathfilename_ship(1)),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Shipping file does not exist: ', trim(pathfilename_ship(1))
        stop
    endif

    ship_index_count=0
    ship_value=0.
    ship_index=0
    ship_array_index=0
    
    !Open the file for reading
    unit_in=20
    open(unit_in,file=pathfilename_ship(1),access='sequential',status='old',readonly)  
    write(unit_logfile,'(a)') ' Opening shipping file '//trim(pathfilename_ship(1))
    
    rewind(unit_in)

    subsource_index=1
    havbase_data_type=.true.
    !Read header old: ddlatitude;ddlongitude;totalnoxemission;totalparticulatematteremission;fk_vessellloydstype;fk_ais_norwegianmainvesselcategory;date;time
    !Read header new: mmsi;date_time_utc;lat;lon;lloydstype;norvesselcategory;sizegroupgrosston;vesselname;imonumber;dist_nextpoint;sec_nextpoint;fuelconsumption;me_fuelquality;co2emission;so2emission;particulatematteremission;noxemission;nmvocemission;ch4emission;n2oemission;coemission;blackcarbonemission;organiccarbonemission
    read(unit_in,'(A)') temp_str
    write(*,*) trim(temp_str)
    count=0
    do while(.not.eof(unit_in))
        read(unit_in,'(A)') temp_str
        !read(unit_in,*) temp_str
        !write(*,*) trim(temp_str)
        ddlatitude=0.;ddlongitude=0.;totalnoxemission=0.;totalparticulatematteremission=0.
        if (havbase_data_type) then
            !Extract the values in the temp_str
            index_val=index(temp_str,';',back=.false.);temp_str=temp_str(index_val+1:) !Skip mmsi
            index_val=index(temp_str,';',back=.false.);temp_str=temp_str(index_val+1:) !Skip date_time_utc
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) ddlatitude !Read an entry
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) ddlongitude !Read an entry
            index_val=index(temp_str,';',back=.false.);temp_str=temp_str(index_val+1:) !Skip lloydstype
            index_val=index(temp_str,';',back=.false.);temp_str=temp_str(index_val+1:) !Skip norvesselcategory
            index_val=index(temp_str,';',back=.false.);temp_str=temp_str(index_val+1:) !Skip sizegroupgrosston
            index_val=index(temp_str,';',back=.false.);temp_str=temp_str(index_val+1:) !Skip vesselname
            index_val=index(temp_str,';',back=.false.);temp_str=temp_str(index_val+1:) !Skip imonumber
            index_val=index(temp_str,';',back=.false.);temp_str=temp_str(index_val+1:) !Skip dist_nextpoint
            index_val=index(temp_str,';',back=.false.);temp_str=temp_str(index_val+1:) !Skip sec_nextpoint
            index_val=index(temp_str,';',back=.false.);temp_str=temp_str(index_val+1:) !Skip fuelconsumption
            index_val=index(temp_str,';',back=.false.);temp_str=temp_str(index_val+1:) !Skip me_fuelquality
            index_val=index(temp_str,';',back=.false.);temp_str=temp_str(index_val+1:) !Skip co2emission
            index_val=index(temp_str,';',back=.false.);temp_str=temp_str(index_val+1:) !Skip so2emission
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) totalparticulatematteremission !Read an entry
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) totalnoxemission !Read an entry
        else
            
            !Extract the values in the temp_str
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
            read(temp_str1,*) ddlatitude
            !write (*,*) ddlatitude
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
            read(temp_str1,*) ddlongitude
            !write (*,*) ddlongitude
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
            !write(*,*) index_val,trim(temp_str1),trim(temp_str)
            if (index_val.gt.1) read(temp_str1,*) totalnoxemission
            !write (*,*) totalnoxemission
            index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
            if (index_val.gt.1) read(temp_str1,*) totalparticulatematteremission
            !index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
            !if (index_val.gt.1) read(temp_str1,*) temp_int
            !index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
            !if (index_val.gt.1) read(temp_str1,*) temp_int
            !index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
            !if (index_val.gt.1) read(temp_str1,*) temp_str2
            !write (*,*) trim(temp_str1)
            !temp_str1=temp_str
            !if (len(temp_str1).gt.0) read(temp_str1,*) temp_str2
            !write (*,*) trim(temp_str1)
        endif
        
        !write(*,*) count,ddlatitude,ddlongitude,totalnoxemission,totalparticulatematteremission
        count=count+1
        if (mod(count,10000).eq.0) write(*,'(2i12,2f12.2,2e12.2)') count,ship_index_count,ddlatitude,ddlongitude,totalnoxemission,totalparticulatematteremission  
        
        if (totalnoxemission.gt.0.or.totalparticulatematteremission.gt.0) then
            
            !Convert lat lon to utm coords
            if (projection_type.eq.UTM_projection_index) then
                call LL2UTM(1,utm_zone,ddlatitude,ddlongitude,y_ship,x_ship)
            elseif (projection_type.eq.LAEA_projection_index) then
                call LL2LAEA(x_ship,y_ship,ddlongitude,ddlatitude,projection_attributes)
            endif
        
            !Find the grid index it belongs to. This assumes a minimum UTM grid at 0 so the index can be negative
            !Not certain if this 0.5 is correct
            i_ship_index=floor((x_ship-0)/ship_delta+0.5)
            j_ship_index=floor((y_ship-0)/ship_delta+0.5)
        
            i_ship_min=min(i_ship_index,i_ship_min)
            i_ship_max=max(i_ship_index,i_ship_max)
            j_ship_min=min(j_ship_index,j_ship_min)
            j_ship_max=max(j_ship_index,j_ship_max)
        
            if (i_ship_index.lt.i_ship_dim_min.or.i_ship_index.gt.i_ship_dim_max.or.j_ship_index.lt.j_ship_dim_min.or.j_ship_index.gt.j_ship_dim_max) then
                write(*,*) i_ship_dim_min,i_ship_dim_max,j_ship_dim_min,j_ship_dim_max
                write(*,*) i_ship_index,j_ship_index
                stop
            endif
        
            i_count=ship_array_index(i_ship_index,j_ship_index)
            if (i_count.eq.0) then
                ship_index_count=ship_index_count+1
                ship_array_index(i_ship_index,j_ship_index)=ship_index_count
                i_count=ship_index_count
            endif
            ship_index(i_count,ship_count_dim_index)=ship_index(i_count,ship_count_dim_index)+1
            if (totalparticulatematteremission.gt.0.and..not.isnan(totalparticulatematteremission)) ship_value(i_count,ship_pm_dim_index)=ship_value(i_count,ship_pm_dim_index)+totalparticulatematteremission
            if (totalnoxemission.gt.0.and..not.isnan(totalnoxemission)) ship_value(i_count,ship_nox_dim_index)=ship_value(i_count,ship_nox_dim_index)+totalnoxemission
            ship_index(i_count,ship_i_dim_index)=i_ship_index
            ship_index(i_count,ship_j_dim_index)=j_ship_index
        
            ship_value(i_count,ship_x_dim_index)=(ship_index(i_count,ship_i_dim_index)+.5)*ship_delta
            ship_value(i_count,ship_y_dim_index)=(ship_index(i_count,ship_j_dim_index)+.5)*ship_delta
            call PROJ2LL(ship_value(i_count,ship_x_dim_index),ship_value(i_count,ship_y_dim_index),ship_value(i_count,ship_lon_dim_index),ship_value(i_count,ship_lat_dim_index),projection_attributes,projection_type)
            call UTM2LL(utm_zone,ship_value(i_count,ship_y_dim_index),ship_value(i_count,ship_x_dim_index),ship_value(i_count,ship_lat_dim_index),ship_value(i_count,ship_lon_dim_index))
           
        !if (mod(count,10000).eq.0) write(*,'(2i,2f,2e)') count,ship_index_count,ddlatitude,ddlongitude,totalnoxemission,totalparticulatematteremission  
        !if (mod(count,10000).eq.0) write(*,'(3i12,4f14.4,2es14.5)') count,i_count,ship_index(i_count,ship_count_dim_index),ship_value(i_count,ship_y_dim_index),ship_value(i_count,ship_x_dim_index),ship_value(i_count,ship_lat_dim_index),ship_value(i_count,ship_lon_dim_index),ship_value(i_count,ship_pm_dim_index),ship_value(i_count,ship_nox_dim_index)
        !write(*,*) ship_index_count,x_ship,y_ship
        !write(*,*) i_count,i_ship_index,j_ship_index
        !write(*,*) ship_index(i_count,ship_i_dim_index),ship_index(i_count,ship_j_dim_index),ship_index(i_count,ship_count_dim_index)
 
        endif
        
        !Determine if the grid has been found before and add to it
        if(1.eq.2) then
        if ((totalnoxemission.gt.0.or.totalparticulatematteremission.gt.0)) then
            found_index=.false.
            do i_count=1,ship_index_count
                if (ship_index(i_count,ship_i_dim_index).eq.i_ship_index.and.ship_index(i_count,ship_j_dim_index).eq.j_ship_index) then
                    ship_index(i_count,ship_count_dim_index)=ship_index(i_count,ship_count_dim_index)+1
                    ship_value(i_count,ship_pm_dim_index)=ship_value(i_count,ship_pm_dim_index)+totalparticulatematteremission
                    ship_value(i_count,ship_nox_dim_index)=ship_value(i_count,ship_nox_dim_index)+totalnoxemission
                    found_index=.true.
                    exit
                endif           
            enddo
            if (.not.found_index) then
                ship_index_count=ship_index_count+1
                i_count=ship_index_count
                ship_index(i_count,ship_i_dim_index)=i_ship_index
                ship_index(i_count,ship_j_dim_index)=j_ship_index
                ship_index(i_count,ship_count_dim_index)=1
                ship_value(i_count,ship_pm_dim_index)=totalparticulatematteremission
                ship_value(i_count,ship_nox_dim_index)=totalnoxemission
            endif
            
            !write(*,*) count,ship_index_count      
        endif
        endif
        
    enddo
    write(unit_logfile,'(A,2I)') 'Shipping counts = ',count,ship_index_count  
 
    close(unit_in)
    


    !Calculate the x,y, lat,lon of the centre of the valid grids
    !do i_count=1,ship_index_count
        !ship_value(i_count,ship_x_dim_index)=(ship_index(i_count,ship_i_dim_index)+.5)*ship_delta
        !ship_value(i_count,ship_y_dim_index)=(ship_index(i_count,ship_j_dim_index)+.5)*ship_delta
        !call UTM2LL(utm_zone,ship_value(i_count,ship_y_dim_index),ship_value(i_count,ship_x_dim_index),ship_value(i_count,ship_lat_dim_index),ship_value(i_count,ship_lon_dim_index))
    !enddo
    
    !Save the data in the same format as previously provided
    pathfilename_ship(1)=trim(pathname_ship(1))//'Aggregated_'//trim(filename_ship(1))
    
    !Open the file for reading
    unit_in=20
    open(unit_in,file=pathfilename_ship(1),access='sequential',status='unknown')  
    write(unit_logfile,'(a)') ' Writing to aggregated shipping file '//trim(pathfilename_ship(1))

    write(unit_in,'(a)') 'ddlatitude;ddlongitude;totalnoxemission;totalparticulatematteremission;fk_vessellloydstype;fk_ais_norwegianmainvesselcategory;date;time'
    do i_count=1,ship_index_count
        if (ship_index(i_count,ship_count_dim_index).gt.1) then
            write(unit_in,'(f14.5,a1,f14.5,a1,e12.2,a1,e12.2,a4)') ship_value(i_count,ship_lat_dim_index),';',ship_value(i_count,ship_lon_dim_index),';',ship_value(i_count,ship_nox_dim_index),';',ship_value(i_count,ship_pm_dim_index),';;;;'
        endif
    enddo
      
    close(unit_in)
    
    !Max and min array dimensions found
    write(*,*) i_ship_min,i_ship_max,j_ship_min,j_ship_max
    stop
    
    end subroutine uEMEP_preaggregate_shipping_asi_data
    