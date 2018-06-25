!uEMEP_read_shipping_asi_data.f90
    
    subroutine uEMEP_read_shipping_asi_data
 
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
    integer, allocatable :: count_subgrid(:,:)
    real, allocatable :: temp1_subgrid(:,:),temp2_subgrid(:,:),temp3_subgrid(:,:)
   
    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading shipping asi data  (uEMEP_read_shipping_asi_data)'
	write(unit_logfile,'(A)') '================================================================'
    
    source_index=shipping_index
    n_subsource(source_index)=1
    proxy_emission_subgrid(:,:,source_index,:)=0.
    t=1

    allocate (count_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index)))
    count_subgrid(:,:)=0
    allocate (temp1_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index)))
    allocate (temp2_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index)))
    allocate (temp3_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index)))
    
    if (read_existing_grid_data(proxy_emission_file_index(source_index))) then
        do subsource_index=1,n_subsource(source_index)
        temp_name=trim(pathname_grid(proxy_emission_file_index(source_index)))//trim(filename_grid(proxy_emission_file_index(source_index)))//trim(subsource_str(subsource_index))//'_'//trim(var_name_nc(conc_nc_index,compound_index,allsource_index))//'_'//trim(file_tag)//'.asc'
            inquire(file=temp_name,exist=exists)
            if (.not.exists) then
                write(unit_logfile,'(A)')'ERROR: '//trim(temp_name)//' does not exist.'
                stop
            endif
            call read_esri_ascii_file(unit_logfile,temp_name,emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index),emission_subgrid_delta(x_dim_index,source_index), &
                    temp1_subgrid,temp2_subgrid,temp3_subgrid)
                
            proxy_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index,subsource_index)=temp1_subgrid
            x_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index)=temp2_subgrid
            y_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index)=temp3_subgrid
                
        enddo
        write(unit_logfile,'(A)')'Reading existing shipping emissions'
        write(unit_logfile,'(A)')'Reading: '//trim(temp_name)
        return
    endif
    
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
            
        !Convert lat lon to utm coords
        call LL2UTM(1,utm_zone,ddlatitude,ddlongitude,y_ship,x_ship)
        
        !Find the grid index it belongs to
        i_ship_index=1+floor((x_ship-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
        j_ship_index=1+floor((y_ship-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))
        !(x_subgrid(i,j)-subgrid_min(1))/+subgrid_delta(1)+1=i
        
        !Add to subgrid
        
        if (i_ship_index.ge.1.and.i_ship_index.le.emission_subgrid_dim(x_dim_index,source_index) &
            .and.j_ship_index.ge.1.and.j_ship_index.le.emission_subgrid_dim(y_dim_index,source_index)) then
            if  (totalnoxemission.gt.0.and.compound_index.eq.nox_nc_index) then
            proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,subsource_index)=proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,subsource_index)+totalnoxemission
            count_subgrid(i_ship_index,j_ship_index)=count_subgrid(i_ship_index,j_ship_index)+1
            elseif (totalparticulatematteremission.gt.0.and.(compound_index.eq.pm25_nc_index.or.compound_index.eq.pm10_nc_index)) then
            proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,subsource_index)=proxy_emission_subgrid(i_ship_index,j_ship_index,source_index,subsource_index)+totalparticulatematteremission
            count_subgrid(i_ship_index,j_ship_index)=count_subgrid(i_ship_index,j_ship_index)+1
            endif
        endif
            
        endif
        
    enddo
    write(unit_logfile,'(A,I)') 'Shipping counts = ',count
    write(unit_logfile,'(A,es12.3)') 'Total emission = ',sum(proxy_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index,subsource_index))
 
    !Remove shipping with 1 or less values in the grid, but not if aggregated emissions are being used
    if (.not.use_aggregated_shipping_emissions_flag) then
    do j=1,emission_subgrid_dim(y_dim_index,source_index)
    do i=1,emission_subgrid_dim(x_dim_index,source_index)
        if (count_subgrid(i,j).le.1) then
            proxy_emission_subgrid(i,j,source_index,subsource_index)=0.
        endif
    enddo
    enddo
    endif

    
    close(unit_in)
    
    
    deallocate (count_subgrid)

    
    if (save_intermediate_files) then
    if (.not.read_existing_grid_data(proxy_emission_file_index(source_index))) then
        temp_name=trim(pathname_grid(proxy_emission_file_index(source_index)))//trim(filename_grid(proxy_emission_file_index(source_index)))//trim(subsource_str(subsource_index))//'_'//trim(var_name_nc(conc_nc_index,compound_index,allsource_index))//'_'//trim(file_tag)//'.asc'
        write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
        !write(*,*) emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index),emission_subgrid_delta(x_dim_index,source_index)
        !write(*,*) size(emission_subgrid,1),size(emission_subgrid,2),size(emission_subgrid,3),size(emission_subgrid,4),size(emission_subgrid,5)
        !write(*,*) size(x_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index),1),size(x_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index),2)
        !write(*,*) size(y_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index),1),size(y_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index),2)
       ! write(*,*) size(emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),t,source_index,subsource_index),1),size(emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),t,source_index,subsource_index),2)
        
        temp1_subgrid=proxy_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index,subsource_index)
        temp2_subgrid=x_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index)
        temp3_subgrid=y_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index)
        !write(*,*) size(temp_subgrid,1),size(temp_subgrid,2)
        call write_esri_ascii_file(unit_logfile,temp_name,emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index),emission_subgrid_delta(x_dim_index,source_index), &
                temp1_subgrid,temp2_subgrid,temp3_subgrid)
    endif
    endif
    
    deallocate (temp1_subgrid,temp2_subgrid,temp3_subgrid)
    
    end subroutine uEMEP_read_shipping_asi_data
    
    
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
            call LL2UTM(1,utm_zone,ddlatitude,ddlongitude,y_ship,x_ship)
        
            !Find the grid index it belongs to. This assumes a minimum UTM grid at 0 so the index can be negative
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
    do i_count=1,ship_index_count
        !ship_value(i_count,ship_x_dim_index)=(ship_index(i_count,ship_i_dim_index)+.5)*ship_delta
        !ship_value(i_count,ship_y_dim_index)=(ship_index(i_count,ship_j_dim_index)+.5)*ship_delta
        !call UTM2LL(utm_zone,ship_value(i_count,ship_y_dim_index),ship_value(i_count,ship_x_dim_index),ship_value(i_count,ship_lat_dim_index),ship_value(i_count,ship_lon_dim_index))
    enddo
    
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
    