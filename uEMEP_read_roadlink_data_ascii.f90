    subroutine uEMEP_read_roadlink_data_ascii
    
    use uEMEP_definitions
    
    implicit none
    
    integer i,j,k
    character(256) search_str,temp_str
    real temp
    integer unit_in
    integer rl_length_short
    integer exists
    logical nxtdat_flag
    real sub_nodes_x(5000),sub_nodes_y(5000)
    integer temp_id,n_subnodes,temp_road_type,temp_nlanes
    real temp_adt,temp_hdv,temp_speed,temp_width
    integer counter
    real size_major(3)
    integer n_loop,loop_step
    real x_grid_min,x_grid_max,y_grid_min,y_grid_max
    integer counter_major,counter_sub
    logical :: show_diagnostics=.false.
    real diagnostic_val(10)
    
    real, allocatable :: inputdata_rl_temp(:,:)
    integer, allocatable :: inputdata_int_rl_temp(:,:)
    
    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading road link data ascii (uEMEP_read_roadlink_data_ascii)'
	write(unit_logfile,'(A)') '================================================================'
    
    min_adt=10.
    min_link_size=1.

    pathfilename_rl(1)=trim(pathname_rl(1))//trim(filename_rl(1))
    !write(*,*) pathname_rl(2),filename_rl(2),pathfilename_rl(2)
    
    !Test existence of the road link filename (2). If does not exist then use default
    inquire(file=trim(pathfilename_rl(1)),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Road link file ascii does not exist: ', trim(pathfilename_rl(1))
        stop
    endif

   !Open the file for reading to test the available links in the region
    if (reduce_roadlink_region_flag) then
        unit_in=20
        open(unit_in,file=pathfilename_rl(1),access='sequential',status='old',readonly)  
        write(unit_logfile,'(a)') ' Opening road link file(ascii) '//trim(pathfilename_rl(1))
    
        rewind(unit_in)
        call NXTDAT(unit_in,nxtdat_flag)
        !read the header to find out how many links there are
        !read(unit_in,'(a)',ERR=20) temp_str
        read(unit_in,*,ERR=20) n_roadlinks_major,n_roadlinks
        write(unit_logfile,'(a,i)') ' Number of major road links= ', n_roadlinks_major
        write(unit_logfile,'(a,i)') ' Number of sub road links= ', n_roadlinks
             
        !Allocate the arrays after reading in the number of roads
        !allocate (inputdata_rl_temp(n_roadlinks,num_var_rl))
        !allocate (inputdata_int_rl_temp(n_roadlinks,num_int_rl))
        allocate (valid_link_flag(n_roadlinks_major))
    
        valid_link_flag=.false.
   
        !Initialise
        !inputdata_rl_temp=0.
        !inputdata_int_rl_temp=0

        counter_major=0
        counter_sub=0
        !Read the data to find roads in the tile
        x_grid_min=emission_subgrid_min(x_dim_index,traffic_index)
        x_grid_max=emission_subgrid_min(x_dim_index,traffic_index)+(emission_subgrid_dim(x_dim_index,traffic_index)+1)*emission_subgrid_delta(x_dim_index,traffic_index)
        y_grid_min=emission_subgrid_min(y_dim_index,traffic_index)
        y_grid_max=emission_subgrid_min(y_dim_index,traffic_index)+(emission_subgrid_dim(y_dim_index,traffic_index)+1)*emission_subgrid_delta(y_dim_index,traffic_index)
    
        do i=1,n_roadlinks_major
            !ID ADT HDV ROAD_TYPE SPEED N_SUBLINKS
            read(unit_in,*,ERR=20) temp_id,temp_adt,temp_hdv,temp_road_type,temp_speed,temp_width,temp_nlanes,n_subnodes
            !read(unit_in,*,ERR=20) !temp_id
            !write(*,*) temp_id,temp_adt,n_subnodes
            read(unit_in,*) sub_nodes_x(1)
            read(unit_in,*) sub_nodes_y(1)
    
            !Test position within emission region
            if (sub_nodes_x(1).ge.x_grid_min.and.sub_nodes_x(1).le.x_grid_max &
                .and.sub_nodes_y(1).ge.y_grid_min.and.sub_nodes_y(1).le.y_grid_max) then
                counter_major=counter_major+1
                counter_sub=counter_sub+n_subnodes-1
                valid_link_flag(i)=.true.
            endif
            
        enddo
        close(unit_in,status='keep')
           
        write(unit_logfile,'(a,i,i)') ' Number of major and sub road links within the region = ', counter_major,counter_sub
    endif

    !Open the file for reading
    unit_in=20
    open(unit_in,file=pathfilename_rl(1),access='sequential',status='old',readonly)  
    write(unit_logfile,'(a)') ' Opening road link file(ascii) '//trim(pathfilename_rl(1))
    
    rewind(unit_in)
    call NXTDAT(unit_in,nxtdat_flag)
    !read the header to find out how many links there are
    !read(unit_in,'(a)',ERR=20) temp_str
    read(unit_in,*,ERR=20) n_roadlinks_major,n_roadlinks
    write(unit_logfile,'(a,i)') ' Number of major road links= ', n_roadlinks_major
    write(unit_logfile,'(a,i)') ' Number of sub road links= ', n_roadlinks
             
    !Allocate the arrays after reading in the number of roads
    if (reduce_roadlink_region_flag) then
        n_roadlinks=counter_sub
        n_roadlinks_major_selected=counter_major
    else
        !Set all road links to valid when not regionally selecting
        allocate (valid_link_flag(n_roadlinks_major))
        valid_link_flag=.true.
    endif
    allocate (inputdata_rl_temp(n_roadlinks,num_var_rl))
    allocate (inputdata_int_rl_temp(n_roadlinks,num_int_rl))
    
    !Initialise
    inputdata_rl_temp=0.
    inputdata_int_rl_temp=0
    
    counter=0
    counter_major=0

    !Read the data
    do i=1,n_roadlinks_major
    if (.not.valid_link_flag(i)) then
        read(unit_in,*,ERR=20) 
        !write(*,*) temp_id,temp_adt,n_subnodes
        read(unit_in,*) 
        read(unit_in,*) 
    else
        !ID ADT HDV ROAD_TYPE SPEED N_SUBLINKS
        read(unit_in,*,ERR=20) temp_id,temp_adt,temp_hdv,temp_road_type,temp_speed,temp_width,temp_nlanes,n_subnodes
        !write(*,*) temp_id,temp_adt,n_subnodes
        read(unit_in,*) sub_nodes_x(1:n_subnodes)
        read(unit_in,*) sub_nodes_y(1:n_subnodes)
        !write(*,*) sub_nodes_x(1:n_subnodes),sub_nodes_y(1:n_subnodes)
        !put in the road link data
        
        !Test size of major link. If less than min_link_size then treat it as a single link from the start sub_node to the end sub_node
        size_major(1)=maxval(sub_nodes_x(1:n_subnodes))-minval(sub_nodes_x(1:n_subnodes))
        size_major(2)=maxval(sub_nodes_y(1:n_subnodes))-minval(sub_nodes_y(1:n_subnodes))
        size_major(3)=sqrt(size_major(1)**2+size_major(2)**2)
        if (size_major(3).gt.min_link_size) then
            n_loop=n_subnodes-1
            loop_step=1
        else
            n_loop=1
            loop_step=n_subnodes-1
            !write(*,*) size_major(3)
        endif
        
        counter_major=counter_major+1

        if (temp_adt.ge.min_adt) then
        do j=1,n_loop
            counter=counter+1          
            inputdata_int_rl_temp(counter,major_index_rl_index)=counter_major
            inputdata_int_rl_temp(counter,id_rl_index)=temp_id
            inputdata_rl_temp(counter,adt_rl_index)=temp_adt
            inputdata_rl_temp(counter,hdv_rl_index)=temp_hdv
            inputdata_int_rl_temp(counter,roadtype_rl_index)=temp_road_type
            inputdata_rl_temp(counter,speed_rl_index)=temp_speed
            inputdata_rl_temp(counter,width_rl_index)=temp_width
            inputdata_rl_temp(counter,x1_rl_index)=sub_nodes_x(j)
            inputdata_rl_temp(counter,x2_rl_index)=sub_nodes_x(j+loop_step)
            inputdata_int_rl_temp(counter,nlanes_rl_index)=temp_nlanes
            inputdata_rl_temp(counter,y1_rl_index)=sub_nodes_y(j)
            inputdata_rl_temp(counter,y2_rl_index)=sub_nodes_y(j+loop_step)
            !write(*,*) inputdata_int_rl(counter,id_rl_index),inputdata_rl(counter,x1_rl_index),inputdata_rl(counter,y2_rl_index)
        enddo
        endif
    
    endif    
    enddo
    n_roadlinks=counter
    write(unit_logfile,'(a,i)') ' Number of road links used = ', n_roadlinks
 
    close(unit_in,status='keep')
 
    !Allocate the arrays after reading in the number of roads
    allocate (inputdata_rl(n_roadlinks,num_var_rl))
    allocate (inputdata_int_rl(n_roadlinks,num_int_rl))

    inputdata_rl=inputdata_rl_temp(1:n_roadlinks,:)
    inputdata_int_rl=inputdata_int_rl_temp(1:n_roadlinks,:)
    
    deallocate (inputdata_rl_temp)
    deallocate (inputdata_int_rl_temp)
    
    !No speed in the files currently. Set all to 50 km/hr. Temporary
    !inputdata_rl(:,speed_rl_index)=50.
    !inputdata_rl(:,width_rl_index)=10.
    !inputdata_int_rl(:,nlanes_rl_index)=2
    
    !Set the road type, normal or tunnel (tunnel or jet). When a tunnel then there is no retention, always dry
    !do i=1,n_roadlinks        
    !    if (inputdata_int_rl(i,roadtype_rl_index).eq.5.or.inputdata_int_rl(i,roadtype_rl_index).eq.6) then           
    !        inputdata_int_rl(i,roadtype_rl_index)=tunnel_roadtype
    !    else
    !        inputdata_int_rl(i,roadtype_rl_index)=normal_roadtype
    !    endif
    !enddo

    !Calculate some additional values
    inputdata_rl(:,x0_rl_index)=(inputdata_rl(:,x1_rl_index)+inputdata_rl(:,x2_rl_index))/2.
    inputdata_rl(:,y0_rl_index)=(inputdata_rl(:,y1_rl_index)+inputdata_rl(:,y2_rl_index))/2.
    inputdata_rl(:,length_rl_index)=sqrt((inputdata_rl(:,x1_rl_index)-inputdata_rl(:,x2_rl_index))**2+(inputdata_rl(:,y1_rl_index)-inputdata_rl(:,y2_rl_index))**2)
    
    !Calculate road orientation and check for range overflows for length as well
    !inputdata_rl(angle_rl_index,:)=180./3.14159*acos((inputdata_rl(y2_rl_index,:)-inputdata_rl(y1_rl_index,:))/inputdata_rl(length_rl_index,:))
    do i=1,n_roadlinks
        call UTM2LL(utm_zone,inputdata_rl(i,y1_rl_index),inputdata_rl(i,x1_rl_index),inputdata_rl(i,lat1_rl_index),inputdata_rl(i,lon1_rl_index))
        call UTM2LL(utm_zone,inputdata_rl(i,y2_rl_index),inputdata_rl(i,x2_rl_index),inputdata_rl(i,lat2_rl_index),inputdata_rl(i,lon2_rl_index))
        call UTM2LL(utm_zone,inputdata_rl(i,y0_rl_index),inputdata_rl(i,x0_rl_index),inputdata_rl(i,lat0_rl_index),inputdata_rl(i,lon0_rl_index))
    enddo
        
    !Check lengths
    rl_length_short=0
    do i=1,n_roadlinks
        if (inputdata_rl(i,length_rl_index).eq.0.0) then
            rl_length_short=rl_length_short+1
            !write(unit_logfile,'(a,2i,f12.5)') ' WARNING: Zero link length, setting to 1 m ',i,inputdata_int_rl(i,id_rl_index),inputdata_rl(i,length_rl_index)
            inputdata_rl(i,length_rl_index)=1.
        endif
    enddo
    !write(*,*) 'Max length: ',maxval(inputdata_rl(:,length_rl_index)),inputdata_rl(maxloc(inputdata_rl(:,length_rl_index)),lat0_rl_index),inputdata_rl(maxloc(inputdata_rl(:,length_rl_index)),lon0_rl_index)
    write(unit_logfile,*) 'Number of road links with 0 length: ',rl_length_short,'  Setting to 1 m' 
    write(unit_logfile,*) 'Max road link x and y: ',maxval(inputdata_rl(:,x0_rl_index)),maxval(inputdata_rl(:,y0_rl_index))
    write(unit_logfile,*) 'Min road link x and y: ',minval(inputdata_rl(:,x0_rl_index)),minval(inputdata_rl(:,y0_rl_index))
    
    if (n_roadlinks.gt.0) then
    write(unit_logfile,'(a14,12a10)') ' LINK ','ID','X1','X2','Y1','Y2','WIDTH','LENGTH','ADT','LON','LAT','N_LANES','TYPE'
    i=1
    write(unit_logfile,'(a14,i10,7f10.1,2f10.4,2i10)') ' First link = ',inputdata_int_rl(i,id_rl_index),inputdata_rl(i,x1_rl_index),inputdata_rl(i,x2_rl_index) &
        ,inputdata_rl(i,y1_rl_index),inputdata_rl(i,y2_rl_index),inputdata_rl(i,width_rl_index) &
        ,inputdata_rl(i,length_rl_index),inputdata_rl(i,adt_rl_index) &
        ,inputdata_rl(i,lon0_rl_index),inputdata_rl(i,lat0_rl_index) &
        ,inputdata_int_rl(i,nlanes_rl_index),inputdata_int_rl(i,roadtype_rl_index)
    i=n_roadlinks
    write(unit_logfile,'(a14,i10,7f10.1,2f10.4,2i10)') ' Last link = ',inputdata_int_rl(i,id_rl_index),inputdata_rl(i,x1_rl_index),inputdata_rl(i,x2_rl_index) &
        ,inputdata_rl(i,y1_rl_index),inputdata_rl(i,y2_rl_index),inputdata_rl(i,width_rl_index) &
        ,inputdata_rl(i,length_rl_index),inputdata_rl(i,adt_rl_index) &
        ,inputdata_rl(i,lon0_rl_index),inputdata_rl(i,lat0_rl_index) &
        ,inputdata_int_rl(i,nlanes_rl_index),inputdata_int_rl(i,roadtype_rl_index)
    else
        write(unit_logfile,'(a)') 'No road links available in this region'
    endif

    !Calculate the veh km totals
    if (show_diagnostics) then
        diagnostic_val=0.
        do i=1,n_roadlinks
            !Total kilometres
            diagnostic_val(1)=diagnostic_val(1)+inputdata_rl(i,length_rl_index)/1000.
            !Total veh kilometres
            diagnostic_val(2)=diagnostic_val(2)+inputdata_rl(i,length_rl_index)/1000.*inputdata_rl(i,adt_rl_index)*365.
            !Light veh kilometres
            diagnostic_val(3)=diagnostic_val(3)+(1.-inputdata_rl(i,hdv_rl_index)/100.)*inputdata_rl(i,length_rl_index)/1000.*inputdata_rl(i,adt_rl_index)*365.
            !Light veh kilometres
            diagnostic_val(4)=diagnostic_val(4)+(inputdata_rl(i,hdv_rl_index)/100.)*inputdata_rl(i,length_rl_index)/1000.*inputdata_rl(i,adt_rl_index)*365.
        enddo
        write(unit_logfile,'(a,es12.4)') 'Total km= ',diagnostic_val(1)
        write(unit_logfile,'(a,es12.4)') 'Total veh.km= ',diagnostic_val(2)
        write(unit_logfile,'(a,es12.4)') 'Total light veh.km= ',diagnostic_val(3)
        write(unit_logfile,'(a,es12.4)') 'Total heavy veh.km= ',diagnostic_val(4)
    endif
    
    return
20  write(unit_logfile,'(2A)') 'ERROR reading road link file: ',trim(pathfilename_rl(2))
    stop
    
    
    end subroutine uEMEP_read_roadlink_data_ascii
    
    subroutine uEMEP_read_roadlink_emission_data
    !Reads in NORTRIP formatted emission data and places it in the correct road links
    
    use uEMEP_definitions
    
    implicit none
    
    integer i,j,k
    character(256) search_str,temp_str
    real temp
    integer unit_in
    integer rl_length_short
    integer exists
    logical nxtdat_flag
    real sub_nodes_x(5000),sub_nodes_y(5000)
    integer temp_id,n_subnodes,temp_road_type,temp_nlanes
    real temp_adt,temp_hdv,temp_speed,temp_width
    integer counter
    real size_major(3)
    integer n_loop,loop_step

    real, allocatable :: inputdata_rl_temp(:)
    integer, allocatable :: inputdata_int_rl_id(:)
    
    integer emission_date_array(6)
    integer n_roadlink_emission_compound
    character(16) n_roadlink_emission_compound_str(10)
    character(256) n_roadlink_emission_unit_str
    character(256) n_roadlink_emission_date_str
    integer n_roadlink_emission,n_roadlink_emission_time
    integer time_index_temp,t_match_index,t
    integer date_array_temp(6)
    integer n_roadlink_emission_selected
    character(256) format_temp
    
    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading road link emission data (uEMEP_read_roadlink_emission_data)'
	write(unit_logfile,'(A)') '================================================================'
    

    pathfilename_rl(2)=trim(pathname_rl(2))//trim(filename_rl(2))
    !write(*,*) pathname_rl(2),filename_rl(2),pathfilename_rl(2)
    
    !Test existence of the road link filename (2). If does not exist then use default
    inquire(file=trim(pathfilename_rl(2)),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Road link emission file does not exist: ', trim(pathfilename_rl(2))
        stop
    endif

    !Open the file for reading
    unit_in=20
    open(unit_in,file=pathfilename_rl(2),access='sequential',status='old',readonly)  
    write(unit_logfile,'(a)') ' Opening road link file(ascii) '//trim(pathfilename_rl(2))
    
    rewind(unit_in)
    call NXTDAT(unit_in,nxtdat_flag)
    !read the header to find out how many links there are
    !read(unit_in,'(a)',ERR=20) temp_str
    read(unit_in,*,ERR=20) n_roadlink_emission_compound
    write(unit_logfile,'(a,i)') ' Number of road link emission compounds= ', n_roadlink_emission_compound
    call NXTDAT(unit_in,nxtdat_flag)
    read(unit_in,*,ERR=20) n_roadlink_emission_compound_str(1:n_roadlink_emission_compound)
    write(unit_logfile,'(a,<n_roadlink_emission_compound>a16)') ' Road link emission compounds= ', n_roadlink_emission_compound_str(1:n_roadlink_emission_compound)
    call NXTDAT(unit_in,nxtdat_flag)
    read(unit_in,*,ERR=20) n_roadlink_emission_unit_str
    write(unit_logfile,'(a,a)') ' Road link emission units= ', trim(n_roadlink_emission_unit_str)
    call NXTDAT(unit_in,nxtdat_flag)
    read(unit_in,*,ERR=20) n_roadlink_emission_date_str
    write(unit_logfile,'(a,a)') ' Road link emission start date= ', trim(n_roadlink_emission_date_str)
    call NXTDAT(unit_in,nxtdat_flag)
    read(unit_in,*,ERR=20) n_roadlink_emission,n_roadlink_emission_time
    write(unit_logfile,'(a,i)') ' Number of road links= ', n_roadlink_emission
    write(unit_logfile,'(a,i)') ' Number of time steps= ', n_roadlink_emission_time

    if (n_roadlink_emission.ne.n_roadlinks_major) then
        write(unit_logfile,'(A,2i12)') 'ERROR: Number of emission road links is not the same as the static road links: ',n_roadlink_emission,n_roadlinks_major
        stop
    endif
    
    !Check that start time and end time are covered in the emission data before progessing further
    !DOES NOT WORK WITH SINGLE TIME LOOP. FIX!!!
    format_temp='yyyymmddHH'
    call datestr_to_date(n_roadlink_emission_date_str,format_temp,emission_date_array)
    if (use_single_time_loop_flag) then
        !Does not check for t_match but assumes it is 1
        time_index_temp=end_time_loop_index
        t_match_index=1
    else
        time_index_temp=subgrid_dim(t_dim_index)
        t_match_index=0
        !write(*,*) shape(val_dim_nc)
        do t=1,time_index_temp
            call number_to_date(val_dim_nc(t,time_dim_nc_index),date_array_temp,ref_year_EMEP)
            if (date_array_temp(1).eq.emission_date_array(1).and.date_array_temp(2).eq.emission_date_array(2).and. &
                date_array_temp(3).eq.emission_date_array(3).and.date_array_temp(4).eq.emission_date_array(4)) then
                t_match_index=t
            endif
        enddo
        if (t_match_index.eq.0) then
            write(unit_logfile,'(A,6i6)') 'ERROR: No starting date found in road emission data: ',emission_date_array
            stop
        else
            write(unit_logfile,'(A,6i6)') ' Road link emission starting date found. Index: ',t_match_index
        endif
        if (n_roadlink_emission_time-t_match_index+1.lt.time_index_temp) then
            write(unit_logfile,'(A,2i6)') 'ERROR: Not enough time data in road link emission files: ',n_roadlink_emission_time-t_match_index+1,time_index_temp
            stop
        endif 
    endif
    
    !Allocate the arrays after reading in the number of roads
    n_roadlink_emission_selected=n_roadlink_emission
    if (reduce_roadlink_region_flag) then
        n_roadlink_emission_selected=n_roadlinks_major_selected
    endif
    allocate (inputdata_rl_emissions(n_roadlink_emission_selected,n_roadlink_emission_time,n_roadlink_emission_compound))
    allocate (inputdata_rl_temp(n_roadlink_emission_time))
    allocate (inputdata_int_rl_id(n_roadlink_emission_selected))
    
    !Initialise
    inputdata_rl_temp=0.

    counter=0
    !Read the data
    call NXTDAT(unit_in,nxtdat_flag)
    do i=1,n_roadlink_emission
        if (valid_link_flag(i)) then
            counter=counter+1
            read(unit_in,*,ERR=20) inputdata_int_rl_id(counter)
            !write(*,*) i,inputdata_int_rl_id(i)
            do j=1,n_roadlink_emission_compound
                read(unit_in,*) inputdata_rl_temp(1:n_roadlink_emission_time)
                inputdata_rl_emissions(counter,1:time_index_temp,j)=inputdata_rl_temp(t_match_index:t_match_index+time_index_temp-1)
            enddo
        else
            read(unit_in,*,ERR=20) 
            do j=1,n_roadlink_emission_compound
                read(unit_in,*) 
            enddo           
        endif
        
    enddo
    write(unit_logfile,'(a,i)') ' Number of road links that should be read = ', n_roadlink_emission_selected
    write(unit_logfile,'(a,i)') ' Number of road links read = ', counter
 
    close(unit_in,status='keep')

    !Check that road link ID's match
    do i=1,n_roadlinks
        if (inputdata_int_rl(i,id_rl_index).ne.inputdata_int_rl_id(inputdata_int_rl(i,major_index_rl_index))) then
            write(unit_logfile,'(A,3i12)') 'ERROR: Mismatch of road link IDs in the emission files: ',i,inputdata_int_rl(i,id_rl_index),inputdata_int_rl_id(i)
            stop
        endif
    enddo
    
    if (allocated(inputdata_rl_temp)) deallocate(inputdata_rl_temp)
    if (allocated(inputdata_int_rl_id)) deallocate(inputdata_int_rl_id)
    
    return
20  write(unit_logfile,'(2A)') 'ERROR reading road link emission file: ',trim(pathfilename_rl(2))
    stop
    
    
    end subroutine uEMEP_read_roadlink_emission_data
    
    
    subroutine uEMEP_change_road_data
    
    use uEMEP_definitions
    
    implicit none
    
    integer i,j,k
    character(256) pathfilename_rl_change
    integer unit_in
    integer change_offset_index,change_scale_index,change_replace_index
    parameter (change_offset_index=1,change_scale_index=2,change_replace_index=3)
    !integer change_adt_index,change_hdv_index,change_speed_index
    !parameter (change_adt_index=1,change_hdv_index=2,change_speed_index=3)
    real change_val(num_var_rl,3)
    character(256) temp_str
    real change_x(2),change_y(2)
    integer count
    integer change_loop(3),change_index
    integer :: n_change_loop=3
    integer exists
    
    change_loop(1)=adt_rl_index
    change_loop(2)=hdv_rl_index
    change_loop(3)=speed_rl_index
    
    !If there is no path or file name for the replacement file then do not calculate
    if (pathname_rl_change.eq.''.or.filename_rl_change.eq.'') return
    
    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Changing road link data (uEMEP_change_road_data)'
	write(unit_logfile,'(A)') '================================================================'
    
    !Read in replacement file
    pathfilename_rl_change=trim(pathname_rl_change)//trim(filename_rl_change)
    
    !Test existence of the road link filename (2). If does not exist then use default
    inquire(file=trim(pathfilename_rl_change),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Road link change file does not exist: ', trim(pathfilename_rl_change)
        stop
    endif
    
     !Open the file for reading
    unit_in=20
    open(unit_in,file=pathfilename_rl_change,access='sequential',status='old',readonly)  
    write(unit_logfile,'(a)') ' Opening road link change file(ascii) '//trim(pathfilename_rl_change)
    rewind(unit_in)
   
    !Skip over coordinates header
    read(unit_in,*) temp_str

    !Read coordinates
    read(unit_in,*) change_x(1),change_y(1),change_x(2),change_y(2)
    !write(*,*) change_x(1),change_y(1),change_x(2),change_y(2)
    
    !Skip over values header
    read(unit_in,*) temp_str
    
    !Read values
    do j=1,n_change_loop
        change_index=change_loop(j)
        !write(*,*) change_index
        read(unit_in,*) temp_str,change_val(change_index,change_offset_index),change_val(change_index,change_scale_index),change_val(change_index,change_replace_index)
        !write(*,*) trim(temp_str),change_val(change_index,change_offset_index),change_val(change_index,change_scale_index),change_val(change_index,change_replace_index)
    enddo
    close(unit_in)
    
    !Search for road links. If found change them 
    count=0
    do i=1,n_roadlinks
        if(inputdata_rl(i,x0_rl_index).ge.change_x(1).and.inputdata_rl(i,x0_rl_index).le.change_x(2).and. &
           inputdata_rl(i,y0_rl_index).ge.change_y(1).and.inputdata_rl(i,y0_rl_index).le.change_y(2)) then
           count=count+1
           
             do j=1,n_change_loop
                change_index=change_loop(j)
                inputdata_rl(i,change_index)=change_val(change_index,change_offset_index)+change_val(change_index,change_scale_index)*inputdata_rl(i,change_index)
                if (change_val(change_index,change_replace_index).ne.0) inputdata_rl(i,change_index)=change_val(change_index,change_replace_index)
             enddo
         
        endif
        
    enddo
    
    write(unit_logfile,'(A,i)') 'Number of road links changed = ',count

    end subroutine uEMEP_change_road_data
