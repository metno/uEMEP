    subroutine uEMEP_read_roadlink_data_ascii
    
    use uEMEP_definitions
    
    implicit none
    
    character(256) search_str,temp_str
    real temp
    integer unit_in
    integer rl_length_short
    integer exists
    logical nxtdat_flag
    integer n_roadlinks_major
    real sub_nodes_x(5000),sub_nodes_y(5000)
    integer temp_id,n_subnodes,temp_road_type,temp_nlanes
    real temp_adt,temp_hdv,temp_speed,temp_width
    integer counter
    real size_major(3)
    integer n_loop,loop_step

    real, allocatable :: inputdata_rl_temp(:,:)
    integer, allocatable :: inputdata_int_rl_temp(:,:)
    
    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading road link data ascii (uEMEP_read_roadlink_data_ascii)'
	write(unit_logfile,'(A)') '================================================================'
    
    if (read_existing_grid_data(proxy_emission_file_index(traffic_index))) then 
        write(unit_logfile,'(A)')'Reading existing traffic emissions'
        return
    endif
    

    min_adt=100.
    ratio_truck_car_emission=10.
    min_link_size=50.

    pathfilename_rl(1)=trim(pathname_rl(1))//trim(filename_rl(1))
    !write(*,*) pathname_rl(2),filename_rl(2),pathfilename_rl(2)
    
    !Test existence of the road link filename (2). If does not exist then use default
    inquire(file=trim(pathfilename_rl(1)),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Eoad link file ascii does not exist: ', trim(pathfilename_rl(1))
        stop
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
    allocate (inputdata_rl_temp(n_roadlinks,num_var_rl))
    allocate (inputdata_int_rl_temp(n_roadlinks,num_int_rl))
    
    !Initialise
    inputdata_rl_temp=0.
    inputdata_int_rl_temp=0

    counter=0
    !Read the data
    do i=1,n_roadlinks_major
        !ID ADT HDV ROAD_TYPE SPEED N_SUBLINKS
        read(unit_in,*,ERR=20) temp_id,temp_adt,temp_hdv,temp_road_type,temp_speed,temp_width,temp_nlanes,n_subnodes
        !write(*,*) i,temp_adt
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
        
        if (temp_adt.ge.min_adt) then
        do j=1,n_loop
            counter=counter+1
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
    
    return
20  write(unit_logfile,'(2A)') 'ERROR reading road link file: ',trim(pathfilename_rl(2))
    stop
    
    
    end subroutine uEMEP_read_roadlink_data_ascii