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
    real sub_nodes_x(5000),sub_nodes_y(5000),sub_nodes_lon(5000),sub_nodes_lat(5000)
    integer temp_id,n_subnodes,temp_road_type,temp_nlanes
    real temp_adt,temp_hdv,temp_speed,temp_width
    integer counter
    real size_major(3)
    integer n_loop,loop_step
    real x_grid_min,x_grid_max,y_grid_min,y_grid_max
    integer counter_major,counter_sub
    logical :: show_diagnostics=.false.
    real diagnostic_val(10)
    integer temp_category,temp_structure_type,temp_region_id,temp_surface_id,temp_route_id
    real temp_length,temp_tunnel_length
    logical :: road_data_in_latlon=.false. 
    real temp_real
    
    real, allocatable :: inputdata_rl_temp(:,:)
    integer, allocatable :: inputdata_int_rl_temp(:,:)
    real, allocatable :: inputdata_rl_multi(:,:)
    integer, allocatable :: inputdata_int_rl_multi(:,:)
    real, allocatable :: inputdata_rl_multi_new(:,:)
    integer, allocatable :: inputdata_int_rl_multi_new(:,:)
    integer n_multi_roadlinks_new,n_multi_roadlinks
    integer m
    integer n_road_link_file_loop
    logical :: first_road_link_file_read=.false. 
    real temp_val
    integer temp_int
    
    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading road link data ascii (uEMEP_read_roadlink_data_ascii)'
	write(unit_logfile,'(A)') '================================================================'
    
    
    min_adt=10.
    min_link_size=1.
    
    if (num_multiple_roadlink_files.eq.0) then
        n_road_link_file_loop=1
    else
        n_road_link_file_loop=num_multiple_roadlink_files
    endif
    
    !Start the file loop here
    !---------------------------
    n_multi_roadlinks_new=0
    n_multi_roadlinks=0
    first_road_link_file_read=.false.
    
    do m=1,n_road_link_file_loop
    
    if (num_multiple_roadlink_files.gt.0) then
        !pathname_rl(1)=pathname_mrl(m)
        filename_rl(1)=filename_mrl(m)
        write(*,*) m,n_road_link_file_loop,trim(filename_rl(1))
    endif
    
    pathfilename_rl(1)=trim(pathname_rl(1))//trim(filename_rl(1))
    !write(*,*) pathname_rl(2),filename_rl(2),pathfilename_rl(2)
    
    !Test existence of the road link filename (2). If does not exist then stop
    inquire(file=trim(pathfilename_rl(1)),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Road link file ascii does not exist: ', trim(pathfilename_rl(1))
        stop
    endif

    !Check filename for the string 'latlon' that specifies if the data is in lat lon coordinates
    if (index(filename_rl(1),'latlon').gt.0) then
        road_data_in_latlon=.true. 
        write(unit_logfile,'(A,A)') ' Reading road data positions as lat lon: ', trim(pathfilename_rl(1))
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
        if (no_header_roadlink_data_flag) then
            write(unit_logfile,'(a)') ' Reading road link file(ascii) without header: '//trim(pathfilename_rl(1))
            i=0
            do while(.not.eof(unit_in))
                if (.not.eof(unit_in)) read(unit_in,*,ERR=20) temp_real,temp_adt,temp_hdv,temp_road_type,temp_speed,temp_width,temp_nlanes,n_subnodes
                if (.not.eof(unit_in)) read(unit_in,*,ERR=20) sub_nodes_x(1) !Read x nodes
                if (.not.eof(unit_in)) read(unit_in,*,ERR=20) sub_nodes_y(1) !Read y nodes
                if (temp_real.gt.0.and..not.eof(unit_in)) then
                    i=i+1
                    n_roadlinks=n_roadlinks+n_subnodes-1
                endif
                if (int(temp_real).ne.temp_real) then
                    write(unit_logfile,'(a,i,f)') ' Problem with record at point with ID: ',i,temp_real
                    stop
                endif
                
            enddo
            n_roadlinks_major=i
            !n_roadlinks=0
            rewind(unit_in)
            call NXTDAT(unit_in,nxtdat_flag)
        else
            read(unit_in,*,ERR=20) n_roadlinks_major,n_roadlinks
        endif
        
        if (n_roadlinks.eq.0) then
            write(unit_logfile,'(a)') ' Reading road link file(ascii) with header but without subnode counts: '//trim(pathfilename_rl(1))
            rewind(unit_in)
            call NXTDAT(unit_in,nxtdat_flag)
            read(unit_in,*,ERR=20) n_roadlinks_major,n_roadlinks
            i=0
            do while(.not.eof(unit_in))
                if (.not.eof(unit_in)) read(unit_in,*,ERR=20) temp_id,temp_adt,temp_hdv,temp_road_type,temp_speed,temp_width,temp_nlanes,n_subnodes !Read attributes up to n_subnodes
                if (.not.eof(unit_in)) read(unit_in,*,ERR=20) sub_nodes_x(1) !Read x nodes
                if (.not.eof(unit_in)) read(unit_in,*,ERR=20) sub_nodes_y(1) !Read y nodes
                !write(*,*) temp_id
                if (temp_id.gt.0.and..not.eof(unit_in)) then
                i=i+1
                n_roadlinks=n_roadlinks+n_subnodes-1
                endif
            enddo
            !Have commented out this in the cases where the number of links are not as written. Can happen with OSM data
            !n_roadlinks_major=i
            rewind(unit_in)
            call NXTDAT(unit_in,nxtdat_flag)
            read(unit_in,*,ERR=20) temp_int,temp_int
        endif

        write(unit_logfile,'(a,i)') ' Number of major road links in header = ', n_roadlinks_major
        write(unit_logfile,'(a,i)') ' Number of sub road links in header= ', n_roadlinks
            
        
        !Allocate the arrays after reading in the number of roads
        !allocate (inputdata_rl_temp(n_roadlinks,num_var_rl))
        !allocate (inputdata_int_rl_temp(n_roadlinks,num_int_rl))
        if (allocated(valid_link_flag)) deallocate (valid_link_flag)
        allocate (valid_link_flag(n_roadlinks_major))
        valid_link_flag=.false.
   
        !Initialise
        !inputdata_rl_temp=0.
        !inputdata_int_rl_temp=0

        counter_major=0
        counter_sub=0
        !Read the data to find roads in the tile
        !x_grid_min=emission_subgrid_min(x_dim_index,traffic_index)
        !x_grid_max=emission_subgrid_min(x_dim_index,traffic_index)+(emission_subgrid_dim(x_dim_index,traffic_index)+1)*emission_subgrid_delta(x_dim_index,traffic_index)
        !y_grid_min=emission_subgrid_min(y_dim_index,traffic_index)
        !y_grid_max=emission_subgrid_min(y_dim_index,traffic_index)+(emission_subgrid_dim(y_dim_index,traffic_index)+1)*emission_subgrid_delta(y_dim_index,traffic_index)
    
        x_grid_min=init_emission_subgrid_min(x_dim_index,traffic_index)
        x_grid_max=init_emission_subgrid_min(x_dim_index,traffic_index)+(init_emission_subgrid_dim(x_dim_index,traffic_index)+1)*init_emission_subgrid_delta(x_dim_index,traffic_index)
        y_grid_min=init_emission_subgrid_min(y_dim_index,traffic_index)
        y_grid_max=init_emission_subgrid_min(y_dim_index,traffic_index)+(init_emission_subgrid_dim(y_dim_index,traffic_index)+1)*init_emission_subgrid_delta(y_dim_index,traffic_index)

        !Under the special case where multiple roads are read from different countries then us the intital grid to determine which to keep
        !This is not implemented yet but needs to select the region when reading multiple files
        !if 
        !init_subgrid_min(x_dim_index)=subgrid_min(x_dim_index)
        !init_subgrid_min(y_dim_index)=subgrid_min(y_dim_index)
        !init_subgrid_max(x_dim_index)=subgrid_max(x_dim_index)
        !init_subgrid_max(y_dim_index)=subgrid_max(y_dim_index)
        !endif
        
        !rewind(unit_in)
        !call NXTDAT(unit_in,nxtdat_flag)
        !read(unit_in,*,ERR=20) temp_int,temp_int
    
        do i=1,n_roadlinks_major
            !ID ADT HDV ROAD_TYPE SPEED N_SUBLINKS
            read(unit_in,*,ERR=20) temp_id,temp_adt,temp_hdv,temp_road_type,temp_speed,temp_width,temp_nlanes,n_subnodes
            !read(unit_in,*,ERR=20) !temp_id
            !write(*,*) temp_id,temp_adt,n_subnodes
    
            if (road_data_in_latlon) then
                !Order is reversed in these files
                read(unit_in,*) sub_nodes_y(1)
                read(unit_in,*) sub_nodes_x(1)
        
                if  (projection_type.eq.UTM_projection_index) then
                    !write(*,*) i,sub_nodes_x(1),sub_nodes_y(1)
                    call LL2UTM(1,utm_zone,sub_nodes_y(1),sub_nodes_x(1),sub_nodes_y(1),sub_nodes_x(1))
                    !write(*,*) i,sub_nodes_x(1),sub_nodes_y(1)
                elseif  (projection_type.eq.LTM_projection_index) then
                    call LL2LTM(1,ltm_lon0,sub_nodes_y(1),sub_nodes_x(1),sub_nodes_y(1),sub_nodes_x(1))
                elseif (projection_type.eq.LCC_projection_index) then
                elseif (projection_type.eq.PS_projection_index) then
                elseif (projection_type.eq.LAEA_projection_index) then
                    call LL2LAEA(sub_nodes_x(1),sub_nodes_y(1),sub_nodes_x(1),sub_nodes_y(1),projection_attributes)
                endif
        
            else
                read(unit_in,*) sub_nodes_x(1)
                read(unit_in,*) sub_nodes_y(1)            
            endif
            !Convert to EMEP coordinates from UTM to lambert. No choices here. Special case but should be changed
            if (save_emissions_for_EMEP(traffic_index)) then
                call PROJ2LL(sub_nodes_x(1),sub_nodes_y(1),sub_nodes_lon(1),sub_nodes_lat(1),projection_attributes,projection_type)
                !call UTM2LL(utm_zone,sub_nodes_y(1),sub_nodes_x(1),sub_nodes_lat(1),sub_nodes_lon(1))
                if (projection_type.eq.LCC_projection_index) then
                    call lb2lambert2_uEMEP(sub_nodes_x(1),sub_nodes_y(1),sub_nodes_lon(1),sub_nodes_lat(1),EMEP_projection_attributes)
                elseif (projection_type.eq.PS_projection_index) then
                    call LL2PS_spherical(sub_nodes_x(1),sub_nodes_y(1),sub_nodes_lon(1),sub_nodes_lat(1),EMEP_projection_attributes)
                else
                    !Remains as lat lon
                endif
                
                !write(*,*) sub_nodes_x(1),sub_nodes_y(1),sub_nodes_lon(1),sub_nodes_lat(1)
            endif
                !write(*,*) sub_nodes_x(1),sub_nodes_y(1)
                !write(*,*) x_grid_min,x_grid_max,y_grid_min,y_grid_max
            
            !Test position within emission region
            if (sub_nodes_x(1).ge.x_grid_min.and.sub_nodes_x(1).le.x_grid_max &
                .and.sub_nodes_y(1).ge.y_grid_min.and.sub_nodes_y(1).le.y_grid_max) then
                counter_major=counter_major+1
                counter_sub=counter_sub+n_subnodes-1
                valid_link_flag(i)=.true.
            else
                valid_link_flag(i)=.false.
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

        if (no_header_roadlink_data_flag) then
            write(unit_logfile,'(a)') ' Reading road link file(ascii) without header: '//trim(pathfilename_rl(1))
            i=0
            do while(.not.eof(unit_in))
                if (.not.eof(unit_in)) read(unit_in,*,ERR=20) temp_id,temp_adt,temp_hdv,temp_road_type,temp_speed,temp_width,temp_nlanes,n_subnodes !Read attributes
                if (.not.eof(unit_in)) read(unit_in,*,ERR=20) sub_nodes_x(1) !Read x nodes
                if (.not.eof(unit_in)) read(unit_in,*,ERR=20) sub_nodes_y(1) !Read y nodes
                if (temp_id.gt.0.and..not.eof(unit_in)) then
                i=i+1
                endif
            enddo
            n_roadlinks_major=i
            n_roadlinks=n_roadlinks+n_subnodes-1
            rewind(unit_in)
            call NXTDAT(unit_in,nxtdat_flag)
        else
            read(unit_in,*,ERR=20) n_roadlinks_major,n_roadlinks
        endif
        
        if (n_roadlinks.eq.0.and..not.reduce_roadlink_region_flag) then
            write(unit_logfile,'(a)') ' Reading road link file(ascii) with header but without subnode counts: '//trim(pathfilename_rl(1))
            rewind(unit_in)
            call NXTDAT(unit_in,nxtdat_flag)
            read(unit_in,*,ERR=20) n_roadlinks_major,n_roadlinks
            i=0
            do while(.not.eof(unit_in))
                if (.not.eof(unit_in)) read(unit_in,*,ERR=20) temp_id,temp_adt,temp_hdv,temp_road_type,temp_speed,temp_width,temp_nlanes,n_subnodes !Read attributes
                if (.not.eof(unit_in)) read(unit_in,*,ERR=20) sub_nodes_x(1) !Read x nodes
                if (.not.eof(unit_in)) read(unit_in,*,ERR=20) sub_nodes_y(1) !Read y nodes
                !write(*,*) temp_id
                if (temp_id.gt.0.and..not.eof(unit_in)) then
                i=i+1
                n_roadlinks=n_roadlinks+n_subnodes-1
               endif
            enddo
            !n_roadlinks_major=i
            rewind(unit_in)
            call NXTDAT(unit_in,nxtdat_flag)
        endif

        write(unit_logfile,'(a,i)') ' Number of major road links= ', n_roadlinks_major
        write(unit_logfile,'(a,i)') ' Number of sub road links= ', n_roadlinks
            
    !Allocate the arrays after reading in the number of roads
    if (reduce_roadlink_region_flag) then
        n_roadlinks=counter_sub
        n_roadlinks_major_selected=counter_major
    else
        !Set all road links to valid when not regionally selecting. THis will fail because already allocated
        if (allocated(valid_link_flag)) deallocate (valid_link_flag)
        allocate (valid_link_flag(n_roadlinks_major))
        valid_link_flag=.true.
    endif
    
    if (allocated(inputdata_rl_temp)) deallocate (inputdata_rl_temp)
    if (allocated(inputdata_int_rl_temp)) deallocate (inputdata_int_rl_temp)
    allocate (inputdata_rl_temp(n_roadlinks,num_var_rl))
    allocate (inputdata_int_rl_temp(n_roadlinks,num_int_rl))
    
    !Initialise
    inputdata_rl_temp=0.
    inputdata_int_rl_temp=0
    
    counter=0
    counter_major=0

    !rewind(unit_in)
    !call NXTDAT(unit_in,nxtdat_flag)
    !read(unit_in,*,ERR=20) temp_int,temp_int
    !Read the data
    do i=1,n_roadlinks_major
!    if (.not.valid_link_flag(i)) then
!        if (read_OSM_roadlink_data_flag) then
            !ID ADT HDV ROAD_TYPE SPEED N_SUBLINKS
!            if (.not.eof(unit_in)) read(unit_in,*,ERR=20) temp_id,temp_adt,temp_hdv,temp_road_type,temp_speed,temp_width,temp_nlanes,n_subnodes
!        else
            !ID ADT HDV ROAD_ACTIVITY_TYPE SPEED ROAD_WIDTH N_LANES N_SUBNODES ROAD_CATEGORY ROAD_LENGTH ROAD_STRUCTURE_TYPE REGION_ID ROAD_SURFACE_ID TUNNEL_LENGTH ROUTE_ID
!            if (.not.eof(unit_in)) read(unit_in,*,ERR=20) temp_id,temp_adt,temp_hdv,temp_road_type,temp_speed,temp_width,temp_nlanes,n_subnodes &
!            ,temp_category,temp_length,temp_structure_type,temp_region_id,temp_surface_id,temp_tunnel_length,temp_route_id
!        endif
        !write(*,*) temp_id,temp_adt,n_subnodes
!        if (.not.eof(unit_in)) read(unit_in,*,ERR=20) temp_val
!        if (.not.eof(unit_in)) read(unit_in,*,ERR=20) temp_val
!    else
        if (read_OSM_roadlink_data_flag) then
            !ID ADT HDV ROAD_TYPE SPEED N_SUBLINKS
            if (.not.eof(unit_in)) read(unit_in,*,ERR=20) temp_id,temp_adt,temp_hdv,temp_road_type,temp_speed,temp_width,temp_nlanes,n_subnodes
            temp_category=0;temp_length=0;temp_structure_type=0;temp_region_id=0;temp_surface_id=0;temp_tunnel_length=0;temp_route_id=0
        else
            !ID ADT HDV ROAD_ACTIVITY_TYPE SPEED ROAD_WIDTH N_LANES N_SUBNODES ROAD_CATEGORY ROAD_LENGTH ROAD_STRUCTURE_TYPE REGION_ID ROAD_SURFACE_ID TUNNEL_LENGTH ROUTE_ID
            if (.not.eof(unit_in)) read(unit_in,*,ERR=20) temp_id,temp_adt,temp_hdv,temp_road_type,temp_speed,temp_width,temp_nlanes,n_subnodes &
            ,temp_category,temp_length,temp_structure_type,temp_region_id,temp_surface_id,temp_tunnel_length,temp_route_id
        endif
        !write(*,*) i,temp_id,temp_adt,n_subnodes
        if (.not.eof(unit_in)) read(unit_in,*) sub_nodes_x(1:n_subnodes)
        if (.not.eof(unit_in)) read(unit_in,*) sub_nodes_y(1:n_subnodes)
        !write(*,*) sub_nodes_x(1:n_subnodes),sub_nodes_y(1:n_subnodes)
        !put in the road link data
        
        !Test size of major link. If less than min_link_size then treat it as a single link from the start sub_node to the end sub_node
    if (valid_link_flag(i)) then
        
        !Do not do this with latlon input data as it does not work
        if (.not.road_data_in_latlon) then
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
        else
            n_loop=n_subnodes-1
            loop_step=1
        endif
        
        counter_major=counter_major+1
        !if (n_subnodes.gt.5000) then
        !    write(*,*) 'Greater than 5000 ',i
        !endif
        
        !Place the data in the road links
        if (temp_adt.ge.min_adt) then
        do j=1,n_loop
            counter=counter+1
            if (counter.le.n_roadlinks) then
            inputdata_int_rl_temp(counter,major_index_rl_index)=counter_major
            inputdata_int_rl_temp(counter,id_rl_index)=temp_id
            inputdata_rl_temp(counter,adt_rl_index)=temp_adt**osm_adt_power_scale
            inputdata_rl_temp(counter,hdv_rl_index)=temp_hdv
            inputdata_int_rl_temp(counter,roadtype_rl_index)=temp_road_type
            inputdata_rl_temp(counter,speed_rl_index)=temp_speed
            inputdata_rl_temp(counter,width_rl_index)=temp_width
            inputdata_int_rl_temp(counter,nlanes_rl_index)=temp_nlanes
            inputdata_rl_temp(counter,x1_rl_index)=sub_nodes_x(j)
            inputdata_rl_temp(counter,x2_rl_index)=sub_nodes_x(j+loop_step)
            inputdata_rl_temp(counter,y1_rl_index)=sub_nodes_y(j)
            inputdata_rl_temp(counter,y2_rl_index)=sub_nodes_y(j+loop_step)
            !write(*,*) inputdata_int_rl(counter,id_rl_index),inputdata_rl(counter,x1_rl_index),inputdata_rl(counter,y2_rl_index)
            inputdata_rl_temp(counter,tunnel_length_rl_index)=temp_tunnel_length
            endif
        enddo
            !write(*,*) counter, inputdata_int_rl_temp(counter,major_index_rl_index)
        endif
    
    endif    
    enddo
    write(unit_logfile,'(a,i)') ' Number of counted road links= ', counter
    write(unit_logfile,'(a,i)') ' Number of road links allocated = ', n_roadlinks
    n_roadlinks=min(counter,n_roadlinks)
    close(unit_in,status='keep')
 
    !Allocate the arrays after reading in the number of roads
    if (allocated(inputdata_rl)) deallocate (inputdata_rl)
    if (allocated(inputdata_int_rl)) deallocate (inputdata_int_rl)
    allocate (inputdata_rl(n_roadlinks,num_var_rl))
    allocate (inputdata_int_rl(n_roadlinks,num_int_rl))

    inputdata_rl=inputdata_rl_temp(1:n_roadlinks,:)
    inputdata_int_rl=inputdata_int_rl_temp(1:n_roadlinks,:)
    
    if (road_data_in_latlon) then
        !Order is reversed in these files so they are reversed in the projection calls as well
        
        do i=1,n_roadlinks
            if  (projection_type.eq.UTM_projection_index) then
                call LL2UTM(1,utm_zone,inputdata_rl_temp(i,x1_rl_index),inputdata_rl_temp(i,y1_rl_index),inputdata_rl(i,y1_rl_index),inputdata_rl(i,x1_rl_index))
                call LL2UTM(1,utm_zone,inputdata_rl_temp(i,x2_rl_index),inputdata_rl_temp(i,y2_rl_index),inputdata_rl(i,y2_rl_index),inputdata_rl(i,x2_rl_index))
            elseif  (projection_type.eq.LTM_projection_index) then
                call LL2LTM(1,ltm_lon0,inputdata_rl_temp(i,x1_rl_index),inputdata_rl_temp(i,y1_rl_index),inputdata_rl(i,y1_rl_index),inputdata_rl(i,x1_rl_index))
                call LL2LTM(1,ltm_lon0,inputdata_rl_temp(i,x2_rl_index),inputdata_rl_temp(i,y2_rl_index),inputdata_rl(i,y2_rl_index),inputdata_rl(i,x2_rl_index))
            elseif (projection_type.eq.LCC_projection_index) then
            elseif (projection_type.eq.PS_projection_index) then
            elseif (projection_type.eq.LAEA_projection_index) then
                call LL2LAEA(inputdata_rl(i,x1_rl_index),inputdata_rl(i,y1_rl_index),inputdata_rl_temp(i,y1_rl_index),inputdata_rl_temp(i,x1_rl_index),projection_attributes)
                call LL2LAEA(inputdata_rl(i,x2_rl_index),inputdata_rl(i,y2_rl_index),inputdata_rl_temp(i,y2_rl_index),inputdata_rl_temp(i,x2_rl_index),projection_attributes)
            endif
        enddo
        
    endif

    if (num_multiple_roadlink_files.gt.1) then
        
        if (n_roadlinks.gt.0) then
        
        n_multi_roadlinks_new=n_multi_roadlinks+n_roadlinks
        
        !For the first roadloop 
        if (.not.first_road_link_file_read) then
            !Allocate the multi road array for the first time
            n_multi_roadlinks_new=n_roadlinks
            n_multi_roadlinks=n_multi_roadlinks_new
            if (.not.allocated(inputdata_rl_multi)) allocate (inputdata_rl_multi(n_multi_roadlinks,num_var_rl))
            if (.not.allocated(inputdata_int_rl_multi)) allocate (inputdata_int_rl_multi(n_multi_roadlinks,num_int_rl))
            inputdata_rl_multi(1:n_multi_roadlinks,:)=inputdata_rl(1:n_multi_roadlinks,:)
            inputdata_int_rl_multi(1:n_multi_roadlinks,:)=inputdata_int_rl(1:n_roadlinks,:)          
            first_road_link_file_read=.true.
       else
            
            !Allocate the new multi array with all roads so far
            allocate (inputdata_rl_multi_new(n_multi_roadlinks_new,num_var_rl))
            allocate (inputdata_int_rl_multi_new(n_multi_roadlinks_new,num_int_rl))

            !Place the current multi roads in the new multiroads
            !NOTE: inputdata_rl_multi is not allocated until later, debugging error but do not use this anyway!
            inputdata_rl_multi_new(1:n_multi_roadlinks,:)=inputdata_rl_multi(1:n_multi_roadlinks,:)
            inputdata_int_rl_multi_new(1:n_multi_roadlinks,:)=inputdata_int_rl_multi(1:n_roadlinks,:)
        
            !Place the last read road links in the new multiroads
            inputdata_rl_multi_new(n_multi_roadlinks+1:n_multi_roadlinks_new,:)=inputdata_rl(1:n_roadlinks,:)
            inputdata_int_rl_multi_new(n_multi_roadlinks+1:n_multi_roadlinks_new,:)=inputdata_int_rl(1:n_roadlinks,:)

            !Deallocate the old multi road arrays
            if (allocated(inputdata_rl_multi)) deallocate(inputdata_rl_multi)
            if (allocated(inputdata_int_rl_multi)) deallocate(inputdata_int_rl_multi)

            n_multi_roadlinks=n_multi_roadlinks_new
        
            !Allocate the multi road array
            allocate (inputdata_rl_multi(n_multi_roadlinks,num_var_rl))
            allocate (inputdata_int_rl_multi(n_multi_roadlinks,num_int_rl))
   
            !Put the new data in the old one
            inputdata_rl_multi=inputdata_rl_multi_new
            inputdata_int_rl_multi=inputdata_int_rl_multi_new

            !Deallocate the new multi road arrays
            if (allocated(inputdata_rl_multi_new)) deallocate(inputdata_rl_multi_new)
            if (allocated(inputdata_int_rl_multi_new)) deallocate(inputdata_int_rl_multi_new)

        endif
        
        endif

        !Deallocate the other arrays as well
        if (allocated(inputdata_rl)) deallocate(inputdata_rl)
        if (allocated(inputdata_int_rl)) deallocate(inputdata_int_rl)   
        if (allocated(inputdata_rl_temp)) deallocate (inputdata_rl_temp)
        if (allocated(inputdata_int_rl_temp)) deallocate (inputdata_int_rl_temp)
        if (allocated(valid_link_flag)) deallocate (valid_link_flag)

        write(unit_logfile,'(a,i)') ' Number of accumulated multi-road links used = ', n_multi_roadlinks_new
        
    endif
    
    !End the multiloop here
    enddo
    
    if (num_multiple_roadlink_files.gt.1) then
        allocate (inputdata_rl(n_multi_roadlinks,num_var_rl))
        allocate (inputdata_int_rl(n_multi_roadlinks,num_int_rl))
        inputdata_rl=inputdata_rl_multi
        inputdata_int_rl=inputdata_int_rl_multi
        if (allocated(inputdata_rl_multi)) deallocate(inputdata_rl_multi)
        if (allocated(inputdata_int_rl_multi)) deallocate(inputdata_int_rl_multi)
        n_roadlinks=n_multi_roadlinks_new
        write(unit_logfile,'(a,i)') ' Number of accumulated multi-road links used = ', n_multi_roadlinks
    endif
    
    

    if (allocated(inputdata_rl_temp)) deallocate (inputdata_rl_temp)
    if (allocated(inputdata_int_rl_temp)) deallocate (inputdata_int_rl_temp)
    
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
        call PROJ2LL(inputdata_rl(i,x1_rl_index),inputdata_rl(i,y1_rl_index),inputdata_rl(i,lon1_rl_index),inputdata_rl(i,lat1_rl_index),projection_attributes,projection_type)
        call PROJ2LL(inputdata_rl(i,x2_rl_index),inputdata_rl(i,y2_rl_index),inputdata_rl(i,lon2_rl_index),inputdata_rl(i,lat2_rl_index),projection_attributes,projection_type)
        call PROJ2LL(inputdata_rl(i,x0_rl_index),inputdata_rl(i,y0_rl_index),inputdata_rl(i,lon0_rl_index),inputdata_rl(i,lat0_rl_index),projection_attributes,projection_type)
        !call UTM2LL(utm_zone,inputdata_rl(i,y1_rl_index),inputdata_rl(i,x1_rl_index),inputdata_rl(i,lat1_rl_index),inputdata_rl(i,lon1_rl_index))
        !call UTM2LL(utm_zone,inputdata_rl(i,y2_rl_index),inputdata_rl(i,x2_rl_index),inputdata_rl(i,lat2_rl_index),inputdata_rl(i,lon2_rl_index))
        !call UTM2LL(utm_zone,inputdata_rl(i,y0_rl_index),inputdata_rl(i,x0_rl_index),inputdata_rl(i,lat0_rl_index),inputdata_rl(i,lon0_rl_index))
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
    
    double precision emission_date_number_start,emission_date_number
    double precision date_to_number
    
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
    
    !Check that start time and end time are covered in the emission data before progressing further
    format_temp='yyyymmddHH'
    call datestr_to_date(n_roadlink_emission_date_str,format_temp,emission_date_array)
    if (use_single_time_loop_flag) then
        time_index_temp=end_time_loop_index
    else
        time_index_temp=subgrid_dim(t_dim_index)
    endif
    
        t_match_index=0
        !write(*,*) shape(val_dim_nc)
        !do t=1,time_index_temp
        !    call number_to_date(val_dim_nc(t,time_dim_nc_index),date_array_temp,ref_year_EMEP)
        !    if (date_array_temp(1).eq.emission_date_array(1).and.date_array_temp(2).eq.emission_date_array(2).and. &
        !        date_array_temp(3).eq.emission_date_array(3).and.date_array_temp(4).eq.emission_date_array(4)) then
        !        t_match_index=t
        !    endif
        !    write(*,'(4i)') date_array_temp(1:4)
        !enddo
        emission_date_number_start=date_to_number(emission_date_array,ref_year_EMEP)
        call number_to_date(val_dim_nc(1,time_dim_nc_index),date_array_temp,ref_year_EMEP)
        do t=1,n_roadlink_emission_time
            emission_date_number=emission_date_number_start+(t-1)/24.
            call number_to_date(emission_date_number,emission_date_array,ref_year_EMEP)
            if (date_array_temp(1).eq.emission_date_array(1).and.date_array_temp(2).eq.emission_date_array(2).and. &
                date_array_temp(3).eq.emission_date_array(3).and.date_array_temp(4).eq.emission_date_array(4)) then
                t_match_index=t
            endif
            !write(*,'(4i)') emission_date_array(1:4)
        enddo
        if (t_match_index.eq.0) then
            write(unit_logfile,'(A,6i6)') 'ERROR: No starting date found in road emission data: ',emission_date_array
            stop
        else
            write(unit_logfile,'(A,6i6)') ' Road link emission starting date found. Index: ',t_match_index
            write(unit_logfile,'(A,6i6)') ' Road link emission starting date found. Index: ',t_match_index
        endif
        if (n_roadlink_emission_time-t_match_index+1.lt.time_index_temp) then
            write(unit_logfile,'(A,2i6)') 'ERROR: Not enough time data in road link emission files: ',n_roadlink_emission_time-t_match_index+1,time_index_temp
            stop
        else
            write(unit_logfile,'(A,6i6)') ' Road link emission ending date found. Index: ',t_match_index+time_index_temp-1
        endif 
    !endif
    
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
                if (j.le.n_pollutant_loop) then
                inputdata_rl_emissions(counter,1:time_index_temp,j)=inputdata_rl_temp(t_match_index:t_match_index+time_index_temp-1)
                endif
            enddo
            !write(*,*) counter,inputdata_rl_emissions(counter,10,1),inputdata_rl_emissions(counter,10,2)
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

    subroutine read_country_bounding_box_data
    !This routine reads in a file that provides information on the country bounding box
    !as well as reference to the filename used in the OSM files
    !Is also useful for other purposes but used here only for OSM file names
    
    use uEMEP_definitions
    
    implicit none

    integer exists
    integer i
    character(256) CNTR_ID
    character(256) OSM_country,Long_name
    real min_lat,min_lon,max_lat,max_lon,min_y_3035,min_x_3035,max_y_3035,max_x_3035
    real lon_grid_min,lat_grid_min,lon_grid_max,lat_grid_max
    real lon_grid_min2,lat_grid_min2,lon_grid_max2,lat_grid_max2
    real lon_new_min,lat_new_min,lon_new_max,lat_new_max
    integer unit_in
    character(256) temp_str
    integer count
    logical found_country
    real x_out(4),y_out(4)
    

    !Will fail at lon=-180
    
    !Read in replacement file
    pathfilename_boundingbox=trim(pathname_boundingbox)//trim(filename_boundingbox)
    
    !Test existence of the road link filename (2). If does not exist then use default
    inquire(file=trim(pathfilename_boundingbox),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Bounding box file does not exist: ', trim(pathfilename_boundingbox)
        stop
    endif
    
    if (trim(select_country_by_name).ne.'') then
        
    !Open the file for reading
    unit_in=20
    open(unit_in,file=pathfilename_boundingbox,access='sequential',status='old',readonly)  
    write(unit_logfile,'(a)') ' Opening Bounding box file  '//trim(filename_boundingbox)
    rewind(unit_in)
   
    !Skip over the header
    !Index	CNTR_ID	OSM_country	min_lat	min_lon	max_lat	max_lon	min_y_3035	min_x_3035	max_y_3035	max_x_3035		Long_name
    read(unit_in,*) temp_str

    !Read coordinates
    found_country=.false.
    do while (.not.eof(unit_in))
        read(unit_in,*) i,CNTR_ID,OSM_country,min_lon,min_lat,max_lon,max_lat,min_x_3035,min_y_3035,max_x_3035,max_y_3035,Long_name
        !write(*,*) i,trim(CNTR_ID),trim(OSM_country),min_lon,min_lat,max_lon,max_lat,min_x_3035,min_y_3035,max_x_3035,max_y_3035,trim(Long_name)
        
        if (index(trim(select_country_by_name),trim(CNTR_ID)).gt.0) then
            
            !Set the min and max lat and lon values for the current grid
            if  (projection_type.eq.UTM_projection_index) then
                call LL2UTM(1,utm_zone,min_lat,min_lon,y_out(1),x_out(1))
                call LL2UTM(1,utm_zone,max_lat,max_lon,y_out(2),x_out(2))
                call LL2UTM(1,utm_zone,max_lat,min_lon,y_out(3),x_out(3))
                call LL2UTM(1,utm_zone,min_lat,max_lon,y_out(4),x_out(4))
            elseif  (projection_type.eq.LTM_projection_index) then
                call LL2LTM(1,ltm_lon0,min_lat,min_lon,y_out(1),x_out(1))
                call LL2LTM(1,ltm_lon0,max_lat,max_lon,y_out(2),x_out(2))
                call LL2LTM(1,ltm_lon0,max_lat,min_lon,y_out(3),x_out(3))
                call LL2LTM(1,ltm_lon0,min_lat,max_lon,y_out(4),x_out(4))
            elseif (projection_type.eq.LAEA_projection_index) then
                call LL2LAEA(x_out(1),y_out(1),min_lon,min_lat,projection_attributes)
                call LL2LAEA(x_out(2),y_out(2),max_lon,max_lat,projection_attributes)
                call LL2LAEA(x_out(3),y_out(3),min_lon,max_lat,projection_attributes)
                call LL2LAEA(x_out(4),y_out(4),max_lon,min_lat,projection_attributes)
            endif
            
            subgrid_min(x_dim_index)=minval(x_out)
            subgrid_max(x_dim_index)=maxval(x_out)
            subgrid_min(y_dim_index)=minval(y_out)
            subgrid_max(y_dim_index)=maxval(y_out)
            
            !Snap to nearest 10 km
            subgrid_min(x_dim_index)=floor(subgrid_min(x_dim_index)/10000.)*10000
            subgrid_min(y_dim_index)=floor(subgrid_min(y_dim_index)/10000.)*10000
            subgrid_max(x_dim_index)=ceiling(subgrid_max(x_dim_index)/10000.)*10000
            subgrid_max(y_dim_index)=ceiling(subgrid_max(y_dim_index)/10000.)*10000
            
            write(unit_logfile,'(i,a)') i,' Setting grid for ID: '//trim(CNTR_ID)//'   OSM name: '//trim(OSM_country)//'   Name: '//trim(Long_name)
            write(unit_logfile,'(a,f12.1)') 'subgrid_min(x_dim_index)=',subgrid_min(x_dim_index)
            write(unit_logfile,'(a,f12.1)') 'subgrid_min(y_dim_index)=',subgrid_min(y_dim_index)
            write(unit_logfile,'(a,f12.1)') 'subgrid_max(x_dim_index)=',subgrid_max(x_dim_index)
            write(unit_logfile,'(a,f12.1)') 'subgrid_max(y_dim_index)=',subgrid_max(y_dim_index)
            found_country=.true.
        endif
    
        !Reset the initial subgrid as well, needed for EMEP and receptor selection
        init_subgrid_min(x_dim_index)=subgrid_min(x_dim_index)
        init_subgrid_min(y_dim_index)=subgrid_min(y_dim_index)
        init_subgrid_max(x_dim_index)=subgrid_max(x_dim_index)
        init_subgrid_max(y_dim_index)=subgrid_max(y_dim_index)

    enddo
    
10  close(unit_in)  
      
    if (.not.found_country) then
            write(unit_logfile,'(a)') ' No country with this ID found: '//trim(select_country_by_name)
    endif

    endif
    
    !Open the file for reading
    
    if (auto_select_OSM_country_flag) then

    !Set the min and max lat and lon values for the current grid
        
    call PROJ2LL(subgrid_min(x_dim_index),subgrid_min(y_dim_index),lon_grid_min,lat_grid_min,projection_attributes,projection_type)
    call PROJ2LL(subgrid_max(x_dim_index),subgrid_max(y_dim_index),lon_grid_max,lat_grid_max,projection_attributes,projection_type)
    call PROJ2LL(subgrid_min(x_dim_index),subgrid_max(y_dim_index),lon_grid_min2,lat_grid_max2,projection_attributes,projection_type)
    call PROJ2LL(subgrid_max(x_dim_index),subgrid_min(y_dim_index),lon_grid_max2,lat_grid_min2,projection_attributes,projection_type)
    lon_grid_max=max(lon_grid_max,lon_grid_max2)
    lon_grid_min=min(lon_grid_min,lon_grid_min2)
    lat_grid_max=max(lat_grid_max,lat_grid_max2)
    lat_grid_min=min(lat_grid_min,lat_grid_min2)
    
    !write(*,*) projection_attributes
    !write(*,*) projection_type
    
    !Test
    !lon_grid_max=5.;lat_grid_max=50.
    !call LL2LAEA(subgrid_max(x_dim_index),subgrid_max(y_dim_index),lon_grid_max,lat_grid_max,projection_attributes,projection_type)
    !call LAEA2LL(subgrid_max(x_dim_index),subgrid_max(y_dim_index),lon_grid_max,lat_grid_max,projection_attributes,projection_type)
    
    write(unit_logfile,'(a)') ' Opening Bounding box file  '//trim(filename_boundingbox)
    write(unit_logfile,'(a,2f12.6)') ' Lon (min,max)  ',lon_grid_min,lon_grid_max
    write(unit_logfile,'(a,2f12.6)') ' Lat (min,max)  ',lat_grid_min,lat_grid_max

    
    unit_in=20
    open(unit_in,file=pathfilename_boundingbox,access='sequential',status='old',readonly)  
    rewind(unit_in)
   
    !Skip over the header
    !Index	CNTR_ID	OSM_country	min_lat	min_lon	max_lat	max_lon	min_y_3035	min_x_3035	max_y_3035	max_x_3035		Long_name
    read(unit_in,*) temp_str

    !Read coordinates
    count=0
    filename_mrl=''
    do while (.not.eof(unit_in))
        read(unit_in,*) i,CNTR_ID,OSM_country,min_lon,min_lat,max_lon,max_lat,min_x_3035,min_y_3035,max_x_3035,max_y_3035,Long_name
        !write(*,*) i,trim(CNTR_ID),trim(OSM_country),min_lon,min_lat,max_lon,max_lat,min_x_3035,min_y_3035,max_x_3035,max_y_3035,trim(Long_name)
        
        !test the bounding box in lat lon coordinates
        lon_new_min=max(min_lon,lon_grid_min)
        lat_new_min=max(min_lat,lat_grid_min)
        lon_new_max=min(max_lon,lon_grid_max)
        lat_new_max=min(max_lat,lat_grid_max)
        
        if (lon_new_min.gt.lon_new_max.or.lat_new_min.gt.lat_new_max) then
            !No intersection so do nothing
        elseif (index('none',trim(OSM_country)).le.0) then
            !update and attribute the filename
            count=count+1
            if (count.gt.50) then
                write(unit_logfile,'(a)') ' Max files are 50. Stopping  '
                stop
            endif
                     
            filename_mrl(count)='Road_data_OSM_'//trim(OSM_country)//'_latlon.txt'
            write(unit_logfile,'(2i,a)') count,i, ' Including OSM file:  '//trim(filename_mrl(count))
        endif
        
    enddo
    
    if (count.eq.0) then
        write(unit_logfile,'(a)') ' No countries overlap this area. Setting OSM to default file. No traffic will be calculated'
        !stop
    else
        write(unit_logfile,'(a,i)') ' Specifying this many OSM road link files to be read',count        
        num_multiple_roadlink_files=count
    endif
    

20  close(unit_in)  
        

    endif
    
    !stop
    

    end subroutine read_country_bounding_box_data
