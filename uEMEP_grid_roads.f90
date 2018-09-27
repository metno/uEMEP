!==========================================================================
!   NORTRIP_grid_roads.f90
!   Places line source proxy emissions in the traffic subgrid
!==========================================================================
    subroutine uEMEP_grid_roads
    
    use uEMEP_definitions

    implicit none
    
    integer i,j,k
    integer ro
    real, allocatable :: f_subgrid(:)
    real, allocatable :: adt_temp(:,:)
    real, allocatable :: adt_car_temp(:)
    real, allocatable :: adt_truck_temp(:)
    real x_subgrid_in(2),y_subgrid_in(2)
    real x_line_in(2),y_line_in(2),lat_line_in(2),lon_line_in(2)
    integer i_traffic_index(2),j_traffic_index(2)
    character(256) temp_name
    logical exists
    integer i_start,i_end,j_start,j_end
    integer source_index,t
    integer :: subsource_index=1
    integer i_roadlink_emission_compound(n_pollutant_loop)
    integer tt,ttt
    integer major_ro
    integer t_start_temp,t_end_temp
        
    integer i_pollutant
    
    !functions
    real line_fraction_in_grid_func
    real sigma0_traffic_func
    real minFF_traffic_func
        
    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Gridding road link proxy data (uEMEP_grid_roads)'
	write(unit_logfile,'(A)') '================================================================'

    allocate (f_subgrid(n_roadlinks))
    allocate (adt_temp(n_roadlinks,n_pollutant_loop))
    allocate (adt_car_temp(n_roadlinks))
    allocate (adt_truck_temp(n_roadlinks))
    !allocate (traffic_emission_subgrid(subgrid_dim(1),subgrid_dim(2)),n_emission_subgrid_index)
    
    source_index=traffic_index
    t=1
   
    proxy_emission_subgrid(:,:,source_index,:)=0.
    
    if (use_traffic_for_sigma0_flag) then
        emission_properties_subgrid(:,:,emission_sigy00_index,source_index)=0.
        emission_properties_subgrid(:,:,emission_sigz00_index,source_index)=0.
    endif
 !   emission_properties_subgrid(:,:,emission_minFF_index,source_index)=0.
 !   if (use_traffic_for_minFF_flag) then
 !       emission_properties_subgrid(:,:,emission_minFF_index,source_index)=0.
 !   endif
    !write(*,*) n_roadlinks
    !write(*,*) shape(inputdata_rl)
    
    !Possible to split the traffic source into different subsources at this point if necessary, e.g. light and heavy traffic
    !Here we weight the adt by the emission ratio and give an emission factor valid for cars
    adt_car_temp=inputdata_rl(1:n_roadlinks,adt_rl_index)*(1.-inputdata_rl(1:n_roadlinks,hdv_rl_index)/100.)
    adt_truck_temp=inputdata_rl(1:n_roadlinks,adt_rl_index)*inputdata_rl(1:n_roadlinks,hdv_rl_index)/100.
    do i_pollutant=1,n_pollutant_loop
    adt_temp(:,i_pollutant)=adt_car_temp+adt_truck_temp*ratio_truck_car_emission(pollutant_loop_index(i_pollutant))
    enddo
    
    !Calculate the pseudo traffic emissions in each grid
    write(unit_logfile,*)'Gridding traffic emission proxy data'
    
    !Convert from uEMEP to NORTRIP
    if (use_NORTRIP_emission_data) then

        do i_pollutant=1,n_pollutant_loop
            if (pollutant_loop_index(i_pollutant).eq.pm10_nc_index) then
                i_roadlink_emission_compound(i_pollutant)=1
            elseif (pollutant_loop_index(i_pollutant).eq.pm25_nc_index) then
                i_roadlink_emission_compound(i_pollutant)=2
            elseif (pollutant_loop_index(i_pollutant).eq.pmex_nc_index) then
                i_roadlink_emission_compound(i_pollutant)=3
            elseif (pollutant_loop_index(i_pollutant).eq.nox_nc_index) then
                i_roadlink_emission_compound(i_pollutant)=4
            else
                write(unit_logfile,*) 'STOPPING: No valid compound chosen for NORTRIP. Stopping uEMEP_grid_roads'
                stop
            endif             
        enddo

        emission_subgrid(:,:,:,traffic_index,:)=0.
        !Set the time to be used based on the time loop flag
        if (use_single_time_loop_flag) then
            t_start_temp=t_loop
            t_end_temp=t_loop
        else
            t_start_temp=1
            t_end_temp=subgrid_dim(t_dim_index)
        endif
        
    endif
    
    do ro=1,n_roadlinks

        x_line_in=inputdata_rl(ro,x1_rl_index:x2_rl_index)
        y_line_in=inputdata_rl(ro,y1_rl_index:y2_rl_index)
        
        !Convert to EMEP coordinates from UTM to lambert. No choices here
        if (save_emissions_for_EMEP(traffic_index)) then
            do i=1,2
                call UTM2LL(utm_zone,y_line_in(i),x_line_in(i),lat_line_in(i),lon_line_in(i))
                call lb2lambert2_uEMEP(x_line_in(i),y_line_in(i),lon_line_in(i),lat_line_in(i),EMEP_projection_attributes)
            enddo
            !write(*,*) x_line_in(1),y_line_in(1),lon_line_in(1),lat_line_in(1)
        endif
        
        
        i_traffic_index=1+floor((x_line_in-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
        j_traffic_index=1+floor((y_line_in-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))
 
        if ((i_traffic_index(1).ge.1.or.i_traffic_index(2).ge.1).and.(j_traffic_index(1).ge.1.or.j_traffic_index(2).ge.1).and. &
            (i_traffic_index(1).le.emission_subgrid_dim(x_dim_index,source_index).or.i_traffic_index(2).le.emission_subgrid_dim(x_dim_index,source_index)).and. &
            (j_traffic_index(1).le.emission_subgrid_dim(y_dim_index,source_index).or.j_traffic_index(2).le.emission_subgrid_dim(y_dim_index,source_index))) then
            
            !write(*,*) ro,i_traffic_index,j_traffic_index
            !Limit the loop if it is near the edge
            i_start=max(1,minval(i_traffic_index))
            i_end=min(emission_subgrid_dim(x_dim_index,source_index),maxval(i_traffic_index))
            j_start=max(1,minval(j_traffic_index))
            j_end=min(emission_subgrid_dim(y_dim_index,source_index),maxval(j_traffic_index))
 
            !if (i_end-i_start.gt.2.or.j_end-j_start.gt.2) write(*,*) ro,i_start,i_end,j_start,j_end
            
            do j=j_start,j_end
            do i=i_start,i_end

                x_subgrid_in(1)=x_emission_subgrid(i,j,source_index)-emission_subgrid_delta(x_dim_index,source_index)/2.
                x_subgrid_in(2)=x_emission_subgrid(i,j,source_index)+emission_subgrid_delta(x_dim_index,source_index)/2.
                y_subgrid_in(1)=y_emission_subgrid(i,j,source_index)-emission_subgrid_delta(y_dim_index,source_index)/2.
                y_subgrid_in(2)=y_emission_subgrid(i,j,source_index)+emission_subgrid_delta(y_dim_index,source_index)/2.
                
                f_subgrid(ro)=line_fraction_in_grid_func(x_subgrid_in,y_subgrid_in,x_line_in,y_line_in)
                
                !do subsource_index=1,n_subsource(source_index)
                    proxy_emission_subgrid(i,j,source_index,:)=proxy_emission_subgrid(i,j,source_index,:) &
                        +inputdata_rl(ro,length_rl_index)*f_subgrid(ro)*adt_temp(ro,:)
                    
                    !Put the temporally changing emissions straight into the emission subgrid
                    !Will not be overwritten in uEMEP_convert_proxy_to_emissions if use_NORTRIP_emission_data=true
                    if (use_NORTRIP_emission_data) then
                        major_ro=inputdata_int_rl(ro,major_index_rl_index)
                        do tt=t_start_temp,t_end_temp
                            !Convert from g/km/hour to ug/s/subgrid
                            if (use_single_time_loop_flag) then
                                ttt=t_loop
                                t=1
                            else
                                t=tt
                                ttt=tt
                            endif
                            !write(*,*) 'ro,major_ro,tt,t,ttt',ro,major_ro,tt,t,ttt
                            !write(*,*) 'emission_grid',shape(emission_subgrid)
                            !write(*,*) 'inputdata_rl_emissions',shape(inputdata_rl_emissions)

                            !if (t_loop.eq.2) stop
                            
                            do i_pollutant=1,n_pollutant_loop
                                emission_subgrid(i,j,t,source_index,i_pollutant)=emission_subgrid(i,j,t,source_index,i_pollutant)+ &
                               +inputdata_rl(ro,length_rl_index)*f_subgrid(ro)*inputdata_rl_emissions(major_ro,ttt,i_roadlink_emission_compound(i_pollutant)) &
                                *1.e6/1.e3/3600.
                            enddo
                            
                        enddo
                    endif

                    !Set the sigma values according to traffic speed and road width using proxy weighting
                    if (use_traffic_for_sigma0_flag.and..not.save_emissions_for_EMEP(traffic_index)) then
                        emission_properties_subgrid(i,j,emission_sigy00_index,source_index)=emission_properties_subgrid(i,j,emission_sigy00_index,source_index) &
                            +inputdata_rl(ro,length_rl_index)*f_subgrid(ro)*adt_temp(ro,1)*sqrt((inputdata_rl(ro,width_rl_index)/2.)**2+sigma0_traffic_func(inputdata_rl(ro,speed_rl_index))**2)
                        emission_properties_subgrid(i,j,emission_sigz00_index,source_index)=emission_properties_subgrid(i,j,emission_sigz00_index,source_index) &
                            +inputdata_rl(ro,length_rl_index)*f_subgrid(ro)*adt_temp(ro,1)*sigma0_traffic_func(inputdata_rl(ro,speed_rl_index))
                    endif
                    !Should not be used
!                    if (use_traffic_for_minFF_flag) then
!                        emission_properties_subgrid(i,j,emission_minFF_index,source_index,subsource_index)=emission_properties_subgrid(i,j,emission_minFF_index,source_index,subsource_index) &
!                            +inputdata_rl(ro,length_rl_index)*f_subgrid(ro)*adt_temp(ro)*minFF_traffic_func(inputdata_rl(ro,speed_rl_index),adt_temp(ro),inputdata_rl(ro,width_rl_index))
!                    endif
                    
                !enddo
                !write(*,*) ro,i,j,f_subgrid(ro)
                !write(*,*) ro,f_subgrid(ro),traffic_emission_subgrid(i,j,x_emission_subgrid_index),traffic_emission_subgrid(i,j,y_emission_subgrid_index),x_line_in,y_line_in
            enddo
            enddo
            
            !write(*,*) 'Gridding traffic emission',ro,' of ',n_roadlinks
            
        endif
           
        !if (mod(ro,10000).eq.0) write(*,*) 'Gridding traffic emission',ro,' of ',n_roadlinks    
    enddo
    
    !Set the road properties based on ADT weighting
    if (use_traffic_for_sigma0_flag) then
        emission_properties_subgrid(:,:,emission_sigy00_index,source_index)=emission_properties_subgrid(:,:,emission_sigy00_index,source_index)/proxy_emission_subgrid(:,:,source_index,1)
        emission_properties_subgrid(:,:,emission_sigz00_index,source_index)=emission_properties_subgrid(:,:,emission_sigz00_index,source_index)/proxy_emission_subgrid(:,:,source_index,1)
    endif
!    if (use_traffic_for_minFF_flag) then
!        emission_properties_subgrid(:,:,emission_minFF_index,source_index,:)=emission_properties_subgrid(:,:,emission_minFF_index,source_index,:)/proxy_emission_subgrid(:,:,source_index,:)
!     endif
    
    deallocate (f_subgrid)
    deallocate (adt_temp)
    deallocate (adt_car_temp)
    deallocate (adt_truck_temp)
    
 
    !Deallocate road link arrays after gridding but not when the external time step is used and not when the multiple receptor grids are used
    !because gridding roads is called again
    if (use_single_time_loop_flag.or.use_multiple_receptor_grids_flag) then
        !Do not deallocate because they will be used again
    else
        if (allocated(inputdata_rl)) deallocate(inputdata_rl)
        if (allocated(inputdata_int_rl)) deallocate(inputdata_int_rl)
        if (allocated(inputdata_rl_emissions)) deallocate(inputdata_rl_emissions)
    endif
    
    end subroutine uEMEP_grid_roads
    
    function sigma0_traffic_func(speed)
    implicit none
    real:: speed
    real:: sigma0_traffic_func
    real :: min_sigma=0.5
    real :: max_sigma=3.
    real :: min_speed=40.
    real :: max_speed=100.
    real :: gradient
    
        gradient=(max_sigma-min_sigma)/(max_speed-min_speed)
        sigma0_traffic_func=min(max(min_sigma+(speed-min_speed)*gradient,min_sigma),max_sigma)
        
    end function sigma0_traffic_func
    
    function minFF_traffic_func(speed,adt,width)
    implicit none
    real:: speed,adt,width
    real:: minFF_traffic_func
    real :: min_FF=0.0
    real :: max_FF=2.
    real :: min_speed=40.
    real :: max_speed=100.
    real :: gradient
    
       ! gradient=(max_FF-min_FF)/((max_speed-min_speed)*100000./25.)
       ! minFF_traffic_func=min(max(speed*adt/width*gradient,min_FF),max_FF)
        gradient=(max_FF-min_FF)/((max_speed-min_speed))
        minFF_traffic_func=min(max(speed*gradient,min_FF),max_FF)
        
    end function minFF_traffic_func
    
!==========================================================================
!   NORTRIP model line_fraction_in_grid_func
!==========================================================================
    function line_fraction_in_grid_func(x_grid,y_grid,x_line,y_line)
    
    implicit none
    
    !real, intent(in) :: x_grid_in(2),y_grid_in(2),x_line_in(2),y_line_in(2)
    real :: line_fraction_in_grid_func
    
    real :: x_grid(2),y_grid(2),x_line(2),y_line(2)
    real :: x_int(2),y_int(2)
    real :: length_line,length_int
    real :: dx,dy
    real :: x_temp,y_temp
    
    integer node,anti_node
    integer node_x_grid,node_y_grid
    integer :: n_intersection
    
    !Set to local variables
    !x_grid=x_grid_in
    !y_grid=y_grid_in
    !x_line=x_line_in
    !y_line=y_line_in
    
    !Set the initial fraction
    line_fraction_in_grid_func=0.
    
    !return

    !Check first for lines that cannot have an intersection. Will return 0
    if (x_line(1).lt.x_grid(1).and.x_line(2).lt.x_grid(1)) return
    if (x_line(1).ge.x_grid(2).and.x_line(2).ge.x_grid(2)) return
    if (y_line(1).lt.y_grid(1).and.y_line(2).lt.y_grid(1)) return
    if (y_line(1).ge.y_grid(2).and.y_line(2).ge.y_grid(2)) return

    !if ((x_line(1).lt.x_grid(1).and.x_line(2).lt.x_grid(1)).or.(x_line(1).ge.x_grid(2).and.x_line(2).ge.x_grid(2)).or.(y_line(1).lt.y_grid(1).and.y_line(2).lt.y_grid(1)).or.(y_line(1).ge.y_grid(2).and.y_line(2).ge.y_grid(2))) return

    !Set length of road link
    length_line=sqrt((x_line(1)-x_line(2))**2+(y_line(1)-y_line(2))**2)
    
    
    !Set the initial intercepts
    x_int(1:2)=x_line
    y_int(1:2)=y_line

    !write(*,*) x_grid(:),y_grid(:),x_line(:),y_line(:),length_line
    
    if (length_line.eq.0) return
  
    dx=MAXVAL(x_grid)-MINVAL(x_grid)
    dy=MAXVAL(y_grid)-MINVAL(y_grid)
       
    !Check for lines that are completely inside the grid
    if (x_line(1).ge.x_grid(1).and.x_line(2).ge.x_grid(1) &
        .and.x_line(1).lt.x_grid(2).and.x_line(2).lt.x_grid(2) &
        .and.y_line(1).ge.y_grid(1).and.y_line(2).ge.y_grid(1) &
        .and.y_line(1).lt.y_grid(2).and.y_line(2).lt.y_grid(2)) then
        line_fraction_in_grid_func=1.
        x_int=x_line
        y_int=y_line
        return
    endif
        
    !Check for lines with the one of the nodes within
    do node=1,2
        
        if (node.eq.1) anti_node=2
        if (node.eq.2) anti_node=1
                    
        if (x_line(node).ge.x_grid(1).and.x_line(node).lt.x_grid(2) &
            .and.y_line(node).ge.y_grid(1).and.y_line(node).lt.y_grid(2)) then
            !This node is in the grid
            !write(*,*) 'One node in grid'
            
            !Shift parallel and equal lines when they are on the grid edge
            if (x_line(node).eq.x_line(anti_node).and.x_line(node).eq.x_grid(1)) then
               x_line=x_line+dx*1e-6
            endif
            if (y_line(node).eq.y_line(anti_node).and.y_line(node).eq.y_grid(1)) then
               y_line=y_line+dy*1e-6
            endif
            
            !Can't intersect since it is parallel to the horizontal grid lines
            if (y_line(node).ne.y_line(anti_node)) then
            
                !Check intersection with the horizontal grid faces
                do node_y_grid=1,2
                    x_temp=x_line(node)+(y_grid(node_y_grid)-y_line(node))*(x_line(anti_node)-x_line(node))/(y_line(anti_node)-y_line(node))
                    y_temp=y_grid(node_y_grid)
                    !write(*,*) node,x_line(node),y_line(node),x_temp,y_temp,MINVAL(y_line),MAXVAL(y_line)
                    if (y_temp.ge.MINVAL(y_line).and.y_temp.le.MAXVAL(y_line).and.y_temp.ne.y_line(node).and.x_temp.ge.MINVAL(x_grid).and.x_temp.le.MAXVAL(x_grid)) then
                        y_int(anti_node)=y_grid(node_y_grid)
                        x_int(anti_node)=x_temp
                        x_int(node)=x_line(node)
                        y_int(node)=y_line(node)
                        length_int=sqrt((x_int(node)-x_int(anti_node))**2+(y_int(node)-y_int(anti_node))**2)
                        line_fraction_in_grid_func=length_int/length_line                    
                        return
                    endif
                enddo
            endif
            
            !Can't intersect since it is parallel with the vertical grid lines
            if (x_line(node).ne.x_line(anti_node)) then
                
                !Check intersection with the vertical grid faces
                do node_x_grid=1,2
                    y_temp=y_line(node)+(x_grid(node_x_grid)-x_line(node))*(y_line(anti_node)-y_line(node))/(x_line(anti_node)-x_line(node))
                    x_temp=x_grid(node_x_grid)
                    !write(*,*) node,x_line(node),y_line(node),x_temp,y_temp,MINVAL(x_line),MAXVAL(x_line)
                    if (x_temp.ge.MINVAL(x_line).and.x_temp.le.MAXVAL(x_line).and.x_temp.ne.x_line(node).and.y_temp.ge.MINVAL(y_grid).and.y_temp.le.MAXVAL(y_grid)) then
                        x_int(anti_node)=x_grid(node_x_grid)
                        y_int(anti_node)=y_temp
                        y_int(node)=y_line(node)
                        x_int(node)=x_line(node)
                        length_int=sqrt((x_int(node)-x_int(anti_node))**2+(y_int(node)-y_int(anti_node))**2)
                        line_fraction_in_grid_func=length_int/length_line                       
                        return
                    endif
                enddo
            endif
        endif
    
    enddo !node
    
    !Only posibility left is that both nodes are outside the grid
    !Find 2 intersections then
    n_intersection=0
    node=1
    anti_node=2
    if (y_line(node).ne.y_line(anti_node)) then !Can't intersect since it is parallel            
        do node_y_grid=1,2           
            !Check intersection with the horizontal grid faces
            x_temp=x_line(node)+(y_grid(node_y_grid)-y_line(node))*(x_line(anti_node)-x_line(node))/(y_line(anti_node)-y_line(node))                   
            y_temp=y_grid(node_y_grid)
            if (y_temp.ge.MINVAL(y_line).and.y_temp.le.MAXVAL(y_line).and.x_temp.ge.MINVAL(x_grid).and.x_temp.le.MAXVAL(x_grid).and.n_intersection.lt.2) then
                n_intersection=n_intersection+1
                y_int(n_intersection)=y_temp
                x_int(n_intersection)=x_temp
            endif           
        enddo
    endif
    if (x_line(node).ne.x_line(anti_node)) then !Can't intersect since it is parallel
        do node_x_grid=1,2
            y_temp=y_line(node)+(x_grid(node_x_grid)-x_line(node))*(y_line(anti_node)-y_line(node))/(x_line(anti_node)-x_line(node))
            x_temp=x_grid(node_x_grid)
            !Use y_temp.lt.MAXVAL(y_grid) incase it is in one of the corners
            if (x_temp.ge.MINVAL(x_line).and.x_temp.le.MAXVAL(x_line).and.y_temp.ge.MINVAL(y_grid).and.y_temp.lt.MAXVAL(y_grid).and.n_intersection.lt.2) then
                n_intersection=n_intersection+1
                x_int(n_intersection)=x_temp
                y_int(n_intersection)=y_temp
            endif               
        enddo
    endif
      
    if (n_intersection.eq.2) then
        length_int=sqrt((x_int(node)-x_int(anti_node))**2+(y_int(node)-y_int(anti_node))**2)
        line_fraction_in_grid_func=length_int/length_line
    endif
    
    end function line_fraction_in_grid_func

!==========================================================================
!   NORTRIP model save_gridded_lines_test_routine
!   THis routine used only for testing of gridding for line source
!   Not used in modelling
!==========================================================================
    subroutine save_gridded_lines_test_routine
    
    use uEMEP_definitions

    implicit none
    
    real :: x_grid(10,2),y_grid(10,2),x_line(50,2),y_line(50,2)
    real :: line(50,4),length_line(50)
    integer n_grid,n_line
    integer l,g
    real :: f(50)
    
    !functions
    real line_fraction_in_grid_func
    
    n_grid=2
    x_grid(1,:)=(/-1,1/)
    y_grid(1,:)=(/-1,1/)
    x_grid(2,:)=(/1,3/)
    y_grid(2,:)=(/1,3/)

    line(1,:)=(/.5,.5,1.,2./)!x1,y1,x2,y2
    line(2,:)=(/.5,0.,-2.,-0./)!x1,y1,x2,y2
    line(3,:)=(/0.,-0.2,-0.,-2./)!x1,y1,x2,y2
    line(4,:)=(/2.,3.,0.5,2./)!x1,y1,x2,y2
    line(5,:)=(/-2.,-3.,1.5,1.5/)!x1,y1,x2,y2
    line(6,:)=(/.7,-.9,.2,.7/)!x1,y1,x2,y2
    line(7,:)=(/-1.,-3.,-1.,+1./)!x1,y1,x2,y2
    line(8,:)=(/-.5,-1.,3.,-1./)!x1,y1,x2,y2
    line(9,:)=(/-.5,1.,3.,1./)!x1,y1,x2,y2
    line(10,:)=(/1.,-3.,1.,+0./)!x1,y1,x2,y2
    line(11,:)=(/-.7,-3.,-.7,+2./)!x1,y1,x2,y2
    line(12,:)=(/.5,1.5,1.5,.6/)!x1,y1,x2,y2
    line(13,:)=(/-1.,1.,1.,-1./)!x1,y1,x2,y2
    line(14,:)=(/-1.,-1.,1.,1./)!x1,y1,x2,y2
    line(15,:)=(/-1.,1.,1.,1./)!x1,y1,x2,y2
    line(16,:)=(/-1.5,.3,1.5,.3/)!x1,y1,x2,y2
    line(17,:)=(/-3.,2.,1.5,-3./)!x1,y1,x2,y2
    line(18,:)=(/-3.,-2.,1.,1./)!x1,y1,x2,y2
    line(19,:)=(/+3.,-2.,-1.,1./)!x1,y1,x2,y2
    line(20,:)=(/-3.,0.,1.,-1./)!x1,y1,x2,y2
    n_line=20

    write(*,*) 'input data'
    g=1
    do l=1,n_line
        x_line(l,1)=line(l,1)
        x_line(l,2)=line(l,3)
        y_line(l,1)=line(l,2)
        y_line(l,2)=line(l,4)
        length_line(l)=sqrt((x_line(l,1)-x_line(l,2))**2+(y_line(l,1)-y_line(l,2))**2)
        !write(*,*) g,l,x_grid(g,:),y_grid(g,:),x_line(l,:),y_line(l,:),length_line(l)
    enddo

    write(*,*) 'starting gridding'
    do g=1,n_grid
    do l=1,n_line
        
        f(l)=line_fraction_in_grid_func(x_grid(g,:),y_grid(g,:),x_line(l,:),y_line(l,:))
        write(*,*) g,l,f(l)
    enddo
    enddo
    
    stop
    end subroutine save_gridded_lines_test_routine
    
