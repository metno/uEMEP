!uEMEP_tiling_routines.f90
    
!Routines for calculating the positions, size and resolution of the tiling regions
    
    subroutine uEMEP_set_tile_grids

    use uEMEP_definitions
     
    implicit none

    integer n_tiles_population_classes,n_tiles_traffic_classes,n_tiles_shipping_classes
    parameter (n_tiles_population_classes=5,n_tiles_traffic_classes=5,n_tiles_shipping_classes=1)
    real limit_val_tile_population(n_tiles_population_classes)
   ! data limit_val_tile_population /0.,100.,1000.,5000.,10000./
    data limit_val_tile_population /0.,100.,1000.,10000.,100000./
    integer :: num_tiles_with_population(n_tiles_population_classes)=0
    real limit_val_tile_traffic(n_tiles_traffic_classes)
    data limit_val_tile_traffic /0.,1000.,10000.,100000.,1000000./
    integer :: num_tiles_with_traffic(n_tiles_traffic_classes)=0
    real limit_val_tile_shipping(n_tiles_shipping_classes)
    data limit_val_tile_shipping /0./
    integer :: num_tiles_with_shipping(n_tiles_shipping_classes)=0

    real tile_subgrid_delta(n_dim_index)
    real tile_subgrid_min(n_dim_index)
    real tile_subgrid_max(n_dim_index)
    integer tile_subgrid_dim(n_dim_index)
    
    real, allocatable :: tile_subgrid(:,:,:)
    real, allocatable :: x_tile_subgrid(:,:)
    real, allocatable :: y_tile_subgrid(:,:)
    real, allocatable :: lon_tile_subgrid(:,:)
    real, allocatable :: lat_tile_subgrid(:,:)
    real, allocatable :: crossreference_emission_to_tile_subgrid(:,:,:,:)
    real, allocatable :: crossreference_population_to_tile_subgrid(:,:,:)
    integer, allocatable :: tile_class_subgrid(:,:)
    integer, allocatable :: tile_municipality_subgrid(:,:)
    real, allocatable :: aggregated_tile_subgrid(:,:,:)
    real, allocatable :: aggregated_x_tile_subgrid(:,:,:)
    real, allocatable :: aggregated_y_tile_subgrid(:,:,:)
    integer, allocatable :: aggregated_tile_class_subgrid(:,:,:)
   
    integer i,j,i_class,i_tile,j_tile,i_source,k,l
    integer n_tile_index
    integer tile_population_index !,tile_municipality_index,tile_class_index
    
    character(256) temp_name,temp_str,temp_str1
    integer*8 ssb_id
    integer municipality_id
    real x_ssb,f_easting,ssb_dx,y_ssb,ssb_dy
    integer num_tiles_with_municipality
    logical exists
    integer unit_in
    integer count, index_val
    integer num_tile_classes(10)
    real resolution_tile_classes(10)
    character(8) count_str
    integer :: unit_tile=21
    real population_tile_scale
    integer reduce_grid_class
    real aggregation_tile_scale(10)
    integer n_aggregated_tiles
    parameter (n_aggregated_tiles=4)
    
    real aggregated_tile_subgrid_delta(n_dim_index,n_aggregated_tiles)
    real aggregated_tile_subgrid_min(n_dim_index,n_aggregated_tiles)
    real aggregated_tile_subgrid_max(n_dim_index,n_aggregated_tiles)
    integer aggregated_tile_subgrid_dim(n_dim_index,n_aggregated_tiles)
    integer count_tile_class(4)
    integer i_tile_class(4)
    integer j_tile_class(4)
    logical :: use_aggregated_tiling=.true.
    logical :: save_as_seperate_files=.true.
    integer sum_count,max_count
    integer count_class(n_aggregated_tiles,10)
    logical OK
    integer max_counter,zero_counter

    write(unit_logfile,'(A)') ''
    write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Calculating tile distribution and resolution (uEMEP_set_tile_grids)'
	write(unit_logfile,'(A)') '================================================================'

    !Set tile index values
    n_tile_index=n_source_index+1
    tile_population_index=n_source_index+1
    !tile_municipality_index=n_source_index+2
    !tile_class_index=n_source_index+3
    
    !Reduce size of grids at and above this value
    reduce_grid_class=6 !Set to 6 does not reduce
    
    !Specify the tiling region to cover all of Norway
    population_tile_scale=1. !For 10 km
    !population_tile_scale=3. !For 20 km
    population_tile_scale=0.625 !For 5 km, to give 250 people /km^2
    population_tile_scale=0.5 !For 5 km, to give 200 people /km^2
	tile_subgrid_delta(x_dim_index)=5000.
	tile_subgrid_delta(y_dim_index)=5000.
    tile_subgrid_min(x_dim_index)=-70000.-40000
    tile_subgrid_min(y_dim_index)=6440000.-40000
    tile_subgrid_max(x_dim_index)=1110000.+40000.
    tile_subgrid_max(y_dim_index)=7950000.+40000.
    
    limit_val_tile_population=limit_val_tile_population*population_tile_scale
    
    !Set all tile subgrids relative to the target subgrid
    tile_subgrid_dim(x_dim_index)=floor((tile_subgrid_max(x_dim_index)-tile_subgrid_min(x_dim_index))/tile_subgrid_delta(x_dim_index)) !New definition
    tile_subgrid_dim(y_dim_index)=floor((tile_subgrid_max(y_dim_index)-tile_subgrid_min(y_dim_index))/tile_subgrid_delta(y_dim_index)) !New definition
    tile_subgrid_dim(t_dim_index)=1 !Not used

    
    !aggregation_tile_scale
    if (n_aggregated_tiles.eq.4) then
        aggregation_tile_scale(1)=8.
        aggregation_tile_scale(2)=4.
        aggregation_tile_scale(3)=2.
        aggregation_tile_scale(4)=1.
    elseif (n_aggregated_tiles.eq.3) then
        aggregation_tile_scale(1)=4.
        aggregation_tile_scale(2)=2.
        aggregation_tile_scale(3)=1.
    endif
    
    do k=1,n_aggregated_tiles
        aggregated_tile_subgrid_min(:,k)=tile_subgrid_min(:)
        aggregated_tile_subgrid_max(:,k)=tile_subgrid_max(:)
        aggregated_tile_subgrid_delta(:,k)=tile_subgrid_delta(:)*aggregation_tile_scale(k)
        aggregated_tile_subgrid_dim(x_dim_index,k)=floor((aggregated_tile_subgrid_max(x_dim_index,k)-aggregated_tile_subgrid_min(x_dim_index,k))/aggregated_tile_subgrid_delta(x_dim_index,k))
        aggregated_tile_subgrid_dim(y_dim_index,k)=floor((aggregated_tile_subgrid_max(y_dim_index,k)-aggregated_tile_subgrid_min(y_dim_index,k))/aggregated_tile_subgrid_delta(y_dim_index,k))
    enddo
    
    
    !Allocate tile subgrids
    if (.not.allocated(tile_subgrid)) allocate (tile_subgrid(tile_subgrid_dim(x_dim_index),tile_subgrid_dim(y_dim_index),n_tile_index)) 
    if (.not.allocated(x_tile_subgrid)) allocate (x_tile_subgrid(tile_subgrid_dim(x_dim_index),tile_subgrid_dim(y_dim_index)))
    if (.not.allocated(y_tile_subgrid)) allocate (y_tile_subgrid(tile_subgrid_dim(x_dim_index),tile_subgrid_dim(y_dim_index)))
    if (.not.allocated(lon_tile_subgrid)) allocate (lon_tile_subgrid(tile_subgrid_dim(x_dim_index),tile_subgrid_dim(y_dim_index)))
    if (.not.allocated(lat_tile_subgrid)) allocate (lat_tile_subgrid(tile_subgrid_dim(x_dim_index),tile_subgrid_dim(y_dim_index)))
    if (.not.allocated(crossreference_emission_to_tile_subgrid)) allocate (crossreference_emission_to_tile_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),2,n_source_index)) 
    if (.not.allocated(crossreference_population_to_tile_subgrid)) allocate (crossreference_population_to_tile_subgrid(population_subgrid_dim(x_dim_index),population_subgrid_dim(y_dim_index),2)) 
    if (.not.allocated(tile_class_subgrid)) allocate (tile_class_subgrid(tile_subgrid_dim(x_dim_index),tile_subgrid_dim(y_dim_index)))
    if (.not.allocated(tile_municipality_subgrid)) allocate (tile_municipality_subgrid(tile_subgrid_dim(x_dim_index),tile_subgrid_dim(y_dim_index)))
    if (.not.allocated(aggregated_tile_subgrid)) allocate (aggregated_tile_subgrid(tile_subgrid_dim(x_dim_index),tile_subgrid_dim(y_dim_index),n_aggregated_tiles)) 
    if (.not.allocated(aggregated_x_tile_subgrid)) allocate (aggregated_x_tile_subgrid(tile_subgrid_dim(x_dim_index),tile_subgrid_dim(y_dim_index),n_aggregated_tiles))
    if (.not.allocated(aggregated_y_tile_subgrid)) allocate (aggregated_y_tile_subgrid(tile_subgrid_dim(x_dim_index),tile_subgrid_dim(y_dim_index),n_aggregated_tiles))
    if (.not.allocated(aggregated_tile_class_subgrid)) allocate (aggregated_tile_class_subgrid(tile_subgrid_dim(x_dim_index),tile_subgrid_dim(y_dim_index),n_aggregated_tiles))

    !Set to 0 all tile values
    tile_subgrid=0.
    tile_municipality_subgrid=0
    tile_class_subgrid=0
    
    !Read in SSB file containing gridded municipality ids
    ssb_dx=1000.
    ssb_dy=1000.
    f_easting=2.e6
    pathfilename_population(municipality_index)=trim(pathname_population(municipality_index))//trim(filename_population(municipality_index))
    !Test existence of the heating filename. If does not exist then use default
    inquire(file=trim(pathfilename_population(municipality_index)),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: SSB file with municipality IDs does not exist: ', trim(pathfilename_population(municipality_index))
        stop
    endif
    temp_name=pathfilename_population(municipality_index)
    !Open the file for reading
    unit_in=20
    open(unit_in,file=temp_name,access='sequential',status='old',readonly)  
    write(unit_logfile,'(a)') ' Opening SSB municipality file '//trim(temp_name)   
    rewind(unit_in)    
    !Read header SSBID0250M;kommunenum
    read(unit_in,'(A)') temp_str
    write(unit_logfile,'(A)') 'Header: '//trim(temp_str)
    count=0
    do while(.not.eof(unit_in))
        ssb_id=0;municipality_id=0
        !Read in file string    
        read(unit_in,'(A)') temp_str
        index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:);if (index_val.gt.1) read(temp_str1,*) ssb_id
        read(temp_str,*) municipality_id
        count=count+1
        !if (mod(count,100000).eq.0) write(*,*) count,ssb_id,municipality_id
        !Convert id to grid centre coordinates that are already in UTM33 for SSB data
        x_ssb=floor(ssb_id/10000000.)-f_easting+ssb_dx/2.
        y_ssb=mod(ssb_id,10000000)+ssb_dy/2.
        !Find the tile this ssb grid is in
        i_tile=1+floor((x_ssb-tile_subgrid_min(x_dim_index))/tile_subgrid_delta(x_dim_index)) !New definition
        j_tile=1+floor((y_ssb-tile_subgrid_min(y_dim_index))/tile_subgrid_delta(y_dim_index)) !New definition
        !Add the municipality id to give the grid a value
        if (municipality_id.gt.0) then
            tile_municipality_subgrid(i_tile,j_tile)=tile_municipality_subgrid(i_tile,j_tile)+municipality_id
        endif
    enddo
    close(unit_in)
        
    !Determine tile subgrid
    do j=1,tile_subgrid_dim(y_dim_index)
    do i=1,tile_subgrid_dim(x_dim_index)                 
        x_tile_subgrid(i,j)=tile_subgrid_min(x_dim_index)+tile_subgrid_delta(x_dim_index)*(i-0.5) !New definition
        y_tile_subgrid(i,j)=tile_subgrid_min(y_dim_index)+tile_subgrid_delta(y_dim_index)*(j-0.5) !New definition
        if (projection_type.eq.RDM_projection_index) then
            call RDM2LL(y_tile_subgrid(i,j),x_tile_subgrid(i,j),lat_tile_subgrid(i,j),lon_tile_subgrid(i,j))
        elseif (projection_type.eq.UTM_projection_index) then
            call UTM2LL(utm_zone,y_tile_subgrid(i,j),x_tile_subgrid(i,j),lat_tile_subgrid(i,j),lon_tile_subgrid(i,j))
        endif   
        !If the EMEP projection is lambert then set the proj coordinates to lambert, otherwise to lat-lon
        !if (EMEP_projection_type.eq.LCC_projection_index) then
        !    call lb2lambert2_uEMEP(xproj_tile_subgrid(i,j),yproj_tile_subgrid(i,j),lon_tile_subgrid(i,j),lat_tile_subgrid(i,j),EMEP_projection_attributes)
        !else
        !    xproj_tile_subgrid(i,j)=lon_tile_subgrid(i,j)
        !    yproj_tile_subgrid(i,j)=lat_tile_subgrid(i,j)            
        !endif
    enddo
    enddo

    !Determine the aggregated tile subgrids, in UTM only   
    do k=1,n_aggregated_tiles
    do j=1,aggregated_tile_subgrid_dim(y_dim_index,k)
    do i=1,aggregated_tile_subgrid_dim(x_dim_index,k)                 
        aggregated_x_tile_subgrid(i,j,k)=aggregated_tile_subgrid_min(x_dim_index,k)+aggregated_tile_subgrid_delta(x_dim_index,k)*(i-0.5)
        aggregated_y_tile_subgrid(i,j,k)=aggregated_tile_subgrid_min(y_dim_index,k)+aggregated_tile_subgrid_delta(y_dim_index,k)*(j-0.5)
    enddo
    enddo
    enddo
    
    !Create a cross reference grid
    do i_source=1,n_source_index
    if (calculate_source(i_source)) then
    do j=1,emission_subgrid_dim(y_dim_index,i_source)
    do i=1,emission_subgrid_dim(x_dim_index,i_source)
        crossreference_emission_to_tile_subgrid(i,j,x_dim_index,i_source)=1+floor((x_emission_subgrid(i,j,i_source)-tile_subgrid_min(x_dim_index))/tile_subgrid_delta(x_dim_index))
        crossreference_emission_to_tile_subgrid(i,j,y_dim_index,i_source)=1+floor((y_emission_subgrid(i,j,i_source)-tile_subgrid_min(y_dim_index))/tile_subgrid_delta(y_dim_index))      
    enddo
    enddo
    endif
    enddo
    
    do j=1,population_subgrid_dim(y_dim_index)
    do i=1,population_subgrid_dim(x_dim_index)
        crossreference_population_to_tile_subgrid(i,j,x_dim_index)=1+floor((x_population_subgrid(i,j)-tile_subgrid_min(x_dim_index))/tile_subgrid_delta(x_dim_index))
        crossreference_population_to_tile_subgrid(i,j,y_dim_index)=1+floor((y_population_subgrid(i,j)-tile_subgrid_min(y_dim_index))/tile_subgrid_delta(y_dim_index))      
    enddo
    enddo
    
    !Calculate the population within each tile
    do j=1,population_subgrid_dim(y_dim_index)
    do i=1,population_subgrid_dim(x_dim_index)
        i_tile=crossreference_population_to_tile_subgrid(i,j,x_dim_index)
        j_tile=crossreference_population_to_tile_subgrid(i,j,y_dim_index)
        tile_subgrid(i_tile,j_tile,tile_population_index)=tile_subgrid(i_tile,j_tile,tile_population_index)+population_subgrid(i,j,population_data_type)
    enddo
    enddo

    !Calculate the veh.km within each tile
    i_source=traffic_index
    do j=1,emission_subgrid_dim(y_dim_index,i_source)
    do i=1,emission_subgrid_dim(x_dim_index,i_source)
        i_tile=crossreference_emission_to_tile_subgrid(i,j,x_dim_index,i_source)
        j_tile=crossreference_emission_to_tile_subgrid(i,j,y_dim_index,i_source)
        tile_subgrid(i_tile,j_tile,traffic_index)=tile_subgrid(i_tile,j_tile,traffic_index)+sum(proxy_emission_subgrid(i,j,traffic_index,:))/1000. !ADT*km
    enddo
    enddo

    !Calculate the shipping emissions within each tile
    i_source=shipping_index
    do j=1,emission_subgrid_dim(y_dim_index,i_source)
    do i=1,emission_subgrid_dim(x_dim_index,i_source)
        i_tile=crossreference_emission_to_tile_subgrid(i,j,x_dim_index,i_source)
        j_tile=crossreference_emission_to_tile_subgrid(i,j,y_dim_index,i_source)
        tile_subgrid(i_tile,j_tile,i_source)=tile_subgrid(i_tile,j_tile,i_source)+sum(proxy_emission_subgrid(i,j,i_source,:)) !emission for the time period
    enddo
    enddo

    !Write summary results
    num_tiles_with_municipality=0
    do j_tile=1,tile_subgrid_dim(y_dim_index)
    do i_tile=1,tile_subgrid_dim(x_dim_index)
        if (tile_municipality_subgrid(i_tile,j_tile).gt.0) num_tiles_with_municipality=num_tiles_with_municipality+1
    enddo
    enddo
    write(unit_logfile,'(a,i)') 'MUNICIPALITY TILE: ',num_tiles_with_municipality

    !Write summary results
    num_tiles_with_population=0

    do j_tile=1,tile_subgrid_dim(y_dim_index)
    do i_tile=1,tile_subgrid_dim(x_dim_index)
        do i_class=1,n_tiles_population_classes-1
        if (tile_subgrid(i_tile,j_tile,tile_population_index).gt.limit_val_tile_population(i_class) &
            .and.tile_subgrid(i_tile,j_tile,tile_population_index).le.limit_val_tile_population(i_class+1)) num_tiles_with_population(i_class)=num_tiles_with_population(i_class)+1
        enddo
        i_class=n_tiles_population_classes
        if (tile_subgrid(i_tile,j_tile,tile_population_index).gt.limit_val_tile_population(i_class)) num_tiles_with_population(i_class)=num_tiles_with_population(i_class)+1
    enddo
    enddo
    do i_class=1,n_tiles_population_classes-1
    write(unit_logfile,'(a,i,f12.1,a,f12.1,i)') 'POPULATION TILE: ',i_class,limit_val_tile_population(i_class),' -',limit_val_tile_population(i_class+1),num_tiles_with_population(i_class)
    enddo
    i_class=n_tiles_population_classes
    write(unit_logfile,'(a,i,f12.1,a,a12,i)') 'POPULATION TILE: ',i_class,limit_val_tile_population(i_class),' <',' ',num_tiles_with_population(i_class)
    
    !Write summary results
    num_tiles_with_traffic=0

    do j_tile=1,tile_subgrid_dim(y_dim_index)
    do i_tile=1,tile_subgrid_dim(x_dim_index)
        do i_class=1,n_tiles_traffic_classes-1
        if (tile_subgrid(i_tile,j_tile,traffic_index).gt.limit_val_tile_traffic(i_class) &
            .and.tile_subgrid(i_tile,j_tile,traffic_index).le.limit_val_tile_traffic(i_class+1)) num_tiles_with_traffic(i_class)=num_tiles_with_traffic(i_class)+1
        enddo
        i_class=n_tiles_traffic_classes
        if (tile_subgrid(i_tile,j_tile,traffic_index).gt.limit_val_tile_traffic(i_class)) num_tiles_with_traffic(i_class)=num_tiles_with_traffic(i_class)+1
    enddo
    enddo
    do i_class=1,n_tiles_traffic_classes-1
    write(unit_logfile,'(a,i,f12.1,a,f12.1,i)') 'TRAFFIC TILE: ',i_class,limit_val_tile_traffic(i_class),' -',limit_val_tile_traffic(i_class+1),num_tiles_with_traffic(i_class)
    enddo
    i_class=n_tiles_traffic_classes
    write(unit_logfile,'(a,i,f12.1,a,a12,i)') 'TRAFFIC TILE: ',i_class,limit_val_tile_traffic(i_class),' <',' ',num_tiles_with_traffic(i_class)

    !Write summary results
    num_tiles_with_shipping=0
    do j_tile=1,tile_subgrid_dim(y_dim_index)
    do i_tile=1,tile_subgrid_dim(x_dim_index)
        if (tile_subgrid(i_tile,j_tile,shipping_index).gt.limit_val_tile_shipping(1)) num_tiles_with_shipping(1)=num_tiles_with_shipping(1)+1
    enddo
    enddo
    write(unit_logfile,'(a,i)') 'SHIPPING TILE: ',num_tiles_with_shipping(1)

    num_tile_classes=0
    !Class 1: 500 m. No emissions at all so interpolated, irrespective of population
    !Class 2: 250 m. Shipping emissions > 0 and Traffic < 1000 (1) or (population < 1000 and Traffic < 10000 (2))
    !Class 3: 125 m. Traffic < 10000 (2) or (population > 1000 (2) population < 5000 (2)
    !Class 4: 50 m. Traffic > 10000 (2) and population > 5000 (3)
    
    tile_class_subgrid=0
    
    do j_tile=1,tile_subgrid_dim(y_dim_index)
    do i_tile=1,tile_subgrid_dim(x_dim_index)
        
        !Only choose tiles that are part of a municipality
        if (tile_municipality_subgrid(i_tile,j_tile).gt.0) then
            tile_class_subgrid(i_tile,j_tile)=1
            if (tile_subgrid(i_tile,j_tile,shipping_index).le.limit_val_tile_shipping(1).and. &
                tile_subgrid(i_tile,j_tile,traffic_index).le.limit_val_tile_traffic(1).and. &
                tile_subgrid(i_tile,j_tile,tile_population_index).le.limit_val_tile_population(2)) then
                    tile_class_subgrid(i_tile,j_tile)=1 !No sources and population less than 100
            elseif (tile_subgrid(i_tile,j_tile,shipping_index).ge.limit_val_tile_shipping(1).and. &
                tile_subgrid(i_tile,j_tile,traffic_index).le.limit_val_tile_traffic(2).and. &
                tile_subgrid(i_tile,j_tile,tile_population_index).ge.limit_val_tile_population(1)) then
                    tile_class_subgrid(i_tile,j_tile)=2 !Little traffic but any shipping and any population
            elseif ((tile_subgrid(i_tile,j_tile,tile_population_index).gt.limit_val_tile_population(2).and. &
                tile_subgrid(i_tile,j_tile,tile_population_index).le.limit_val_tile_population(4)).or. &
                (tile_subgrid(i_tile,j_tile,traffic_index).ge.limit_val_tile_traffic(2).and. &
                tile_subgrid(i_tile,j_tile,tile_population_index).le.limit_val_tile_population(4))) then
                    tile_class_subgrid(i_tile,j_tile)=3 !Population from 1000 to 5000 or road traffic
            elseif (tile_subgrid(i_tile,j_tile,tile_population_index).gt.limit_val_tile_population(4).and. &
                tile_subgrid(i_tile,j_tile,tile_population_index).le.limit_val_tile_population(5)) then
                    tile_class_subgrid(i_tile,j_tile)=4 !Population > 10000
            elseif (tile_subgrid(i_tile,j_tile,tile_population_index).gt.limit_val_tile_population(5)) then
                    tile_class_subgrid(i_tile,j_tile)=4 !Population > 100000
            endif
            !if (tile_class_subgrid(i_tile,j_tile).eq.2) tile_class_subgrid(i_tile,j_tile)=3
            num_tile_classes(tile_class_subgrid(i_tile,j_tile))=num_tile_classes(tile_class_subgrid(i_tile,j_tile))+1
            
        endif
    enddo
    enddo
   
    write(unit_logfile,'(a,i)') 'TILES OF CLASS 1 (500m): ',num_tile_classes(1)
    write(unit_logfile,'(a,i)') 'TILES OF CLASS 2 (250m): ',num_tile_classes(2)
    write(unit_logfile,'(a,i)') 'TILES OF CLASS 3 (125m): ',num_tile_classes(3)
    if (reduce_grid_class.eq.4) then
        write(unit_logfile,'(a,i)') 'TILES OF CLASS 4 ( 50m): ',num_tile_classes(4)*4
    else
        write(unit_logfile,'(a,i)') 'TILES OF CLASS 4 ( 50m): ',num_tile_classes(4)
    endif
    if (reduce_grid_class.eq.5) then
        write(unit_logfile,'(a,i)') 'TILES OF CLASS 5 ( 25m): ',num_tile_classes(5)*4
    else
        write(unit_logfile,'(a,i)') 'TILES OF CLASS 5 ( 25m): ',num_tile_classes(5)
    endif

    resolution_tile_classes(1)=500.
    resolution_tile_classes(2)=250.
    resolution_tile_classes(3)=125.
    resolution_tile_classes(4)=50.
    resolution_tile_classes(5)=25.
   
    !Preallocate the smallest aggregated tiles with the calculated tile value
    aggregated_tile_class_subgrid=0
    aggregated_tile_class_subgrid(:,:,n_aggregated_tiles)=tile_class_subgrid
    !Aggregate tiles. Using classes of 2, 3 and 4 only
    do k=n_aggregated_tiles-1,1,-1
    do j_tile=1,aggregated_tile_subgrid_dim(y_dim_index,k)
    do i_tile=1,aggregated_tile_subgrid_dim(x_dim_index,k)
        
        count=0
        do j=1,aggregated_tile_subgrid_dim(y_dim_index,k+1)
        do i=1,aggregated_tile_subgrid_dim(x_dim_index,k+1)
            if (aggregated_x_tile_subgrid(i,j,k+1).ge.aggregated_x_tile_subgrid(i_tile,j_tile,k)-aggregated_tile_subgrid_delta(x_dim_index,k)/2. &
                .and.aggregated_x_tile_subgrid(i,j,k+1).lt.aggregated_x_tile_subgrid(i_tile,j_tile,k)+aggregated_tile_subgrid_delta(x_dim_index,k)/2. &
                .and.aggregated_y_tile_subgrid(i,j,k+1).ge.aggregated_y_tile_subgrid(i_tile,j_tile,k)-aggregated_tile_subgrid_delta(y_dim_index,k)/2. &
                .and.aggregated_y_tile_subgrid(i,j,k+1).lt.aggregated_y_tile_subgrid(i_tile,j_tile,k)+aggregated_tile_subgrid_delta(y_dim_index,k)/2.) then
                !.and.aggregated_tile_class_subgrid(i,j,k+1).ne.-1) then
                !This k+1 aggregated tile is within the k tile. Give it a class
                count=count+1
                count_tile_class(count)=aggregated_tile_class_subgrid(i,j,k+1)
                i_tile_class(count)=i
                j_tile_class(count)=j
                endif
        enddo
        enddo
             
        if (count.eq.4) then
            !4 x k+1 tiles have been found. Check if they are the same and not equal to the last tile value
            sum_count=sum(count_tile_class(1:count))
            max_count=maxval(count_tile_class(1:count))
            !if (count_tile_class(1)*count.le.sum_count.and.count_tile_class(2)*count.le.sum_count &
            !    .and.count_tile_class(3)*count.le.sum_count.and.count_tile_class(4)*count.le.sum_count.and.maxval(count_tile_class).lt.4) then
            OK=.true.
            max_counter=0
            zero_counter=0
            do l=1,count
                if (count_tile_class(l).eq.max_count) max_counter=max_counter+1
                if (count_tile_class(l).eq.0) zero_counter=zero_counter+1
            enddo
            
            do l=1,count
                if (((count_tile_class(l).le.max_count.and.count_tile_class(l).ge.0)).and.OK.and.max_count.lt.4 &
                    .and.((max_counter+zero_counter.ge.1.and.k.eq.3).or.(max_counter+zero_counter.ge.1.and.k.eq.2).or.(max_counter+zero_counter.ge.3.and.k.eq.1)) &
                    .and..not.(k.eq.1.and.max_count.eq.3)) then
                    OK=.true.
                else 
                    OK=.false.
                endif
            enddo
            
            if (OK) then
                
            !if (count_tile_class(1).le.max_count.and.count_tile_class(2).le.max_count &
            !    .and.count_tile_class(3).le.max_count.and.count_tile_class(4).le.max_count.and.max_count.lt.4) then
                
                !Allocate the class to the k tile
                aggregated_tile_class_subgrid(i_tile,j_tile,k)=max_count !sum_count/count
                !Remove the class from the k+1 tiles
                do l=1,count
                    aggregated_tile_class_subgrid(i_tile_class(l),j_tile_class(l),k+1)=-1
                enddo
                !write(*,*)  'Count is 4 and using'
            else
                !write(*,*)  'Count is 4 and not using'
                aggregated_tile_class_subgrid(i_tile,j_tile,k)=-1
            endif
            
        else
            !write(*,*) 'Can not find k+1 tiles',count
        endif
               
    enddo
    enddo
    enddo
    
    !Save results in a single file
    if (.not.use_aggregated_tiling) then
    temp_name=trim(pathname_tiles)//trim(filename_tiles)
    write(unit_logfile,'(a,2a)') 'Saving to: ', trim(temp_name)
    open(unit_tile,file=temp_name,access='sequential',status='unknown')
    count=0
    do j_tile=1,tile_subgrid_dim(y_dim_index)
    do i_tile=1,tile_subgrid_dim(x_dim_index)
        if (num_tile_classes(tile_class_subgrid(i_tile,j_tile)).gt.0) then
            if (tile_class_subgrid(i_tile,j_tile).lt.reduce_grid_class) then
                count=count+1
                write(unit_tile,'(a,i0.5)') 'tile_tag= ',count
	            write(unit_tile,'(a,f12.2)') 'subgrid_delta(x_dim_index)=',resolution_tile_classes(tile_class_subgrid(i_tile,j_tile))
	            write(unit_tile,'(a,f12.2)') 'subgrid_delta(y_dim_index)=',resolution_tile_classes(tile_class_subgrid(i_tile,j_tile))
	            write(unit_tile,'(a,f12.2)') 'subgrid_min(x_dim_index)=',x_tile_subgrid(i_tile,j_tile)-tile_subgrid_delta(x_dim_index)/2.
	            write(unit_tile,'(a,f12.2)') 'subgrid_min(y_dim_index)=',y_tile_subgrid(i_tile,j_tile)-tile_subgrid_delta(y_dim_index)/2.
	            write(unit_tile,'(a,f12.2)') 'subgrid_max(x_dim_index)=',x_tile_subgrid(i_tile,j_tile)+tile_subgrid_delta(x_dim_index)/2.
	            write(unit_tile,'(a,f12.2)') 'subgrid_max(y_dim_index)=',y_tile_subgrid(i_tile,j_tile)+tile_subgrid_delta(y_dim_index)/2.
            else
                !Divide into 4 lesser grids
                do j=0,1
                do i=0,1
                    count=count+1
                    write(unit_tile,'(a,i0.5)') 'tile_tag= ',count
	                write(unit_tile,'(a,f12.2)') 'subgrid_delta(x_dim_index)=',resolution_tile_classes(tile_class_subgrid(i_tile,j_tile))
	                write(unit_tile,'(a,f12.2)') 'subgrid_delta(y_dim_index)=',resolution_tile_classes(tile_class_subgrid(i_tile,j_tile))
	                write(unit_tile,'(a,f12.2)') 'subgrid_min(x_dim_index)=',x_tile_subgrid(i_tile,j_tile)-tile_subgrid_delta(x_dim_index)/2.+tile_subgrid_delta(x_dim_index)/2.*i
	                write(unit_tile,'(a,f12.2)') 'subgrid_min(y_dim_index)=',y_tile_subgrid(i_tile,j_tile)-tile_subgrid_delta(y_dim_index)/2.+tile_subgrid_delta(y_dim_index)/2.*j
	                write(unit_tile,'(a,f12.2)') 'subgrid_max(x_dim_index)=',x_tile_subgrid(i_tile,j_tile)+tile_subgrid_delta(x_dim_index)/2.+tile_subgrid_delta(x_dim_index)/2.*(i-1.)
	                write(unit_tile,'(a,f12.2)') 'subgrid_max(y_dim_index)=',y_tile_subgrid(i_tile,j_tile)+tile_subgrid_delta(y_dim_index)/2.+tile_subgrid_delta(y_dim_index)/2.*(j-1.)                   
                enddo
                enddo
            endif
        endif
    enddo
    enddo
    close(unit_tile)
 
    !Save results in multiple files
    if (save_as_seperate_files) then

    count=0
    do j_tile=1,tile_subgrid_dim(y_dim_index)
    do i_tile=1,tile_subgrid_dim(x_dim_index)
        if (num_tile_classes(tile_class_subgrid(i_tile,j_tile)).gt.0) then
            if (tile_class_subgrid(i_tile,j_tile).lt.reduce_grid_class) then
                count=count+1
                write(count_str,'(i8)') count
                temp_name=trim(pathname_tiles)//trim(ADJUSTL(count_str))//'_'//trim(filename_tiles)
                open(unit_tile,file=temp_name,access='sequential',status='unknown')
                write(unit_tile,'(a,i0.5)') 'tile_tag= ',count
	            write(unit_tile,'(a,f12.2)') 'subgrid_delta(x_dim_index)=',resolution_tile_classes(tile_class_subgrid(i_tile,j_tile))
	            write(unit_tile,'(a,f12.2)') 'subgrid_delta(y_dim_index)=',resolution_tile_classes(tile_class_subgrid(i_tile,j_tile))
	            write(unit_tile,'(a,f12.2)') 'subgrid_min(x_dim_index)=',x_tile_subgrid(i_tile,j_tile)-tile_subgrid_delta(x_dim_index)/2.
	            write(unit_tile,'(a,f12.2)') 'subgrid_min(y_dim_index)=',y_tile_subgrid(i_tile,j_tile)-tile_subgrid_delta(y_dim_index)/2.
	            write(unit_tile,'(a,f12.2)') 'subgrid_max(x_dim_index)=',x_tile_subgrid(i_tile,j_tile)+tile_subgrid_delta(x_dim_index)/2.
	            write(unit_tile,'(a,f12.2)') 'subgrid_max(y_dim_index)=',y_tile_subgrid(i_tile,j_tile)+tile_subgrid_delta(y_dim_index)/2.
                close(unit_tile)
            else
                !Divide into 4 lesser grids
                do j=0,1
                do i=0,1
                    count=count+1
                    write(count_str,'(i8)') count
                    temp_name=trim(pathname_tiles)//trim(ADJUSTL(count_str))//'_'//trim(filename_tiles)
                    open(unit_tile,file=temp_name,access='sequential',status='unknown')
                    write(unit_tile,'(a,i0.5)') 'tile_tag= ',count
	                write(unit_tile,'(a,f12.2)') 'subgrid_delta(x_dim_index)=',resolution_tile_classes(tile_class_subgrid(i_tile,j_tile))
	                write(unit_tile,'(a,f12.2)') 'subgrid_delta(y_dim_index)=',resolution_tile_classes(tile_class_subgrid(i_tile,j_tile))
	                write(unit_tile,'(a,f12.2)') 'subgrid_min(x_dim_index)=',x_tile_subgrid(i_tile,j_tile)-tile_subgrid_delta(x_dim_index)/2.+tile_subgrid_delta(x_dim_index)/2.*i
	                write(unit_tile,'(a,f12.2)') 'subgrid_min(y_dim_index)=',y_tile_subgrid(i_tile,j_tile)-tile_subgrid_delta(y_dim_index)/2.+tile_subgrid_delta(y_dim_index)/2.*j
	                write(unit_tile,'(a,f12.2)') 'subgrid_max(x_dim_index)=',x_tile_subgrid(i_tile,j_tile)+tile_subgrid_delta(x_dim_index)/2.+tile_subgrid_delta(x_dim_index)/2.*(i-1.)
	                write(unit_tile,'(a,f12.2)') 'subgrid_max(y_dim_index)=',y_tile_subgrid(i_tile,j_tile)+tile_subgrid_delta(y_dim_index)/2.+tile_subgrid_delta(y_dim_index)/2.*(j-1.)                   
                    close(unit_tile)
                enddo
                enddo
            endif

        endif
    enddo
    enddo
    
    endif
    
    write(unit_logfile,'(a,i)') 'Tiles before aggregation: ',count 

    endif
    
    if (use_aggregated_tiling) then
    temp_name=trim(pathname_tiles)//trim(filename_tiles)
    write(unit_logfile,'(a,2a)') 'Saving to: ', trim(temp_name)
    open(unit_tile,file=temp_name,access='sequential',status='unknown')
    count=0
    do k=1,n_aggregated_tiles
        count_class(k,:)=0
    do j_tile=1,aggregated_tile_subgrid_dim(y_dim_index,k)
    do i_tile=1,aggregated_tile_subgrid_dim(x_dim_index,k)
        !if (num_tile_classes(tile_class_subgrid(i_tile,j_tile)).gt.0) then
            if (aggregated_tile_class_subgrid(i_tile,j_tile,k).gt.0) then
                count=count+1
                write(unit_tile,'(a,i0.5)') 'tile_tag= ',count
	            write(unit_tile,'(a,f12.2)') 'subgrid_delta(x_dim_index)=',resolution_tile_classes(aggregated_tile_class_subgrid(i_tile,j_tile,k))
	            write(unit_tile,'(a,f12.2)') 'subgrid_delta(y_dim_index)=',resolution_tile_classes(aggregated_tile_class_subgrid(i_tile,j_tile,k))
	            write(unit_tile,'(a,f12.2)') 'subgrid_min(x_dim_index)=',aggregated_x_tile_subgrid(i_tile,j_tile,k)-aggregated_tile_subgrid_delta(x_dim_index,k)/2.
	            write(unit_tile,'(a,f12.2)') 'subgrid_min(y_dim_index)=',aggregated_y_tile_subgrid(i_tile,j_tile,k)-aggregated_tile_subgrid_delta(y_dim_index,k)/2.
	            write(unit_tile,'(a,f12.2)') 'subgrid_max(x_dim_index)=',aggregated_x_tile_subgrid(i_tile,j_tile,k)+aggregated_tile_subgrid_delta(x_dim_index,k)/2.
	            write(unit_tile,'(a,f12.2)') 'subgrid_max(y_dim_index)=',aggregated_y_tile_subgrid(i_tile,j_tile,k)+aggregated_tile_subgrid_delta(y_dim_index,k)/2.
                count_class(k,aggregated_tile_class_subgrid(i_tile,j_tile,k))=count_class(k,aggregated_tile_class_subgrid(i_tile,j_tile,k))+1
            endif
        !endif
    enddo
    enddo
        write(unit_logfile,'(a,i,f12.2)') 'TILE SIZE: ',k,aggregated_tile_subgrid_delta(x_dim_index,k)/1000.
        write(unit_logfile,'(a,i)') 'TILES OF CLASS 1 (500m): ',count_class(k,1)
        write(unit_logfile,'(a,i)') 'TILES OF CLASS 2 (250m): ',count_class(k,2)
        write(unit_logfile,'(a,i)') 'TILES OF CLASS 3 (125m): ',count_class(k,3)
        write(unit_logfile,'(a,i)') 'TILES OF CLASS 4 ( 50m): ',count_class(k,4)
        write(unit_logfile,'(a,i)') 'TILES OF CLASS 5 ( 25m): ',count_class(k,5)
    enddo
    close(unit_tile)
    
    !Save results in multiple files
    if (save_as_seperate_files) then
    
    count=0
    do k=1,n_aggregated_tiles
    do j_tile=1,tile_subgrid_dim(y_dim_index)
    do i_tile=1,tile_subgrid_dim(x_dim_index)
            if (aggregated_tile_class_subgrid(i_tile,j_tile,k).gt.0) then
                count=count+1
                write(count_str,'(i8)') count
                temp_name=trim(pathname_tiles)//trim(ADJUSTL(count_str))//'_'//trim(filename_tiles)
                open(unit_tile,file=temp_name,access='sequential',status='unknown')
                write(unit_tile,'(a,i0.5)') 'tile_tag= ',count
	            write(unit_tile,'(a,f12.2)') 'subgrid_delta(x_dim_index)=',resolution_tile_classes(aggregated_tile_class_subgrid(i_tile,j_tile,k))
	            write(unit_tile,'(a,f12.2)') 'subgrid_delta(y_dim_index)=',resolution_tile_classes(aggregated_tile_class_subgrid(i_tile,j_tile,k))
	            write(unit_tile,'(a,f12.2)') 'subgrid_min(x_dim_index)=',aggregated_x_tile_subgrid(i_tile,j_tile,k)-aggregated_tile_subgrid_delta(x_dim_index,k)/2.
	            write(unit_tile,'(a,f12.2)') 'subgrid_min(y_dim_index)=',aggregated_y_tile_subgrid(i_tile,j_tile,k)-aggregated_tile_subgrid_delta(y_dim_index,k)/2.
	            write(unit_tile,'(a,f12.2)') 'subgrid_max(x_dim_index)=',aggregated_x_tile_subgrid(i_tile,j_tile,k)+aggregated_tile_subgrid_delta(x_dim_index,k)/2.
	            write(unit_tile,'(a,f12.2)') 'subgrid_max(y_dim_index)=',aggregated_y_tile_subgrid(i_tile,j_tile,k)+aggregated_tile_subgrid_delta(y_dim_index,k)/2.
                close(unit_tile)
             endif
    enddo
    enddo
    enddo

    endif
    
    endif
    
    write(unit_logfile,'(a,i)') 'Tiles saved: ',count 
   
    write(unit_logfile,'(a)') ' Stopping after calculating tiles'
    stop
    
    end subroutine uEMEP_set_tile_grids