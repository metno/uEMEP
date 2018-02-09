!uEMEP_read_receptor_data
!Reads in receptor positions and names
    
    subroutine uEMEP_read_receptor_data
    
    use uEMEP_definitions
 
    implicit none
    
    logical exists
    character(256) temp_str
    integer unit_in
    integer count
    integer :: use_region=1
    
    if (.not.use_receptor_positions_for_auto_subgrid_flag) then 
        return
    endif

    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading receptor positions (uEMEP_read_receptor_data)'
	write(unit_logfile,'(A)') '================================================================'
    
    pathfilename_receptor=trim(pathname_receptor)//trim(filename_receptor)
    !write(*,*) pathname_rl(2),filename_rl(2),pathfilename_rl(2)
    
    !Test existence of the road link filename (2). If does not exist then use default
    inquire(file=trim(pathfilename_receptor),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Receptor file does not exist: ', trim(pathfilename_receptor)
        stop
    endif

    !Open the file for reading
    unit_in=20
    open(unit_in,file=pathfilename_receptor,access='sequential',status='old',readonly)  
    write(unit_logfile,'(a)') ' Opening receptor file '//trim(pathfilename_receptor)
    
    rewind(unit_in)
    !call NXTDAT(unit_in,nxtdat_flag)
    !read the header to find out how many links there are
    read(unit_in,'(a)',ERR=20) temp_str
    k=0
    do while(.not.eof(unit_in))
        k=k+1
        read(unit_in,*,ERR=20) name_receptor(k,1),lon_receptor(k),lat_receptor(k)!,name_receptor(k,2)
        !write(*,*) trim(name_receptor(k,1)),lon_receptor(k),lat_receptor(k),trim(name_receptor(k,2))
    enddo
    
20  close(unit_in)
    
    n_receptor=k
    write(unit_logfile,'(a,i)') ' Number of receptor points = ', n_receptor
    
    !Convert to x,y positions
    do k=1,n_receptor
        if (projection_type.eq.RDM_projection_index) then
            !No conversion exists for RDM
        elseif (projection_type.eq.UTM_projection_index) then
            call LL2UTM(1,utm_zone,lat_receptor(k),lon_receptor(k),y_receptor(k),x_receptor(k))
        endif
    enddo
    
    !Find the target grid positions of the receptor points
    use_subgrid=.false.
    count=0
    do k=1,n_receptor
        
        i_receptor_subgrid(k)=1+floor((x_receptor(k)-subgrid_min(x_dim_index))/subgrid_delta(x_dim_index)+0.5)
        j_receptor_subgrid(k)=1+floor((y_receptor(k)-subgrid_min(y_dim_index))/subgrid_delta(y_dim_index)+0.5)

        !Set subgrid use or not. At grid and surrounding grids in case of interpolation later
        if (i_receptor_subgrid(k).gt.1.and.i_receptor_subgrid(k).lt.subgrid_dim(x_dim_index).and.j_receptor_subgrid(k).gt.1.and.j_receptor_subgrid(k).lt.subgrid_dim(y_dim_index)) then
            use_subgrid(i_receptor_subgrid(k)-use_region:i_receptor_subgrid(k)+use_region,j_receptor_subgrid(k)-use_region:j_receptor_subgrid(k)+use_region,:)=.true.
            !write(*,*) trim(name_receptor(k,1)),i_receptor_subgrid(k),j_receptor_subgrid(k)
            count=count+1
        endif
        
     enddo
     write(unit_logfile,'(a,i)') ' Number of receptor points available in region = ', count
     
     count=0
     do j=1,subgrid_dim(y_dim_index)
     do i=1,subgrid_dim(x_dim_index)
         if (use_subgrid(i,j,allsource_index)) count=count+1
     enddo
     enddo
     write(unit_logfile,'(a,i)') ' Number of subgrids to be calculated = ', count
        
    end subroutine uEMEP_read_receptor_data