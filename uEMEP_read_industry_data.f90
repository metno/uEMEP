!uEMEP_read_industry_data.f90
!Reads in and places in subgrid industry emissions
!These emissions come from norgeutslipp.no but are reformatted for use in uEMEP

    subroutine uEMEP_read_industry_data

    use uEMEP_definitions
    
    implicit none
    
    integer i,j,k
    integer unit_in
    character(256) temp_str
    integer count
    integer n_industries
    logical exists
    integer subsource_index
    
    integer, allocatable :: industry_ref(:)
    character(256), allocatable :: industry_num(:)
    real, allocatable :: industry_lb_pos(:,:)
    real, allocatable :: industry_xy_pos(:,:)
    real, allocatable :: industry_height(:)
    character(256), allocatable :: industry_code(:)
    integer industry_emission_year
    character(256) industry_emission_num,industry_emission_comp_str,industry_emission_unit
    real industry_emission_comp_val
    integer industry_number 
    integer i_industry_index,j_industry_index,source_index,i_pollutant
    real ratio_industry_pm25_to_pm10
    real x_industry,y_industry
    integer, allocatable :: count_subgrid(:,:,:)

    
    subsource_index=1
    source_index=industry_index

    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading industry data  (uEMEP_read_industry_data)'
	write(unit_logfile,'(A)') '================================================================'

    pathfilename_industry(1)=trim(pathname_industry(1))//trim(filename_industry(1)) !Metadata
    pathfilename_industry(2)=trim(pathname_industry(2))//trim(filename_industry(2)) !Emission data
    
    !Test existence of the shipping filename. If does not exist then use default
    inquire(file=trim(pathfilename_industry(1)),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Industry metadata file does not exist: ', trim(pathfilename_industry(1))
        stop
    endif
    inquire(file=trim(pathfilename_industry(2)),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: Industry emission file does not exist: ', trim(pathfilename_industry(2))
        stop
    endif

    !Open the metadata file for reading
    unit_in=20
    open(unit_in,file=pathfilename_industry(1),access='sequential',status='old',readonly)  
    write(unit_logfile,'(a)') ' Opening industry metadata file '//trim(pathfilename_industry(1))
    
    rewind(unit_in)
    
    !Read header: AnleggIdRef	AnleggNummer	GeografiskLatitude	GeografiskLongitude	Utm33Nord	Utm33Ost	Height	NACEKode
    read(unit_in,'(A)') temp_str

    !write(*,*) trim(temp_str)
    
    !Count how many lines for allocation of arrays
    count=0
    do while(.not.eof(unit_in))
        read(unit_in,'(A)') temp_str
        count=count+1
    enddo
    
    !Allocate arrays
    n_industries=count
    allocate(industry_ref(n_industries))
    allocate(industry_num(n_industries))
    allocate(industry_lb_pos(n_industries,2))
    allocate(industry_xy_pos(n_industries,2))
    allocate(industry_height(n_industries))
    allocate(industry_code(n_industries))
    
    rewind(unit_in)
    !Read header again
    read(unit_in,'(A)') temp_str

    do i=1,n_industries

        read(unit_in,*) industry_ref(i),industry_num(i),industry_lb_pos(i,2),industry_lb_pos(i,1),industry_xy_pos(i,2),industry_xy_pos(i,1),industry_height(i),industry_code(i)
        !write(unit_logfile,'(i12,a16,2f12.4,3f12.1,a16)' ) industry_ref(i),trim(industry_num(i)),industry_lb_pos(i,2),industry_lb_pos(i,1),industry_xy_pos(i,2),industry_xy_pos(i,1),industry_height(i),trim(industry_code(i))
        
    enddo
    
    close(unit_in)
    
    write(unit_logfile,'(a,i)') ' Number of industries: ',n_industries
    
    !Read in emission file
    allocate (count_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index),n_pollutant_loop))
    count_subgrid=0

    !Open the metadata file for reading
    unit_in=20
    open(unit_in,file=pathfilename_industry(2),access='sequential',status='old',readonly)  
    write(unit_logfile,'(a)') ' Opening industry metadata file '//trim(pathfilename_industry(2))
    
    rewind(unit_in)
    
    !Read header: År	AnleggNummer	Komponent	Samlet_mengde	Enhet
    read(unit_in,'(A)') temp_str

    !write(*,*) trim(temp_str)
    
    proxy_emission_subgrid(:,:,source_index,:)=0.
    emission_properties_subgrid(:,:,emission_h_index,source_index)=0.
    ratio_industry_pm25_to_pm10=1.0
    
    !Count how many lines for allocation of arrays
    count=0
    do while(.not.eof(unit_in))
        read(unit_in,*) industry_emission_year,industry_emission_num,industry_emission_comp_str,industry_emission_comp_val,industry_emission_unit
        !write(unit_logfile,'(i12,2a16,f12.2,a16)' ) industry_emission_year,industry_emission_num,trim(industry_emission_comp_str),industry_emission_comp_val,industry_emission_unit
        !count=count+1

        !Find index for the industry
        industry_number=0
        do i=1,n_industries
            if (trim(industry_emission_num).eq.trim(industry_num(i))) then
                industry_number=i
                exit
            endif
        enddo
        
        if (industry_number.eq.0) then
            write(unit_logfile,*) 'No matching industry ID for the emissions: '//trim(industry_emission_num)
        else
            !write(unit_logfile,*) 'Matching industry ID found: '//trim(industry_emission_num)
        endif

        if (industry_number.gt.0) then
        
        !Convert lat lon to utm coords
        call LL2UTM(1,utm_zone,industry_lb_pos(industry_number,2),industry_lb_pos(industry_number,1),y_industry,x_industry)
        
        !Find the grid index it belongs to
        i_industry_index=1+floor((x_industry-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
        j_industry_index=1+floor((y_industry-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))
        
        !Add to subgrid       
        if (i_industry_index.ge.1.and.i_industry_index.le.emission_subgrid_dim(x_dim_index,source_index) &
            .and.j_industry_index.ge.1.and.j_industry_index.le.emission_subgrid_dim(y_dim_index,source_index)) then
            do i_pollutant=1,n_pollutant_loop
                if  (trim(industry_emission_comp_str).eq.'nox'.and.pollutant_loop_index(i_pollutant).eq.nox_nc_index) then
                    proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant)=proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant) &
                        +industry_emission_comp_val
                    count_subgrid(i_industry_index,j_industry_index,i_pollutant)=count_subgrid(i_industry_index,j_industry_index,i_pollutant)+1
                endif
                if (trim(industry_emission_comp_str).eq.'pm10'.and.(pollutant_loop_index(i_pollutant).eq.pm25_nc_index)) then
                    proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant)=proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant) &
                        +industry_emission_comp_val*ratio_industry_pm25_to_pm10
                    count_subgrid(i_industry_index,j_industry_index,i_pollutant)=count_subgrid(i_industry_index,j_industry_index,i_pollutant)+1    
                endif
                if (trim(industry_emission_comp_str).eq.'pm10'.and.pollutant_loop_index(i_pollutant).eq.pm10_nc_index) then
                    proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant)=proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant) &
                        +industry_emission_comp_val
                    count_subgrid(i_industry_index,j_industry_index,i_pollutant)=count_subgrid(i_industry_index,j_industry_index,i_pollutant)+1
                 endif
                if (trim(industry_emission_comp_str).eq.'pm25'.and.(pollutant_loop_index(i_pollutant).eq.pm25_nc_index)) then
                    proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant)=proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant) &
                        +industry_emission_comp_val
                    count_subgrid(i_industry_index,j_industry_index,i_pollutant)=count_subgrid(i_industry_index,j_industry_index,i_pollutant)+1
                endif
                if (trim(industry_emission_comp_str).eq.'pm25'.and.(pollutant_loop_index(i_pollutant).eq.pm10_nc_index)) then
                    proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant)=proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant) &
                        +industry_emission_comp_val
                    count_subgrid(i_industry_index,j_industry_index,i_pollutant)=count_subgrid(i_industry_index,j_industry_index,i_pollutant)+1
                endif
                
                !Needs to be changed to be an average if there are more in a grid, or a weighted average based on emission strengths from one of the compounds?
                !If this is done then emission heights can be different for different compounds
                !Chane the second index of h_emis to be pollutant, add an emission property subgrid dimension that is pollutant
                !Find out how high the chimneys need to be for any particular emission. Depends how many there are I guess
                !emission_properties_subgrid(i_industry_index,j_industry_index,emission_h_index,source_index)=industry_height(industry_number)
                emission_properties_subgrid(i_industry_index,j_industry_index,emission_h_index,source_index)=h_emis(industry_index,1)

            enddo
            count=count+1
            endif
        
        endif
        
    enddo

    write(unit_logfile,'(A,I)') 'Industry counts = ',count
    do i_pollutant=1,n_pollutant_loop   
    write(unit_logfile,'(A,es12.3)') 'Total emission '//trim(pollutant_file_str(pollutant_loop_index(i_pollutant)))//' = ',sum(proxy_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index,i_pollutant))
    enddo    
    
    
    close(unit_in)
    
    deallocate(industry_ref)
    deallocate(industry_num)
    deallocate(industry_lb_pos)
    deallocate(industry_xy_pos)
    deallocate(industry_height)
    deallocate(industry_code)
    deallocate (count_subgrid)

    
    end subroutine uEMEP_read_industry_data
    