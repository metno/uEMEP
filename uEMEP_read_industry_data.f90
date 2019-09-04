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
    real, allocatable :: industry_emission(:,:)
    character(256), allocatable :: industry_code(:)
    integer industry_emission_year
    character(256) industry_emission_num,industry_emission_comp_str,industry_emission_unit
    real industry_emission_comp_val
    integer industry_number 
    integer i_industry_index,j_industry_index,source_index,i_pollutant
    real ratio_industry_pm25_to_pm10
    real x_industry,y_industry
    integer, allocatable :: count_subgrid(:,:,:)
    real, allocatable :: emission_height_subgrid(:,:,:)
    integer pollutant_count

    
    subsource_index=1
    source_index=industry_index
    
    !Set default pm25/pm10 ratio for when pm10 is given but not pm25. Value of 1 appropriate for exhaust but not for fugitive dust emissions
    ratio_industry_pm25_to_pm10=1.0
    


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
    allocate(industry_emission(n_industries,n_compound_index))
    
    
    rewind(unit_in)
    !Read header again
    read(unit_in,'(A)') temp_str

    do i=1,n_industries

        read(unit_in,*) industry_ref(i),industry_num(i),industry_lb_pos(i,2),industry_lb_pos(i,1),industry_xy_pos(i,2),industry_xy_pos(i,1),industry_height(i),industry_code(i)
        !write(unit_logfile,'(i12,a16,2f12.4,3f12.1,a16)' ) industry_ref(i),trim(industry_num(i)),industry_lb_pos(i,2),industry_lb_pos(i,1),industry_xy_pos(i,2),industry_xy_pos(i,1),industry_height(i),trim(industry_code(i))
        
    enddo
    
    close(unit_in)
    
    write(unit_logfile,'(a,i)') ' Number of industries: ',n_industries
    
    !Read in emission file and write emisions to industry array
    allocate (count_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index),n_pollutant_loop))
    allocate (emission_height_subgrid(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index),n_pollutant_loop))
    count_subgrid=0
    emission_height_subgrid=0.

    
    !Open the emission file for reading
    industry_emission=0.
    
    unit_in=20
    open(unit_in,file=pathfilename_industry(2),access='sequential',status='old',readonly)  
        write(unit_logfile,'(a)') ' Opening industry emission file '//trim(pathfilename_industry(2))
    
        rewind(unit_in)
    
        !Read header: År	AnleggNummer	Komponent	Samlet_mengde	Enhet
        read(unit_in,'(A)') temp_str

        do while(.not.eof(unit_in))
            read(unit_in,*) industry_emission_year,industry_emission_num,industry_emission_comp_str,industry_emission_comp_val,industry_emission_unit
            !write(unit_logfile,'(i12,2a16,f12.2,a16)' ) industry_emission_year,trim(industry_emission_num),trim(industry_emission_comp_str),industry_emission_comp_val,trim(industry_emission_unit)

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

            !If an industry is found then put the emissions in the industry_emission array
            if (industry_number.gt.0) then
                if  (trim(industry_emission_comp_str).eq.'nox') then
                    industry_emission(industry_number,nox_index)=industry_emission(industry_number,nox_index)+industry_emission_comp_val
                endif
                if  (trim(industry_emission_comp_str).eq.'pm10') then
                    industry_emission(industry_number,pm10_index)=industry_emission(industry_number,pm10_index)+industry_emission_comp_val
                endif
                if  (trim(industry_emission_comp_str).eq.'pm25') then
                    industry_emission(industry_number,pm25_index)=industry_emission(industry_number,pm25_index)+industry_emission_comp_val
                endif
        
            endif         
        enddo

    close(unit_in)    
 
    !Adjust pm10 and pm25 dependent on if they exist or not
    do i=1,n_industries
        !Set pm10 to pm25 in the case when it is less than pm25.
        if (industry_emission(i,pm10_index).lt.industry_emission(i,pm25_index)) then
            industry_emission(i,pm10_index)=industry_emission(i,pm25_index)
        endif
        !Set pm25 to pm10*ratio in cases where pm10 exists but pm25 does not
        if (industry_emission(i,pm25_index).eq.0.and.industry_emission(i,pm10_index).gt.0) then
            industry_emission(i,pm25_index)=industry_emission(i,pm10_index)*ratio_industry_pm25_to_pm10
        endif
    enddo

    !Initialise the industry emission arrays
    proxy_emission_subgrid(:,:,source_index,:)=0.
    emission_properties_subgrid(:,:,emission_h_index,source_index)=0.
    
    !Count the number of industry emission grid placements
    count=0
    
    do industry_number=1,n_industries
        
    !Convert lat lon to utm coords
    call LL2UTM(1,utm_zone,industry_lb_pos(industry_number,2),industry_lb_pos(industry_number,1),y_industry,x_industry)
    
    !Special case when saving emissions, convert to either latlon or lambert
    if (save_emissions_for_EMEP(industry_index)) then
        if (projection_type.eq.LL_projection_index) then
            x_industry=industry_lb_pos(industry_number,1)
            y_industry=industry_lb_pos(industry_number,2)        
        elseif (projection_type.eq.LCC_projection_index) then
            call lb2lambert2_uEMEP(x_industry,y_industry,industry_lb_pos(industry_number,1),industry_lb_pos(industry_number,2),EMEP_projection_attributes)
        endif
    endif
    
        
        !Find the grid index it belongs to
        i_industry_index=1+floor((x_industry-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
        j_industry_index=1+floor((y_industry-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))
        
        !Add to subgrid if it is within the subgrid range
        if (i_industry_index.ge.1.and.i_industry_index.le.emission_subgrid_dim(x_dim_index,source_index) &
            .and.j_industry_index.ge.1.and.j_industry_index.le.emission_subgrid_dim(y_dim_index,source_index)) then
            do i_pollutant=1,n_pollutant_loop
                if (pollutant_loop_index(i_pollutant).eq.nox_nc_index.and.industry_emission(industry_number,nox_index).gt.0) then
                    proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant)=proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant) &
                        +industry_emission(industry_number,nox_index)
                    emission_height_subgrid(i_industry_index,j_industry_index,i_pollutant)=emission_height_subgrid(i_industry_index,j_industry_index,i_pollutant)+industry_height(industry_number)*industry_emission(industry_number,nox_index)
                    count_subgrid(i_industry_index,j_industry_index,i_pollutant)=count_subgrid(i_industry_index,j_industry_index,i_pollutant)+1
                    count=count+1
                    !write(*,'(a,3i,f12.2)') 'Industry height nox: ',industry_number,i_industry_index,j_industry_index,industry_height(industry_number)
                endif
                if (pollutant_loop_index(i_pollutant).eq.pm10_nc_index.and.industry_emission(industry_number,pm10_index).gt.0) then
                    proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant)=proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant) &
                        +industry_emission(industry_number,pm10_index)
                    emission_height_subgrid(i_industry_index,j_industry_index,i_pollutant)=emission_height_subgrid(i_industry_index,j_industry_index,i_pollutant)+industry_height(industry_number)*industry_emission(industry_number,pm10_index)
                    count_subgrid(i_industry_index,j_industry_index,i_pollutant)=count_subgrid(i_industry_index,j_industry_index,i_pollutant)+1
                    count=count+1
                    !write(*,'(a,3i,f12.2)') 'Industry height pm10: ',industry_number,i_industry_index,j_industry_index,industry_height(industry_number)
               endif
                if (pollutant_loop_index(i_pollutant).eq.pm25_nc_index.and.industry_emission(industry_number,pm25_index).gt.0) then
                    proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant)=proxy_emission_subgrid(i_industry_index,j_industry_index,source_index,i_pollutant) &
                        +industry_emission(industry_number,pm25_index)
                    emission_height_subgrid(i_industry_index,j_industry_index,i_pollutant)=emission_height_subgrid(i_industry_index,j_industry_index,i_pollutant)+industry_height(industry_number)*industry_emission(industry_number,pm25_index)
                    count_subgrid(i_industry_index,j_industry_index,i_pollutant)=count_subgrid(i_industry_index,j_industry_index,i_pollutant)+1
                    count=count+1
                    !write(*,'(a,3i,f12.2)') 'Industry height pm25: ',industry_number,i_industry_index,j_industry_index,industry_height(industry_number)
                endif
                
            enddo                
                    
        endif
                
    enddo

    !Loop through the emission subgrid and take the average emission height
    !This is not weighted as it is averaged over all pollutants as well. Could do this but will not.
    !Probably need a pollutant dependent emission property
    do j=1,emission_subgrid_dim(y_dim_index,source_index)
    do i=1,emission_subgrid_dim(x_dim_index,source_index)
        pollutant_count=0
        do i_pollutant=1,n_pollutant_loop
        if (proxy_emission_subgrid(i,j,source_index,i_pollutant).gt.0) then
            emission_height_subgrid(i,j,i_pollutant)=emission_height_subgrid(i,j,i_pollutant)/proxy_emission_subgrid(i,j,source_index,i_pollutant)
            write(unit_logfile,'(2a,2i6,f12.2)') 'Emission height: ',trim(pollutant_file_str(pollutant_loop_index(i_pollutant))),i,j,emission_height_subgrid(i,j,i_pollutant)
            !Take the average of the pollutants
            emission_properties_subgrid(i,j,emission_h_index,source_index)=emission_properties_subgrid(i,j,emission_h_index,source_index)+emission_height_subgrid(i,j,i_pollutant)
            pollutant_count=pollutant_count+1
        endif
        enddo
        if (pollutant_count.gt.0) then
            emission_properties_subgrid(i,j,emission_h_index,source_index)=emission_properties_subgrid(i,j,emission_h_index,source_index)/pollutant_count
            write(unit_logfile,'(2a,2i6,f12.2)') 'Final emission height: ',trim('mean'),i,j,emission_properties_subgrid(i,j,emission_h_index,source_index)
        endif
        !Set the industry emission heights to that given in the config file if it is > 0
        if (h_emis(industry_index,1).ge.0) then
            emission_properties_subgrid(i,j,emission_h_index,source_index)=h_emis(industry_index,1)
        endif
    enddo
    enddo
    
    write(unit_logfile,'(A,I)') 'Industry counts = ',count
    do i_pollutant=1,n_pollutant_loop   
    write(unit_logfile,'(A,es12.3)') 'Total emission '//trim(pollutant_file_str(pollutant_loop_index(i_pollutant)))//' = ',sum(proxy_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index,i_pollutant))
    enddo    
    write(unit_logfile,'(A,f)') 'Average industry emission height for all industries = ',sum(industry_height(1:n_industries))/n_industries
    
    close(unit_in)
    
    deallocate(industry_ref)
    deallocate(industry_num)
    deallocate(industry_lb_pos)
    deallocate(industry_xy_pos)
    deallocate(industry_height)
    deallocate(industry_code)
    deallocate (count_subgrid)
    deallocate(industry_emission)

    
    end subroutine uEMEP_read_industry_data
    