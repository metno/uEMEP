!uEMEP_read_agriculture_asi_data.f90
    
    subroutine uEMEP_read_agriculture_rivm_data
 
    use uEMEP_definitions
    
    implicit none
    
    integer i,j,k
    character(256) temp_name
    character(256) temp_str,temp_str1,temp_str2
    real temp_val
    integer unit_in
    integer exists
    integer count,index_val
    integer temp_int
    real ddlatitude,ddlongitude,totalnh3emission
    real y_agriculture,x_agriculture
    integer i_agriculture_index,j_agriculture_index
    real nh3emission_scale,nh3_gridsize
    integer i_range,j_range,i_file
    integer i_start,i_end,j_start,j_end
    logical, allocatable :: agriculture_emission_data_available(:,:)
    integer, allocatable :: agriculture_emission_emep_subgrid_count(:,:)
    integer ii,jj,iii,jjj
    integer source_index,subsource_index
    integer t
   
    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading agriculture rivm nh3 data  (uEMEP_read_agriculture_rivm_data)'
	write(unit_logfile,'(A)') '================================================================'
    
    projection_type=RDM_projection_index
    source_index=agriculture_index
    n_subsource(source_index)=2
    !unit_conversion(source_index)=1. !Converts kg to mg
    t=1
    
    !Internal check
    allocate (agriculture_emission_data_available(emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index)))
    
    proxy_emission_subgrid(:,:,source_index,:)=0.
    agriculture_emission_data_available=.false.
    
    pathfilename_agriculture(1)=trim(pathname_agriculture(1))//trim(filename_agriculture(1))
    pathfilename_agriculture(2)=trim(pathname_agriculture(2))//trim(filename_agriculture(2))
    
    !Test existence of the agricultural files. If does not exist then stop
    inquire(file=trim(pathfilename_agriculture(1)),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: agriculture file does not exist: ', trim(pathfilename_agriculture(1))
        stop
    endif
    inquire(file=trim(pathfilename_agriculture(2)),exist=exists)
    if (.not.exists) then
        write(unit_logfile,'(A,A)') ' ERROR: agriculture file does not exist: ', trim(pathfilename_agriculture(2))
        stop
    endif
    
    if (read_existing_grid_data(proxy_emission_file_index(source_index))) then
        do subsource_index=1,n_subsource(source_index)
        temp_name=trim(pathname_grid(proxy_emission_file_index(source_index)))//trim(filename_grid(proxy_emission_file_index(source_index)))//trim(subsource_str(subsource_index))//'_'//trim(file_tag)//'.asc'
        inquire(file=temp_name,exist=exists)
        if (.not.exists) then
            write(unit_logfile,*)'ERROR: '//trim(temp_name)//' does not exist.'
            stop
        endif
        call read_esri_ascii_file(unit_logfile,temp_name,emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index),emission_subgrid_delta(x_dim_index,source_index), &
                proxy_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index,subsource_index), &
                x_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index), &
                y_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index))
        enddo
        return
    endif
    
    
    !Open the file for reading
    do i_file=1,2
    unit_in=20
    open(unit_in,file=pathfilename_agriculture(i_file),access='sequential',status='old',readonly)  
    write(unit_logfile,'(a)') ' Opening agriculture file '//trim(pathfilename_agriculture(i_file))
    
    rewind(unit_in)

    !Read header x,y,emission
    read(unit_in,'(A)') temp_str
    write(unit_logfile,'(A)') trim(temp_str)
    count=0
    do while(.not.eof(unit_in))
        read(unit_in,'(A)') temp_str
        
        ddlatitude=0.;ddlongitude=0.;totalnh3emission=0.;
        !Extract the values in the temp_str
        index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
        read(temp_str1,*) x_agriculture
        !write (*,*) ddlatitude
        index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
        read(temp_str1,*) y_agriculture
        !write (*,*) ddlongitude
        index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
        !write(*,*) index_val,trim(temp_str1),trim(temp_str)
        !write(*,*) index_val,trim(temp_str)
        read(temp_str,*) totalnh3emission
        !write (*,*) trim(temp_str1),totalnh3emission
        !index_val=index(temp_str,';',back=.false.);temp_str1=temp_str(1:index_val-1);temp_str=temp_str(index_val+1:)
        !if (index_val.gt.1) read(temp_str1,*) totalparticulatematteremission
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
        !Convert lat lon to utm coords
 
        call RDM2LL(y_agriculture,x_agriculture,ddlatitude,ddlongitude)

        if (mod(count,10000).eq.0) write(unit_logfile,'(i,5f12.4)') count,y_agriculture,x_agriculture,ddlatitude,ddlongitude,totalnh3emission
        
        
        !Find the grid index it belongs to
        i_agriculture_index=1+floor((x_agriculture-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
        j_agriculture_index=1+floor((y_agriculture-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))
        !(x_subgrid(i,j)-subgrid_min(1))/+subgrid_delta(1)+1=i
        
        !Add to subgrid
        if (i_file.eq.1) then
            !Point sources
            subsource_index=1
            nh3emission_scale=1.
            i_start=i_agriculture_index;i_end=i_agriculture_index;
            j_start=j_agriculture_index;j_end=j_agriculture_index;
            h_emis(source_index,subsource_index)=5
        else
            !Area sources on 1 km grid
            if (n_subsource(source_index).eq.2) then
                subsource_index=2
            else
                subsource_index=1
            endif
            h_emis(source_index,subsource_index)=1.
           
            nh3_gridsize=1000.
            nh3emission_scale=emission_subgrid_delta(x_dim_index,source_index)*emission_subgrid_delta(y_dim_index,source_index)/(nh3_gridsize)**2
            i_range=floor(nh3_gridsize/emission_subgrid_delta(x_dim_index,source_index)/2.)
            j_range=floor(nh3_gridsize/emission_subgrid_delta(y_dim_index,source_index)/2.)
            if (mod(i_range,2).ne.0) then
                i_start=i_agriculture_index-i_range;i_end=i_agriculture_index+i_range-1
            else
                i_start=i_agriculture_index-i_range;i_end=i_agriculture_index+i_range         
            endif
            if (mod(j_range,2).ne.0) then
                j_start=j_agriculture_index-j_range;j_end=j_agriculture_index+j_range-1
            else
                j_start=j_agriculture_index-j_range;j_end=j_agriculture_index+j_range
            endif
            
         endif
        
        !Put emmissions into grids
        do i=i_start,i_end
        do j=j_start,j_end
        if (i.gt.0.and.i.le.emission_subgrid_dim(x_dim_index,source_index).and.j.gt.0.and.j.le.emission_subgrid_dim(y_dim_index,source_index)) then
            proxy_emission_subgrid(i,j,source_index,subsource_index)=proxy_emission_subgrid(i,j,source_index,subsource_index)+totalnh3emission*nh3emission_scale
            agriculture_emission_data_available(i,j)=.true.
        endif
        enddo
        enddo
             
    enddo
    write(unit_logfile,'(A,I)') ' Agriculture counts = ',count
      
    close(unit_in)
    enddo !file loop
    
    !After populating the grid with emission data (kg/yr) then fill all the other untouched grids with EMEP emission data (mg/m^2/yr).
    !Doesn't work. Would need to have a Nederlands mask instead of checking if emissions have been written or not.
    
    !First find out how many agriculture emission subgrids there are in each EMEP grid
    allocate (agriculture_emission_emep_subgrid_count(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))
    
    !This loop not necessary
    agriculture_emission_emep_subgrid_count=0
    do j=1,emission_subgrid_dim(y_dim_index,source_index)
    do i=1,emission_subgrid_dim(x_dim_index,source_index)
        
            !ii=crossreference_emission_to_target_subgrid(i,j,x_dim_index,source_index)
            !jj=crossreference_emission_to_target_subgrid(i,j,y_dim_index,source_index)
            !iii=crossreference_target_to_emep_subgrid(ii,jj,x_dim_index)
            !jjj=crossreference_target_to_emep_subgrid(ii,jj,y_dim_index)
            iii=crossreference_emission_to_emep_subgrid(i,j,x_dim_index,source_index)
            jjj=crossreference_emission_to_emep_subgrid(i,j,y_dim_index,source_index)
            agriculture_emission_emep_subgrid_count(iii,jjj)=agriculture_emission_emep_subgrid_count(iii,jjj)+1
    
    enddo
    enddo
    
    do j=1,emission_subgrid_dim(y_dim_index,source_index)
    do i=1,emission_subgrid_dim(x_dim_index,source_index)
        
        if (.not.agriculture_emission_data_available(i,j)) then
            !ii=crossreference_emission_to_target_subgrid(i,j,x_dim_index,source_index)
            !jj=crossreference_emission_to_target_subgrid(i,j,y_dim_index,source_index)
            !iii=crossreference_target_to_emep_subgrid(ii,jj,x_dim_index)
            !jjj=crossreference_target_to_emep_subgrid(ii,jj,y_dim_index)
            iii=crossreference_emission_to_emep_subgrid(i,j,x_dim_index,source_index)
            jjj=crossreference_emission_to_emep_subgrid(i,j,y_dim_index,source_index)
            if (agriculture_emission_emep_subgrid_count(iii,jjj).ne.0) then
                proxy_emission_subgrid(i,j,source_index,subsource_index)=var3d_nc(iii,jjj,1,emis_nc_index,agriculture_nc_index)/1.0e6 &
                    *emission_subgrid_delta(x_dim_index,source_index)*emission_subgrid_delta(y_dim_index,source_index)
                
            endif
        endif
    !write(*,*) agriculture_emission_emep_subgrid_count(iii,jjj)
    enddo
    enddo

    if (save_intermediate_files) then
    if (.not.read_existing_grid_data(proxy_emission_file_index(source_index))) then
        do subsource_index=1,n_subsource(source_index)
        temp_name=trim(pathname_grid(proxy_emission_file_index(source_index)))//trim(filename_grid(proxy_emission_file_index(source_index)))//trim(subsource_str(subsource_index))//'_'//trim(file_tag)//'.asc'
        write(unit_logfile,'(a)')'Writing to: '//trim(temp_name)
        call write_esri_ascii_file(unit_logfile,temp_name,emission_subgrid_dim(x_dim_index,source_index),emission_subgrid_dim(y_dim_index,source_index),emission_subgrid_delta(x_dim_index,source_index), &
                proxy_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index,subsource_index), &
                x_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index), &
                y_emission_subgrid(1:emission_subgrid_dim(x_dim_index,source_index),1:emission_subgrid_dim(y_dim_index,source_index),source_index))
        enddo
    endif
    endif
    
    deallocate (agriculture_emission_data_available)
    
    end subroutine uEMEP_read_agriculture_rivm_data
    
