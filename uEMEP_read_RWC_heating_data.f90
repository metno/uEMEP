!uEMEP_read_RWC_heating_data.f90
!Reads in MetVed data in SSB format at 250 m
!Reads in HDD cdf file for all of Norway at 2.5km
!Allocates emissions to 250 m grid
    
    subroutine uEMEP_read_RWC_heating_data
    
    use uEMEP_definitions
    use netcdf

    implicit none
    
    logical exists
    integer source_index
    
    write(unit_logfile,'(A)') ''
	write(unit_logfile,'(A)') '================================================================'
	write(unit_logfile,'(A)') 'Reading RWC heating data  (uEMEP_read_RWC_heating_data)'
	write(unit_logfile,'(A)') '================================================================'
    
    source_index=heating_index
    n_subsource(source_index)=1

    proxy_emission_subgrid(:,:,source_index,:)=0.

    pathfilename_heating(RWC_heating_index)=trim(pathname_heating(RWC_heating_index))//trim(filename_heating(RWC_heating_index))
 
        !Test existence of the heating filename. If does not exist then stop
        inquire(file=trim(pathfilename_heating(RWC_heating_index)),exist=exists)
        if (.not.exists) then
            write(unit_logfile,'(A,A)') ' ERROR: SSB file does not exist: ', trim(pathfilename_heating(RWC_heating_index))
            stop
        endif

    
    end subroutine uEMEP_read_RWC_heating_data