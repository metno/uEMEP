module read_agriculture_asi_data

    use mod_lambert_projection, only: lb2lambert2_uEMEP, LL2PS_spherical
    use mod_area_interpolation, only: area_weighted_extended_interpolation_function

    implicit none
    private

    public :: uEMEP_read_agriculture_rivm_data, uEMEP_read_emission_rivm_data

contains

!uEMEP_read_agriculture_asi_data.f90

    subroutine uEMEP_read_agriculture_rivm_data

        use uEMEP_definitions
        use mod_rdm2ll, only: RDM2LL

        implicit none

        integer i,j
        character(256) temp_str,temp_str1
        integer unit_in
        integer exists
        integer count,index_val
        real ddlatitude,ddlongitude,totalnh3emission
        real y_agriculture,x_agriculture
        integer i_agriculture_index,j_agriculture_index
        real nh3emission_scale,nh3_gridsize(2)
        integer i_range,j_range,i_file
        integer i_start,i_end,j_start,j_end
        logical, allocatable :: agriculture_emission_data_available(:,:)
        integer, allocatable :: agriculture_emission_emep_subgrid_count(:,:)
        integer iii,jjj
        integer source_index,subsource_index
        integer t

        real x_temp(3,3),y_temp(3,3),z_temp(3,3)
        real temp_emission

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Reading agriculture rivm nh3 data  (uEMEP_read_agriculture_rivm_data)'
        write(unit_logfile,'(A)') '================================================================'

        if (.not.save_emissions_for_EMEP(agriculture_index)) then
            projection_type=RDM_projection_index
        endif

        source_index=agriculture_index
        n_subsource(source_index)=1
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
                if (save_emissions_for_EMEP(agriculture_index)) then
                    if (projection_type.eq.LL_projection_index) then
                        y_agriculture=ddlatitude
                        x_agriculture=ddlongitude
                    elseif (projection_type.eq.LCC_projection_index) then
                        call lb2lambert2_uEMEP(x_agriculture,y_agriculture,ddlongitude,ddlatitude,EMEP_projection_attributes)
                    elseif (projection_type.eq.PS_projection_index) then
                        call LL2PS_spherical(x_agriculture,y_agriculture,ddlongitude,ddlatitude,EMEP_projection_attributes)
                    endif
                endif

                if (mod(count,10000).eq.0) write(unit_logfile,'(i,5f12.4)') count,y_agriculture,x_agriculture,ddlatitude,ddlongitude,totalnh3emission


                !Find the grid index it belongs to
                i_agriculture_index=1+floor((x_agriculture-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
                j_agriculture_index=1+floor((y_agriculture-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))

                !Add to subgrid
                if (i_file.eq.1) then
                    !Point sources go straight into the emission grid
                    subsource_index=1
                    nh3emission_scale=1.
                    i_start=i_agriculture_index;i_end=i_agriculture_index;
                    j_start=j_agriculture_index;j_end=j_agriculture_index;
                    proxy_emission_subgrid(i_agriculture_index,j_agriculture_index,source_index,pollutant_loop_back_index(nh3_nc_index))= &
                        proxy_emission_subgrid(i_agriculture_index,j_agriculture_index,source_index,pollutant_loop_back_index(nh3_nc_index)) &
                        +totalnh3emission
                    agriculture_emission_data_available(i_agriculture_index,j_agriculture_index)=.true.
                else
                    !Area sources on 1 km grid
                    if (n_subsource(source_index).eq.2) then
                        subsource_index=2
                    else
                        subsource_index=1
                    endif
                    !h_emis(source_index,subsource_index)=3.

                    nh3_gridsize(:)=1000.
                    nh3emission_scale=emission_subgrid_delta(x_dim_index,source_index)*emission_subgrid_delta(y_dim_index,source_index)/(nh3_gridsize(x_dim_index)*nh3_gridsize(y_dim_index))

                    if (save_emissions_for_EMEP(agriculture_index)) then
                        if (projection_type.eq.LL_projection_index) then
                            !Approximate grid size in degrees
                            nh3_gridsize(x_dim_index)=nh3_gridsize(x_dim_index)/110570./cos(3.14159/180.*y_agriculture)
                            nh3_gridsize(y_dim_index)=nh3_gridsize(y_dim_index)/110570.
                            nh3emission_scale=emission_subgrid_delta(x_dim_index,source_index)*emission_subgrid_delta(y_dim_index,source_index) &
                                /(nh3_gridsize(x_dim_index)*nh3_gridsize(y_dim_index))
                        elseif (projection_type.eq.LCC_projection_index) then
                            !Do nothing
                        elseif (projection_type.eq.PS_projection_index) then
                            !Do nothing
                        endif
                    endif

                    !Set the loop range to be big enough
                    i_range=floor(nh3_gridsize(x_dim_index)/emission_subgrid_delta(x_dim_index,source_index)/2.)+1
                    j_range=floor(nh3_gridsize(y_dim_index)/emission_subgrid_delta(y_dim_index,source_index)/2.)+1
                    i_start=i_agriculture_index-i_range;i_end=i_agriculture_index+i_range
                    j_start=j_agriculture_index-j_range;j_end=j_agriculture_index+j_range
                    !write(*,*) i_range,j_range,nh3_gridsize(x_dim_index),emission_subgrid_delta(x_dim_index,source_index)
                    !write(*,*) i_start,i_end,j_start,j_end


                    !Put emmissions into emission grids
                    !Set up a 3 x 3 grid to interpolate from
                    !write(*,*) count
                    !write(*,*) x_agriculture,y_agriculture,nh3_gridsize(x_dim_index),nh3_gridsize(y_dim_index)
                    !write(*,*) x_emission_subgrid(i_start,j_start,source_index),y_emission_subgrid(i_start,j_start,source_index),emission_subgrid_delta(:,source_index)
                    x_temp(1,:)=x_agriculture-nh3_gridsize(x_dim_index);x_temp(2,:)=x_agriculture;x_temp(3,:)=x_agriculture+nh3_gridsize(x_dim_index);
                    y_temp(:,1)=y_agriculture-nh3_gridsize(y_dim_index);y_temp(:,2)=y_agriculture;y_temp(:,3)=y_agriculture+nh3_gridsize(y_dim_index);
                    z_temp(:,:)=0;z_temp(2,2)=totalnh3emission
                    do i=i_start,i_end
                        do j=j_start,j_end
                            if (i.gt.1.and.i.lt.emission_subgrid_dim(x_dim_index,source_index).and.j.gt.1.and.j.lt.emission_subgrid_dim(y_dim_index,source_index)) then
                                temp_emission=area_weighted_extended_interpolation_function(x_temp,y_temp,z_temp,3,3, &
                                    nh3_gridsize,x_emission_subgrid(i,j,source_index),y_emission_subgrid(i,j,source_index),emission_subgrid_delta(:,source_index))
                                proxy_emission_subgrid(i,j,source_index,pollutant_loop_back_index(nh3_nc_index))= &
                                    proxy_emission_subgrid(i,j,source_index,pollutant_loop_back_index(nh3_nc_index))+temp_emission*nh3emission_scale
                                if (proxy_emission_subgrid(i,j,source_index,pollutant_loop_back_index(nh3_nc_index)).gt.0.) agriculture_emission_data_available(i,j)=.true.

                                !write(*,*)i,j,proxy_emission_subgrid(i,j,source_index,pollutant_loop_back_index(nh3_nc_index)),agriculture_emission_data_available(i,j)!temp_emission*nh3emission_scale/totalnh3emission
                            endif
                        enddo
                    enddo

                endif
            enddo
            write(unit_logfile,'(A,I)') ' Agriculture counts = ',count

            close(unit_in)
        enddo !file loop

        !After populating the grid with emission data (kg/yr) then fill all the other untouched grids with EMEP emission data (mg/m^2/yr) or (mg/m^2/hr).
        !Doesn't work. Would need to have a Nederlands mask instead of checking if emissions have been written or not.
        !!Do not use!!
        if (1.eq.2) then

            if (.not.save_emissions_for_EMEP(agriculture_index)) then

                !First find out how many agriculture emission subgrids there are in each EMEP grid
                allocate (agriculture_emission_emep_subgrid_count(dim_length_nc(x_dim_nc_index),dim_length_nc(y_dim_nc_index)))

                !This loop determines if there is an emep grid associated with an emission grid which should always be true
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

                !Note, in current reading we have average emissions from EMEP not total
                do j=1,emission_subgrid_dim(y_dim_index,source_index)
                    do i=1,emission_subgrid_dim(x_dim_index,source_index)
                        !write(*,*) i,j,agriculture_emission_data_available(i,j)
                        if (.not.agriculture_emission_data_available(i,j)) then
                            !ii=crossreference_emission_to_target_subgrid(i,j,x_dim_index,source_index)
                            !jj=crossreference_emission_to_target_subgrid(i,j,y_dim_index,source_index)
                            !iii=crossreference_target_to_emep_subgrid(ii,jj,x_dim_index)
                            !jjj=crossreference_target_to_emep_subgrid(ii,jj,y_dim_index)
                            iii=crossreference_emission_to_emep_subgrid(i,j,x_dim_index,source_index)
                            jjj=crossreference_emission_to_emep_subgrid(i,j,y_dim_index,source_index)
                            if (agriculture_emission_emep_subgrid_count(iii,jjj).ne.0) then
                                !Convert mg/m2/hr to kg/emission subgrid. Need to put in time here for the hourly runs. Need to fix
                                proxy_emission_subgrid(i,j,source_index,pollutant_loop_back_index(nh3_index))=var3d_nc(iii,jjj,1,emis_nc_index,agriculture_nc_index,pollutant_loop_back_index(nh3_nc_index))/1.0e6 &
                                    *emission_subgrid_delta(x_dim_index,source_index)*emission_subgrid_delta(y_dim_index,source_index)*24.*365.
                                !If hourly data then convert to mean hour emission emissions
                                !if (hourly_calculations) proxy_emission_subgrid(i,j,source_index,pollutant_loop_back_index(nh3_index)) &
                                !    =proxy_emission_subgrid(i,j,source_index,pollutant_loop_back_index(nh3_index))/24./365.
                                !write(*,*) i,j,iii,jjj
                                !write(*,*) pollutant_loop_back_index(nh3_nc_index),proxy_emission_subgrid(i,j,source_index,pollutant_loop_back_index(nh3_index)),var3d_nc(iii,jjj,1,emis_nc_index,agriculture_nc_index,pollutant_loop_back_index(nh3_nc_index))
                            endif
                        endif
                        !write(*,*) agriculture_emission_emep_subgrid_count(iii,jjj)
                    enddo
                enddo

            endif

            deallocate (agriculture_emission_data_available)

        endif

    end subroutine uEMEP_read_agriculture_rivm_data


    subroutine uEMEP_read_emission_rivm_data

        use uEMEP_definitions

        implicit none

        integer i,j
        character(256) temp_str
        integer unit_in
        integer exists
        integer count
        real ddlatitude,ddlongitude,totalemission
        real y_emission,x_emission
        integer i_emission_index,j_emission_index
        real emission_scale
        integer i_file
        integer source_index

        character(256) component_str
        integer i_source
        real height
        integer snap, compound_nc_index
        integer i_pollutant
        integer, allocatable :: count_subgrid(:,:,:)
        real :: height_mean(n_source_index)=0
        integer :: count_mean(n_source_index)=0

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Reading emission rivm data  (uEMEP_read_emission_rivm_data)'
        write(unit_logfile,'(A)') '================================================================'

        !Set the projection to the dutch one
        !Should already have been specified. Not generic enough
        !projection_type=RDM_projection_index

        !Set the sources to be downscaled to 0
        do i_source=1,n_source_index
            if (calculate_source(i_source)) then
                proxy_emission_subgrid(:,:,i_source,:)=0.
                emission_properties_subgrid(:,:,emission_h_index,i_source)=0.
            endif
        enddo

        allocate (count_subgrid(emission_max_subgrid_dim(x_dim_index),emission_max_subgrid_dim(y_dim_index),n_source_index))
        count_subgrid=0

        !Read in as g/s and convert to ug/s which is used in the model
        emission_scale=1.0e6

        !Set filename, only 1 file used
        pathfilename_emission_rivm(1)=trim(pathname_emission_rivm(1))//trim(filename_emission_rivm(1))
        !pathfilename_emission_rivm(2)=trim(pathname_emission_rivm(2))//trim(filename_emission_rivm(2))

        !Test existence of the emission files. Only one file used. If does not exist then stop
        inquire(file=trim(pathfilename_emission_rivm(1)),exist=exists)
        if (.not.exists) then
            write(unit_logfile,'(A,A)') ' ERROR: emission RIVM file does not exist: ', trim(pathfilename_emission_rivm(1))
            stop
        endif
        !inquire(file=trim(pathfilename_emission_rivm(2)),exist=exists)
        !if (.not.exists) then
        !    write(unit_logfile,'(A,A)') ' ERROR: emission RIVM file does not exist: ', trim(pathfilename_emission_rivm(2))
        !    stop
        !endif


        !Open the files for reading
        i_file=1
        unit_in=20
        open(unit_in,file=pathfilename_emission_rivm(i_file),access='sequential',status='old',readonly)
        write(unit_logfile,'(a)') ' Opening emmision RIVM file '//trim(pathfilename_emission_rivm(i_file))

        rewind(unit_in)

        !Read header x,y,emission
        read(unit_in,'(A)') temp_str
        write(unit_logfile,'(A)') trim(temp_str)
        count=0
        do while(.not.eof(unit_in))
            !read(unit_in,'(A)') temp_str

            ddlatitude=0.;ddlongitude=0.;totalemission=0.;height=0;snap=0;component_str='';
            !read(unit_in,'(2f,e,f,i,a)') x_emission,y_emission,totalemission,height,snap,component_str
            read(unit_in,*) x_emission,y_emission,totalemission,height,snap,component_str
            !write(*,'(2f,e,f,i,a)') x_emission,y_emission,totalemission,height,snap,trim(component_str)

            compound_nc_index=0
            if (index(component_str,'NOx').ne.0) compound_nc_index=nox_nc_index
            if (index(component_str,'PM10').ne.0) compound_nc_index=pm10_nc_index
            if (index(component_str,'NH3').ne.0) compound_nc_index=nh3_nc_index
            if (index(component_str,'PM25').ne.0) compound_nc_index=pm25_nc_index

            !write(*,'(i,a)') compound_nc_index,trim(component_str)


            !Find the source sector in SNAP
            source_index=0
            do i_source=1,n_source_index

                if (calculate_source(i_source).and.uEMEP_to_EMEP_sector(i_source).eq.snap) then
                    source_index=i_source
                endif
            enddo


            !write(*,*) count
            count=count+1
            !Convert lat lon to utm coords

            !call RDM2LL(y_emission,x_emission,ddlatitude,ddlongitude)

            !if (mod(count,100000).eq.0) write(unit_logfile,'(i,2f12.2,e12.2,f12.2,i,2a)') count,x_emission,y_emission,totalemission,height,snap,'  ',trim(component_str)
            !write(*,*) source_index,compound_nc_index
            !Assumes the projection for the subgrid and the emissions are the same
            if (source_index.gt.0.and.compound_nc_index.gt.0) then
                !Find the subgrid index it belongs to
                i_emission_index=1+floor((x_emission-emission_subgrid_min(x_dim_index,source_index))/emission_subgrid_delta(x_dim_index,source_index))
                j_emission_index=1+floor((y_emission-emission_subgrid_min(y_dim_index,source_index))/emission_subgrid_delta(y_dim_index,source_index))
                !write(*,*) i_emission_index,j_emission_index,emission_subgrid_min(x_dim_index,source_index),emission_subgrid_min(y_dim_index,source_index)
                !Check that it is valid and add it to the subgrid
                if (i_emission_index.ge.1.and.i_emission_index.le.emission_subgrid_dim(x_dim_index,source_index).and.j_emission_index.ge.1.and.j_emission_index.le.emission_subgrid_dim(y_dim_index,source_index)) then
                    !Add to subgrid
                    proxy_emission_subgrid(i_emission_index,j_emission_index,source_index,pollutant_loop_back_index(compound_nc_index))= &
                        proxy_emission_subgrid(i_emission_index,j_emission_index,source_index,pollutant_loop_back_index(compound_nc_index)) &
                        +totalemission*emission_scale

                    emission_properties_subgrid(i_emission_index,j_emission_index,emission_h_index,source_index)= &
                        emission_properties_subgrid(i_emission_index,j_emission_index,emission_h_index,source_index)+height
                    !emission_properties_subgrid(i_emission_index,j_emission_index,emission_h_index,source_index)= height

                    count_subgrid(i_emission_index,j_emission_index,source_index)=count_subgrid(i_emission_index,j_emission_index,source_index)+1
                    !write(*,*) count_subgrid(i_emission_index,j_emission_index,source_index)

                endif

            endif

        enddo
        write(unit_logfile,'(A,I)') ' Emission counts = ',count

        close(unit_in)
        !enddo !file loop

        !Average the emission height sum and check output by looking at means
        do i_source=1,n_source_index
            if (calculate_source(i_source)) then
                do j=1,emission_subgrid_dim(y_dim_index,i_source)
                    do i=1,emission_subgrid_dim(x_dim_index,i_source)
                        if (count_subgrid(i,j,i_source).gt.0) then
                            emission_properties_subgrid(i,j,emission_h_index,i_source)=emission_properties_subgrid(i,j,emission_h_index,i_source)/count_subgrid(i,j,i_source)
                            height_mean(i_source)=emission_properties_subgrid(i,j,emission_h_index,i_source)+height_mean(i_source)
                            count_mean(i_source)=count_mean(i_source)+1
                        endif
                    enddo
                enddo

            endif
        enddo

        where (count_mean.gt.0) height_mean=height_mean/count_mean

        !Show results
        do i_source=1,n_source_index
            if (calculate_source(i_source)) then
                do i_pollutant=1,n_pollutant_loop
                    write(unit_logfile,'(A,A,A,ES10.2,A,f6.2)') 'Emission source ',trim(source_file_str(i_source))//' '//trim(pollutant_file_str(pollutant_loop_index(i_pollutant))),': Total RIVM emissions (ug/s)=', &
                        sum(proxy_emission_subgrid(1:emission_subgrid_dim(x_dim_index,i_source),1:emission_subgrid_dim(y_dim_index,i_source),i_source,i_pollutant)),'  Mean emission height (m)=',height_mean(i_source)
                enddo
            endif
        enddo

        deallocate (count_subgrid)


    end subroutine uEMEP_read_emission_rivm_data

end module read_agriculture_asi_data

