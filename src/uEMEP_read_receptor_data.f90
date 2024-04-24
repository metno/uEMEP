module read_receptor_data

    use utility_functions, only: ll2utm, ll2ltm
    use mod_lambert_projection, only: LL2LAEA

    implicit none
    private

    public :: uEMEP_read_receptor_data, uEMEP_set_loop_receptor_grid, uEMEP_grid_receptor_data

contains

!uEMEP_read_receptor_data
!Reads in receptor positions and names

    subroutine uEMEP_read_receptor_data

        use uEMEP_definitions

        implicit none

        integer k
        logical exists
        character(256) temp_str
        integer unit_in
        integer count
        logical unique_receptor(n_receptor_max)
        integer kk

        use_receptor=.true.

        if (use_receptor_positions_for_auto_subgrid_flag.or.use_multiple_receptor_grids_flag.or.save_netcdf_receptor_flag) then
        else
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
            if (use_multiple_receptor_grids_flag) then
                write(unit_logfile,'(A,A)') ' ERROR: Receptor file does not exist. Cannot calculate: ', trim(pathfilename_receptor)
                stop
            else
                write(unit_logfile,'(A,A)') ' WARNING: Receptor file does not exist. Will not provide receptor output: ', trim(pathfilename_receptor)
                n_receptor=0
                n_receptor_in=n_receptor
                n_valid_receptor=0
                n_valid_receptor_in=n_valid_receptor
                return
            endif
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
        do while(.not.eof(unit_in).and.k.lt.n_receptor_max)
            k=k+1
            read(unit_in,*,ERR=20) name_receptor(k,1),lon_receptor(k),lat_receptor(k),height_receptor(k)!,name_receptor(k,2)
            !write(*,*) trim(name_receptor(k,1)),lon_receptor(k),lat_receptor(k),trim(name_receptor(k,2))
        enddo

20      close(unit_in)

        n_receptor=k
        write(unit_logfile,'(a,2i)') ' Number of receptor points and max allowable = ', n_receptor,n_receptor_max

        !Convert to x,y positions
        do k=1,n_receptor
            if (projection_type.eq.RDM_projection_index) then
                !No conversion exists for RDM
            elseif (projection_type.eq.UTM_projection_index) then
                call ll2utm(1,utm_zone,lat_receptor(k),lon_receptor(k),y_receptor(k),x_receptor(k))
            elseif (projection_type.eq.LTM_projection_index) then
                call ll2ltm(1,ltm_lon0,lat_receptor(k),lon_receptor(k),y_receptor(k),x_receptor(k))
            elseif (projection_type.eq.LAEA_projection_index) then
                call LL2LAEA(x_receptor(k),y_receptor(k),lon_receptor(k),lat_receptor(k),projection_attributes)
            endif
        enddo

        !Save the receptor data in the 'in' array as the other arrays can be changed
        n_receptor_in=n_receptor
        if (use_multiple_receptor_grids_flag) then
            n_receptor_in=n_receptor
            name_receptor_in=name_receptor
            lon_receptor_in=lon_receptor
            lat_receptor_in=lat_receptor
            x_receptor_in=x_receptor
            y_receptor_in=y_receptor
            height_receptor_in=height_receptor
        endif

        !Identify receptors within the initial subgrid region and only calculate these
        init_subgrid_dim(x_dim_index)=floor((init_subgrid_max(x_dim_index)-init_subgrid_min(x_dim_index))/init_subgrid_delta(x_dim_index))+1
        init_subgrid_dim(y_dim_index)=floor((init_subgrid_max(y_dim_index)-init_subgrid_min(y_dim_index))/init_subgrid_delta(y_dim_index))+1

        write(unit_logfile,'(a,2i)') ' Number of initial subgrid = ', init_subgrid_dim(x_dim_index),init_subgrid_dim(y_dim_index)
        write(unit_logfile,'(a,2f12.1)') ' Size of initial subgrid = ', init_subgrid_max(x_dim_index)-init_subgrid_min(x_dim_index),init_subgrid_max(y_dim_index)-init_subgrid_min(y_dim_index)

        !Remove identically named receptors
        count=0
        unique_receptor=.true.
        do k=1,n_receptor
            do kk=1,n_receptor
                if (trim(name_receptor(k,1)).eq.trim(name_receptor(kk,1)).and.unique_receptor(k).and.k.ne.kk) then
                    unique_receptor(kk)=.false.
                endif
            enddo
        enddo

        !Select receptors within the initial grid and with unique names
        do k=1,n_receptor

            i_receptor_subgrid(k)=1+floor((x_receptor(k)-init_subgrid_min(x_dim_index))/init_subgrid_delta(x_dim_index))
            j_receptor_subgrid(k)=1+floor((y_receptor(k)-init_subgrid_min(y_dim_index))/init_subgrid_delta(y_dim_index))
            !write(*,*) trim(name_receptor(k,1)),i_receptor_subgrid(k),j_receptor_subgrid(k)
            !Set subgrid use or not. At grid and surrounding grids in case of interpolation later
            if (i_receptor_subgrid(k).gt.1.and.i_receptor_subgrid(k).lt.init_subgrid_dim(x_dim_index).and.j_receptor_subgrid(k).gt.1.and.j_receptor_subgrid(k).lt.init_subgrid_dim(y_dim_index).and.unique_receptor(k)) then
                use_receptor(k)=.true.
                !write(*,*) trim(name_receptor(k,1)),i_receptor_subgrid(k),j_receptor_subgrid(k)
                count=count+1
                valid_receptor_index(count)=k
                valid_receptor_inverse_index(k)=count
                write(unit_logfile,'(a,a12,3f12.1)') ' Receptor and grid positions (x,y,h) = ',trim(name_receptor(k,1)), x_receptor(k)-init_subgrid_min(x_dim_index),y_receptor(k)-init_subgrid_min(y_dim_index),height_receptor(k)
            else
                use_receptor(k)=.false.
                valid_receptor_inverse_index(k)=0
            endif
        enddo
        n_valid_receptor=count
        n_valid_receptor_in=n_valid_receptor


        write(unit_logfile,'(a,i)') ' Total number of receptors to be calculated = ', n_valid_receptor

    end subroutine uEMEP_read_receptor_data


    subroutine uEMEP_grid_receptor_data

        use uEMEP_definitions

        implicit none

        integer i,j,k
        integer count
        logical use_receptor_temp
        !integer :: use_region=2 ! +/- number of grids to loop around so that receptor positions can be interpolated linearly

        if (use_receptor_positions_for_auto_subgrid_flag.or.use_multiple_receptor_grids_flag) then
        else
            return
        endif

        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Gridding receptor positions (uEMEP_grid_receptor_data)'
        write(unit_logfile,'(A)') '================================================================'


        !Find the target grid positions of the receptor points
        use_subgrid=.false.
        count=0

        do k=1,n_receptor
            !Always true when using use_multiple_receptor_grids_flag as this is inside the use_receptor loop
            if (use_multiple_receptor_grids_flag) then
                use_receptor_temp=.true.
            else
                use_receptor_temp=use_receptor(k)
            endif
            if (use_receptor_temp) then
                i_receptor_subgrid(k)=1+floor((x_receptor(k)-subgrid_min(x_dim_index))/subgrid_delta(x_dim_index))
                j_receptor_subgrid(k)=1+floor((y_receptor(k)-subgrid_min(y_dim_index))/subgrid_delta(y_dim_index))

                !write(*,*) 'HERE2: ',i_receptor_subgrid(k),j_receptor_subgrid(k)
                !Set subgrid use or not. At grid and surrounding grids in case of interpolation later
                if (i_receptor_subgrid(k).gt.use_receptor_region.and.i_receptor_subgrid(k).lt.subgrid_dim(x_dim_index)-use_receptor_region+1.and.j_receptor_subgrid(k).gt.use_receptor_region.and.j_receptor_subgrid(k).lt.subgrid_dim(y_dim_index)-use_receptor_region+1) then
                    use_subgrid(i_receptor_subgrid(k)-use_receptor_region:i_receptor_subgrid(k)+use_receptor_region,j_receptor_subgrid(k)-use_receptor_region:j_receptor_subgrid(k)+use_receptor_region,:)=.true.
                    !write(*,*) trim(name_receptor(k,1)),i_receptor_subgrid(k),j_receptor_subgrid(k)
                    count=count+1
                endif

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

    end subroutine uEMEP_grid_receptor_data


    subroutine uEMEP_set_loop_receptor_grid

        use uEMEP_definitions

        implicit none

        integer k
        integer count
        !integer :: use_region=2 ! +/- number of grids to loop around so that receptor positions can be interpolated linearly
        real x_ref,y_ref

        if (.not.use_multiple_receptor_grids_flag) then
            return
        endif

        !if (g_loop.eq.start_grid_loop_index) then
        write(unit_logfile,'(A)') ''
        write(unit_logfile,'(A)') '================================================================'
        write(unit_logfile,'(A)') 'Starting receptor loop (uEMEP_set_loop_receptor_grid)'
        write(unit_logfile,'(A)') '================================================================'
        ! endif

        k=1

        if (use_multiple_receptor_grids_flag) then
            name_receptor(k,:)=name_receptor_in(g_loop,:)
            lon_receptor(k)=lon_receptor_in(g_loop)
            lat_receptor(k)=lat_receptor_in(g_loop)
            x_receptor(k)=x_receptor_in(g_loop)
            y_receptor(k)=y_receptor_in(g_loop)
            height_receptor(k)=height_receptor_in(g_loop)
        endif

        !Set lowest left edge of subgrid receptor position would be in
        x_ref=(floor((x_receptor(k)-subgrid_receptor_offset(x_dim_index))/subgrid_delta(x_dim_index)+0.0))*subgrid_delta(x_dim_index)+subgrid_receptor_offset(x_dim_index)
        y_ref=(floor((y_receptor(k)-subgrid_receptor_offset(y_dim_index))/subgrid_delta(y_dim_index)+0.0))*subgrid_delta(y_dim_index)+subgrid_receptor_offset(y_dim_index)
        !Set limits
        subgrid_min(x_dim_index)=x_ref-subgrid_delta(x_dim_index)*(use_receptor_region)*1.0
        subgrid_min(y_dim_index)=y_ref-subgrid_delta(y_dim_index)*(use_receptor_region)*1.0
        subgrid_max(x_dim_index)=x_ref+subgrid_delta(x_dim_index)*(use_receptor_region+1)*1.0
        subgrid_max(y_dim_index)=y_ref+subgrid_delta(y_dim_index)*(use_receptor_region+1)*1.0

        subgrid_dim(x_dim_index)=floor((subgrid_max(x_dim_index)-subgrid_min(x_dim_index))/subgrid_delta(x_dim_index))+1
        subgrid_dim(y_dim_index)=floor((subgrid_max(y_dim_index)-subgrid_min(y_dim_index))/subgrid_delta(y_dim_index))+1

        z_rec=height_receptor(k)

        write(unit_logfile,'(a,i12,a)') ' Receptor loop number = ', g_loop,' '//trim(name_receptor(k,1))
        write(unit_logfile,'(a,4f12.1)') ' Receptor and grid positions (x,y) = ', x_receptor(k),x_ref,y_receptor(k),y_ref
        write(unit_logfile,'(a,2i)') ' Number of receptor subgrids = ', subgrid_dim(x_dim_index),subgrid_dim(y_dim_index)
        write(unit_logfile,'(a,f12.1)') ' Receptor height (m) = ', z_rec(allsource_index,1)

        !Find the target grid positions of the receptor points
        count=0
        !do k=1,n_receptor

        i_receptor_subgrid(k)=1+floor((x_receptor(k)-subgrid_min(x_dim_index))/subgrid_delta(x_dim_index))
        j_receptor_subgrid(k)=1+floor((y_receptor(k)-subgrid_min(y_dim_index))/subgrid_delta(y_dim_index))
        write(unit_logfile,'(a,2i)') ' Receptor subgrid index = ', i_receptor_subgrid(k),j_receptor_subgrid(k)

        !Set subgrid use or not. At grid and surrounding grids in case of interpolation later
        if (i_receptor_subgrid(k).gt.use_receptor_region.and.i_receptor_subgrid(k).lt.subgrid_dim(x_dim_index)-use_receptor_region+1.and.j_receptor_subgrid(k).gt.use_receptor_region.and.j_receptor_subgrid(k).lt.subgrid_dim(y_dim_index)-use_receptor_region+1) then
            use_subgrid(i_receptor_subgrid(k)-use_receptor_region:i_receptor_subgrid(k)+use_receptor_region,j_receptor_subgrid(k)-use_receptor_region:j_receptor_subgrid(k)+use_receptor_region,:)=.true.
            !write(*,*) trim(name_receptor(k,1)),i_receptor_subgrid(k),j_receptor_subgrid(k)
            count=count+1
        endif

        !enddo
        write(unit_logfile,'(a,i)') ' Number of receptor points available in region = ', count

        !count=0
        !do j=1,subgrid_dim(y_dim_index)
        !do i=1,subgrid_dim(x_dim_index)
        !    if (use_subgrid(i,j,allsource_index)) count=count+1
        !enddo
        !enddo
        !write(unit_logfile,'(a,i)') ' Number of subgrids to be calculated = ', count

    end subroutine uEMEP_set_loop_receptor_grid

end module read_receptor_data

